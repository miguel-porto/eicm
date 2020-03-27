#' Generate EICM model following a Beta frequency distribution
#' 
#' Generates a randomly parameterized EICM model (predictors, environmental coefficients and species interactions),
#' ensuring that communities simulated with this model have a frequency distribution that matches the given Beta
#' distribution of frequencies as much as possible.
#' 
#' This function is useful for generating a realistic random model for simulation, i.e. a model that, when simulated, will
#' yield species communities with a distribution of frequencies akin of real communities: with most species being rare.
#' The generated coefficients are not assumed to follow any distribution, but are optionally shrinked so that
#' their values will remain "decent". The values of the environmental predictors are drawn from a gaussian distribution.
#' 
#' @param nspecies the number of species to generate.
#' @param nsamples the number of samples to generate.
#' @param nenv the number of environmental predictors to generate.
#' @param ninteractions the number of species interactions to generate.
#' @param shape1 the shape1 parameter of the Beta distribution.
#' @param shape2 the shape2 parameter of the Beta distribution.
#' @param nbins the number of histogram bins for comparing distributions.
#' @param nrepetitions the number of times to realize a model candidate to average their distributions.
#' @param shrinkage the shrinkage factor for generated coefficients, when computing fitness criterion.
#'        Ensures that the generated coefficients remain at plausible values.
#' @param bounds the allowed range for the coefficients \code{c(-bounds, +bounds)}.
#' @param swarm.size the swarm size of the particle swarm optimization.
#' @param maxit.stagnate the number of iterations without improvement until the optimization stops.
#' @return A EICM model of class \code{eicm}
#' @examples
#' \donttest{
#' # Generate model with 32 species, 30 species interactions and 2 environmental predictors
#' # for 500 samples with a frequency distribution following a Beta(1.5, 3)
#' model <- generateEICM(nspecies=32, nsamples=500, nenv=2, ninteractions=30,
#'   shape1=1.5, shape2=3)
#' 
#' # make one realization
#' data <- predict(model, nrepetitions=1)
#'
#' # plot frequency histogram: should follow a Beta distribution.
#' hist(apply(data, 2, sum) / nrow(data), breaks=seq(0, 1, by=0.1), xlim=c(0, 1),
#'   main="Frequency distribution of one realization", xlab="Frequency in samples",
#'   ylab="Number of species")
#' }
#' @export
generateEICM <- function(nspecies, nsamples, nenv, ninteractions, shape1, shape2,
	nbins=10, nrepetitions=5, shrinkage=2, bounds=10, swarm.size=floor((ninteractions + nspecies * (nenv + 1)) * 0.5),
	maxit.stagnate=150) {
	
	# create a target distribution of frequencies
	target <- lapply(1:100, function(i) stats::rbeta(nspecies, shape1, shape2))
	t1 <- sapply(target, function(t) table(findInterval(t, seq(0, 1, len=nbins + 1), rightmost.closed=TRUE))[as.character(seq_len(nbins))])
	t1[is.na(t1)] <- 0
	rownames(t1) <- seq_len(nbins)
	target <- apply(t1, 1, mean)

	# number of parameters to optimize
	npars <- (nenv + 1) * nspecies + ninteractions
	
	message(sprintf("Number of parameters to be generated: %d\n", npars))

	# Generate random environmental variables for the samples
	ENV <- matrix(stats::rnorm(nsamples * nenv), ncol=nenv, nrow=nsamples)
	if(nenv > 0)
		colnames(ENV) <- sprintf("env%02d", 1:nenv)

	# Create a random DAG of interactions
	repeat {
		spmask <- matrix(0, ncol=nspecies, nrow=nspecies)
		spmask[sample(setdiff(seq_len(nspecies * nspecies), seq(1, by=nspecies + 1, len=nspecies)), ninteractions)] <- 1
		if(!isCyclic(spmask)) break
	}

	# fitness function
	fitness <- function(par, spmask, nenv, env, targetdistr) {
		# convert parameter vector to model coefficients
		coefs <- par.to.coefs(par, spmask, nenv, env)
		true <- as.eicm(env=env, env.coefs=coefs[[1]], sp.coefs=coefs[[2]])

		# simulate 10 communities
		allpres <- lapply(seq_len(nrepetitions), function(i) predict.eicm(true, 1))
		# average their frequency distributions
		t <- sapply(allpres, function(pres) {
			t <- table(findInterval(apply(pres, 2, sum) / nsamples, seq(0, 1, len=nbins + 1), rightmost.closed=TRUE))[as.character(seq_len(nbins))]
			t[is.na(t)] <- 0
			return(t)
		})
		t2 <- apply(t, 1, mean)
	
		# the fitness is the difference between the two distributions
		return(sum(abs(t2 - targetdistr)) + shrinkage * (mean(abs(coefs[[1]]))  ))	#+ sum(abs(coefs[[2]])) / sum(coefs[[2]] != 0)
	}
		
	f <- pso::psoptim(par=rep(NA, npars), fn=fitness, gr=NULL, spmask=spmask, nenv=nenv, env=ENV, target, lower=-bounds, upper=bounds
		, control=list(trace=1, s=swarm.size, hybrid=FALSE, maxit=10000, maxit.stagnate=maxit.stagnate, vectorize=F))
	coefs <- par.to.coefs(f$par, spmask, nenv, ENV)
	model <- as.eicm(env=ENV, env.coefs=coefs[[1]], sp.coefs=coefs[[2]])
	attr(model, "target.distribution") <- target
	return(model)
}

par.to.coefs <- function(par, spmask, nenv, env, factor=3) {
	nspecies <- ncol(spmask)
	envcoefs <- matrix(par[1:(nspecies * (nenv + 1))], nrow=nspecies, ncol=nenv+1)
	# here we divide env betas by two because, typically, the intercept betas are larger-valued than env betas
	envcoefs[, -1] <- envcoefs[, -1] / factor
	
	spcoefs <- matrix(0, ncol=nspecies, nrow=nspecies)
	# the same / 2 for sp betas
	spcoefs[spmask != 0] <- par[(nspecies * (nenv + 1) + 1):length(par)] / factor

	if(nenv > 0)
		dimnames(envcoefs) <- list(sprintf("sp%03d", seq_len(nspecies)), c("(Intercept)", colnames(env)[1:nenv]))
	else
		dimnames(envcoefs) <- list(sprintf("sp%03d", seq_len(nspecies)), "(Intercept)")
	dimnames(spcoefs) <- list(rownames(envcoefs), rownames(envcoefs))
	return(list(envcoefs, spcoefs))
}

