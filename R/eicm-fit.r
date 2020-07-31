#' Estimate a EICM model
#'
#' Estimates the parameter values of a EICM model from the provided observation data.
#' This is the low-level estimation function. Users should use \code{\link{eicm}} instead, particularly
#' if estimating latent variables and species interactions.
#'
#' By default, all species interactions are estimated. Uers can control which species interactions
#' are to be estimated with the arguments \code{forbidden}, \code{mask.sp} and \code{exclude.prevalence},
#' which place cumulative restrictions on which interactions to estimate. See \code{vignette("eicm")}
#' for commented examples.
#'
#' @inheritParams eicm.data
#' @param n.latent the number of latent variables to estimate.
#' @param forbidden a formula (or list of) defining which species interactions are not to be estimated. See details.
#'        This constraint is cumulative with other constraints (\code{mask.sp} and \code{exclude.prevalence}).
#' @param allowed a formula (or list of) defining which species interactions are to be estimated. See details.
#'        This constraint is cumulative with other constraints (\code{mask.sp} and \code{exclude.prevalence}).
#' @param mask.sp a scalar or a binary square species x species matrix defining which species interactions to exclude
#'        (0) or include (1) \emph{a priori}. If a scalar (0 or 1), 0 excludes all interactions, 1 allows all interactions.
#'        If a matrix, species in the columns affect species in the rows, so, setting \code{mask.sp[3, 8] <- 0}
#'        means that species #8 is assumed \emph{a priori} to not affect species #3.
#'        This constraint is cumulative with other constraints (\code{forbidden} and \code{exclude.prevalence}).
#' @param exclude.prevalence exclude species interactions which are caused by species
#'        with prevalence equal or lower than this value. This constraint is cumulative with
#'        other constraints (\code{forbidden} and \code{mask.sp})
#' @param options a \code{eicm.options} object defining options for fitting. Usually not needed, use
#'        \code{forbidden}, \code{mask.sp} and \code{exclude.prevalence} instead.
#' @param initial.values the starting values for all parameters. Used only for speeding up
#'        fitting when there are previous estimates available.
#' @param regularization a two-element numeric vector defining the regularization lambdas used for
#'        environmental coefficients and for species interactions respectively. See details.
#' @param regularization.type one of "lasso", "ridge" or "hybrid", defining the type of penalty to apply.
#'        Type "hybrid" applies ridge penalty to environmental coefficients and LASSO to interaction coefficients.
#' @param fast a logical defining whether to do a fast - but less accurate - estimation, or a normal estimation.
#' @param n.cores the number of CPU cores to use in the L-BFGS-B optimization.
#' @param optim.method the optimization function to use. Should be set to the default.
#' @param optim.control the optimization parameters to use. Should be set to the defaults.
#'
#' @note If estimating latent variables and species interactions, use \code{\link{eicm}} instead.
#'
#' @return A fitted \code{eicm} object.
#'
#' @seealso \code{\link{eicm}}, \code{\link{confint.eicm}}, \code{\link{plot.eicm}}
#'
#' @examples
#' # Simulate some random occurrence data
#' nenv <- 2
#' nsp <- 10
#' nsamples <- 200
#' 
#' env <- matrix(rnorm(nenv * nsamples), ncol=nenv, nrow=nsamples)
#' env.coefs <- matrix(runif((nenv + 1) * nsp, -4, 4), nrow=nsp)
#' sp.coefs <- matrix(0, nrow=nsp, ncol=nsp)
#' sp.coefs[3, 5] <- 3
#' sp.coefs[4, 8] <- 2
#' 
#' # Define a true model
#' truemodel <- as.eicm(env=env, env.coefs=env.coefs, sp.coefs=sp.coefs)
#' 
#' # realize the model
#' simulated.data <- predict(truemodel, nrepetitions=1)
#'
#' \donttest{
#' # fit the model without species interactions
#' fittedNoInt <- eicm.fit(simulated.data, env, mask.sp=0)
#'
#' # fit the model with all species interactions
#' fittedInt <- eicm.fit(simulated.data, env, mask.sp=1)
#' 
#' # compute confidence intervals for all parameters
#' fittedInt <- confint(fittedInt, ncores=2)
#' 
#' # plot estimated parameters and confidence intervals
#' plot(fittedInt, type="confint", truemodel=truemodel)
#' }
#' @export
eicm.fit <- function(occurrences, env=NULL, traits=NULL, intercept=TRUE,
	n.latent=0, forbidden=NULL, allowed=NULL, mask.sp=NULL, exclude.prevalence=0, options=NULL, initial.values=NULL,
	regularization=c(ifelse(n.latent > 0, 0.5, 0), 1), regularization.type="hybrid",
	fast=FALSE, n.cores=parallel::detectCores(),
	optim.method="L-BFGS-B",
	optim.control=if(fast) list(trace=1, maxit=10000, ndeps=0.0001, factr=1e10) else
		list(trace=1, maxit=10000, ndeps=0.0001)
	) {
	
	if(!(regularization.type %in% c("ridge", "lasso", "hybrid")))
		stop("Regularization type must be one of: ridge, lasso, hybrid")
		
	if(is.null(env))
		env <- matrix(0, nrow=nrow(occurrences), ncol=0)
		
	# add to the environmental matrix latents and intercept
	if(intercept && all(apply(env, 2, stats::var) > 0.00001)) {
		if(is.null(colnames(env))) {
			env <- cbind(1, env)
			colnames(env) <- c("(Intercept)", sprintf("env%02d", seq_len(ncol(env) - 1)))
		} else {
			env <- cbind("(Intercept)"=1, env)
		}

		message("Added a column for the intercept")
	}


	if(length(regularization) != 2)
		stop("Regularization must be a 2-element vector: the lambda for environmental component and for the interaction component")
		
	attr(regularization, "type") <- regularization.type

	call <- match.call()

	if(is.null(options))
		options <- eicm.options(mask=list(sp=mask.sp))
	else {
		if(!is.null(mask.sp)) stop("Either specify 'mask.sp' or 'options', not both.")
		options <- eicm.options(options)
	}

	if(n.latent > 0 && sum(options$mask$sp) > 0)
		warning("Estimating latent variables and interactions simultaneously leads to inconsistent results. Use 'eicm' for the conducting the complete fitting workflow.")
	
	nenv <- ncol(env) + n.latent
	nsamples <- nrow(occurrences)
	nspecies <- ncol(occurrences)
	
	# Interaction exclusions
	if(!is.null(exclude.prevalence) && exclude.prevalence > 0) {
		options <- excludePrevalence(options, exclude.prevalence, occurrences)
#		message(sprintf("Excluded from estimation interactions involving species with %d or less presences, or with %d or more presences.", exclude.prevalence, nsamples - exclude.prevalence))
		message(sprintf("Excluded from estimation interactions departing from (caused by) species with %d or less presences.", exclude.prevalence))
	}

	if(is.null(traits))
		traits <- matrix(ncol=0, nrow=nspecies, dimnames=list(colnames(occurrences), NULL))

	if(!is.null(forbidden)) {
		tmpmask1 <- getMaskFromForbidden(forbidden, traits, data=occurrences, invert=FALSE)
		# ensure the mask conforms to the occurrences
		tmpmask1 <- tmpmask1[colnames(occurrences), colnames(occurrences)]
	}

	if(!is.null(allowed)) {
		tmpmask2 <- getMaskFromForbidden(allowed, traits, data=occurrences, invert=TRUE)
		# ensure the mask conforms to the occurrences
		tmpmask2 <- tmpmask2[colnames(occurrences), colnames(occurrences)]
	}

	if(exists("tmpmask1"))
		options <- options + list(mask=list(sp=tmpmask1))

	if(exists("tmpmask2"))
		options <- options + list(mask=list(sp=tmpmask2))

	if(is.null(options$offset))
		options$offset <- list(
			env=matrix(0, ncol=nenv, nrow=nspecies)
			, sp=matrix(0, ncol=nspecies, nrow=nspecies))
	
	npars <- getNumberOfParameters(nspecies, nenv, options) + n.latent * nsamples
	spnames <- colnames(occurrences)
	envnames <- colnames(env)
	if(n.latent > 0) {
		envnames <- c(envnames, sprintf("latent%02d", 1:n.latent))
	}	
		
	message(sprintf("Number of parameters to be estimated: %d; %d environmental vars", npars, nenv))
	
	fixed.pars <- list(env=env, occurrences=occurrences
		, nlatent=n.latent, regularization=regularization
		, options=options, fast=fast)

	if(is.null(initial.values)) {
		# message("Estimating with initial values set to random and intercepts set to mean")
			
		prev <- apply(occurrences, 2, sum, na.rm=TRUE)
		nr.notna <- apply(!is.na(occurrences), 2, sum)
		# handle species with only 0s or 1s
		prev[prev == nr.notna] <- nr.notna[prev == nr.notna] - 1
		prev[prev == 0] <- 1
		freq <- prev / nr.notna
		initial.values <- list(
			env=matrix(c(
				stats::binomial()$linkfun(ifelse(freq == Inf, 0.5, ifelse(freq > 0.99999, 0.99999, freq)))
				#binomial()$linkfun(prev / nr.notna)#runif(nspecies, -2.5, 2.5)
				, stats::runif(nspecies * (nenv - 1), -0.5, 0.5) # rep(0, nspecies * (nenv - 1))
			), nrow=nspecies)
			, sp=matrix(stats::runif(nspecies * nspecies, -0.5, 0.5), ncol=nspecies) # matrix(0, ncol=nspecies, nrow=nspecies)
			, samples=matrix(stats::runif(n.latent * nsamples, -0.5, 0.5)	#rep(0.1, n.latent * nsamples)
				, ncol=n.latent, nrow=nsamples)
		)
		
	}
	init <- makeParsFromModelMatrices(initial.values, options$mask)
				

	if(!is.null(optim.control$ndeps) && length(optim.control$ndeps) == 1)
		optim.control$ndeps <- rep(optim.control$ndeps, npars)

	# Estimate parameters
	fitted <- switch(optim.method,
		"L-BFGS-B"={
			if(n.cores > 1) {
				message(sprintf("Optimizing with parallel L-BFGS-B%s", ifelse(fast, ", with approximate likelihood", "")))
				if(npars > 10000 & n.cores > 10) {
					n.cores <- 10
					message("Decreased number of cores to 10.")
				}
				cls <- parallel::makeCluster(n.cores, outfile="")
#				readline("Ready ")
#				save(fixed.pars, file="fp")
#				tmp <- optimParallel::optimParallel(init, single.objective.function, fixed.pars=fixed.pars,
				tmp <- optimParallelMP2(init, single.objective.function, fixed.pars=fixed.pars,
					control=optim.control, parallel=list(forward=TRUE, loginfo=FALSE, cl=cls))
				parallel::stopCluster(cls)
				tmp
			} else {
				message(sprintf("Optimizing with L-BFGS-B%s", ifelse(fast, ", with approximate likelihood", "")))
				stats::optim(init, single.objective.function, fixed.pars=fixed.pars, method="L-BFGS-B",
					control=optim.control)
			}
		}, "ucminf"={
			message(sprintf("Optimizing with ucminf%s", ifelse(fast, ", with approximate likelihood", "")))
			ucminf::ucminf(init, single.objective.function, fixed.pars=fixed.pars, hessian=0, control=optim.control)
		}
	)
	fitted$lastvalue <- fitted$value
	
	optim.mat <- makeModelMatricesFromPars(fitted$par, nspecies, nenv, spnames=spnames, envnames=envnames
		, options=options, n.latent, nsamples
	)
	
	if(n.latent > 0)
		colnames(optim.mat$samples) <- sprintf("latent%02d", 1:n.latent)

	if(n.latent > 0) {
		# Scale the estimated latents to have unit variance.
		# Accordingly counter-scale the corresponding betas.
#		for(c in seq_len(ncol(optim.mat$samples))) {
#			f <- sd(optim.mat$samples[, c])
#			n <- colnames(optim.mat$samples)[c]
#			optim.mat$samples[, c] <- optim.mat$samples[, c] / f
#			optim.mat$env[, n] <- optim.mat$env[, n] * f
#		}
	}
	
	# note that we add the estimated latents to the environmental data matrix	
	out <- as.eicm(
		env.coefs=optim.mat$env, sp.coefs=optim.mat$sp, latent=optim.mat$samples, options=options,
		occurrences=occurrences, env=cbind(env, optim.mat$samples), traits=traits, regularization=regularization
	)
	out[["optim"]] <- fitted
	out[["call"]] <- call
	return(out)
}

single.objective.function <- function(par, fixed.pars) {
	nlatent <- fixed.pars$nlatent
	nenv <- ncol(fixed.pars$env) + nlatent
	nspecies <- ncol(fixed.pars$occurrences)
	spnames <- colnames(fixed.pars$occurrences)
	nsamples <- nrow(fixed.pars$env)

#	matrices <- makeModelMatricesFromPars(par, nspecies, nenv, spnames=spnames
#		, mask=fixed.pars$mask
#	)

	if(is.null(fixed.pars$options$mask) && is.null(fixed.pars$options$offset)) {
		matrices <- .Call(SR__makeModelMatricesFromPars, par, as.integer(nspecies), as.integer(nenv), NULL, NULL)
		rownames(matrices$env) <- spnames
		rownames(matrices$sp) <- spnames
		colnames(matrices$sp) <- spnames
		
		if(nlatent > 0)
			matrices$samples <- matrix(par[seq(length(par) - nlatent * nsamples + 1, length(par))], ncol=nlatent, nrow=nsamples)
		else
			matrices$samples <- matrix(NA, ncol=0, nrow=nsamples)
	} else {
		matrices <- makeModelMatricesFromPars(par, nspecies, nenv, spnames=spnames
			, options=fixed.pars$options, nlatent=nlatent, nsamples=nsamples
		)
	}
	
	if(nlatent > 0)
		colnames(matrices$samples) <- sprintf("latent%02d", 1:nlatent)
	
	newenv <- cbind(fixed.pars$env, matrices$samples)
	
	model <- as.eicm(
		env=newenv, intercept=FALSE,
		env.coefs=matrices$env, sp.coefs=matrices$sp, latent=matrices$samples,
		regularization=fixed.pars$regularization)
	
#	model <- eicm(newenv, env.coefs=matrices$env, sp.coefs=matrices$sp, intercept=FALSE, # we've got intercept already
#		regularization=fixed.pars$regularization, regularization.type=attr(fixed.pars$regularization, "type"))
	
	if(fixed.pars$fast)
		llh <- logLikValue.eicm(model, fixed.pars$occurrences)
	else
		llh <- logLik.eicm(model, fixed.pars$occurrences)

	rm("model", "matrices")
	return(-llh)
}

