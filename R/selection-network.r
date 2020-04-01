# This is intended to be called with a fitted model with no interactions.
# It selects the plausible interactions assuming that all interactions are possible,
# or assuming a reduced model if masksp is provided.
# popsize.sel can accept and expression where n is number of needed individuals to have one individual per bit.
# This function is not yet "export quality" but it may be exported some day.
model.selection.network <- function(object, regularization, regularization.type, penalty,
	masksp=NULL, forbidden=NULL, exclude.prevalence=0,
	estimate.latents=FALSE,	select.direction=TRUE, 
	criterion=function(object) {object$optim$lastvalue},
	popsize.factor=1, maxit.stagnated=100,
	parallel=TRUE, pmutation=0.01, crossover="singlepoint",
	initial=NULL, mutation.mean=0.1,
	fast=FALSE, optim.method = "L-BFGS-B", optim.control = list(trace=1, maxit=10000, ndeps=0.001), #, reltol=1e-11
	control = list(), monitor=NULL, max.cached=15000, ...) {

#	set.seed(NULL)
	attr(regularization, "type") <- regularization.type
	
	if(!estimate.latents) object$model$samples <- NULL
	
	nsp <- ncol(object$data$occurrences)
	if(is.null(masksp) && is.null(forbidden) && exclude.prevalence == 0) {
		# Select the full network
		nbits <- (nsp * (nsp - 1)) / 2	# only one bit per species pair.
		if(select.direction) {
			bit.type <- rep(2, nbits)
		} else {
			bit.type <- rep(1, nbits)
		}
	} else {	# Select a subnetwork
		if(!select.direction) stop("If masksp is provided, select.direction must be TRUE") # because it's easier that the mask determines it
		
		if(!is.null(forbidden)) {
			maskForbidden <- getMaskFromForbidden(forbidden, object$data$traits)
			# ensure the mask conforms to the presences
			maskForbidden <- maskForbidden[colnames(object$data$occurrences), colnames(object$data$occurrences)]
		} else
			maskForbidden <- matrix(1L, ncol=nsp, nrow=nsp)
		
		if(is.null(masksp))
			masksp <- maskForbidden
		else
			masksp <- masksp & maskForbidden	# combine both masks with AND

		# Exclude interactions involving species with <= the given prevalence		
		if(exclude.prevalence > 0) {
			exclude <- apply(object$data$occurrences, 2, sum) <= exclude.prevalence
			masksp[exclude, ] <- 0
			masksp[, exclude] <- 0
		}
		lt1 <- masksp[lower.tri(masksp)] == 1
		ut1 <- t(masksp)[lower.tri(masksp)] == 1
		bit.type <- ifelse(lt1 & ut1, 2, ifelse(lt1, 1, ifelse(ut1, -1, 0)))
		nbits <- sum(bit.type != 0)		
	}
		
	corrmat <- matrix(seq_len(nsp * nsp), ncol=nsp)
	# correspondence between each bit and the respective species pair (the index into the mask matrix)
	corresp <- rep(NA, length(bit.type))
	corresp[bit.type == 2] <- t(corrmat)[lower.tri(corrmat)][which(bit.type == 2)]
	corresp[bit.type == 1] <- corrmat[lower.tri(corrmat)][which(bit.type == 1)]
	corresp[bit.type == -1] <- t(corrmat)[lower.tri(corrmat)][which(bit.type == -1)]
	
	bit.type <- abs(bit.type[bit.type != 0])
	corresp <- corresp[!is.na(corresp)]

	needed <- sum(bit.type)
	# Create the initial population.
	# This ensures that all the bits are represented in the initial population, in both states.
	initial <- matrix(0, ncol=nbits, nrow=needed)
	cols <- rep.int(seq_len(nbits), times=bit.type)
	initial[cbind(seq_len(needed), cols)] <- (1 - c(1, diff(cols))) + 1
	
	popsize.sel <- popsize.factor * needed

	if(popsize.sel < needed) {
		initial <- initial[sort(sample(needed, popsize.sel)), ]
		message(sprintf("Population size is lower than %d, which is not recommended.", needed))
		warning(sprintf("Population size is lower than %d, which is not recommended.", needed))
	} else if(popsize.sel > needed) {
		initial <- rbind(initial, initial[sort(sample(needed, popsize.sel - needed, replace=TRUE)), ])
	}
	
	message(sprintf("Optimizing with %d bits of which %d are bidirectional, and a population of %d.", nbits, sum(bit.type == 2), popsize.sel))

	crossover <- switch(pmatch(crossover, c("singlepoint", "uniform"), nomatch=1), {
		message("Using single point crossover")
		GA::gabin_spCrossover
	}, {
		message("Using uniform crossover")
		GA::gabin_uCrossover
	})

	# run the genetic algorithm to search the best model, ranked by criterion
	results <- my.ga(
		#ifelse(select.direction, "ternary", "binary"), fitness=bitStringToModel.network, criterion=criterion
		fitness=bitStringToModel.network, criterion=criterion, modelobject=object, bit.types=bit.type, correspondence=corresp
		, regularization=regularization, penalty=penalty, estimate.latents=estimate.latents, optim.method=optim.method, optim.control=optim.control
		, select.direction=select.direction, fast=fast
		, nBits = nbits, popSize = popsize.sel, maxiter = 100000, run=maxit.stagnated, parallel=parallel
		, pmutation=pmutation, mutation=make.mutator(mutation.mean, bit.type)
		, crossover=crossover, max.cached=max.cached
		, monitor=monitor, suggestions=initial
		)
	return(results)
}

# Builds a suitable model mask according the provided mixed binary/ternary bit string
bitStringToMBO.network.with.direction <- function(string, modelobject, bit.types, correspondence) {
	nsp <- ncol(modelobject$data$occurrences)

	if(sum(bit.types > 0) != length(string)) {
		stop(sprintf("Bit length %d is not the same as the number of estimated parameters %d", length(string), sum(bit.types > 0)))
	}
		
	# construct a model mask based on the bits
	masksp <- matrix(0, nrow=nsp, ncol=nsp)
	masksp[correspondence] <- ifelse(bit.types == 1, string == 1, ifelse(string == 2, 1, ifelse(string == 1, NA, 0)))
	# the NAs are those that should be set to 1 transposed
	masksp[is.na(t(masksp))] <- 1
	masksp[is.na(masksp)] <- 0

	maskenv <- 1
		
	storage.mode(maskenv) <- "integer"
	storage.mode(masksp) <- "integer"

	mask <- list(env=maskenv, sp=masksp)

	# merge mask with model builder options
	if(is.null(modelobject$model$options))
		options <- list(mask=mask)
	else {
		options <- modelobject$model$options + list(mask=mask)
	}

	rm(mask)
	return(options)
}

# Builds a suitable model mask according the provided bit string
# NOTE: should we deprecate this feature?
bitStringToMBO.network <- function(string, modelobject) {
	nsp <- ncol(modelobject$data$occurrences)
	# how many bits do we need to encode model selection
	nbits <- (nsp * (nsp - 1)) / 2

	if(nbits != length(string)) {
		stop(sprintf("Bit length %d is not the same as the number of estimated parameters %d", length(string), nbits))
	}
		
	# construct a model mask based on the bits
	# i.e. mask out all coefficients with the 0 bit
	masksp <- matrix(0, nrow=nsp, ncol=nsp)
	masksp[lower.tri(masksp)] <- string
	# model selection does not account for direction. Direction is estimated each time
	masksp <- masksp + t(masksp)
	maskenv <- 1
	
	if(any(masksp > 1)) stop("Unexpected error in model mask")
	
	storage.mode(maskenv) <- "integer"
	storage.mode(masksp) <- "integer"
	
	mask <- list(env=maskenv, sp=masksp)

	# merge mask with model builder options
	if(is.null(modelobject$model$options))
		options <- list(mask=mask)
	else {
	# TODO do we need this?
		options <- list(mask=mask)
		# combine the calculated mask AND the given mask (to account for forbidden interactions)
		# NOTE THIS IS NEW, UNTESTED
	#	mask$sp <- (object$model$model.builder.options$mask$sp != 0) & (mask$sp != 0)
	#	options <- object$model$model.builder.options + list(mask=mask)
	}
	# fit model
	
	rm(mask)
	return(options)
}

bitStringToModel.network <- function(string, modelobject, bit.types, correspondence, regularization, estimate.latents=FALSE
	, similar.model=NULL, optim.method=NULL, optim.control=NULL, select.direction=TRUE, fast=FALSE, ...) {

	if(select.direction)
		options <- bitStringToMBO.network.with.direction(string, modelobject, bit.types, correspondence)
	else
		options <- bitStringToMBO.network(string, modelobject)

	init.values <- if(is.null(similar.model)) modelobject$model else similar.model
	
	# We need this to force random initial values in all suggestions
	# TODO: WHY? Because if we start those at 0, we'll probably be stuck in a local optimum?
#	init.values$sp[init.values$sp == 0] <- rnorm(sum(init.values$sp == 0))
#	init.values$env[init.values$env == 0] <- rnorm(sum(init.values$env == 0))

#	init.values$sp[init.values$sp == 0] <- sample(c(0.1, -0.1), sum(init.values$sp == 0), replace=TRUE)

#	init.values$sp[lower.tri(init.values$sp)][init.values$sp[lower.tri(init.values$sp)] == 0] <-
#		runif(sum(init.values$sp[lower.tri(init.values$sp)] == 0), -2, 2)
#		sample(c(0.1, -0.1), sum(init.values$sp[lower.tri(init.values$sp)] == 0), replace=TRUE)
	
#	init.values$env[init.values$env == 0] <- sample(c(0.1, -0.1), sum(init.values$env == 0), replace=TRUE)

	# Symmetrize because we don't want starting values to account for direction
	# we want to start with random direction to minimize propagating biases along generations
	init.values$sp <- init.values$sp + t(init.values$sp)
	init.values$sp <- init.values$sp * options$mask$sp
	
	if(estimate.latents) {
	# remove estimated latents from model data
		modelobject$data$env <- modelobject$data$env[, 1:(ncol(modelobject$data$env) - ncol(init.values$samples))]
	}

	npars <- getNumberOfParameters(ncol(modelobject$data$occurrences), ncol(modelobject$data$env), options)
	out <- eicm.fit(modelobject$data$occurrences, modelobject$data$env, intercept=FALSE, n.latent=ifelse(estimate.latents, ncol(init.values$samples), 0)
		, regularization=regularization, regularization.type=attr(regularization, "type")
		# if we have a similar model, use it for starting values
		, initial.values=init.values
		, fast=fast, optim.method = optim.method, optim.control = optim.control
		, options=options, ...)
	
	rm(options)
	gc()
	return(out)
}


gaMonitor.eicm <- function(object, bestmodel, worstmodel) {
	fitness.stats <- stats::quantile(stats::na.exclude(object@fitness), probs=c(0, 0.5, 1))
	nterms.stats <- stats::quantile(apply(object@population, 1, sum), probs=c(0, 0.5, 1))
	message(sprintf("\rIt %d (%.0fs, ciT %d+%d/%d) | Fit %.1f %.1f %.1f | #term %d %d %d",
		object@iter, object@time.took, object@cached, object@informed, object@cache.size
		, fitness.stats[3], fitness.stats[2], fitness.stats[1], nterms.stats[1], as.integer(nterms.stats[2]), nterms.stats[3]))
	utils::flush.console()
}


