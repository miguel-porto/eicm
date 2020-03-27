is.scalar <- function(x) {
	!is.null(x) && is.atomic(x) && length(x) == 1L && !is.character(x)
}

# When an offset is non-zero, the parameter is fixed, and will not be estimated
# i.e. it will be removed from the mask
combineMaskAndOffset <- function(mask, offset=NULL) {
	if(!is.null(offset)) {
		if(is.scalar(mask$env))	# we need to expand the all-zero matrix to combine
			mask$env <- matrix(mask$env, ncol=ncol(offset$env), nrow=nrow(offset$env))

		mask$env <- (mask$env != 0) & (offset$env == 0)

		if(is.scalar(mask$sp))	# we need to expand the all-zero matrix to combine
			mask$sp <- matrix(mask$sp, ncol=ncol(offset$sp), nrow=nrow(offset$sp))

		mask$sp <- (mask$sp != 0) & (offset$sp == 0)
		mask$sp <- mask$sp & (t(offset$sp) == 0)
	}
	
	if(!is.scalar(mask$env)) {	# we've got a matrix, so convert to binary (0-not to be estimated)
		mask$env <- mask$env != 0
		mode(mask$env) <- "integer"
	}
	if(!is.scalar(mask$sp)) {	# we've got a matrix, so convert to binary (0-not to be estimated)
		mask$sp <- mask$sp != 0
		mode(mask$sp) <- "integer"
	}

	return(mask)
}

# used for GA optimization
makeModelMatricesFromPars <- function(par, nspecies, nenv, spnames=NULL
	, envnames=NULL, options=eicm.options(), nlatent=0, nsamples=0) {

	if(!inherits(options, "eicm.options"))
		options <- eicm.options(options)
		
	if(is.null(options$offset)) {
		options$offset <- list(
			env=matrix(0, ncol=nenv, nrow=nspecies)
			, sp=matrix(0, ncol=nspecies, nrow=nspecies))
	}
		
#return(list(env=matrix(as.integer(mask$env), ncol=nenv), sp=matrix(as.integer(mask$sp), ncol=nspecies)))
	out <- .Call(SR__makeModelMatricesFromPars, par, as.integer(nspecies), as.integer(nenv),
		options$mask, options$offset)
			
	dimnames(out$env) <- list(spnames, envnames)
	dimnames(out$sp) <- list(spnames, spnames)

	if(nlatent > 0)
		out$samples <- matrix(par[seq(length(par) - nlatent * nsamples + 1, length(par))], ncol=nlatent, nrow=nsamples)
	else
		out$samples <- matrix(NA, ncol=0, nrow=nsamples)
		
	return(out)
}

# used for GA optimization
makeParsFromModelMatrices <- function(matrices, mask=NULL) {
	# TODO ensure order
	if(!is.null(mask)) {
		mask$env <- as.integer(mask$env)
		mask$sp <- as.integer(mask$sp)
	}
	
	out <- .Call(SR__makeParsFromModelMatrices, matrices, mask)
	
	if(!is.null(matrices$samples))
		out <- c(out, as.vector(matrices$samples))
	
	return(out)
}

# used for GA optimization
getNumberOfParameters <- function(nspecies, nenv, options=eicm.options()) {
	return(.Call(SR__getNumberOfParameters, as.integer(nspecies), as.integer(nenv), options$mask))
}

# used for GA optimization
getMostSimilarModel <- function(popToEval, cachedModelList) {
	return(.Call(SR__getMostSimilarModel, popToEval, cachedModelList))
}

# Is the adjacency matrix a cyclic graph?
isCyclic <- function(coefs) {
	return(.Call(SR__isCyclic, coefs))
}

# mask out from interactions species that affect others and that have a prevalence <= than given value
# combine with given mask
excludePrevalence <- function(options, prevalence.threshold, occurrences) {
	prev <- apply(occurrences, 2, sum)
	nsamples <- nrow(occurrences)
	nspecies <- ncol(occurrences)
	
	exclude <- prev <= prevalence.threshold
#	exclude <- exclude | (prev >= (nsamples - prevalence.threshold))

	spmask <- matrix(1, ncol=nspecies, nrow=nspecies)
#	spmask[exclude, ] <- 0
	spmask[, exclude] <- 0	# only exclude when species is in the columns (i.e. it's the "affector")
	options <- options + list(mask=list(sp=spmask))
	return(options)
}

# Converts a list of formulae that depict forbidden (or allowed) interactions, into a model mask
# Formulae must be of the form A ~ B + C + ..., in which all variables are either species names or trait categories.
# The meaning is A must not depend on B nor C
getMaskFromForbidden <- function(forbidden, traits, invert=FALSE, data=NULL) {
 	if(inherits(forbidden, "formula"))
 		forbidden <- list(forbidden)
	spnames <- rownames(traits)
	nspecies <- nrow(traits)
	trait.categories <- lapply(traits, levels)
	all.traits <- unlist(trait.categories)
	if(any(duplicated(all.traits))) stop("There must not be duplicated trait categories between traits.")
	if(any(duplicated(c(spnames, all.traits)))) stop("There must not be duplicated trait categories and species names.")
	
	lterms <- lapply(forbidden, stats::terms, data=data)
	f.vars <- lapply(lterms, function(a) as.character(attr(a, "variables"))[-1])
	resp <- mapply(function(a, b) a[attr(b, "response")], f.vars, lterms)
	rhs <- lapply(lterms, function(a) colnames(attr(a, "factors")))
#	rhs <- mapply(function(a, b) a[colnames(attr(b, "factors"))], f.vars, lterms, SIMPLIFY=FALSE)

	mask <- matrix(ifelse(invert, 0, 1), ncol=nspecies, nrow=nspecies, dimnames=list(spnames, spnames))
	diag(mask) <- 0
	for(i in seq_along(forbidden)) {
		if(resp[i] %in% spnames) {
			rowindex <- resp[i]
		} else {
			which.trait <- which(sapply(trait.categories, function(t) resp[i] %in% t))
			rowindex <- traits[, which.trait] == resp[i]
		}

		colindex <- rep(FALSE, nspecies)
		for(vars in rhs[[i]]) {
			if(vars %in% spnames) {
				colindex[vars == spnames] <- TRUE
			} else {
				which.trait <- which(sapply(trait.categories, function(t) vars %in% t))
				colindex[traits[, which.trait] == vars] <- TRUE
			}
		}
		
		mask[rowindex, colindex] <- ifelse(invert, 1, 0)
	}
	return(mask)
}

