#' Predict method for EICM fits
#'
#' Obtains probability predictions (conditional or unconditional) or a model realization from a
#' parameterized EICM model.
#'
#' The interaction network of the model must be an \strong{acyclic graph}.
#' Predictions are obtained by realizing the model multiple times and averaging realizations,
#' because there is not a closed-form expression for their calculation.
#' 
#' To obtain conditional predictions, include presence/absence columns with species names in \code{newdata}.
#' Named columns for all the environmental predictors must always be included.
#'
#' @param object the fitted EICM model.
#' @param nrepetitions the number of realizations to conduct for computing probabilities. Set to 1
#'        if you only need simulated community data.
#' @param newdata optionally, a matrix in which to look for variables with which to predict.
#'        If omitted, the original data (used to fit the model) is used.
#' @param ... not used.
#' 
#' @note If the eicm was fit without regularization, unconditional predictions are numerically equal
#' to those of simple binomial GLMs.
#' 
#' @return A species x sample matrix with predictions. If nrepetitions=1, predictions correspond to one
#' realization, otherwise they are probabilities.
#' 
#' @examples
#' # Load the included parameterized model
#' data(truemodel)
#' 
#' # for reference, plot the species interaction network
#' plot(truemodel, type="network")
#' 
#' # Unconditional predictions
#' # let's fix environmental predictors at 0, for simplicity.
#' predict(truemodel, newdata=cbind(env01=0, env02=0))
#' 
#' # Conditional predictions
#' # predict probabilities for all species conditional on the
#' # known presence of sp011 (compare sp014 and sp004 with the above)
#' predict(truemodel, newdata=cbind(env01=0, env02=0, sp011=1))
#' 
#' # Propagation of indirect species effects
#' # predict probabilities for all species conditional on the known 
#' # absence (first line) and known presence (second line) of sp005 and sp023
#' predict(truemodel, newdata=cbind(env01=0, env02=0, sp012=c(0, 1), sp018=c(0, 1)), nrep=100000)
#' 
#' # Notice the effects on sp026 and the effect propagation to those
#' # species which depend on it (sp013, sp008)
#' # Also compare with unconditional predictions
#' @export
predict.eicm <- function(object, nrepetitions=10000, newdata=NULL, ...) {
	if(!inherits(object, "eicm")) stop("Object must be of class 'eicm'")
	
	nspecies <- nrow(object$model$env)
	nrepetitions <- as.integer(nrepetitions)
	
	# build a species dependency tree
	# TODO this could be in C?
	cycles <- calculateCycles(object$model$sp)
	
	if(inherits(cycles, "character")) {
		stop("Cyclic species dependency graph with species ", paste(cycles, collapse=", "))
	}

	if(is.null(newdata)) {
		env.table <- object$data$env
		sp.table <- NULL
	} else {
#		if(!inherits(newdata, "eicm.data")) stop("'newdata' must be of class 'eicm.data'")
		# remove intercept column
		object$data$env <- object$data$env[, colnames(object$data$env) != "(Intercept)", drop=FALSE]
		newdata <- newdata[, colnames(newdata) != "(Intercept)", drop=FALSE]
		if(!all(colnames(object$data$env) %in% colnames(newdata)))	# TODO this doesn't need to be so
			stop("Newdata matrix must have all the environment variables of the fitted model")

		# extract and put new data vars in the same order as the fitted model
		env.table <- newdata[, colnames(object$data$env), drop=FALSE]
		# remove those from newdata
		newdata <- newdata[, which(!(colnames(newdata) %in% colnames(object$data$env))), drop=FALSE]

		missingspecies <- !(colnames(newdata) %in% rownames(object$model$env))
		if(any(missingspecies)) {
			warning("Not all provided species belong to the model, those were ignored")
			newdata <- newdata[, !missingspecies, drop=FALSE]
		}

		sp.table <- matrix(NA, nrow=nrow(env.table), ncol=nspecies, dimnames=list(NULL, rownames(object$model$env)))
		sp.table[, colnames(newdata)] <- newdata
		# so now we have a table with all model species, filled with NAs for those that were not provided as newdata
	}
	# now add again the intercept
	#if(all(apply(env.table, 2, var) > 0.00001)) {
	if(!("(Intercept)" %in% colnames(env.table))) {
		env.table <- cbind("(Intercept)"=1, env.table)
		#message("Added a column for the intercept")
	}
			
	# fill in presence matrix along the dependency tree, start with those that only depend on environment
	randomSeed <- as.integer(stats::runif(1, 0, 32000))
	fitted <- .Call(SR__simulate_community_probability, nrepetitions
		, env.table, if(is.null(sp.table)) NULL else as.integer(sp.table)
		, object$model$env, object$model$sp, cycles, randomSeed)
	
	colnames(fitted) <- rownames(object$model$env)
	return(fitted)
}


# from a square pairwise coefficient matrix, compute the species dependency tree and check for cyclic graphs
calculateCycles <- function(coefs) {
	cycles <- list()
	tmpspcoefs <- coefs
	repeat {
		leaf.nodes <- apply(tmpspcoefs != 0, 1, sum, na.rm=TRUE) == 0	# NOTE: small coefs have been zeroed elsewhere
		if(all(!leaf.nodes))	# cyclic graph!
			return(if(is.null(colnames(tmpspcoefs))) NA else colnames(tmpspcoefs) )

		#cycles <- c(cycles, list(names(leaf.nodes)[leaf.nodes]))
		cycles <- c(cycles, list(as.integer(match(names(leaf.nodes)[leaf.nodes], colnames(coefs)))))
	
		if(sum(!leaf.nodes) == 0) break;
		tmpspcoefs <- tmpspcoefs[!leaf.nodes, !leaf.nodes, drop=FALSE]
	}
	return(cycles)
}

