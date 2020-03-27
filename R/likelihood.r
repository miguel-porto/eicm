# TODO port to C?
#' EICM (penalized) log-likelihood
#' 
#' Compute the (penalized) log-likelihood of the data matrix included in the EICM model,
#' or the log-likelihood of a new occurrence data matrix given the model.
#' 
#' @param object a EICM model
#' @param occurrences the occurrence data matrix. If omitted, the data matrix used to fit the model is used.
#' @param allow.na logical. Allow NAs in the occurrence matrix? If no NAs exist, it's faster to set to FALSE.
#' @param ... additional argument(s) for methods.
#' @return A \code{logLik} object.
#' @export
logLik.eicm <- function(object, occurrences=NULL, allow.na=TRUE, ...) {
	if(!inherits(object, "eicm")) stop("Object must be of class 'eicm'")
	penalized <- FALSE
#	if(!is.null(eicm$optim$possible) && eicm$optim$possible == FALSE) {
#		out <- NA
#	} else {
		if(is.null(occurrences))
			occurrences <- object$data$occurrences
		
		if(is.null(occurrences)) stop("Real data must be given")
	
		if(allow.na)
			llh <- .Call(SR__likelihood_NAallowed, object$data$env, object$model$env, object$model$sp, as.integer(occurrences))
		else
			llh <- .Call(SR__likelihood, object$data$env, object$model$env, object$model$sp, as.integer(occurrences))

		out <- sum(llh)
		attr(out, "lik.samples") <- llh
		attr(out, "df") <- sum(object$model$sp != 0) + sum(object$model$env != 0)
		
		regularization <- attr(object, "regularization")
		if(!is.null(regularization) && any(regularization > 0)) {
			switch(attr(regularization, "type"), lasso={
				out <- out - 
					regularization[1] * sum(abs(object$model$samples)) -
					regularization[1] * sum(abs(object$model$env[, -1])) -	# exclude intercept from penalty
					regularization[2] * sum(abs(object$model$sp))
			}, ridge={
				out <- out - 
					regularization[1] * sum(object$model$samples ^ 2) -
					regularization[1] * sum(object$model$env[, -1] ^ 2) -	# exclude intercept from penalty
					regularization[2] * sum(object$model$sp ^ 2)
			}, hybrid={
				out <- out - 
					regularization[1] * sum(object$model$samples ^ 2) -
					regularization[1] * sum(object$model$env[, -1] ^ 2) -	# exclude intercept from penalty
					regularization[2] * sum(abs(object$model$sp))
			}, stop(sprintf("Invalid regularization type: %s", attr(regularization, "type"))))
			penalized <- TRUE
		}
#	}
	attr(out, "nobs") <- nrow(object$data$env)
	attr(out, "regularization") <- regularization
	if(penalized)
		class(out) <- c("logLik", "penalized")
	else
		class(out) <- "logLik"
	return(out)
}

# Super fast likelihood value with no fancy attributes and no checks. Just the value for internal use.
logLikValue.eicm <- function(eicm, presences, allow.na=TRUE) {
	if(allow.na)
		llh <- .Call(SR__likelihood_superfast_NAallowed, eicm$data$env, eicm$model$env, eicm$model$sp, as.integer(presences))
	else
		llh <- .Call(SR__likelihood_superfast, eicm$data$env, eicm$model$env, eicm$model$sp, as.integer(presences))
		
	regularization <- attr(eicm, "regularization")
	if(!is.null(regularization) && any(regularization > 0)) {
		switch(attr(regularization, "type"), lasso={
			llh <- llh - 
				regularization[1] * sum(abs(eicm$model$samples)) -
				regularization[1] * sum(abs(eicm$model$env[, -1])) -	# exclude intercept from penalty
				regularization[2] * sum(abs(eicm$model$sp))
		}, ridge={
			llh <- llh - 
				regularization[1] * sum(eicm$model$samples ^ 2) -
				regularization[1] * sum(eicm$model$env[, -1] ^ 2) -	# exclude intercept from penalty
				regularization[2] * sum(eicm$model$sp ^ 2)
		}, hybrid={
			llh <- llh - 
				regularization[1] * sum(eicm$model$samples ^ 2) -
				regularization[1] * sum(eicm$model$env[, -1] ^ 2) -	# exclude intercept from penalty
				regularization[2] * sum(abs(eicm$model$sp))
		}, stop(sprintf("Invalid regularization type: %s", attr(regularization, "type"))))
	}
	return(llh)
}

