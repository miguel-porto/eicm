#' Set EICM fitting options
#'
#' Construct a EICM options object to inform the fitting engine. There is usually no need to use this function directly.
#'
#' Possible options are currently: \code{mask} and \code{offset}. Both are lists having the same structure of an \code{eicm.matrix} object:
#' \describe{
#'   \item{$mask}{Binary matrices defining which coefficients are to be estimated in model fitting
#'                OR a scalar, constant for all coefficients. 0 or FALSE exclude the given term from estimation,
#'                i.e., fixes it at 0.
#'     \describe{
#'       \item{$mask$env}{environmental coefficient mask}
#'       \item{$mask$sp}{species interaction mask}
#'     }
#'   }
#'   \item{$offset}{Numeric matrices defining constant terms, to be fixed and not estimated
#'     \describe{
#'       \item{$offset$env}{environmental coefficient offset}
#'       \item{$offset$sp}{species interaction offset}
#'     }
#'   }
#' }
#' When an offset for a term is nonzero, the respective mask element will be automatically zeroed (so it is not estimated).
#'
#' @param ... a named list of options. See details.
#'
#' @return A \code{eicm.options} object with options for model fitting, currently a model mask and model offsets.
#' @export
eicm.options <- function(...) {
	opts <- list(...)
	if(length(opts) == 1 && is.null(names(opts)))	# user provided a plain list of options
		opts <- opts[[1]]
	
	opt.names <- names(opts)
	excess <- opt.names[is.na(pmatch(opt.names, c("mask", "offset")))]
	
	if(length(excess) > 0)
		warning("These options were not recognised and were ignored: ", excess)
	
	mask <- opts$mask
	offset <- opts$offset
	
	if(is.null(mask))	# default mask is estimate all
		mask <- list(env=1L, sp=1L)
	else {
		if(is.null(mask$env))	# default is to estimate all
			mask$env <- 1L
		if(is.null(mask$sp))
			mask$sp <- 1L
			
		# ensure mask is binary and integer
		mask$sp[mask$sp == 0] <- 0
		mask$sp[mask$sp != 0] <- 1
		mode(mask$sp) <- "integer"
		mask$env[mask$env == 0] <- 0
		mask$env[mask$env != 0] <- 1
		mode(mask$env) <- "integer"
		mask <- list(env=mask$env, sp=mask$sp)	# this is just to place it in the correct order cause C code does not check
	}
	
	mask <- combineMaskAndOffset(mask, offset)

	out <- list(
		mask=mask
		, offset=offset
	)
		
	class(out) <- "eicm.options"
	return(out)
}

#' Define a model object for a EICM model
#'
#' Constructs a EICM model object for prediction. The model object contains all coefficient matrices that may be needed for prediction.
#' Usually, you don't need to invoke this function directly, use \code{\link{as.eicm}} instead.
#'
#' The EICM model is a list composed of three matrices plus the fitting options:
#' \enumerate{
#'   \item env
#'   \item sp
#'   \item samples
#' }
#'
#' @param env.coefs the environmental coefficient matrix: a species x variable matrix (including intercept).
#' @param sp.coefs the species interaction coefficient matrix: a species x species matrix, with zero diagonal.
#' @param latent the values for the latent variables in each sample: a sample x latent variable matrix.
#' @param options options for the model fitting.
#' 
#' @return A \code{eicm.matrix} object that can be used for defining a model.
#' @export
eicm.matrix <- function(env.coefs, sp.coefs=NULL, latent=NULL, options=NULL) {
	if(is.null(sp.coefs))
		sp.coefs <- matrix(0, ncol=nrow(env.coefs), nrow=nrow(env.coefs))
		
	if(is.null(colnames(sp.coefs)) && is.null(rownames(env.coefs))) {
		colnames(sp.coefs) <- sprintf("sp%03d", seq_len(ncol(sp.coefs)))
		rownames(sp.coefs) <- colnames(sp.coefs)
		rownames(env.coefs) <- colnames(sp.coefs)
	}
	if(is.null(colnames(env.coefs))) {
		colnames(env.coefs) <- c("(Intercept)", sprintf("env%02d", seq_len(ncol(env.coefs) - 1)))
	}

	out <- list(env=env.coefs, sp=sp.coefs, samples=latent, options=options)
	class(out) <- "eicm.matrix"
	return(out)
}

#' Define a data object for a EICM model
#'
#' Constructs a EICM data object for prediction. The data object contains all data matrices that may be needed for prediction.
#' Usually, you don't need to invoke this function directly, use \code{\link{as.eicm}} instead.

#' @param occurrences a binary (0/1) sample x species matrix, possibly including NAs.
#' @param env an optional sample x environmental variable matrix, for the known environmental predictors.
#' @param traits an optional species x trait matrix. Currently, it is only used for excluding
#'        species interactions \emph{a priori}.
#' @param intercept logical specifying whether to add a column for the species-level intercepts.
#' 
#' @return A \code{eicm.data} object that can be used for defining a model.
#' @export
eicm.data <- function(occurrences=NULL, env=NULL, traits=NULL, intercept=TRUE) {
	if(!is.null(occurrences) && !all(sort(unique(as.vector(occurrences))) == c(0, 1))) {
		stop("'occurrences' must be a presence-absence matrix, with only 0s and 1s.")
	}

	if(is.null(env) && !is.null(occurrences))
		env <- matrix(nrow=nrow(occurrences), ncol=0)

	if(is.null(traits) && !is.null(occurrences))
		traits <- matrix(ncol=0, nrow=ncol(occurrences), dimnames=list(colnames(occurrences), NULL))

	if(intercept && all(apply(env, 2, stats::var) > 0.00001)) {
		if(is.null(colnames(env))) {
			env <- cbind(1, env)
			colnames(env) <- c("(Intercept)", sprintf("env%02d", seq_len(ncol(env) - 1)))
		} else {
			env <- cbind("(Intercept)"=1, env)
		}
		#message("Added a column for the intercept")
	}

	out <- list(
		env=env,
		occurrences=occurrences,
		traits=traits
	)
	class(out) <- "eicm.data"
	return(out)	
}

#' Define a parameterized EICM model
#'
#' Constructs a EICM model object from given coefficients and data. Useful for simulating "true" models,
#' otherwise only used internally.
#' 
#' \code{regularization} is only used for storing the regularization lambdas used in model fitting.
#' It is ignored in simulation.
#'
#' @inheritParams eicm.fit
#' @inheritParams eicm.data
#' @inheritParams eicm.matrix
#'
#' @return A \code{eicm} object that can be used for prediction.
#' @note This function is only useful for simulation purposes. If you want to predict values from a fitted model,
#' a \code{eicm} object is already provided for the fitted model.
#' 
#' @seealso \code{\link{predict.eicm}}
#' 
#' @examples
#' # Generate some coefficients
#' nenv <- 2
#' nsp <- 20
#' nsamples <- 200
#' 
#' env <- matrix(rnorm(nenv * nsamples), ncol=nenv, nrow=nsamples)
#' env.coefs <- matrix(runif((nenv + 1) * nsp, -4, 4), nrow=nsp)
#' sp.coefs <- matrix(0, nrow=nsp, ncol=nsp)
#' sp.coefs[3, 5] <- 3
#' sp.coefs[4, 8] <- 2
#' 
#' # Define a true model (including environmental data)
#' truemodel <- as.eicm(env=env, env.coefs=env.coefs, sp.coefs=sp.coefs)
#'
#' # We can now realize it
#' predict(truemodel)
#' @export
as.eicm <- function(env.coefs, sp.coefs=NULL, latent=NULL, options=NULL, occurrences=NULL, env=NULL,
	traits=NULL, intercept=TRUE, regularization=NULL) {
	data <- eicm.data(occurrences=occurrences, env=env, traits=traits, intercept=intercept)
	modelmatrix <- eicm.matrix(env.coefs=env.coefs, sp.coefs=sp.coefs, latent=latent, options=options)
	
	out <- list(
		data=data,
		model=modelmatrix
	)

	attr(out, "regularization") <- regularization
	class(out) <- "eicm"
	return(out)
}

