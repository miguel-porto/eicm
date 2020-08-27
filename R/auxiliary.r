#' Print EICM model matrix
#' 
#' Prints an excerpt of the EICM model coefficients.
#' 
#' @param x a EICM model matrix.
#' @param ... additional argument(s) for methods.
#' 
#' @return NULL.
#' @export
print.eicm.matrix <- function(x, ...) {
	cat(
		sprintf("EICM model with:\nSpecies: %d\nEnvironmental variables: %d\nSamples: %d\nLatent variables: %d\n",
			ncol(x$sp), ncol(x$env), nrow(x$samples), ncol(x$samples)))
	cat("Environmental coefficients (head)\n")
	print(utils::head(x$env))
	cat("\nInteraction coefficients (head)\n")
	print(utils::head(x$sp))
	cat("\nValues of estimated latents (head)\n")
	print(utils::head(x$samples))
	return(invisible(NULL))
}

#' Extract EICM model coefficients
#' 
#' Extract the EICM model coefficients, organized in three separate matrices.
#' 
#' @param object a EICM model.
#' @param ... additional argument(s) for methods.
#' 
#' @return A list with three coefficient matrices:
#'   \describe{
#'     \item{$env}{environmental coefficients}
#'     \item{$sp}{species interaction coefficients. It reads as: species C (column) affects species R (row) with the coefficient sp[R, C]}
#'     \item{$samples}{estimated latent variable values}
#'   }
#' @export
coef.eicm <- function(object, ...) {
	return(object$model)
}

