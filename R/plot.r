#' Quick plot EICM results
#' 
#' Multiple types of plots for EICM models: coefficients, network topology, confidence intervals and likelihood profiles.
#' Allows to plot a single model or the comparison between an estimated and a true model.
#'
#' If no \code{true.model} is provided, \code{type} must be one of \code{confint}, \code{profile}, \code{network}.
#'
#' If \code{true.model} is provided, \code{type} must be one of \code{network} or \code{coefficients}.
#' In the latter case, see \code{\link{coefficientComparisonPlot}} for possible options.
#'
#' If \code{x} is of type \code{eicm.list} (as returned by \code{\link{eicm}}), this function first
#' tries to plot the model after network selection, then, if it was not computed, the fitted model with the full network.
#' 
#' @param x a EICM model.
#' @param type character. The type of plot, one of \code{confint}, \code{profile}, \code{network} or \code{coefficients}.
#'        See details.
#' @param true.model the true model to compare with (usually, the one used for simulating the data).
#' @param ... further arguments to pass to \code{\link{coefficientComparisonPlot}} or other plotting functions.
#'
#' @seealso \code{\link{coefficientComparisonPlot}}, \code{\link{confint.eicm}}
#' @return NULL.
#' @export
plot.eicm <- function(x, type=ifelse(is.null(true.model), "network", "coefficients"), true.model=NULL, ...) {
	if(is.na(pmatch(type, c("confint", "profile", "network", "coefficients"))))
		stop("Invalid plot 'type'.")
	if(is.null(true.model)) {
		switch(pmatch(type, c("confint", "profile", "network"))
		, {
			if(is.null(x$confint))
				stop("Confidence intervals are not computed yet.")
			plot.confint.eicm(x$confint, ...)
		}, {
			if(is.null(x$profile))
				stop("Likelihood profiles are not computed yet.")
			for(i in seq_len(6)) {
				oldpar <- graphics::par(no.readonly = TRUE)
				on.exit(graphics::par(oldpar))
				graphics::par(mfrow=c(2, 3))
				plot.profile.eicm(x$profile[[i]], ...)
			}
		}, {
			plot.eicm.matrix(x$model, type="network", ...)
		})
	} else {
		if(!inherits(true.model, "eicm")) stop("true.model must be of class 'eicm'")
		plot.eicm.matrix(x$model, true.model, type=type, ...)
	}
	invisible()
}

#' @inherit plot.eicm
#' @export
plot.eicm.list <- function(x, type=ifelse(is.null(true.model), "network", "coefficients"), true.model=NULL, ...) {
	if(!is.null(x$selected.model))
		plot.eicm(x$selected.model, type, true.model, ...)
	else if(!is.null(x$fitted.model)) {
		plot.eicm(x$fitted.model, type, true.model, ...)
	} else stop("No selected model and no fitted model in object.")
}

#' Plot EICM likelihood profile
#' 
#' Plot one likelihood profile with the line representing a confidence level threshold.
#' 
#' @param x a profile.eicm object, created with \code{\link{profile.eicm}}.
#' @param level the significance level desired.
#' @param ... other arguments passed to other functions.
#' @return NULL.
#' @export
plot.profile.eicm <- function(x, level=0.99, ...) {
	percentile <- stats::qchisq(level, 1) / 2
	graphics::plot.default(x, type="l", lwd=1.5, xlab="Parameter value", ylab="LogLik", main=paste(attr(x, "parameter"), collapse="~"), bty="l")
	graphics::points(x)
#	abline(h=max(x[, 2]) - percentile, lty=2)
	graphics::abline(h=attr(x, "maxLogLik") - percentile, lty=2)
	graphics::abline(v=attr(x, "parEstimate"), col="red", lty=2)
	graphics::abline(v=0)
	invisible()
}

#' Plot EICM estimates and confidence intervals
#'
#' Plot EICM estimates and confidence intervals, in a dot-and-whisker plot.
#'
#' @param x a eicm.confint object
#' @param truemodel for validation purposes only. The true model used to simulate data.
#' @param ... other arguments passed to other functions.
#' @return NULL.
#' @export
plot.confint.eicm <- function(x, truemodel=NULL, ...) {
	ci <- x
	signif <- ci[, 2] * ci[, 3] > 0
	y <- seq_len(nrow(ci))
	
	oldpar <- graphics::par(no.readonly = TRUE)
	on.exit(graphics::par(oldpar))
	graphics::par(mar=c(2, 8, 0.1, 0.1), cex=0.7)
	graphics::plot.new()
	graphics::plot.window(xlim=range(ci, finite=TRUE), ylim=range(y))
	graphics::abline(v=0, lwd=2)

	if(!is.null(truemodel)) {
		w <- which(truemodel$model$sp != 0, arr.ind=TRUE)
		true <- data.frame(from=rownames(truemodel$model$sp)[w[, 1]], to=colnames(truemodel$model$sp)[w[, 2]], coef=truemodel$model$sp[w])
		dir <- paste(true[, 1], "~", true[, 2])
		inv <- paste(true[, 2], "~", true[, 1])
		true[, "wasestimated"] <- (dir %in% rownames(ci)[signif] | inv %in% rownames(ci)[signif])
	
		m1 <- match(rownames(ci), dir)
		m1[is.na(m1)] <- 0
		m2 <- match(rownames(ci), inv)
		m2[is.na(m2)] <- 0
		truematch <- m1 + m2
		truematch[truematch == 0] <- NA
		truecoefs <- true[truematch, 3]
		
		graphics::abline(h=y, lty=3, col=ifelse(rownames(ci) %in% dir | rownames(ci) %in% inv, "blue", "gray"))
	} else
		graphics::abline(h=y, lty=3, col="gray")
		
	whi.finite <- is.finite(ci[, 2]) & is.finite(ci[, 3])
	whi.lowinf <- !is.finite(ci[, 2]) & is.finite(ci[, 3])
	whi.higinf <- is.finite(ci[, 2]) & !is.finite(ci[, 3])
	finite.ci <- ci[whi.finite, ]
	lowinf.ci <- ci[whi.lowinf, , drop=FALSE]
	higinf.ci <- ci[whi.higinf, , drop=FALSE]
	graphics::arrows(finite.ci[, 2], y[whi.finite], finite.ci[, 3], y[whi.finite]
		, angle=90, code=3, length=0.05, lwd=2, col=ifelse(signif[whi.finite], "black", "gray"))
		
	if(sum(whi.lowinf) > 0) {
		graphics::arrows(graphics::par("usr")[1], y[whi.lowinf], lowinf.ci[, 3], y[whi.lowinf]
			, angle=45, code=1, length=0.075, lwd=2, col=ifelse(signif[whi.lowinf], "black", "gray"))
		graphics::arrows(graphics::par("usr")[1], y[whi.lowinf], lowinf.ci[, 3], y[whi.lowinf]
			, angle=90, code=2, length=0.05, lwd=2, col=ifelse(signif[whi.lowinf], "black", "gray"))
	}

	if(sum(whi.higinf) > 0) {
		graphics::arrows(higinf.ci[, 2], y[whi.higinf], graphics::par("usr")[2], y[whi.higinf]
			, angle=45, code=2, length=0.075, lwd=2, col=ifelse(signif[whi.higinf], "black", "gray"))
		graphics::arrows(higinf.ci[, 2], y[whi.higinf], graphics::par("usr")[2], y[whi.higinf]
			, angle=90, code=1, length=0.075, lwd=2, col=ifelse(signif[whi.higinf], "black", "gray"))
	}
	# estimated sp coefs
	graphics::points(ci[, 1], y, pch=19, col=ifelse(signif, "black", "gray"), cex=1.5)
	
	if(!is.null(truemodel)) {
		# true sp coefs
		graphics::points(truecoefs, y, pch=19, col="blue", cex=2)
	}
	graphics::axis(1)
	graphics::axis(2, at=y, labels=rownames(ci), las=2)
	invisible()
}

