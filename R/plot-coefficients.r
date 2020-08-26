#' Diagnostic monitor plots for EICM
#'
#' Visually compare the true model with estimation results (final results or during model fitting)
#' and compute accuracy statistics.
#'
#' @param model the EICM model of interest.
#' @param true.model the true model to compare with (usually, the one used for simulating the data).
#' @param nenv.to.plot the number of environmental variables to plot.
#' @param nlatent.to.plot the number of latent variables to plot.
#' @param plot.intercept logical. plot the species-level intercepts?
#' @param plot.interactions logical. Plot interaction coefficient scatterplot?
#' @param excluded.interactions a binary species x species matrix telling which interactions were excluded \emph{a priori}.
#' @param layout logical. Do multi-panel layout?
#' @param noplot logical. Do plots? If TRUE, it will return the accuracy statistics only.
#' @param env.stats logical. Compute accuracy for environmental predictors?
#' @param legend logical. Plot legend?
#' @return A vector with accuracy statistics.
#' @export
coefficientComparisonPlot <- function(model, true.model, nenv.to.plot=0, nlatent.to.plot=0,
	plot.interactions=any(model$sp != 0), plot.intercept=FALSE, excluded.interactions=NULL,
	layout=TRUE, noplot=FALSE, env.stats=TRUE, legend=TRUE) {
	if(!inherits(model, "eicm.matrix")) stop("Object must be of class 'eicm.matrix'")
	if(!inherits(true.model, "eicm")) stop("true.model must be of class 'eicm'")
	mat <- model
	
	if(layout && !noplot) {
		n.plots <- nenv.to.plot + nlatent.to.plot * (ncol(true.model$data$env) - 1) + plot.interactions * 1 + plot.intercept
		if(n.plots < 1) stop("No plots to do!")
		oldpar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(oldpar))
		graphics::par(mfrow=c(1, n.plots))
	}
	
	# reorder latents to give them the same order of the true model, based on their correlation
	crosscor <- stats::cor(mat$env, true.model$model$env)
	nvars <- nrow(crosscor)
	# this is to ensure only one match is selected in each row-col
	ord <- order(abs(crosscor), decreasing=TRUE)
	matches <- integer(nvars)
	matches[(ord[1] - 1) %/% nvars + 1] <- ord[1]
	b <- matrix(nrow=0, ncol=2)
	b <- rbind(b, c((ord[1] - 1) %/% nvars, (ord[1] - 1) %% nvars))
	for(i in seq_along(ord)) {
		if(!( (((ord[i] - 1) %/% nvars) %in% b[, 1] ) || (((ord[i] - 1) %% nvars) %in% b[, 2]))) {
			matches[(ord[i] - 1) %/% nvars + 1] <- ord[i]
			b <- rbind(b, c((ord[i] - 1) %/% nvars, (ord[i] - 1) %% nvars))
		}
	}
	matches <- (matches - 1) %% nvars + 1
	#matches <- apply(crosscor, 1, function(a) which.max(abs(a)))
	mat$env <- mat$env[, matches, drop=FALSE]
	if(nlatent.to.plot > 0) {
		matches <- matches[-1] - 1
		mat$samples <- mat$samples[, matches]
	}

	# INTERCEPTS PLOT
	if(!noplot) {
		if(plot.intercept) {
			graphics::plot(I(mat$env[, 1] ) ~ true.model$model$env[, 1], asp=1, bty="l", xlab="True", ylab="Estimated", pch=19, cex=0.75, main="Intercepts")
			graphics::abline(0, 1, col="gray")
			graphics::abline(0, -1, col="gray")
			graphics::abline(v=0, h=0, lty=2)
		}
	}

	# ENVIRONMENTAL BETAS PLOT
	qual.fit <- list()
	if(nenv.to.plot > 0)
		envbetas <- true.model$model$env[, (ncol(true.model$model$env) - nenv.to.plot + 1):ncol(true.model$model$env), drop=FALSE]

	f <- rep(1, nenv.to.plot)
	for(i in seq_len(nenv.to.plot)) {
		x1 <- envbetas[, i]
		y1 <- mat$env[, ncol(mat$env) - nenv.to.plot + i] * f[i]
		if(!noplot) {
			graphics::plot(I(y1) ~ x1, asp=1, bty="l", xlab="True", ylab="Estimated", pch=19, cex=0.75, main="Environmental betas")#, xlim=c(-8, 8), ylim=c(-8, 8))
#		abline(lm(y1 ~ x1), col="blue")
			graphics::abline(0, 1, col="gray")
			graphics::abline(0, -1, col="gray")
			graphics::abline(v=0, h=0, lty=2)
		}
		if(env.stats)
			qual.fit[[i]] <- stats::glm(y1 ~ envbetas)
	}
		
	if(!noplot) {
		for(i in seq_len(nlatent.to.plot)) {
			for(j in seq_len(ncol(true.model$data$env) - 1)) {
				x2 <- true.model$data$env[, j + 1]
				graphics::plot(I(mat$samples[, i] / f[i]) ~ x2, asp=1, bty="l", xlab=sprintf("True #%d", j), ylab=sprintf("Latent #%d", i), pch=19, cex=0.75)#, xlim=c(-8, 8), ylim=c(-8, 8))
				graphics::abline(0, 1, col="gray")
				graphics::abline(0, -1, col="gray")
				graphics::abline(v=0, h=0, lty=2)
			}
		}
	}
	if(FALSE) {
		if(nlatent.to.plot > 0) {
			r.sq <- sapply(seq_len(ncol(mat$samples)), function(j) {
				summary(lm(mat$samples[, j] ~ true.model$data$env[, -1, drop=FALSE]))$r.squared
			})
			out <- c(r.sq=r.sq)
		} else out <- numeric(0)
	}
	
	out <- numeric(0)
	if(length(qual.fit) > 0) {
		coefs <- sapply(qual.fit, function(m) {
			stats::coef(m)[-1]
		})
		out <- c(out, coefs)

		r2 <- sapply(qual.fit, function(m) {
			(m$null.deviance - m$deviance) / m$null.deviance
		})
		out <- c(out, r2)
#		names(out) <- c(paste0("env", seq_len(nenv.to.plot), ".beta"), paste0("env", seq_len(nenv.to.plot), ".r2"))
		names(out) <- c(rep("env.betas", nenv.to.plot ^ 2), paste0("env", seq_len(nenv.to.plot), ".r2"))
	}
	
	if(plot.interactions)
		out <- c(out, plot.true.estimated(mat, true.model=true.model, time.took=NULL, excluded.interactions=excluded.interactions, noplot=noplot, layout=layout, legend=legend))
	
	return(out)
}
	
plot.true.estimated <- function(estimated.model.coefs, true.model, time.took=NULL,
	estimated.model=NULL, excluded.interactions=NULL, noplot=FALSE, layout=TRUE, legend=TRUE) {

	interactions <- find.correct.spurious(estimated.model.coefs, true.model, excluded.interactions)

	colors <- matrix("#ffffff00", ncol=ncol(true.model$model$sp), nrow=nrow(true.model$model$sp))
	colors[interactions$correctEstimation] <- "#00aa00"
	colors[interactions$spuriousEstimation] <- "red"
	colors[interactions$missedEstimation] <- "gray"
	colors[interactions$nodirection.wrongDirection] <- "blue"
	if(!is.null(interactions$excludedInteractions))
		colors[interactions$excludedInteractions] <- "black"
#	colors[interactions$correctDirection] <- "#00aa00"
	colors[upper.tri(colors)] <- NA
	diag(colors) <- NA

	stats.spurious <- sum(interactions$spuriousEstimation)
	stats.missed <- sum(interactions$missedEstimation)
	stats.severeMissed <- sum(interactions$severeMissedEstimation)
	if(!is.null(interactions$excludedInteractions))
		stats.excluded <- sum(interactions$excludedInteractions)
	else
		stats.excluded <- NULL
	if(!is.null(interactions$excludedTrueNegatives))
		stats.excludedTN <- sum(interactions$excludedTrueNegatives)
	else
		stats.excludedTN <- NULL
	stats.correct <- sum(interactions$correctEstimation) / 2
	stats.wrongdir <- sum(colors=="blue", na.rm=TRUE)
	stats.discarded <- sum(colors=="#ffffff00", na.rm=TRUE)

	if(!noplot) {
		graphics::plot(interactions$nodirection.estimated ~ interactions$nodirection.true, asp=1, type="n", bty="l", col=colors,
		xlab=ifelse(layout, "True", NA), ylab=ifelse(layout, "Estimated", NA),
		main=ifelse(layout, "Species Interaction coefficients", NA),
		pch=19, cex=0.75)

	#rect(-THRESHOLD, usr[3], THRESHOLD, usr[4], col="gray", border=NA)
		graphics::abline(0, 1)
		graphics::abline(v=0, h=0, lty=2)
		graphics::points(interactions$nodirection.estimated ~ interactions$nodirection.true, col=colors, pch=19)
		if(legend) {
			legend("topleft", c(
				sprintf("True positives (%d)", stats.correct),
				sprintf("of which wrong direction (%d)", stats.wrongdir),
				ifelse(is.null(stats.excluded),
					sprintf("False negatives (%d)", stats.missed) ,
					sprintf("False negatives (%d) (%d hard)", stats.missed, stats.excluded)),
				sprintf("False positives (%d)", stats.spurious),
				ifelse(is.null(stats.excludedTN),
					sprintf("True negatives (%d)", stats.discarded),
					sprintf("True negatives (%d) (%d hard)", stats.discarded, stats.excludedTN))
				)
			, bty="n", pch=19, col=c("#00aa00", "blue", "gray", "red", NA), cex=0.8)
		}
		if(!is.null(estimated.model))
			graphics::mtext(sprintf("logLik: %.2f", logLik.eicm(estimated.model)), cex=0.7, line=-0.1)
	}

	FP <- stats.spurious
	FN <- stats.missed
	FN.severe <- stats.severeMissed
	TP <- stats.correct
	TN <- stats.discarded

	return(c(
		spurious=stats.spurious,
		missed=stats.missed,
		severeMissed=stats.severeMissed,
		excluded=stats.excluded,
		correct=stats.correct,
		wrongdir=stats.wrongdir,
		discarded=stats.discarded,
		error.rate=round((FP + FN) / (TP + TN + FN + FP), 3),
		accuracy=round((TP + TN) / (TP + TN + FN + FP), 3),
		accuracyStrong=round((TP + TN) / (TP + TN + FN.severe + FP), 3),
		sensitivity=round(TP / (TP + FN), 3),
		sensitivityStrong=round(TP / (TP + FN.severe), 3),
		specificity=round(TN / (TN + FP), 3),
		precision=round(TP / (TP + FP), 3),
		MCC=round((TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)), 3),
		MCCStrong=round((TP * TN - FP * FN.severe) / sqrt((TP + FP) * (TP + FN.severe) * (TN + FP) * (TN + FN.severe)), 3),
		interaction.betafit=interactions$int.beta,
		interaction.r2fit=interactions$int.r2
	))
}

# Fold a square matrix by its diagonal and merge everythin in LT. RT and diagonal are zeroed.
# If two elements overlap, the maximum (in absolute value) is retained.
# Returns a lower triangle matrix, whose elements are merged with the upper triangle.
foldMatrix <- function(mat, op="max") {
	diag(mat) <- 0
	tmp1 <- mat
	tmp2 <- mat
	tmp1[lower.tri(tmp1)] <- 0
	tmp1 <- t(tmp1)
	tmp2[upper.tri(tmp2)] <- 0
	switch(op
		,"max"=return(ifelse(abs(tmp1) > abs(tmp2), tmp1, tmp2))
		,"and"=return(ifelse(tmp1 != 0 & tmp2 != 0, 1, 0))
		)
#	stopifnot(!any(tmp1 != 0 & tmp2 != 0))
#	return(tmp1 + tmp2)
}

keep.only.maximum <- function(mat) {
	diag(mat) <- 0
	tmp1 <- mat
	tmp2 <- mat
	tmp1[lower.tri(tmp1)] <- 0
	tmp1 <- t(tmp1)
	tmp2[upper.tri(tmp2)] <- 0
	w <- ifelse(abs(tmp1) > abs(tmp2), 1, 2)
	w[upper.tri(w)] <- t(w)[upper.tri(w)]	# symmetrize
	
	ifelse(w == 1, t(tmp1), tmp2)
}

find.correct.spurious <- function(estimated.model.coefs, true.model, excluded.interactions=NULL) {
	opt.mat <- estimated.model.coefs
	opt.mat$sp <- keep.only.maximum(opt.mat$sp)
#	mar <- par("mar")
#	par(mar=rep(0, 4))
#	image(opt.mat$sp, asp=1, zlim=c(-2, 2), col=c(rev(heat.colors(4)), "white", heat.colors(4)))
#	par(mar=mar)

	# copy upper triangle to lower, of the true & estim interaction matrices
	estim.folded <- foldMatrix(opt.mat$sp)
	true.folded <- foldMatrix(true.model$model$sp)

	spuriousEstimation <- (estim.folded != 0) & (true.folded == 0)
	missedEstimation <- (estim.folded == 0) & (true.folded != 0)
	severeMissedEstimation <- (estim.folded == 0) & (abs(true.folded) >= 0.5)	# missed that were higher than 0.5
	correctEstimation <- (estim.folded != 0) & (true.folded != 0)	# this includes both correct and wrong directions
	correctEstimation[upper.tri(correctEstimation)] <- t(correctEstimation)[upper.tri(correctEstimation)]	# symmetrize
	wrongDirection <- (opt.mat$sp != 0) & correctEstimation & (true.model$model$sp == 0)
	correctDirection <- (opt.mat$sp != 0) & correctEstimation & (true.model$model$sp != 0)
	
	if(!is.null(excluded.interactions)) {
		excludedInteractions <- (estim.folded == 0) & (true.folded != 0) & foldMatrix(excluded.interactions, op="and") != 0
		excludedTrueNegatives <- (estim.folded == 0) & (true.folded == 0) & foldMatrix(excluded.interactions, op="and") != 0
	} else {
		excludedInteractions <- NULL
		excludedTrueNegatives <- NULL
	}

	#estimated <- opt.mat$sp != 0
	# copy upper triangle to lower, of the wrong direction matrix
	wrongdir.folded <- ifelse(foldMatrix(wrongDirection) != 0, TRUE, FALSE)

	
	out <- list(spuriousEstimation=spuriousEstimation, missedEstimation=missedEstimation, severeMissedEstimation=severeMissedEstimation
		, correctEstimation=correctEstimation, wrongDirection=wrongDirection, correctDirection=correctDirection
		, nodirection.estimated=estim.folded, nodirection.true=true.folded, nodirection.wrongDirection=wrongdir.folded
		, excludedInteractions=excludedInteractions, excludedTrueNegatives=excludedTrueNegatives)

	were.estimated <- (estim.folded != 0) & (true.folded != 0)
	if(sum(were.estimated) >= 2) {
		match.coefs <- stats::lm(estim.folded[were.estimated] ~ true.folded[were.estimated])
		out$int.beta <- as.numeric(stats::coef(match.coefs)[2])
		out$int.r2 <- summary(match.coefs)$r.squared
	} else {
		out$int.beta <- NA
		out$int.r2 <- NA
	}
	return(out)
}

