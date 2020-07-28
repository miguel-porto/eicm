#' Fit and select an Explicit Interaction Community Model (EICM)
#'
#' Given species occurrence data and (optionally) measured environmental predictors,
#' fits and selects an EICM that models species occurrence probability as a function of
#' measured predictors, unmeasured predictors (latent variables) and direct species interactions.
#'
#' An Explicit Interaction Community Model (EICM) is a simultaneous equation linear model in which each
#' species model integrates all the other species as predictors, along with measured and latent variables.
#'
#' This is the main function for fitting EICM models, and is preferred over using \code{\link{eicm.fit}} directly.
#'
#' This function conducts the fitting and network topology selection workflow, which includes three stages:
#' 1) estimate latent variable values; 2) make preliminary estimates for species interactions;
#' 3) conduct network topology selection over a reduced model (based on the preliminary estimates).
#' 
#' The selection stage is optional. If not conducted, the species interactions are estimated
#' (all or a subset according to the user-provided constraints), but not selected.
#' See \code{vignette("eicm")} for commented examples on a priori excluding interactions.
#' 
#' Missing data in the response matrix is allowed.
#' 
#' @inheritParams eicm.fit
#' @param penalty the penalty applied to the number of species interactions to include, during variable selection.
#' @param theta.threshold exclude species interactions (from network selection) whose preliminary coefficient (in absolute value)
#'        is lower than this value. This exclusion criterion is cumulative with the other user-defined exclusions.
#' @param latent.lambda the regularization applied to latent variables and respective coefficients
#'        when estimating their values in samples.
#' @param fit.all.with.latents logical. Whether to use the previously estimated latent variables
#'        when estimating the preliminary species interactions.
#' @param true.model for validation purposes only: the true model that has generated the data, to which
#'        the estimated coefficients will be compared in each selection algorithm iteration.
#' @param popsize.sel the population size for the genetic algorithm, expressed as the factor to multiply
#'        by the recommended minimum. Ignored if \code{do.selection=FALSE}.
#' @param n.cores the number of CPU cores to use in the variable selection stage. Ignored if \code{do.selection=FALSE}.
#' @param do.selection logical. Conduct the variable selection stage, over species interaction network topology?
#' @param do.plots logical. Plot diagnostic and trace plots?
#'
#' @return A \code{eicm.list} with the following components:
#' \describe{
#'   \item{true.model:}{a copy of the \code{true.model} argument.}
#'   \item{latents.only:}{the model with only the latent variables estimated.}
#'   \item{fitted.model}{the model with only the species interactions estimated.}
#'   \item{selected.model:}{the final model with all coefficients estimated, after network topology selection.
#'                          This is the "best" model given the selection criterion (which depends on
#'                          \code{regularization} and \code{penalty}.}
#' }
#' When accessing the results, remember to pick the model you want (usually, \code{selected.model}).
#' \code{\link{plot}} automatically picks \code{selected.model} or, if NULL, \code{fitted.model}.
#' @seealso \code{\link{eicm-package}}, \code{\link{eicm.fit}}, \code{\link{plot.eicm}}
#'
#' @examples
#' # refer to the vignette for a more detailed explanation
#' \donttest{
#' # This can take some time to run
#' 
#' # Load the included parameterized model
#' data(truemodel)
#' 
#' # make one realization of the model
#' occurrences <- predict(truemodel, nrepetitions=1)
#' 
#' # Fit and select a model with 2 latent variables to be estimated and all
#' # interactions possible
#' m <- eicm(occurrences, n.latent=2, penalty=4, theta.threshold=0.5, n.cores=2)
#' 
#' plot(m)
#' }
#' @export
#' @importClassesFrom GA ga
#' @useDynLib eicm, .registration = TRUE, .fixes = "SR_"
eicm <- function(occurrences, env=NULL, traits=NULL, intercept=TRUE,	# data
	n.latent=0, forbidden=NULL, allowed=NULL, mask.sp=NULL, exclude.prevalence=0,		# formulation
	regularization=c(ifelse(n.latent > 0, 6, 0.5), 1), regularization.type="hybrid",				# regularization
	penalty=4, theta.threshold=0.5, latent.lambda=1, fit.all.with.latents=TRUE,
	popsize.sel=2, n.cores=parallel::detectCores(),
	true.model=NULL, do.selection=TRUE, do.plots=TRUE) {

	if(!(regularization.type %in% c("ridge", "lasso", "hybrid")))
		stop("Regularization type must be one of: ridge, lasso, hybrid")
		
# TODO adjust regularization defaults, which apparently depend of the # of samples
	attr(regularization, "type") <- regularization.type
	
	options <- eicm.options(mask=list(sp=mask.sp))
	
	if(theta.threshold < 0) theta.threshold <- 0
	if(n.latent > 0) {
		message("Fit latents, no interactions")

		time0 <- system.time(latents.only <- fit.latents(env=env, occurrences=occurrences,
			regularization=c(latent.lambda, 0), regularization.type="ridge",
			nlat=n.latent, fast=FALSE))

		# standardize latents and counter-standardize their coefficients TODO do we need this?
		f <- apply(latents.only$model$samples, 2, stats::sd)
		latents.only$data$env[, -1] <- sweep(latents.only$data$env[, -1, drop=FALSE], 2, f, "/")
		latents.only$model$env[, -1] <- sweep(latents.only$model$env[, -1, drop=FALSE], 2, f, "*")
		latents.only$model$samples <- sweep(latents.only$model$samples, 2, f, "/")
		
		if(!is.null(true.model) && do.plots) {
			grDevices::dev.new(width=12, height=4)

			coefficientComparisonPlot(latents.only$model, true.model, nenv.to.plot=n.latent, nlatent.to.plot=n.latent,
				plot.intercept=TRUE, plot.interactions=FALSE)

		}
	}			

	# fit the full network, so that weak interactions may be discarded
	# TODO: if theta.threshold==0, we don't need to fit all interactions
	if(n.latent > 0) {
		model <- latents.only
		model$model$samples <- NULL
		
		if(fit.all.with.latents) {
			# fix latent values
			message("Fit the model with all interactions, with latent values fixed")
			model$model$sp <- matrix(stats::runif(ncol(model$model$sp) * ncol(model$model$sp), -0.5, 0.5), ncol=ncol(model$model$sp))
			# here we use model$data$env, which has the estimated latent values appended to it
			time1 <- system.time(fitted.model <- eicm.fit(
				model$data$occurrences, env=model$data$env, traits=traits, intercept=FALSE,
				n.latent=0, initial.values=model$model,
				options=options, forbidden=forbidden, allowed=allowed, exclude.prevalence=exclude.prevalence,
				fast=TRUE, regularization=regularization, regularization.type=regularization.type))
		} else {
			message("Fit the model with all interactions, with no latents")
			# remove estimated latent coefficients from environmental matrix
			model$model$env <- model$model$env[, -((ncol(model$model$env) - n.latent + 1):ncol(model$model$env)), drop=FALSE]
			model$model$sp <- matrix(stats::runif(ncol(model$model$sp) * ncol(model$model$sp), -0.5, 0.5), ncol=ncol(model$model$sp))
			# here we use the original "env"
			time1 <- system.time(fitted.model <- eicm.fit(
				occurrences, env=env, traits=traits, intercept=TRUE,
				n.latent=0, initial.values=model$model,
				options=options, forbidden=forbidden, allowed=allowed, exclude.prevalence=exclude.prevalence,
				fast=TRUE, regularization=regularization, regularization.type=regularization.type))
		}

		if(!is.null(true.model) && nrow(fitted.model$model$env) != nrow(true.model$model$env)) warning("Some species are missing from the estimated model")
	} else {	# No latents
		message("Fit the model with all interactions, with no latents")
		# fit the model with everything
		time1 <- system.time(fitted.model <- eicm.fit(occurrences, env=env, traits=traits, intercept=TRUE
			, n.latent=0,
			options=options, forbidden=forbidden, allowed=allowed, exclude.prevalence=exclude.prevalence,
			fast=TRUE, regularization=regularization, regularization.type=regularization.type))
			
		latents.only <- fitted.model
		latents.only$model$sp[, ] <- 0
	}
	
	if(!do.selection) {
		# don't select interactions, just return full network estimates
		out <- list(true.model=true.model, latents.only=latents.only, fitted.model=fitted.model)
		if(exists("time0")) attr(out, "time.fit.latents") <- time0
		if(exists("time1")) attr(out, "time.fit.all") <- time1
		attr(out, "regularization") <- regularization
		class(out) <- "eicm.list"
		return(out)
	}

	if(!is.null(true.model)) {
		if(do.plots) {
			monitor.function <- function(model, plot.interactions) {
				return(
					coefficientComparisonPlot(model, true.model, nenv.to.plot=ncol(fitted.model$model$env) - 1,
						nlatent.to.plot=0, plot.intercept=TRUE,
						excluded.interactions=abs(fitted.model$model$sp) < theta.threshold,
						plot.interactions=plot.interactions)
				)
			}
			selection.monitor <- function(object, bestmodel, worstmodel) {
				stats <- monitor.function(bestmodel, plot.interactions=TRUE)
				fitness.stats <- stats::quantile(stats::na.exclude(object@fitness), probs=c(0, 0.5, 1))
				message(sprintf("\rIt %d (%.0fs, ciT %d+%d/%d) | Fit %.1f %.1f %.1f | %d TP %d FN %d FP",
					object@iter, object@time.took, object@cached, object@informed, object@cache.size,
					fitness.stats[3], fitness.stats[2], fitness.stats[1],
					stats["correct"], stats["missed"], stats["spurious"]))
				utils::flush.console()
			}
		} else {
			selection.monitor <- function(object, bestmodel, worstmodel) {
				stats <- coefficientComparisonPlot(bestmodel, true.model, nenv.to.plot=ncol(fitted.model$model$env) - 1,
					nlatent.to.plot=0, plot.intercept=TRUE,
					excluded.interactions=abs(fitted.model$model$sp) < theta.threshold,
					plot.interactions=TRUE, noplot=TRUE)

				fitness.stats <- stats::quantile(stats::na.exclude(object@fitness), probs=c(0, 0.5, 1))
				message(sprintf("\rIt %d (%.0fs, ciT %d+%d/%d) | Fit %.1f %.1f %.1f | %d TP %d FN %d FP",
					object@iter, object@time.took, object@cached, object@informed, object@cache.size,
					fitness.stats[3], fitness.stats[2], fitness.stats[1],
					stats["correct"], stats["missed"], stats["spurious"]))
				utils::flush.console()
			}
		}
	
		if(do.plots && !is.null(latents.only)) {
			grDevices::dev.new(width=12, height=4)
			monitor.function(fitted.model$model, plot.interactions=TRUE)
			graphics::abline(h=c(-theta.threshold, theta.threshold), lty=2)
		}
		message(sprintf("Discarded %d true interactions (out of %d) from the start."
			, sum(true.model$model$sp != 0 & (abs(fitted.model$model$sp) <= theta.threshold | !fitted.model$model$options$mask$sp))
			, sum(true.model$model$sp != 0)
		))
	} else {
		if(do.plots) {
			selection.monitor <- function(object, bestmodel, worstmodel) {
				if(is.null(bestmodel)) return
				plot.eicm.matrix(bestmodel, type="network")
				gaMonitor.eicm(object, bestmodel, worstmodel)
			}
		} else selection.monitor <- NULL
	}
		
	# now hard threshold over interactions
	# NOTE this mask already includes user-specified exclusions (options, prevalence and forbidden)
	# because they were not estimated in the fitted.model
	masksp <- abs(fitted.model$model$sp) > theta.threshold
	tmp <- !masksp & fitted.model$model$options$mask$sp
	diag(tmp) <- FALSE
	message(sprintf("Discarded %d parameters (out of %d) by hard thresholding with |theta| < %.2f.",
		sum(tmp), sum(fitted.model$model$sp != 0), theta.threshold))
#		sum(tmp), sum(fitted.model$model$options$mask$sp), theta.threshold))

#	if(do.plots)
#		grDevices::dev.new(width=12, height=4)
	
	latents.only.copy <- latents.only
	latents.only.copy$model$samples <- NULL
	latents.only.copy$model$sp <- fitted.model$model$sp
	latents.only.copy$model$options$mask$sp <- fitted.model$model$options$mask$sp

	time2 <- system.time(var.selection <- model.selection.network(latents.only.copy, regularization=regularization
		, regularization.type=regularization.type
		, penalty=penalty, masksp=masksp, select.direction=TRUE, exclude.prevalence=exclude.prevalence
		, fast=TRUE, optim.method = "ucminf", optim.control = list(trace=0, maxeval=10000, gradstep=c(0.001, 0.001), grtol=0.1)
		#, fast=FALSE, optim.method = "L-BFGS-B", optim.control = list(trace=0, maxit=10000, ndeps=0.001)
		, parallel=n.cores
		, maxit.stagnated=50, pmutation = 0.01, popsize.factor=popsize.sel
		, monitor=selection.monitor
	))
	
	bits <- var.selection@solution[1, ]
	selected.model <- var.selection@allEvaluations[[paste(bits, collapse="")]]

	selected.model <- list(
		data=eicm.data(occurrences=occurrences, env=env)
		, model=eicm.matrix(selected.model$env, sp.coefs=selected.model$sp, latent=selected.model$samples, options=selected.model$options)
	)
	attr(selected.model, "regularization") <- regularization
	class(selected.model) <- "eicm"

	if(is.function(selection.monitor))
		selection.monitor(NULL, selected.model$model, NULL)

	out <- list(true.model=true.model, fitted.model=fitted.model, latents.only=latents.only, selected.model=selected.model)
	if(exists("time0")) attr(out, "time.fit.latents") <- time0
	if(exists("time1")) attr(out, "time.fit.all") <- time1
	if(exists("time2")) attr(out, "time.selection") <- time2
	attr(out, "regularization") <- regularization
	attr(out, "penalty") <- penalty
	attr(out, "theta.threshold") <- theta.threshold
	attr(out, "latent.lambda") <- latent.lambda
	class(out) <- "eicm.list"
	return(out)
}

# Utility function to estimate latents possibly with given species coefficients
fit.latents <- function(env, occurrences, regularization, regularization.type, nlat, envcoefs=NULL, spcoefs=NULL, fast=FALSE) {
	if(is.null(env))
		env <- matrix(0, nrow=nrow(occurrences), ncol=0)

	# no species interactions
	options <- eicm.options(mask=list(env=1L, sp=0L))

	# This is to estimate latents in new data, from sparse data! We here fix the species coefficients.
	if(!is.null(spcoefs) || !is.null(envcoefs)) {
		if(is.null(spcoefs))
			spcoefs <- matrix(0, ncol=ncol(occurrences), nrow=ncol(occurrences))
		if(is.null(envcoefs))
			envcoefs <- matrix(0, ncol=ncol(env) + nlat, nrow=ncol(occurrences))
			
		options <- options + list(offset=list(env=envcoefs, sp=spcoefs))
	}

	model <- eicm.fit(occurrences, env=env, intercept=TRUE
		, options=options
		, fast=fast
		, n.latent=nlat		
		, regularization=regularization, regularization.type=regularization.type)
	return(model)
}

gaMonitor.eicm <- function(object, bestmodel, worstmodel) {
	fitness.stats <- stats::quantile(stats::na.exclude(object@fitness), probs=c(0, 0.5, 1))
	nterms.stats <- stats::quantile(apply(object@population, 1, sum), probs=c(0, 0.5, 1))
	message(sprintf("\rIt %d (%.0fs, ciT %d+%d/%d) | Fit %.1f %.1f %.1f | #term %d %d %d",
		object@iter, object@time.took, object@cached, object@informed, object@cache.size
		, fitness.stats[3], fitness.stats[2], fitness.stats[1], nterms.stats[1], as.integer(nterms.stats[2]), nterms.stats[3]))
	utils::flush.console()
}

