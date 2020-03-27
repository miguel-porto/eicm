#' Likelihood profiles for EICMs
#'
#' Computes the profile (penalized) likelihood for all (or only one) estimated parameters in a EICM model.
#'
#' Likelihod profiles will use the same regularization settings that were used in model fitting.
#' 
#' @param fitted the fitted EICM model.
#' @param all.pars logical. Compute for all model parameters?
#' @param parmatrix if all.pars=FALSE, in which matrix is the parameter of interest, "env" or "sp"?
#' @param species if all.pars=FALSE, in which row of \code{parmatrix} is the parameter of interest?
#' @param parameter if all.pars=FALSE, in which column of \code{parmatrix} is the parameter of interest?
#' @param step the step increments/decrements at which to compute the likelihood profile points.
#' @param ncores the number of CPU cores to use when computing profiles for all parameters.
#' @param alpha highest significance level that will be guaranteed for this profile.
#' @param ... additional argument(s) for methods
#' 
#' @note Confidence intervals should \strong{not} be computed on a model whose terms have been selected.
#' 
#' This function is optimized for computing profiles of multiple parameters simultaneously (in parallel).
#' 
#' @return The same model object updated with a new \code{profile} component.
#' 
#' @examples
#' \donttest{
#' # load the included parameterized model
#' data(truemodel)
#' 
#' # realize the model
#' occurrences <- predict(truemodel, nrepetitions=1)
#'
#' # fit the model without species interactions
#' fitted <- eicm(occurrences, n.latent=2, mask.sp=0, do.selection=FALSE)$fitted.model
#' 
#' # compute likelihood profiles for all parameters
#' fitted <- profile(fitted, ncores=2)
#' 
#' # plot the first 9 profiles
#' par(mfrow=c(3, 3))
#' dummy <- lapply(fitted$profile[1:9], plot)
#' }
#' @export
profile.eicm <- function(fitted, all.pars=TRUE, parmatrix, species, parameter, 
	step=0.3, ncores=parallel::detectCores(), alpha=0.01, ...) {

	i <- NULL
	`%DO%` <- foreach::`%dopar%`
	
	if(all.pars) {
		testsEnv <- data.frame(which(fitted$model$env != 0, arr.ind=TRUE), parmatrix=rep("env", sum(fitted$model$env != 0)), row.names=NULL, stringsAsFactors=FALSE)
		testsSp <- data.frame(which(fitted$model$sp != 0, arr.ind=TRUE), parmatrix=rep("sp", sum(fitted$model$sp != 0)), row.names=NULL, stringsAsFactors=FALSE)

		tests <- as.matrix(rbind(testsEnv, testsSp))
		
		if(!is.null(ncores)) {
			cls <- parallel::makeCluster(ncores)
			progress <- function(n) message(sprintf("[Likelihood profiles] %.1f%% complete (%d of %d)", round(n / nrow(tests) * 100), n, nrow(tests)))
			doSNOW::registerDoSNOW(cls)
			opts <- list(progress=progress)
			mcoptions <- list(preschedule = FALSE)
		}
		message("Computing ", nrow(tests), " likelihood profiles with ", foreach::getDoParWorkers()," cores...")
		
		w <- iterators::icount(nrow(tests))
		all.ll <- foreach::foreach(i=w, .options.multicore = mcoptions, .options.snow=opts) %DO% {
			pars <- tests[i, ]
			profile.eicm(fitted, all.pars=FALSE, species=as.integer(pars[1]), parameter=as.integer(pars[2])
				, parmatrix=as.character(pars[3]), step=step, ncores=NULL, alpha=alpha)
		}
		
		if(!is.null(ncores))
			parallel::stopCluster(cls)
		
		names(all.ll) <- paste(
			c(rownames(fitted$model$env)[testsEnv[, 1]], rownames(fitted$model$sp)[testsSp[, 1]])
			,"~", c(colnames(fitted$model$env)[testsEnv[, 2]], colnames(fitted$model$sp)[testsSp[, 2]])
		)
		
		fitted$profile <- all.ll
		return(fitted)
	}
	
	# We NEED the profile to cross the significance threshold, so we add a tolerance just in case...
	percentile <- (stats::qchisq(1 - alpha, 1) / 2) + 0.5
	nenv <- ncol(fitted$model$env)
	nspecies <- ncol(fitted$model$sp)
	# compute confidence intervals for environmental coefficients
	mask <- fitted$model
	zerooffset <- list(
		env=matrix(0, ncol=nenv, nrow=nspecies)
		, sp=matrix(0, ncol=nspecies, nrow=nspecies)
	)
	
	parval <- switch(parmatrix,
		"env"=mask$env[species, parameter],
		"sp"=mask$sp[species, parameter])

	# the max loglik	
	ll <- compute.one.logLik(parval, fitted, parmatrix, species, parameter
		, zerooffset=zerooffset)
	values <- parval
	
	# stopping criteria for assessing profile flatness (at the head or tail, separately)
	wnd <- 3	# the number of points (minus one) to fit a LM (with the heading or trailing likelihood values)
	slope.threshold <- 0.05		# the minimum slope to consider as non-flat
	
	# compute profile points to the left, step by step
	nrep <- 0
	repeat {
		values <- c(values[1] - step, values)
		
		newll <- compute.one.logLik(values[1], fitted, parmatrix, species, parameter
			, zerooffset=zerooffset)

		ll <- c(newll, ll)

		#if((max(ll) - ll[1]) >= percentile || diff(range(values)) >= stop.threshold) break
		if((max(ll) - ll[1]) >= percentile) break
		
		if(length(ll) > wnd && nrep > wnd + 2) {
			# did it stabilize?
			if(abs(stats::coef(stats::lm(ll[1:(1 + wnd)] ~ I(seq_len(wnd + 1))))[2]) < slope.threshold) break
		}
		nrep <- nrep + 1
	}

	nrep <- 0
	repeat {
		values <- c(values, values[length(values)] + step)
			
		newll <- compute.one.logLik(values[length(values)], fitted, parmatrix, species, parameter
			, zerooffset=zerooffset)

		ll <- c(ll, newll)

		#if((max(ll) - ll[length(ll)]) >= percentile || diff(range(values)) >= stop.threshold * 2) break
		if((max(ll) - ll[length(ll)]) >= percentile) break

		if(nrep > wnd + 2) {
			# did it stabilize?
			if(abs(stats::coef(stats::lm(ll[(length(ll) - wnd):length(ll)] ~ I(seq_len(wnd + 1))))[2]) < slope.threshold) break
		}
		nrep <- nrep + 1
	}
	
	out <- matrix(c(values, ll), ncol=2, dimnames=list(NULL, c("values", "logLik")))
	attr(out, "parmatrix") <- parmatrix
	attr(out, "parameter") <- c(rownames(fitted$model[[parmatrix]])[species], colnames(fitted$model[[parmatrix]])[parameter])
	attr(out, "maxLogLik") <- logLik.eicm(fitted)
	attr(out, "parEstimate") <- parval
	class(out) <- "profile.eicm"
	return(out)
}


compute.one.logLik <- function(parval, fitted, parmatrix, species, parameter, zerooffset=NULL) {
	if(is.null(zerooffset)) {
		nenv <- ncol(fitted$model$env)
		nspecies <- ncol(fitted$model$sp)

		zerooffset <- list(
			env=matrix(0, ncol=nenv, nrow=nspecies)
			, sp=matrix(0, ncol=nspecies, nrow=nspecies))
	}

	# fix the target parameter value
	zerooffset[[parmatrix]][species, parameter] <- parval
	#mbOpts <- modifyList(fitted$model$mbo, list(mask=fitted$model, offset=zerooffset))
	mbOpts <- eicm.options(mask=fitted$model, offset=zerooffset)

	ps1 <- eicm.fit(fitted$data$occurrences, fitted$data$env, intercept=TRUE, options=mbOpts,
		regularization=attr(fitted, "regularization"), initial.values=fitted$model)
		
	gc()
	return(as.numeric(logLik.eicm(ps1)))
}

#' Confidence intervals for EICM parameters
#'
#' Computes the profile (penalized) likelihood confidence intervals for all estimated parameters in a EICM model.
#' If the likelihood profiles are not computed yet, they will be computed first.
#'
#' @param object the fitted EICM model.
#' @param parm currently unused.
#' @param level the confidence level required.
#' @param step the step increments/decrements at which to compute the likelihood profile points.
#' @param ncores the number of CPU cores to use when computing profiles for all parameters.
#' @param ... additional argument(s) for methods
#'
#' @return The same model object with a new \code{confint} component.
#'
#' @examples
#' \donttest{
#' # load the included parameterized model
#' data(truemodel)
#' 
#' # realize the model
#' occurrences <- predict(truemodel, nrepetitions=1)
#'
#' # fit the model without species interactions
#' fitted <- eicm(occurrences, n.latent=2, mask.sp=0, do.selection=FALSE)$fitted.model
#' 
#' # compute confidence intervals for all parameters
#' # this updates the fitted model with the confints
#' fitted <- confint(fitted, ncores=2)
#' 
#' # plot the confidence intervals
#' plot(fitted, type="confint")
#' }
#' @export
confint.eicm <- function(object, parm, level=0.99, step=0.3, ncores=parallel::detectCores(), ...) {
	if(is.null(object$profile)) {
		message("Computing profile likelihood...")
		prof <- profile.eicm(object, all.pars=TRUE, step=step, alpha=1 - level, ncores=ncores)$profile
	} else
		prof <- object$profile
		
	percentile <- stats::qchisq(level, 1) / 2
	# TODO do we assume the new estimate and new likelihood, or keep the original???
	# here we make both: update the estimate if the likelihood is higher than original
	ori.loglik <- logLik.eicm(object)
	ori.estimate <- sapply(prof, function(tab) attr(tab, "parEstimate"))
	new.estimate <- sapply(prof, function(tab) tab[which.max(tab[, 2]), 1])
	new.ll <- sapply(prof, function(tab) max(tab[, 2]))
	estimate <- ifelse(ori.loglik < new.ll, new.estimate, ori.estimate)

	new.maxll <- max(new.ll)
	# TODO use global max logLik or use max logLik for each parameter??
	# cut.thr <- new.maxll - percentile
	
	out <- matrix(ncol=3, nrow=length(prof), dimnames=list(names(prof), c("Estimate", "Lower CI", "Upper CI")))
	#out <- matrix(ncol=3, nrow=length(prof), dimnames=list(paste0("[", seq_along(prof), "] ", names(prof)), c("Estimate", "Lower CI", "Upper CI")))
	out[, 1] <- estimate
	for(i in seq_along(prof)) {
		# update estimates with new values
		attr(prof[[i]], "parEstimate") <- estimate[i]
		attr(prof[[i]], "maxLogLik") <- new.maxll
		
		tab <- prof[[i]]
		cut.thr <- max(tab[, 2]) - percentile
		ll.c <- tab[, 2] - cut.thr
		roots <- which(ll.c[-1] * ll.c[-length(ll.c)] < 0)
		if(length(roots) == 0) {	# if apparently flat likelihood try to increase grid span
			# no, what happens is that the profile is all below the threshold, so we need more resolution around the estimate
			step <- tab[2, 1] - tab[1, 1]
			oritab <- tab
			parmatrix <- attr(oritab, "parmatrix")
			parR <- which(rownames(object$model[[parmatrix]]) == attr(oritab, "parameter")[1])
			parC <- which(colnames(object$model[[parmatrix]]) == attr(oritab, "parameter")[2])
			
			message(sprintf("No roots for parameter %s %s[%d, %d] (%s); extending profile...", parmatrix, parR, parC, attr(oritab, "parameter")))
			
			nrepeats <- 0
			if(max(tab[, 2] - cut.thr) < 0) {
				stop("Increase grid points")
			# TODO shouldn't stop here, but compute needed values
				wmax <- which.max(tab[, 2])
				w2max <- which.max(tab[-which.max(tab[, 2]), 2])
				newrange <- c(max(tab[, 2]), max(tab[-which.max(tab[, 2]), 2]))
			} else {
				# stopping criteria for assessing profile flatness (at the head or tail, separately)
				wnd <- 3	# the number of points (minus one) to fit a LM (with the heading or trailing likelihood values)
				slope.threshold <- 0.05		# the minimum slope to consider as non-flat
# TODO  a kind of "binary search" would be better, start with a long shot
				repeat {
					if(tab[1, 2] > tab[nrow(tab), 2]) {	# increase at the end
						add.ll <- compute.one.logLik(tab[nrow(tab), 1] + step, object, parmatrix, parR, parC)
						prof[[i]] <- rbind(prof[[i]], c(tab[nrow(tab), 1] + step, add.ll))
						inc.tail <- TRUE
					} else {	# increase at the beginning
						add.ll <- compute.one.logLik(tab[1, 1] - step, object, parmatrix, parR, parC)
						prof[[i]] <- rbind(c(tab[1, 1] - step, add.ll), prof[[i]])
						inc.tail <- FALSE
					}
					tab <- prof[[i]]
					ll.c <- tab[, 2] - cut.thr
					roots <- which(ll.c[-1] * ll.c[-length(ll.c)] < 0)

					if(length(roots) > 0) break;
					
					if(nrepeats > 10) {
						# did it stabilize?
						if(!inc.tail && abs(stats::coef(stats::lm(tab[1:(1 + wnd), 2] ~ I(seq_len(wnd + 1))))[2]) < slope.threshold) break
						if(inc.tail && abs(stats::coef(stats::lm(tab[(length(tab[, 2]) - wnd):length(tab[, 2]), 2] ~ I(seq_len(wnd + 1))))[2]) < slope.threshold) break
					}

					# TODO a better stopping criterion. This sometimes doesn't go far enough
					nrepeats <- nrepeats + 1
				}
			}
			
			# update the profile in return object
			mostattributes(prof[[i]]) <- attributes(oritab)
			dim(prof[[i]]) <- c(length(prof[[i]]) / 2, 2)
			colnames(prof[[i]]) <- colnames(oritab)
			class(prof[[i]]) <- "profile.eicm"
			tab <- prof[[i]]

#			if(nrepeats > 6) next	# if we still couldn't find an intersect, then it is likely flat, return NA confint
			if(length(roots) == 0) next	# if we still couldn't find an intersect, then it is likely flat, return NA confint
		}
#		if(length(roots) > 2) {

		# we found a root, go on
		# find the tightest interval where the estimate is located 
		# (there might be some unstability in the curve, nothing to worry about)
		where.is.estimate <- findInterval(estimate[i], tab[, 1])
		where.is.estimate <- findInterval(where.is.estimate, roots, rightmost.closed=TRUE)
		if(where.is.estimate == length(roots))	# estimate is past all roots, so +Inf
			roots <- roots[length(roots)]
		else if(where.is.estimate == 0)			# estimate is lower than all roots, so -Inf
			roots <- roots[1]
		else				# estimate is within two roots, fetch those (the nearest)
			roots <- roots[c(where.is.estimate, where.is.estimate + 1)]
#		}

		if(any(is.na(roots))) {
			stop("Some unexpected error")
			# print(tab)
			# print(ll.c[-1] * ll.c[-length(ll.c)])
			# print(which(ll.c[-1] * ll.c[-length(ll.c)] < 0))
		}
		# TODO instead of stopping, calculate the needed value as above
		if(max(roots) + 1 > nrow(tab)) stop("I need one more value in parameter", names(prof)[i])
		
		# interpolate linearly between to grid points, at the intersection
		from <- tab[roots, , drop=FALSE]
		to <- tab[roots + 1, , drop=FALSE]
		dx <- to[, 1] - from[, 1]
		dy <- to[, 2] - from[, 2]
		m <- dy / dx
		b <- from[, 2] - m * from[, 1]
		x <- (cut.thr - b) / m
	
		if(length(x) == 1) {
			if(ll.c[roots] < 0)		# note that we keep the original estimate, cause this is a sill likelihood
				out[i, ] <- c(ori.estimate[i], x, Inf)
			else
				out[i, ] <- c(ori.estimate[i], -Inf, x)
		} else
			out[i, 2:3] <- x
	}

	object$profile <- prof
	object$confint <- out
	attr(object$confint, "level") <- level
	class(object$confint) <- "confint.eicm"
	return(object)
}

