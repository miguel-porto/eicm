##############################################################################
#                                                                            #
#                        GENETIC ALGORITHMS in R                             #
#                       modified version for eicm                            #
#     adapted from https://cran.r-project.org/package=GA version 3.1.1.      #
##############################################################################

my.ga <- function(
               fitness, criterion, ...,
               lower, upper, nBits,
               population = "gabin_Population",
               selection = "gabin_lrSelection",
               crossover = "gabin_spCrossover",
               mutation,
               popSize = 50, 
               pcrossover = 0.8, 
               pmutation = 0.1, 
               elitism = base::max(1, round(popSize*0.05)), 
               updatePop = FALSE,
               postFitness = NULL,
               maxiter = 100,
               run = maxiter,
               maxFitness = Inf,
               names = NULL,
               suggestions = NULL,
               optim = FALSE,
               optimArgs = list(method = "L-BFGS-B", 
                                poptim = 0.05,
                                pressel = 0.5,
                                control = list(fnscale = -1, maxit = 100)),
               keepBest = FALSE,
               parallel = FALSE,
               monitor = NULL,
############## START CHANGED CODE #########################
               max.cached = 15000,
############### END CHANGED CODE ##########################
               seed = NULL) 
{
  call <- match.call()
  
  
  if(!is.function(population)) population <- get(population, envir=environment(GA::ga))
  if(!is.function(selection))  selection  <- get(selection, envir=environment(GA::ga))
  if(!is.function(crossover))  crossover  <- get(crossover, envir=environment(GA::ga))
  if(!is.function(mutation))   mutation   <- get(mutation, envir=environment(GA::ga))
  
  if(missing(fitness))
    { stop("A fitness function must be provided") }
  if(!is.function(fitness)) 
    { stop("A fitness function must be provided") }
  if(popSize < 10) 
    { warning("The population size is less than 10.") }
  if(maxiter < 1) 
    { stop("The maximum number of iterations must be at least 1.") }
  if(elitism > popSize) 
    { stop("The elitism cannot be larger that population size.") }
  if(pcrossover < 0 | pcrossover > 1)
    { stop("Probability of crossover must be between 0 and 1.") }
  if(is.numeric(pmutation))
  { 
    if(pmutation < 0 | pmutation > 1)
      { stop("If numeric probability of mutation must be between 0 and 1.") }
    else if(!is.function(population))
           { stop("pmutation must be a numeric value in (0,1) or a function.") }
  }

  # check for min and max arguments instead of lower and upper
  callArgs <- list(...)
  penalty <- callArgs$penalty
  callArgs$penalty <- NULL
  
  modelobject <- callArgs$modelobject
  callArgs$modelobject <- NULL

  if(any("min" %in% names(callArgs)))
  {
    lower <- callArgs$min
    callArgs$min <- NULL
    warning("'min' arg is deprecated. Use 'lower' instead.")
  }
  if(any("max" %in% names(callArgs)))
  {
    upper <- callArgs$max
    callArgs$max <- NULL
    warning("'max' arg is deprecated. Use 'upper' instead.")
  }

  if(missing(lower) & missing(upper) & missing(nBits))
    { stop("A lower and upper range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }

	nBits <- as.vector(nBits)[1]
	lower <- upper <- NA
	nvars <- nBits 
	if(is.null(names))
		names <- paste0("x", 1:nvars)

  # check suggestions
  if(is.null(suggestions))
    { suggestions <- matrix(nrow = 0, ncol = nvars) }
  else
    { if(is.vector(suggestions)) 
        { if(nvars > 1) suggestions <- matrix(suggestions, nrow = 1)
          else          suggestions <- matrix(suggestions, ncol = 1) }
      else
        { suggestions <- as.matrix(suggestions) }
      if(nvars != ncol(suggestions))
        stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
    }
  
###########################################################
############## START CHANGED CODE #########################
###########################################################
	if(is.numeric(parallel) && parallel > 1) {
		cls <- parallel::makeCluster(parallel)
#		doParallel::registerDoParallel(cls)	# doParallel does not support progress bar!
		doSNOW::registerDoSNOW(cls)
		`%DO%` <- foreach::`%dopar%`
	} else {
		if(!is.numeric(parallel)) stop("parallel must be the number of cores.")
		cls <- NULL
		`%DO%` <- foreach::`%do%`
	}
###########################################################
############### END CHANGED CODE ##########################
###########################################################
  # set seed for reproducibility  
  if(!is.null(seed)) set.seed(seed)
  i. <- NULL # dummy to trick R CMD check 

  fitnessSummary <- matrix(as.double(NA), nrow = maxiter, ncol = 6)
  colnames(fitnessSummary) <- names(GA::gaSummary(stats::rnorm(10)))
  bestSol <- if(keepBest) vector(mode = "list", length = maxiter)
             else         list()
  Fitness <- rep(NA, popSize)

  object <- methods::new("ga.eicm", 
                call = call, 
                lower = lower, 
                upper = upper, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1, 
                maxiter = maxiter,
                suggestions = suggestions,
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = if(is.numeric(pmutation)) pmutation else NA,
                fitness = Fitness, 
                summary = fitnessSummary,
                bestSol = bestSol)
  
  if(maxiter == 0)
    return(object)
  
  # generate beginning population
  Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  ng <- min(nrow(suggestions), popSize)
  if(ng > 0) # use suggestion if provided
    { Pop[1:ng,] <- suggestions }
  # fill the rest with a random population
  if(popSize > ng)
    { Pop[(ng+1):popSize,] <- population(object)[1:(popSize-ng),] }
  object@population <- Pop

############## START CHANGED CODE #########################
	evaluation.history <- list()	# that's the fitted model cache
	message(sprintf("Percentage of mutation: %.3f; poisson mean: %.3f", pmutation, attr(mutation, "mean")))
############### END CHANGED CODE ##########################

  # start iterations
  for(iter in seq_len(maxiter))
     {
      object@iter <- iter

############ CHANGED LINE
		start.time <- Sys.time()
      # evalute fitness function (when needed) 
      if(!parallel) 	# TODO DEPRECATED
        { for(i in seq_len(popSize))
             if(is.na(Fitness[i]))
               { fit <- do.call(fitness, c(list(Pop[i,]), callArgs)) 
                 if(updatePop)
                   Pop[i,] <- attributes(fit)[[1]]
                 Fitness[i] <- fit
               }
      }
###########################################################
############## START CHANGED CODE #########################
###########################################################
      else {
#		if(gc()["Vcells", 2] > 1000) {
		if((iter %% 25) == 0 && !is.null(cls)) {
			message("Restarting cluster to free memory.")
			parallel::stopCluster(cls)
			cls <- parallel::makeCluster(parallel)
#			doParallel::registerDoParallel(cls)
			doSNOW::registerDoSNOW(cls)
		}

      	current.pop.strings <- apply(Pop, 1, paste, collapse="")

 		# returns fitness if it has survived or if it's cached, NA otherwise.
 		tmp.fit <- vapply(seq_len(popSize), FUN=function(i.) {
			if(is.na(Fitness[i.])) {	# individual has changed, recalculate
				if(!is.null(evaluation.history[[current.pop.strings[i.]]])	# result is cached
					) {	
					return(attr(evaluation.history[[current.pop.strings[i.]]], "criterion"))	# there's an equal model, return fitness
				}
					return(NA)
			} else 
				return(Fitness[i.])
 		}, FUN.VALUE=1)

		object@survived <- 0
		object@cached <- sum(!sapply(tmp.fit, is.na))

 		individuals.to.eval <- sapply(tmp.fit, is.na)

 		# the unique bit strings to evaluate fitness
 		bitstrings.to.eval <- unique(Pop[individuals.to.eval, , drop=FALSE])

		# returns a similar model or NA is there is none 		
		if(length(evaluation.history) == 0)
			similar.models <- rep(NA, nrow(bitstrings.to.eval))
		else {
			similar.models <- .Call(SR__getMostSimilarModel, bitstrings.to.eval, evaluation.history)
			similar.models <- lapply(similar.models, function(i) {
				if(is.na(i))
					return(NA)
				else
					return(evaluation.history[[i]])
			})
		}

		object@informed <- sum(!sapply(similar.models, identical, NA))
		object@evaluated <- length(similar.models)
		
 		if(interactive()) {
	 		message(sprintf("\rFitting %d models...             ", length(similar.models)), appendLF=is.null(cls))
#	 		cat(sprintf("\rFitting %d models...             ", length(similar.models)))
			progress <- function(n) {
				message(sprintf("\rFitting %d models... %d complete", length(similar.models), n), appendLF=FALSE)
#				cat(sprintf("\rFitting %d models... %d complete", length(similar.models), n))
#				writeLines(sprintf("\rFitting %d models... %d complete", length(similar.models), n), conn, sep="")
#				flush(conn)
			}
		} else
			progress <- function(n) {}
#		conn <- base::stdout()
#		parallel::clusterExport(cls, c("progress", "conn"), envir=environment())
 		# no duplicates are evaluated

		# NOTE we need the following loop (a workaround) because, for some reason, sometimes a
		# corrupted occurence matrix is passed to the parallel workers.
		# TODO investigate whether this is due to a memory leak, or an external bug of doParallel.
		tmpEval <- rep(list(NA), length(similar.models))
		while(any(sapply(tmpEval, identical, NA))) {
		
			NAs <- sapply(tmpEval, identical, NA)
			if(sum(NAs) < length(tmpEval))
				message(sprintf("Repeating %d models...", sum(NAs)))
			
			tmpEval[NAs] <- foreach::foreach(i. = seq_along(similar.models[NAs]), .options.multicore=list(preschedule = FALSE)
				, .options.snow=list(progress=progress)) %DO% {
	#		tmpEval <- lapply(seq_along(similar.models), function(i.) {
					if(!all(sort(unique(as.vector(modelobject$data$occurrences))) == c(0, 1))) {
						# corrupted matrix! TODO: why does it happen?
						# saveRDS(modelobject$data$occurrences, file="debug-GA2error.rds")
						return(NA)
					}

					if(identical(similar.models[[i.]], NA)) {
						suppressMessages(
							tmp <- do.call(fitness, args=c(string=list(bitstrings.to.eval[i., ]),
								modelobject=list(modelobject),
								callArgs))	# no similar models, that's for the first evaluation
						)
					} else {
						suppressMessages(
							tmp <- do.call(fitness, args=c(string=list(bitstrings.to.eval[i., ]),
								modelobject=list(modelobject),
								similar.model=list(similar.models[[i.]]),
								callArgs)) # similar model informs fitting
						)
					}

					if(inherits(tmp, "CV.models")) {
					# CROSS VALIDATION TODO
						rmse <- sapply(tmp, function(m) {
							if(!attr(m, "possible")) return(NA)
							
							testdata <- attr(m, "testingdata")
							m$data$env <- testdata$env

							return(logLik.eicm(m, testdata$presences))
							#preds <- predict(m, 1000)
							#sum((preds - testdata$presences) ^ 2)
						})
						criterion.value <- mean(rmse)
						if(is.na(criterion.value))
							out <- NA
						else
							out <- tmp[[1]]$model
					} else {
					# AIC-like
						#if(!attr(tmp, "possible")) {
						#	stop("Fitting not possible or not improved", list(i., bitstrings.to.eval, similar.models, tmp))
						#	out <- NA
						#	criterion.value <- NA
						#} else {
							out <- tmp$model
							criterion.value <- -criterion(tmp) -
								penalty * sum(tmp$model$options$mask$sp != 0)
						#}
					}
					
					if(identical(out, NA) || identical(criterion.value, NA))
						stop("This should never happen")
					else
						attr(out, "criterion") <- criterion.value

					rm(tmp)
					gc()
	#				progress(i., conn)
					return(out)
			}
		}
#	)
		
		# now merge with fitness vector
		string.eval <- apply(bitstrings.to.eval, 1, paste, collapse="")

		if(any(individuals.to.eval))
			tmp.fit[individuals.to.eval] <- sapply(tmpEval, attr, "criterion")[match(current.pop.strings[individuals.to.eval], string.eval)]

		Fitness <- tmp.fit
		
		# update cached
		for(i in seq_along(tmpEval)) {
			cached.copy <- evaluation.history[[string.eval[i]]]
			if(is.null(cached.copy)) {
					# clean up cache to avoid out of memory
					if(length(evaluation.history) > max.cached) {
						tokeep <- which(names(evaluation.history) %in% current.pop.strings)
						evaluation.history[setdiff(1:5000, tokeep)] <- NULL
						gc()
					}
		
					evaluation.history[[string.eval[i]]] <- tmpEval[[i]]
			}
		}
#		splitted.bits <- strsplit(names(evaluation.history), NULL)
	}
###########################################################
############### END CHANGED CODE ##########################
###########################################################
      
      # update object
      object@population <- Pop
      object@fitness <- Fitness
      
      # update iterations summary
      fitnessSummary[iter,] <- GA::gaSummary(object@fitness)
      object@summary <- fitnessSummary
############ CHANGED CODE
      object@time.took <- as.numeric(difftime(Sys.time(), start.time, units="secs")) 
      object@cache.size <- length(evaluation.history)
  
      if(is.function(monitor)) {
        monitor(object, evaluation.history[[current.pop.strings[which.max(Fitness)]]], evaluation.history[[current.pop.strings[which.min(Fitness)]]])
      } else {
	      gaMonitor.eicm(object, evaluation.history[[current.pop.strings[which.max(Fitness)]]], evaluation.history[[current.pop.strings[which.min(Fitness)]]])
      }
############ END CHANGED CODE
      
      if(keepBest) 
        { object@bestSol[[iter]] <- unique(Pop[Fitness == max(Fitness, na.rm = TRUE),,drop=FALSE]) }
      
      # apply a user's defined function to update the GA object
      if(is.function(postFitness))
        { 
          # object <- postFitness(object, ...) 
          object <- do.call(postFitness, c(object, callArgs)) 
          Fitness <- object@fitness
          Pop <- object@population
      }

      # check stopping criteria
      if(iter > 1)
        object@run <- GA::garun(fitnessSummary[seq(iter),1])
      if(object@run >= run) break  
      if(max(Fitness, na.rm = TRUE) >= maxFitness) break
      if(object@iter == maxiter) break  

      ord <- order(Fitness, decreasing = TRUE)
      PopSorted <- Pop[ord,,drop=FALSE]
      FitnessSorted <- Fitness[ord]
        
      # selection
      if(is.function(selection))
        { sel <- selection(object)
          Pop <- sel$population
          Fitness <- sel$fitness
        }
      else
        { sel <- sample(1:popSize, size = popSize, replace = TRUE)
          Pop <- object@population[sel,]
          Fitness <- object@fitness[sel]
        }
      object@population <- Pop
      object@fitness <- Fitness

      # crossover
      if(is.function(crossover) & pcrossover > 0)
        { nmating <- floor(popSize/2)
          mating <- matrix(sample(1:(2*nmating), size = (2*nmating)), ncol = 2)
          for(i in seq_len(nmating))
             { if(pcrossover > stats::runif(1))
                 { parents <- mating[i,]
                   Crossover <- crossover(object, parents)
                   Pop[parents,] <- Crossover$children
                   Fitness[parents] <- Crossover$fitness
                 }
             }             
          object@population <- Pop
          object@fitness <- Fitness
        }

      # mutation
      pm <- if(is.function(pmutation)) pmutation(object) else pmutation
      if(is.function(mutation) & pm > 0)
        { for(i in seq_len(popSize)) 
             { if(pm > stats::runif(1)) 
                 { Mutation <- mutation(object, i)
                   Pop[i,] <- Mutation
                   Fitness[i] <- NA
                 }
             }
          object@population <- Pop
          object@fitness <- Fitness
        }
      # elitism
      if(elitism > 0) # (elitism > 0 & iter > 1) 
        { ord <- order(object@fitness, na.last = TRUE)
          u <- which(!duplicated(PopSorted, margin = 1))
		  possibleElitism <- min(elitism, length(u))
          Pop[ord[1:possibleElitism],] <- PopSorted[u[1:possibleElitism],]
          Fitness[ord[1:possibleElitism]] <- FitnessSorted[u[1:possibleElitism]]
          object@population <- Pop
          object@fitness <- Fitness
      } 

      if(is.function(monitor)) 
############ CHANGED LINE
        { utils::flush.console() }
  }
  
  if(is.function(monitor)) 
    { utils::flush.console() }

  # in case of premature convergence remove NAs from summary 
  # fitness evalutations
  object@summary <- stats::na.exclude(object@summary)
  attr(object@summary, "na.action") <- NULL

  # get solution(s)
  object@fitnessValue <- max(object@fitness, na.rm = TRUE)
  valueAt <- which(object@fitness == object@fitnessValue)
  solution <- object@population[valueAt,,drop=FALSE]
  if(nrow(solution) > 1)
    { # find unique solutions to precision given by default tolerance
      eps <- GA::gaControl("eps")
      solution <- unique(round(solution/eps)*eps, margin = 1)
    }
  #colnames(solution) <- GA:::parNames(object)
  object@solution <- solution
  if(keepBest)
    object@bestSol <- object@bestSol[!sapply(object@bestSol, is.null)]  

object@allEvaluations <- evaluation.history
  
  if(!is.null(cls)) parallel::stopCluster(cls)
  # return an object of class 'ga'
  return(object)
}

# This is a mutator for circular integer variables, that only mutates for the contiguous class,
# and wraps if overflow.
make.mutator <- function(poissonmean=0.8, num.classes) {
	out <- function (object, parent, ...) {
		mutate <- parent <- as.vector(object@population[parent, ])
		n <- length(parent)
		nbitstomutate <- stats::rpois(1, poissonmean) + 1
		j <- sample(1:n, size = nbitstomutate)
		mutate[j] <- (mutate[j] + sample(c(-1,1), length(j), replace=TRUE)) %% (num.classes[j] + 1)
		return(mutate)
	}

	attr(out, "mean") <- poissonmean
	return(out)
}

setClass(Class = "ga.eicm", 
         representation(allEvaluations = "list"
	         , time.took="numeric"
	         , cached="numeric"
	         , informed="numeric"
	         , survived="numeric"
	         , evaluated="numeric"
	         , cache.size="numeric"
	         )
         , contains="ga", package = "eicm")

