## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(timeit = local({
  now = NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res = difftime(Sys.time(), now, units="secs")
      now <<- NULL
      # use options$label if you want the chunk label as well
      if(res > 60)
	      sprintf("<p style=\"text-align:right\">-> running time of the above code chunk: %.1f minutes</p>", res / 60)
      else
	      sprintf("<p style=\"text-align:right\">-> running time of the above chunk: %.1f seconds</p>", res)
    }
  }})
)

library(eicm)
set.seed(2)		# for vignette reproducibility

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # fit & select an EICM with 2 latent variables
#  m <- eicm(occurrences, n.latent=2)
#  
#  # display estimated coefficients (note they are organized in matrices)
#  coef(m$selected.model)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  m <- eicm(occurrences, n.latent=2, regularization=c(6, 0.5), penalty=4, theta.threshold=0.5)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # excluding interactions with the formula syntax
#  m <- eicm(occurrences, forbidden=list(
#      sp3 ~ sp4 + sp5,		# sp3 must not be affected by sp4 nor sp5
#      sp4 ~ .,				# sp4 must not be affected by any other
#      sp1 ~ . - sp8			# sp1 must not be affected by any other except sp8
#      ))
#  
#  # display species interaction coefficients
#  # note the zeroed coefficients are those that were excluded
#  coef(m$fitted.model)$sp

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # Excluding interactions with the matrix syntax
#  
#  # create a square matrix species x species, all zeroes
#  mask <- matrix(0, nrow=ncol(occurrences), ncol=ncol(occurrences))
#  
#  # set to 1 those interactions we want to include
#  mask[4, 2] <- 1		# species #2 may affect species #4
#  mask[6, 1] <- 1		# species #1 may affect species #6
#  mask[, 7] <- 1		# species #7 may affect any other
#  mask[2, ] <- 1		# species #2 may be affected by any other
#  
#  m <- eicm(occurrences, mask.sp=mask)
#  
#  # display species interaction coefficients
#  # note the zeroed coefficients are those that were excluded
#  coef(m$fitted.model)$sp

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # exclude interactions departing from species with 20 or less occurrences
#  m <- eicm(occurrences, mask.sp=mask, exclude.prevalence=20)

