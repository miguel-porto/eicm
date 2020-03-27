#' eicm: Explicit Interaction Community Models
#'
#' Explicit Interaction Community Models 
#' 
#' Model fitting and species biotic interaction network topology selection for explicit interaction community models.
#' Explicit interaction community models are an extension of binomial linear models for joint modelling of species
#' communities, that incorporate both the effects of species biotic interactions and the effects of missing covariates.
#' Species interactions are modelled as direct effects of each species on each of the others, and are estimated
#' alongside the effects of missing covariates, modelled as latent factors.
#' The package includes a penalized maximum likelihood fitting function, and a genetic algorithm for selecting
#' the most parsimonious species interaction network topology.
#' 
#' @section Main functions:
#' \code{\link{eicm}}, \code{\link{eicm.fit}}
#'
#' @docType package
#' @name eicm-package
NULL
