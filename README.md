# eicm (R package) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

**Explicit Interaction Community Models**

Available on [CRAN](https://cran.r-project.org/package=eicm), can be directly installed from within R.

Model fitting and species biotic interaction network topology selection for EICMs.
EICMs are an extension of binomial GLMs for joint modelling of species communities, that incorporate
both the effects of species biotic interactions and the effects of missing covariates. Species
interactions are modelled as direct effects of each species on each of the others, and are estimated
alongside the effects of missing covariates, modelled as latent factors. The package includes a
Penalized Maximum Likelihood fitting function, and a Genetic Algorithm for selecting the most
parsimonious species interaction network topology.

