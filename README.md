
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PointedSDMs

<!-- badges: start -->

[![R-CMD-check](https://github.com/PhilipMostert/PointedSDMs/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/PhilipMostert/PointedSDMs/actions/workflows/R-CMD-check.yaml)[![Codecov
test
coverage](https://codecov.io/gh/PhilipMostert/PointedSDMs/branch/main/graph/badge.svg)](https://app.codecov.io/gh/PhilipMostert/PointedSDMs?branch=ChangingToR6)
[![DOI](https://zenodo.org/badge/368823136.svg)](https://zenodo.org/badge/latestdoi/368823136)

<!-- badges: end -->

The goal of *PointedSDMs* is to simplify the construction of integrated
species distribution models (ISDMs) for large collections of
heterogeneous data. It does so by building wrapper functions around
[inlabru](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13168),
which uses the [INLA
methodology](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2008.00700.x)
to estimate a class of latent Gaussian models.

## Installation

You can install the development version of PointedSDMs from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("PhilipMostert/PointedSDMs")
```

or directly through CRAN using:

``` r
install.packages('PointedSDMs')
```

## Package functionality

*PointedSDMs* includes a selection of functions used to streamline the
construction of ISDMs as well and perform model cross-validation. The
core functions of the package are:

| Function name    | Function description                                                                                          |
|------------------|---------------------------------------------------------------------------------------------------------------|
| `startISDM()`    | Initialize and specify the components used in the integrated model.                                           |
| `startSpecies()` | Initialize and specify the components used in the multi-species integrated model.                             |
| `blockedCV()`    | Perform spatial blocked cross-validation.                                                                     |
| `fitISDM()`      | Estimate and preform inference on the integrated model.                                                       |
| `datasetOut()`   | Perform dataset-out cross-validation, which calculates the impact individual datasets have on the full model. |

The function `intModel()` produces an [R6](https://github.com/r-lib/R6)
object, and as a result there are various *slot functions* available to
further specify the components of the model. These *slot functions*
include:

| `intModel()` slot function   | Function description                                                                                                                                             |
|------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `` `.$help()` ``             | Show documentation for each of the slot functions.                                                                                                               |
| `` `.$plot()` ``             | Used to create a plot of the available data. The output of this function is an object of class [`gg`](https://github.com/tidyverse/ggplot2).                     |
| `` `.$addBias()` ``          | Add an additional spatial field to a dataset to account for sampling bias in unstructured datasets.                                                              |
| `` `.$updateFormula()` ``    | Used to update a formula for a process. The idea is to start specify the full model with `startISDM()`, and then thin components per dataset with this function. |
| `` `.$updateComponents()` `` | Change or add new components used by [inlabru](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13168) in the integrated model.             |
| `` `.$priorsFixed()` ``      | Change the specification of the prior distribution for the fixed effects in the model.                                                                           |
| `` `.$specifySpatial()` ``   | Specify the spatial field in the model using penalizing complexity (PC) priors.                                                                                  |
| `` `.$spatialBlock()` ``     | Used to specify how the points are spatially blocked. Spatial cross-validation is subsequently performed using `blockedCV()`.                                    |
| `` `.$addSamplers()` ``      | Function to add an integration domain for the PO datasets.                                                                                                       |
| `` `.$specifyRandom()` ``    | Specify the priors for the random effects in the model.                                                                                                          |
| `` `.$changeLink()` ``       | Change the link function of a process.                                                                                                                           |

## Example

This is a basic example which shows you how to specify and run an
integrated model, using three disparate datasets containing locations of
the solitary tinamou (*Tinamus solitarius)*.

``` r

library(PointedSDMs)
library(ggplot2)
library(terra)
```

``` r

bru_options_set(inla.mode = "experimental")

#Load data in

data("SolitaryTinamou")

projection <- "+proj=longlat +ellps=WGS84"

species <- SolitaryTinamou$datasets

covariates <- terra::rast(system.file('extdata/SolitaryTinamouCovariates.tif', 
                                      package = "PointedSDMs"))

mesh <- SolitaryTinamou$mesh
```

Setting up the model is done easily with `startISDM()`, where we specify
the required components of the model:

``` r

#Specify model -- here we run a model with one spatial covariate and a shared spatial field

model <- startISDM(species, spatialCovariates = covariates,
                 Projection = projection, Mesh = mesh, responsePA = 'Present')
```

We can also make a quick plot of where the species are located using
`` `.$plot()` ``:

``` r

region <- SolitaryTinamou$region

model$plot(Boundary = FALSE) + 
  geom_sf(data = st_boundary(region))
```

<img src="man/figures/README-plot-1.png" width="100%" />

To improve stability, we specify priors for the intercepts of the model
using `` `.$priorsFixed()` ``

``` r

model$priorsFixed(Effect = 'Intercept',
                  mean.linear = 0, 
                  prec.linear = 0.11)
```

And *PC* priors for the spatial field using `` `.$specifySpatial()` ``:

``` r

model$specifySpatial(sharedSpatial = TRUE,
                     prior.range = c(50, 0.1),
                     prior.sigma = c(0.1, 0.1))
```

We can then estimate the parameters in the model using the `fitISDM()`
function:

``` r

modelRun <- fitISDM(model, options = list(control.inla = list(int.strategy = 'eb',
                                                              diagonal = 0.1), 
                                          safe = TRUE))
summary(modelRun)
#> inlabru version: 2.10.1
#> INLA version: 24.03.20
#> Components:
#> eBird_spatial: main = spde(geometry), group = exchangeable(1L), replicate = iid(1L)
#> Parks_spatial(=eBird_spatial): main = unknown(geometry), group = exchangeable(1L), replicate = iid(1L)
#> Gbif_spatial(=eBird_spatial): main = unknown(geometry), group = exchangeable(1L), replicate = iid(1L)
#> Forest: main = linear(Forest), group = exchangeable(1L), replicate = iid(1L)
#> NPP: main = linear(NPP), group = exchangeable(1L), replicate = iid(1L)
#> Altitude: main = linear(Altitude), group = exchangeable(1L), replicate = iid(1L)
#> eBird_intercept: main = linear(1), group = exchangeable(1L), replicate = iid(1L)
#> Parks_intercept: main = linear(1), group = exchangeable(1L), replicate = iid(1L)
#> Gbif_intercept: main = linear(1), group = exchangeable(1L), replicate = iid(1L)
#> Likelihoods:
#>   Family: 'cp'
#>     Data class: 'sf', 'data.frame'
#>     Predictor: geometry ~ .
#>   Family: 'binomial'
#>     Data class: 'sf', 'data.frame'
#>     Predictor: Present ~ .
#>   Family: 'cp'
#>     Data class: 'sf', 'data.frame'
#>     Predictor: geometry ~ .
#> Time used:
#>     Pre = 2.41, Running = 17.7, Post = 0.319, Total = 20.5 
#> Fixed effects:
#>                   mean    sd 0.025quant 0.5quant 0.975quant   mode kld
#> Forest           0.027 0.004      0.019    0.027      0.034  0.027   0
#> NPP              0.000 0.000      0.000    0.000      0.000  0.000   0
#> Altitude        -0.001 0.000     -0.001   -0.001      0.000 -0.001   0
#> eBird_intercept -6.857 0.283     -7.412   -6.857     -6.301 -6.857   0
#> Parks_intercept -6.089 0.401     -6.875   -6.089     -5.302 -6.089   0
#> Gbif_intercept  -7.415 0.317     -8.037   -7.415     -6.793 -7.415   0
#> 
#> Random effects:
#>   Name     Model
#>     eBird_spatial SPDE2 model
#>    Parks_spatial Copy
#>    Gbif_spatial Copy
#> 
#> Model hyperparameters:
#>                          mean    sd 0.025quant 0.5quant 0.975quant  mode
#> Range for eBird_spatial 2.915 0.278      2.417    2.898      3.511 2.858
#> Stdev for eBird_spatial 0.923 0.058      0.814    0.921      1.042 0.918
#> Beta for Parks_spatial  0.020 0.062     -0.103    0.021      0.141 0.023
#> Beta for Gbif_spatial   0.441 0.054      0.335    0.440      0.549 0.438
#> 
#> Deviance Information Criterion (DIC) ...............: -650.55
#> Deviance Information Criterion (DIC, saturated) ....: NA
#> Effective number of parameters .....................: -1472.94
#> 
#> Watanabe-Akaike information criterion (WAIC) ...: 1580.61
#> Effective number of parameters .................: 331.43
#> 
#> Marginal log-Likelihood:  -1607.03 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```

*PointedSDMs* also includes generic predict and plot functions:

``` r

predictions <- predict(modelRun, mesh = mesh,
                       mask = region, 
                       spatial = TRUE,
                       fun = 'linear')

plot(predictions, variable = c('mean', 'sd'))
```

<img src="man/figures/README-predict_and_plot-1.png" width="100%" />
