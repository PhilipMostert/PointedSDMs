
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PointedSDMs

<!-- badges: start -->

[![R-CMD-check](https://github.com/PhilipMostert/PointedSDMs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PhilipMostert/PointedSDMs/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/PhilipMostert/PointedSDMs/branch/ChangingToR6/graph/badge.svg)](https://app.codecov.io/gh/PhilipMostert/PointedSDMs?branch=ChangingToR6)

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

| Function name  | Function description                                                                                          |
|----------------|---------------------------------------------------------------------------------------------------------------|
| `intModel()`   | Initialize and specify the components used in the integrated model.                                           |
| `blockedCV()`  | Perform spatial blocked cross-validation.                                                                     |
| `runModel()`   | Estimate and preform inference on the integrated model.                                                       |
| `datasetOut()` | Perform dataset-out cross-validation, which calculates the impact individual datasets have on the full model. |

The function `intModel()` produces an [R6](https://github.com/r-lib/R6)
object, and as a result there are various *slot functions* available to
further specify the components of the model. These *slot functions*
include:

| `intModel()` slot function   | Function description                                                                                                                                            |
|------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `` `.$plot()` ``             | Used to create a plot of the available data. The output of this function is an object of class [`gg`](https://github.com/tidyverse/ggplot2).                    |
| `` `.$addBias()` ``          | Add an additional spatial field to a dataset to account for sampling bias in unstructured datasets.                                                             |
| `` `.$updateFormula()` ``    | Used to update a formula for a process. The idea is to start specify the full model with `intModel()`, and then thin components per dataset with this function. |
| `` `.$updateComponents()` `` | Change or add new components used by [inlabru](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13168) in the integrated model.            |
| `` `.$priorsFixed()` ``      | Change the specification of the prior distribution for the fixed effects in the model.                                                                          |
| `` `.$specifySpatial()` ``   | Specify the spatial field in the model using penalizing complexity (PC) priors.                                                                                 |
| `` `.$spatialBlock()` ``     | Used to specify how the points are spatially blocked. Spatial cross-validation is subsequently performed using `blockedCV()`.                                   |
| `` `.$addSamplers()` ``      | Function to add an integration domain for the PO datasets.                                                                                                      |

## Example

This is a basic example which shows you how to specify and run an
integrated model, using three disparate datasets containing locations of
the solitary tinamou (*Tinamus solitarius)*.

``` r
library(PointedSDMs)
library(ggplot2)
library(raster)
```

``` r
#Load data in

data("SolitaryTinamou")

projection <- CRS("+proj=longlat +ellps=WGS84")

species <- SolitaryTinamou$datasets

Forest <- SolitaryTinamou$covariates$Forest

crs(Forest) <- projection

mesh <- SolitaryTinamou$mesh
mesh$crs <- projection
```

Setting up the model is done easily with `intModel()`, where we specify
the required components of the model:

``` r
#Specify model -- here we run a model with one spatial covariate and a shared spatial field

model <- intModel(species, spatialCovariates = Forest, Coordinates = c('X', 'Y'),
                 Projection = projection, Mesh = mesh, responsePA = 'Present')
```

We can also make a quick plot of where the species are located using
`` `.$plot()` ``:

``` r
region <- SolitaryTinamou$region

model$plot(Boundary = FALSE) + gg(region)
```

<img src="man/figures/README-plot-1.png" width="100%" />

We can estimate the parameters in the model using the `runModel()`
function:

``` r
#Run the integrated model

modelRun <- runModel(model, options = list(control.inla = list(int.strategy = 'eb')))
summary(modelRun)
#> Summary of 'bruSDM' object:
#> 
#> inlabru version: 2.5.2
#> INLA version: 22.07.01
#> 
#> Types of data modelled:
#>                                     
#> eBird                   Present only
#> Parks                Present absence
#> Gbif                    Present only
#> Time used:
#>     Pre = 1.43, Running = 16.5, Post = 0.02, Total = 17.9 
#> Fixed effects:
#>                   mean    sd 0.025quant 0.5quant 0.975quant mode   kld
#> Forest          -0.002 0.001     -0.005   -0.002      0.000   NA 0.099
#> eBird_intercept -0.244 0.047     -0.336   -0.244     -0.152   NA 0.505
#> Parks_intercept -0.536 0.178     -0.892   -0.535     -0.191   NA 0.000
#> Gbif_intercept  -0.553 0.048     -0.647   -0.553     -0.459   NA 0.285
#> 
#> Random effects:
#>   Name     Model
#>     shared_spatial SPDE2 model
#> 
#> Model hyperparameters:
#>                            mean    sd 0.025quant 0.5quant 0.975quant  mode
#> Theta1 for shared_spatial -3.39 0.004      -3.39    -3.39      -3.38 -3.39
#> Theta2 for shared_spatial -8.00 0.004      -8.01    -8.00      -7.99 -8.00
#> 
#> Deviance Information Criterion (DIC) ...............: 4204.54
#> Deviance Information Criterion (DIC, saturated) ....: -23256.48
#> Effective number of parameters .....................: 229.82
#> 
#> Watanabe-Akaike information criterion (WAIC) ...: 3088.50
#> Effective number of parameters .................: 653.03
#> 
#> Marginal log-Likelihood:  -3891.83 
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

plot(predictions)
```

<img src="man/figures/README-predict_and_plot-1.png" width="100%" />
