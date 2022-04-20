
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PointedSDMs

<!-- badges: start -->
[![R-CMD-check](https://github.com/PhilipMostert/cSDMs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PhilipMostert/cSDMs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of PointedSDMs is to simplify the construction of integrated
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

## Example

This is a basic example which shows you how to specify and run an
integrated model, using three disparate datasets containing locations of
the solitary tinamou.

``` r
library(PointedSDMs)
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

``` r
#Specify model -- here we run a model with one spatial covariate and a shared spatial field

model <- intModel(species, spatialCovariates = Forest, Coordinates = c('X', 'Y'),
                 Projection = projection, Mesh = mesh, responsePA = 'Present')
```

``` r
#Run the integrated model

modelRun <- runModel(model, options = list(control.inla = list(int.strategy = 'eb')))
summary(modelRun)
#> Summary of 'bruSDM' object:
#> 
#> inlabru version: 2.5.2
#> INLA version: 22.03.16
#> 
#> Types of data modelled:
#>                                     
#> eBird                   Present only
#> Parks                Present absence
#> Gbif                    Present only
#> Time used:
#>     Pre = 3.11, Running = 7.37, Post = 0.0308, Total = 10.5 
#> Fixed effects:
#>                   mean     sd 0.025quant 0.5quant 0.975quant   mode kld
#> Forest           1.101  0.029      1.044    1.101      1.160  1.101   0
#> eBird_intercept  4.186 18.184    -31.516    4.185     39.858  4.186   0
#> Parks_intercept -8.790 18.187    -44.497   -8.790     26.888 -8.790   0
#> Gbif_intercept   2.570 18.185    -33.132    2.570     38.243  2.570   0
#> 
#> Random effects:
#>   Name     Model
#>     shared_spatial SPDE2 model
#> 
#> Model hyperparameters:
#>                            mean   sd 0.025quant 0.5quant 0.975quant  mode
#> Theta1 for shared_spatial -2.08 0.00      -2.08    -2.08      -2.05 -2.08
#> Theta2 for shared_spatial -3.24 0.00      -3.24    -3.24      -3.18 -3.24
#> 
#> Deviance Information Criterion (DIC) ...............: -2835.44
#> Deviance Information Criterion (DIC, saturated) ....: -30379.28
#> Effective number of parameters .....................: -986.73
#> 
#> Watanabe-Akaike information criterion (WAIC) ...: -Inf
#> Effective number of parameters .................: 4.46e+74
#> 
#> Marginal log-Likelihood:  -571.83 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```
