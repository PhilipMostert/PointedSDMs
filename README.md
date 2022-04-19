
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inlabruSDMs

<!-- badges: start -->

[![R-CMD-check](https://github.com/PhilipMostert/cSDMs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PhilipMostert/cSDMs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of inlabruSDMs is to simplify the construction of integrated
species distribution models (ISDMs) for large collections of
heterogeneous data. It does so by building wrapper functions around
[inlabru](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13168),
which uses the [INLA
methodology](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2008.00700.x)
to estimate a class of latent Gaussian models.

## Installation

You can install the development version of inlabruSDMs from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("PhilipMostert/inlabruSDMs")
```

## Example

This is a basic example which shows you how to specify and run an
integrated model, using three disparate datasets containing locations of
the solitary tinamou.

``` r
library(inlabruSDMs)
#> Loading required package: INLA
#> Loading required package: Matrix
#> Loading required package: foreach
#> Loading required package: parallel
#> Loading required package: sp
#> This is INLA_22.03.16 built 2022-03-16 13:24:07 UTC.
#>  - See www.r-inla.org/contact-us for how to get help.
#> Loading required package: inlabru
library(raster)
## basic example code
```

``` r
#Load data in
data("SolitaryTinamou")

projection <- CRS("+proj=longlat +ellps=WGS84")
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded datum Unknown based on WGS84 ellipsoid in Proj4
#> definition

species <- SolitaryTinamou$datasets

covariates <- SolitaryTinamou$covariates
covariates <- scale(stack(covariates))
crs(covariates) <- projection

mesh <- SolitaryTinamou$mesh
mesh$crs <- projection
```

``` r
#Specify model -- here we ryb a model with three spatial covariates and a shared spatial field
model <- intModel(species, spatialCovariates = covariates, Coordinates = c('X', 'Y'),
                 Projection = projection, Mesh = mesh, pointsSpatial = TRUE,
                 pointsIntercept = FALSE)
```

``` r
#Run the integrated model
modelRun <- runModel(model, options = list(control.inla = list(int.strategy = 'eb')))
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'Forest' for 'Forest' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'Forest' for 'Forest' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'Forest' for 'Forest' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'NPP' for 'NPP' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'NPP' for 'NPP' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'NPP' for 'NPP' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'Altitude' for 'Altitude' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'Altitude' for 'Altitude' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(subcomp$input, data = lh$data, env = env, : Model input 'Altitude' for 'Altitude' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat +R=1
#> +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat
#> +R=6378137 +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat
#> +R=6378137 +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=geocent +R=1
#> +units=m +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'Forest' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'NPP' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'Altitude' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat +R=1
#> +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat
#> +R=6378137 +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat
#> +R=6378137 +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=geocent +R=1
#> +units=m +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'Forest' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'NPP' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'Altitude' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat +R=1
#> +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat
#> +R=6378137 +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=longlat
#> +R=6378137 +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj
#> = prefer_proj): Discarded ellps unknown in Proj4 definition: +proj=geocent +R=1
#> +units=m +no_defs +type=crs
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum unknown in Proj4 definition
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'Forest' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'NPP' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
#> Warning in input_eval.bru_input(component[[part]]$input, data, env = component$env, : Model input 'Altitude' for 'main' returned some NA values.
#> Attempting to fill in spatially by nearest available value.
#> To avoid this basic covariate imputation, supply complete data.
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 1 is not projected; GEOS expects planar coordinates
#> Warning in RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance"): Spatial
#> object 2 is not projected; GEOS expects planar coordinates
summary(modelRun)
#> Summary of 'bruSDM' object:
#> 
#> inlabru version: 2.5.2
#> INLA version: 22.03.16
#> 
#> Types of data modelled:
#>                                     
#> eBird                   Present only
#> Parks                   Present only
#> Gbif                    Present only
#> Time used:
#>     Pre = 3.23, Running = 18.3, Post = 0.0341, Total = 21.6 
#> Fixed effects:
#>           mean    sd 0.025quant 0.5quant 0.975quant  mode kld
#> Forest   0.130 0.030      0.072    0.130      0.188 0.130   0
#> NPP      0.098 0.029      0.041    0.098      0.154 0.098   0
#> Altitude 0.113 0.030      0.054    0.113      0.173 0.113   0
#> 
#> Random effects:
#>   Name     Model
#>     shared_spatial SPDE2 model
#> 
#> Model hyperparameters:
#>                            mean   sd 0.025quant 0.5quant 0.975quant  mode
#> Theta1 for shared_spatial -2.72 0.00      -2.72    -2.72      -2.72 -2.72
#> Theta2 for shared_spatial -3.96 0.00      -3.97    -3.96      -3.96 -3.96
#> 
#> Deviance Information Criterion (DIC) ...............: 4666.47
#> Deviance Information Criterion (DIC, saturated) ....: -36354.75
#> Effective number of parameters .....................: -75.42
#> 
#> Watanabe-Akaike information criterion (WAIC) ...: 5299.54
#> Effective number of parameters .................: 546.02
#> 
#> Marginal log-Likelihood:  -4221.12 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```
