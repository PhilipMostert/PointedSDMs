# inlabruSDMs

<!-- badges: start -->
[![R-CMD-check](https://github.com/PhilipMostert/inlabruSDMs/workflows/R-CMD-check/badge.svg)](https://github.com/PhilipMostert/inlabruSDMs/actions)
[![Travis build status](https://travis-ci.com/PhilipMostert/inlabruSDMs.svg?branch=ChangingToR6)](https://travis-ci.com/PhilipMostert/inlabruSDMs)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/PhilipMostert/inlabruSDMs?branch=ChangingToR6&svg=true)](https://ci.appveyor.com/project/PhilipMostert/inlabruSDMs)
<!-- badges: end -->

inlabruSDMs is an R package used to construct species distribution models from a variety of different data sources. It does so by building wrapper functions around [inlabru](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13168), which uses the [INLA methodology](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2008.00700.x) to estimate a class of latent Gaussian models.

## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("PhilipMostert/inlabruSDMs")
