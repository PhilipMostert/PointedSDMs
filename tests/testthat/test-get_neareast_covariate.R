

testthat::test_that('Test that get_nearest_covariate returns non-null covariates at each data point', {
  
  ##Use data from PointedSDMs
 library(PointedSDMs)
  
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  data("SolTin_covariates")
  covariates <- sp::SpatialPointsDataFrame(coords = SolTin_covariates[,c('X','Y')],
                                           data = SolTin_covariates[,c('Forest', 'NPP','Altitude')],
                                           proj4string = Projection) 
  data("SolTin_ebird")
  ebird <- sp::SpatialPoints(SolTin_ebird[,c("X","Y")], proj4string = Projection)
  ##Add POresp to show components are kept
  ebird$POresp <- rep(1,nrow(ebird@coords))
  ##Choose to keep only Forest and NPP
  ##
  nearest_cov <- get_nearest_covariate(points = ebird,
                                       spatialcovariates = covariates,
                                       covariatestokeep = c('Forest','NPP'),
                                       coords = c('X','Y'),
                                       proj = Projection,
                                       componentstokeep = c('POresp','PAresp','weights'))
  
  expect_equal(class(nearest_cov)[1], 'SpatialPointsDataFrame')
  expect_equal(names(nearest_cov), c('POresp','Forest','NPP'))
  expect_equal(any(sapply(nearest_cov@data, is.na)), FALSE)
  
  
})
