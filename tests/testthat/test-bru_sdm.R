testthat::test_that('Test that bru_sdm creates desired outputs based on different arguments', {
  
  ##Use data from PointedSDMs
  library(PointedSDMs)
  library(spatstat)
  library(sp)
  
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  Meshpars <- list(cutoff=5, max.edge=c(5, 10), offset=c(1,1))
  
  ebird <- PointedSDMs::SolTin_ebird
  parks <- PointedSDMs::SolTin_parks
  
  #Add species
  ebird$species <- 'ebird_species'
  parks$species <- 'parks_species'
  
  region.mask=as.owin(cbind(SolTin_covariates[,c("X","Y")], In=rep(TRUE,nrow(SolTin_covariates))),
                      step=c(0.25, 0.3))
  region.mask$m[is.na(region.mask$m)] <- FALSE
  Region.poly <- simplify.owin(as.polygonal(region.mask), dmin=0.5)
  PolyPoints <- cbind(Region.poly$bdry[[1]]$x[c(1,length(Region.poly$bdry[[1]]$x):1)],
                      Region.poly$bdry[[1]]$y[c(1,length(Region.poly$bdry[[1]]$y):1)])
  Pgon <- Polygons(list(region=Polygon(coords=PolyPoints)), ID="region")
  region.polygon=SpatialPolygons(list(Pgon), proj4string = Projection)
  
  data_to_use <- organize_data(ebird, parks, poresp = 'poresp', paresp = 'Present',
                               coords = c('X','Y'), proj = Projection,
                               boundary = region.polygon, meshpars = Meshpars,
                               speciesname = 'species')
  
  integrated_model <- bru_sdm(data_to_use, 
                      spatialcovariates = SolTin_covariates, 
                      tolerance = 0.2)
  
  #Test the class is correct
  expect_equal(class(integrated_model), c("bru_sdm", "bru", "iinla", "inla"))
  
  #Test that each dataset has own spde
  expect_true(all(c('ebird_spde','parks_spde')%in%names(integrated_model$summary.random)))
  
  #Test that link function is picked up
  expect_setequal(unlist(integrated_model$bru_sdm_options$control.family), c('default', 'cloglog'))
  
  #Test sources of information correct
  expect_setequal(unlist(integrated_model$sources_of_information), c('ebird','parks'))
  
  #rerun model with shared spatial
  integrated_model_shared <-  bru_sdm(data_to_use, 
                              spatialcovariates = SolTin_covariates,
                              sharedspatial = TRUE,
                              tolerance = 0.2)
  
  #Test that each dataset has own spde
  expect_true('shared_spatial'%in%names(integrated_model_shared$summary.random))
  
  #Run different covs per dataset
  integrated_model_diff_covs <- bru_sdm(data_to_use, 
                                spatialcovariates = SolTin_covariates,
                                covariatesbydataset = c(ebird = 'Forest', parks = 'NPP'),
                                sharedspatial = TRUE,
                                tolerance = 0.2)
  
  expect_true('Forest'%in%integrated_model_diff_covs$bru_info$lhoods$ebird$include_components)
  expect_false('NPP'%in%integrated_model_diff_covs$bru_info$lhoods$ebird$include_components)
  expect_false('Forest'%in%integrated_model_diff_covs$bru_info$lhoods$parks$include_components)
  expect_true('NPP'%in%integrated_model_diff_covs$bru_info$lhoods$parks$include_components)
  
  #Test species effects are included
  
  integrated_model_species <- bru_sdm(data_to_use, 
                              spatialcovariates = SolTin_covariates,
                              sharedspatial = TRUE,
                              specieseffects = TRUE,
                              tolerance = 0.2)
  
  expect_true('species_spde'%in%names(integrated_model_species$summary.random))
  expect_true(all(c('ebird_species_Forest','parks_species_Forest')%in%row.names(integrated_model_species$summary.fixed)))
  expect_true(all(c('ebird_species_NPP','parks_species_NPP')%in%row.names(integrated_model_species$summary.fixed)))
  expect_true(all(c('ebird_species_Altitude','parks_species_Altitude')%in%row.names(integrated_model_species$summary.fixed)))
  
  
  })
