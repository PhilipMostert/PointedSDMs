testthat::test_that('Test that leave_one out works properly.', {
  
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
                              tolerance = 0.2,
                              options = list(control.inla = list(int.strategy = 'eb')))

  ##Test leaving out ebird dataset
  ebird_out <- leave_one_out(model = integrated_model, dataset = 'ebird',
                             predictions = TRUE)
  
  ##
  expect_false('ebird'%in%names(ebird_out$Leaving_out_ebird$dataset_names))
  expect_false('ebird'%in%names(ebird_out$Leaving_out_ebird$bru_info$lhoods))
  expect_length(ebird_out$Leaving_out_ebird$.args$control.family, 1)
  
  ##Run model with species
  integrated_model_2 <- bru_sdm(data_to_use, 
                        spatialcovariates = SolTin_covariates, 
                        tolerance = 0.2, specieseffects = TRUE,
                        options = list(control.inla = list(int.strategy = 'eb')))
  
  ebird_out_2 <- leave_one_out(model = integrated_model, dataset = 'ebird')
  
  ##Should leave out all the components related to ebird.
  expect_false('ebird_species_intercept'%in%row.names(ebird_out_2$Leaving_out_ebird$summary.fixed))
  expect_false('ebird_species_Altitude'%in%row.names(ebird_out_2$Leaving_out_ebird$summary.fixed))
  expect_false('ebird_species_NPP'%in%row.names(ebird_out_2$Leaving_out_ebird$summary.fixed))
  expect_false('ebird_species_Forest'%in%row.names(ebird_out_2$Leaving_out_ebird$summary.fixed))
  
  
    
})