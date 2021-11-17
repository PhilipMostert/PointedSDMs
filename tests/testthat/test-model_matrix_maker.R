testthat::test_that('Test model_matrix_maker makes a matrix of covariate effects for each species', {
  
  ##Use PointedSDMs
  library(PointedSDMs)
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  data("SolTin_covariates")
  covariates <- sp::SpatialPointsDataFrame(coords = SolTin_covariates[,c('X','Y')],
                                           data = SolTin_covariates[,c('Forest', 'NPP','Altitude')],
                                           proj4string = Projection) 
  data("SolTin_ebird")
  ebird <- sp::SpatialPoints(SolTin_ebird[,c("X","Y")], proj4string = Projection)
  data('SolTin_gbif')
  gbif <- sp::SpatialPoints(SolTin_gbif[, c('X','Y')], proj4string = Projection)
  
  ##Make species
  ebird$species <- sample(c('species_a','species_b'), size = nrow(ebird@coords), replace = TRUE)
  gbif$species <- sample(c('species_c'), size = nrow(gbif@coords), replace = TRUE)
  
  all_species <- c('species_a','species_b','species_c')
  #Make model_matrix
   #Input is a named list
  species_matrix <- model_matrix_maker(datasets = list(ebird = ebird, gbif = gbif), species = 'species',
                                       allspecies = all_species,
                                       coords = c('X','Y'), covariates = covariates,
                                       componentstokeep = c('species'),
                                       proj = Projection)
  
  expect_setequal(names(species_matrix$ebird@data), c("species", "Forest", "NPP" ,"Altitude" , "species_a_intercept" ,"species_b_intercept",
                                                       "species_a_Forest", "species_b_Forest", "species_a_NPP", "species_b_NPP", "species_a_Altitude", "species_b_Altitude"))
  ##Need to fix this:: probably remove all
   #the if only one species then print covariate name
   #Issue regarding is many species in Data 1 but only 1 in 2.
  expect_setequal(names(species_matrix$gbif@data), c("species", "species_c_Forest", "species_c_NPP", "species_c_Altitude" ,
                                                      "species_c_intercept", "Forest", "NPP", "Altitude"))
  
  expect_false(any(sapply(species_matrix$ebirddata, is.na)))
  
  
  })
