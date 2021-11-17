testthat::test_that('Test that formula_maker makes formulas given differnt arguments', {
  
  #Arbitrary projection
  projection <- CRS('+proj=tmerc')
  
  #Make random shape to generate points on
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  
  Poly = Polygon(xy)
  Poly = Polygons(list(Poly),1)
  SpatialPoly = SpatialPolygons(list(Poly), proj4string = projection)
  
  #Make random points
  #Random presence only dataset
  PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  #Add species name
  PO$species <- sample(x = c('a','b'), size = nrow(PO@coords), replace = TRUE)
  
  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  #Add species name
  PA$species <- sample(x = c('c','d'), size = nrow(PA@coords), replace = TRUE)
  #Add response name
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  
  response_variables <- c('POresp', 'PAresp')
  dataset_names <- c('PO', 'PA')
  
  #Arbitrary covariate names
  covariates <- c('Elevation', 'Temperature', 'Precipitation')
  
  #Make formula with: spatial fields per dataset;
  #                   intercepts per dataset.
  
  formula_1 <- formula_maker(response = response_variables, dataset = dataset_names,
               covariates = covariates, pointsspatial = TRUE,
               pointsintercept = TRUE, sharedspatial = FALSE,
               spatialdatasets = NULL, covariatesbydataset = NULL,
               species = NULL, pointcovariates = NULL)
  
  expect_equal(as.character(formula_1[[1]])[2], 'POresp')
  expect_equal(as.character(formula_1[[2]])[2], 'PAresp')
  expect_true(grepl('PO_intercept', as.character(formula_1[[1]])[3]))
  expect_true(grepl('PO_spde', as.character(formula_1[[1]])[3]))
  expect_true(grepl('PA_intercept', as.character(formula_1[[2]])[3]))
  expect_true(grepl('PA_spde', as.character(formula_1[[2]])[3]))
  
  #Make formula with: shared spatial field;
  #                   no intercepts.
  
  formula_2 <- formula_maker(response = response_variables, dataset = dataset_names,
                             covariates = covariates, pointsspatial = TRUE,
                             pointsintercept = FALSE, sharedspatial = TRUE,
                             spatialdatasets = NULL, covariatesbydataset = NULL,
                             species = NULL, pointcovariates = NULL)
  
  expect_false(grepl('PO_intercept', as.character(formula_2[[1]])[3]))
  expect_false(grepl('PA_intercept', as.character(formula_2[[2]])[3]))
  expect_true(grepl('shared_spatial', as.character(formula_2[[1]])[3]))
  expect_true(grepl('shared_spatial', as.character(formula_2[[2]])[3]))
  
  #Make formula with: PO spatial, PA no spatial;
  #                   PO has covariate Elevation;
  #                   PA has covariate Temperature.
  
  formula_3 <- formula_maker(response = response_variables, dataset = dataset_names,
                             covariates = covariates, pointsspatial = TRUE,
                             pointsintercept = TRUE, sharedspatial = FALSE,
                             spatialdatasets = c('PO'), covariatesbydataset = c(PO = 'Elevation', PA = 'Temperature'),
                             species = NULL, pointcovariates = NULL)
  
  expect_true(grepl('Elevation', as.character(formula_3[[1]])[3]))
  expect_false(grepl('Temperature', as.character(formula_3[[1]])[3]))
  expect_true(grepl('PO_spde', as.character(formula_3[[1]])[3]))
  expect_false(grepl('Elevation', as.character(formula_3[[2]])[3]))
  expect_true(grepl('Temperature', as.character(formula_3[[2]])[3]))
  expect_false(grepl('PA_spde', as.character(formula_3[[2]])[3]))
  
  #Make formula with: Species effects;
  #                   Species intercepts;
  #                   No other spatial effects.
  
  speciesdataset <- list(PO = PO$species, PA = PA$species)
  
  formula_4 <- formula_maker(response = response_variables, dataset = dataset_names,
                             covariates = covariates, pointsspatial = FALSE,
                             pointsintercept = TRUE, sharedspatial = FALSE,
                             spatialdatasets = NULL, covariatesbydataset = NULL,
                             species = 'species',
                             speciesdataset = speciesdataset,
                             pointcovariates = NULL)
  
  expect_true(grepl('species_spde', as.character(formula_4[[1]])[3]))
  expect_true(grepl('species_spde', as.character(formula_4[[2]])[3]))
  expect_false(grepl('PO_spde', as.character(formula_4[[1]])[3]))
  expect_false(grepl('PA_spde', as.character(formula_4[[2]])[3]))
  
  expect_true(grepl('a_intercept', as.character(formula_4[[1]])[3]))
  expect_true(grepl('b_intercept', as.character(formula_4[[1]])[3]))
  expect_false(grepl('c_intercept', as.character(formula_4[[1]])[3]))
  expect_false(grepl('d_intercept', as.character(formula_4[[1]])[3]))
  expect_false(grepl('a_intercept', as.character(formula_4[[2]])[3]))
  expect_false(grepl('b_intercept', as.character(formula_4[[2]])[3]))
  expect_true(grepl('c_intercept', as.character(formula_4[[2]])[3]))
  expect_true(grepl('d_intercept', as.character(formula_4[[2]])[3]))
  
  expect_true(grepl('a_Elevation', as.character(formula_4[[1]])[3]))
  expect_true(grepl('b_Elevation', as.character(formula_4[[1]])[3]))
  expect_false(grepl('c_Elevation', as.character(formula_4[[1]])[3]))
  expect_false(grepl('d_Elevation', as.character(formula_4[[1]])[3]))
  expect_false(grepl('a_Elevation', as.character(formula_4[[2]])[3]))
  expect_false(grepl('b_Elevation', as.character(formula_4[[2]])[3]))
  expect_true(grepl('c_Elevation', as.character(formula_4[[2]])[3]))
  expect_true(grepl('d_Elevation', as.character(formula_4[[2]])[3]))
  
  expect_true(grepl('a_Temperature', as.character(formula_4[[1]])[3]))
  expect_true(grepl('b_Temperature', as.character(formula_4[[1]])[3]))
  expect_false(grepl('c_Temperature', as.character(formula_4[[1]])[3]))
  expect_false(grepl('d_Temperature', as.character(formula_4[[1]])[3]))
  expect_false(grepl('a_Temperature', as.character(formula_4[[2]])[3]))
  expect_false(grepl('b_Temperature', as.character(formula_4[[2]])[3]))
  expect_true(grepl('c_Temperature', as.character(formula_4[[2]])[3]))
  expect_true(grepl('d_Temperature', as.character(formula_4[[2]])[3]))
  
  expect_true(grepl('a_Precipitation', as.character(formula_4[[1]])[3]))
  expect_true(grepl('b_Precipitation', as.character(formula_4[[1]])[3]))
  expect_false(grepl('c_Precipitation', as.character(formula_4[[1]])[3]))
  expect_false(grepl('d_Precipitation', as.character(formula_4[[1]])[3]))
  expect_false(grepl('a_Precipitation', as.character(formula_4[[2]])[3]))
  expect_false(grepl('b_Precipitation', as.character(formula_4[[2]])[3]))
  expect_true(grepl('c_Precipitation', as.character(formula_4[[2]])[3]))
  expect_true(grepl('d_Precipitation', as.character(formula_4[[2]])[3]))
  
  
  })