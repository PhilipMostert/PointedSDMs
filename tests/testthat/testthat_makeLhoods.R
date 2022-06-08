test_that('makeLhoods makes a list of likelihoods', {
  skip_on_cran()
  
  ##Set up a model
  ##Set up arbitrary data
  projection <- CRS('+proj=tmerc')
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  
  Poly = Polygon(xy)
  Poly = Polygons(list(Poly),1)
  SpatialPoly = SpatialPolygons(list(Poly), proj4string = projection)
  
  ##Old coordinate names
  #Make random points
  #Random presence only dataset
  PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  ##Add random variable
  PO$numvar <- runif(n = nrow(PO@coords))
  PO$factvar <- sample(x = c('a','b'), size = nrow(PO@coords), replace = TRUE)
  PO$species <- sample(x = c('fish'), size = nrow(PO@coords), replace = TRUE)
  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
  PA$pointcov <- runif(n = nrow(PA@coords))
  PA$binommark <- sample(x = 2:3, size = nrow(PA@data), replace = TRUE)
  PA$marktrial <- sample(x = 3:5, size = nrow(PA@data), replace = TRUE)
  PA$species <- sample(x = c('bird'), nrow(PA@data), replace = TRUE)
  mesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(SpatialPoly), 
                             max.edge = 2)
  iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  ##Make PA a data.frame object
  PA <- data.frame(PA)
  
  coordnames <- colnames(PO@coords)
  responseCounts <- 'count'
  responsePA <- 'PAresp'
  trialName <- 'trial'
  markNames <- c('numvar', 'factvar', 'binommark')
  marksFamily <- c('gaussian', 'multinomial', 'binomial')
  markTrial = 'marktrial'
  pointCovs <- 'pointcov'
  speciesName <- 'species'
  
  cov <- sp::spsample(x = SpatialPoly, n = 1000, type = 'stratified')
  cov$covariate <- rgamma(n = nrow(cov@coords), shape = 2)
  cov <- sp::SpatialPixelsDataFrame(points = cov@coords,
                                    data = data.frame(covariate = cov$covariate),
                                    proj4string = projection,
                                    tolerance = 0.585235)
  cov <- raster::raster(cov)
  
  
  obj <- intModel(PO, PA, Coordinates = coordnames, Projection = projection, Mesh = mesh,
                  IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, trialsMarks = markTrial,
                  responsePA = responsePA, markNames = c('numvar', 'binommark'), markFamily = c('gaussian', 'binomial'),
                  speciesName = speciesName, spatialCovariates = cov)
  
  Lhoods <- makeLhoods(data = obj$.__enclos_env__$private$modelData,
                       formula = obj$.__enclos_env__$private$Formulas,
                       family = obj$.__enclos_env__$private$Family,
                       mesh = obj$.__enclos_env__$private$INLAmesh,
                       ips = obj$.__enclos_env__$private$IPS,
                       paresp = obj$.__enclos_env__$private$responsePA,
                       ntrialsvar = trialName,
                       markstrialsvar = markTrial,
                       speciesname = speciesName,
                       speciesindex = obj$.__enclos_env__$private$speciesIndex)
  
  #Expect a list of Lhoods
  expect_type(Lhoods, 'list')
  expect_true(all(unlist(lapply(Lhoods, class)) == c("bru_like", "list")))
  
  #Names should be in format: dataset_species_response
  expect_setequal(names(Lhoods), c("PO_fish_coordinates", "PO_fish_numvar", "PA_bird_PAresp", "PA_bird_binommark"))

  #Expect familus
  expect_setequal(sapply(Lhoods, function(x) x$family), c( "cp", "gaussian", "binomial", "binomial"))
  
  #Trials are correct for PA response
  expect_identical(Lhoods$PA_bird_PAresp$Ntrials, PA$trial)
  
  #Trials are correct for marks
  expect_identical(Lhoods$PA_bird_binommark$Ntrials, PA$marktrial)
  
  #Correct formulas
  expect_equal(deparse1(Lhoods$PO_fish_coordinates$formula), 'coordinates ~ .')
  expect_equal(deparse1(Lhoods$PO_fish_numvar$formula), 'numvar ~ .')
  expect_equal(deparse1(Lhoods$PA_bird_PAresp$formula), 'PAresp ~ .')
  expect_equal(deparse1(Lhoods$PA_bird_binommark$formula), 'binommark ~ .')
  
  
  
  })
