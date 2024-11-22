test_that('makeLhoods makes a list of likelihoods', {
  skip_on_cran()
  
  ##Set up a model
  ##Set up arbitrary data
  projection <- '+proj=tmerc'
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  xy <- cbind(x, y)
  SpatialPoly <- st_sfc(st_polygon(list(xy)), crs = projection)
  
  ##Old coordinate names
  #Make random points
  #Random presence only dataset
  PO <- st_as_sf(st_sample(SpatialPoly, size = 100, crs = projection))
  st_geometry(PO) <- 'geometry'
  ##Add random variable
  PO$numvar <- runif(n = nrow(PO))
  PO$factvar <- sample(x = c('a','b'), size = nrow(PO), replace = TRUE)
  PO$species <- sample(x = c('fish'), size = nrow(PO), replace = TRUE)
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, size = 100, crs = projection))
  st_geometry(PA) <- 'geometry'
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  PA$pointcov <- runif(n = nrow(PA))
  PA$binommark <- sample(x = 2:3, size = nrow(PA), replace = TRUE)
  PA$marktrial <- sample(x = 3:5, size = nrow(PA), replace = TRUE)
  PA$species <- sample(x = c('bird'), nrow(PA), replace = TRUE)
  mesh <- fmesher::fm_mesh_2d_inla(boundary = fmesher::fm_as_segm(SpatialPoly), 
                                   max.edge = 2, crs = fmesher::fm_crs(projection))
  iPoints <- fmesher::fm_int(samplers = SpatialPoly, domain = mesh)

  coordnames <- c('x.1', 'x.2')
  responseCounts <- 'count'
  responsePA <- 'PAresp'
  trialName <- 'trial'
  markNames <- c('numvar', 'factvar', 'binommark')
  marksFamily <- c('gaussian', 'multinomial', 'binomial')
  markTrial = 'marktrial'
  pointCovs <- 'pointcov'
  speciesName <- 'species'
  
  cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
  terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
  
  
  obj <- startSpecies(PO, PA, Projection = projection, Mesh = mesh,
                  IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                  responsePA = responsePA, speciesName = speciesName, spatialCovariates = cov)
  
  Lhoods <- makeLhoods(data = obj$.__enclos_env__$private$modelData,
                       formula = obj$.__enclos_env__$private$Formulas,
                       family = obj$.__enclos_env__$private$Family,
                       mesh = obj$.__enclos_env__$private$INLAmesh,
                       ips = obj$.__enclos_env__$private$IPS,
                       samplers = obj$.__enclos_env__$private$Samplers,
                       paresp = obj$.__enclos_env__$private$responsePA,
                       ntrialsvar = trialName,
                       markstrialsvar = markTrial,
                       speciesname = speciesName,
                       speciesindex = obj$.__enclos_env__$private$speciesIndex)
  
  #Expect a list of Lhoods
  expect_type(Lhoods, 'list')
  expect_true(all(vapply(Lhoods, function(x) inherits(x, "bru_like"), logical(1))))
  
  #Names should be in format: dataset_species_response
  expect_setequal(names(Lhoods), c("PO_fish_geometry", "PA_bird_PAresp"))

  #Expect familus
  expect_setequal(sapply(Lhoods, function(x) x$family), c( "cp", "binomial"))
  
  #Trials are correct for PA response
  expect_identical(Lhoods$PA_bird_PAresp$Ntrials, PA$trial)
  
  #Correct formulas
  expect_equal(deparse1(Lhoods$PO_fish_geometry$formula), 'geometry ~ .')
  expect_equal(deparse1(Lhoods$PA_bird_PAresp$formula), 'PAresp ~ .')
  
  
  
  })
