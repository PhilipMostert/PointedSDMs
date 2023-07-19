test_that('intModel is able to initialize a dataSDM object as well as correctly add the data to the integrated model', {
  skip_on_cran()
  
  ##Set up arbitrary data
  projection <- '+proj=tmerc'
  
  #Make random shape to generate points on
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  SpatialPoly <- st_sfc(st_polygon(list(xy)), crs = projection)
  
  ##Old coordinate names
  #Make random points
  #Random presence only dataset
  PO <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PO) <- 'geometry'
  ##Add random variable
  PO$numvar <- runif(n = nrow(PO))
  PO$factvar <- sample(x = c('a','b'), size = nrow(PO), replace = TRUE)
  PO$species <- sample(x = c('fish'), size = nrow(PO), replace = TRUE)
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PA) <- 'geometry'
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  PA$pointcov <- runif(n = nrow(PA))
  PA$binommark <- sample(x = 2:5, size = nrow(PA), replace = TRUE)
  PA$marktrial <- sample(x = 0:1, size = nrow(PA), replace = TRUE)
  PA$species <- sample(x = c('bird'), nrow(PA), replace = TRUE)
  mesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(SpatialPoly), 
                             max.edge = 2, crs = inlabru::fm_crs(projection))
  #iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  iPoints <- inlabru::fm_int(samplers = SpatialPoly, domain = mesh)
  ##Make PA a data.frame object
  PA$long <- st_coordinates(PA)[,1]
  PA$lat <- st_coordinates(PA)[,2]
  st_geometry(PA) <- NULL
  PA <- data.frame(PA)
  
  coordnames <- c('long', 'lat')
  responseCounts <- 'count'
  responsePA <- 'PAresp'
  trialName <- 'trial'
  markNames <- c('numvar', 'factvar', 'binommark')
  marksFamily <- c('gaussian', 'multinomial', 'binomial')
  markTrial = 'marktrial'
  pointCovs <- 'pointcov'
  speciesName <- 'species'

  
  obj <- intModel(PO, PA, Coordinates = coordnames, Projection = projection, Mesh = mesh,
                IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                responsePA = responsePA, markNames = markNames, markFamily = marksFamily,
                speciesName = speciesName)
  
  expect_true(all(class(obj) == c('dataSDM', 'R6')))
  expect_setequal(names(obj$.__enclos_env__$private$modelData), c("PO", "PA"))
  
  ##Test warnings: No data added
  expect_warning(intModel(Coordinates = coordnames, Projection = projection, Mesh = mesh,
                        IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                        responsePA = responsePA, markNames = markNames, markFamily = marksFamily,
                        speciesName = speciesName))
  
  ##Test error: Coordinates length > 2
  expect_error(intModel(PO, PA, Coordinates = c('x','y','z'), Projection = projection, Mesh = mesh,
                        IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                        responsePA = responsePA, markNames = markNames, markFamily = marksFamily,
                        speciesName = speciesName))
  
  ##Test error: Coordinates[1] == Coordinates[2]
  expect_error(intModel(PO, PA, Coordinates = c('x','x'), Projection = projection, Mesh = mesh,
                      IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                      responsePA = responsePA, markNames = markNames, markFamily = marksFamily,
                      speciesName = speciesName))
  
  ##Test error: proj not CRS
  expect_error(intModel(PO, PA, Coordinates = coordnames, Projection = list(), Mesh = mesh,
                      IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                      responsePA = responsePA, markNames = markNames, markFamily = marksFamily,
                      speciesName = speciesName))
  
  ##Test error: INLAmesh not an inla.mesh object
  expect_error(intModel(PO, PA, Coordinates = coordnames, Projection = Projection, Mesh = list(),
                      IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                      responsePA = responsePA, markNames = markNames, markFamily = marksFamily,
                      speciesName = speciesName))
  
})
