test_that('intModel is able to initialize a dataSDM object as well as correctly add the data to the integrated model', {

  testthat::skip('Function is outdated')  
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
  mesh <- fmesher::fm_mesh_2d_inla(boundary = fmesher::fm_as_segm(SpatialPoly), 
                             max.edge = 2, crs = fmesher::fm_crs(projection))
  #iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  iPoints <- fmesher::fm_int(samplers = SpatialPoly, domain = mesh)

  
  coordnames <- c('long', 'lat')
  responseCounts <- 'count'
  responsePA <- 'PAresp'
  trialName <- 'trial'
  markNames <- c('numvar', 'factvar', 'binommark')
  marksFamily <- c('gaussian', 'multinomial', 'binomial')
  markTrial = 'marktrial'
  pointCovs <- 'pointcov'
  speciesName <- 'species'

  
  expect_warning(intModel(PO, PA, Projection = projection, Mesh = mesh,
                IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                responsePA = responsePA), 'This function has been depreciated for startISDM. If you wish to create a multi-species model, please use startSpecies.')
  
  obj <- suppressWarnings(intModel(PO, PA, Projection = projection, Mesh = mesh,
           IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
           responsePA = responsePA))
  expect_true(inherits(obj, 'specifyISDM'))
  
})
