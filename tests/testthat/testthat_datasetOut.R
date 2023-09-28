test_that('datasetOut is able to correctly remove the correct datasets and metadata from the integrated model.', {
  
  ##First generate a model:
   ##Just import a model??
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
  PA$binommark <- sample(x = 0:1, size = nrow(PA), replace = TRUE)
  PA$marktrial <- sample(x = 2:5, size = nrow(PA), replace = TRUE)
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
  
  cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
  terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
  names(cov) <- 'covariate'
  
  obj <- intModel(PO, PA, Coordinates = coordnames, Projection = projection, Mesh = mesh,
                 IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, speciesSpatial = 'individual',
                 responsePA = responsePA, markNames = c('factvar'), markFamily = c('multinomial'),
                 speciesName = speciesName, spatialCovariates = cov, pointsSpatial = NULL)
  
  spatMod <- fitISDM(data = obj,
                      options  = list(control.inla=list(int.strategy='eb')))
  
  ##Remove PO from the model
  outPO <- datasetOut(model = spatMod, dataset = 'PO', predictions = T)
  
  expect_setequal(class(outPO), c("datasetOut", "list"))
  expect_true(names(outPO) == 'Leaving_out_PO')
  expect_output(print(outPO),'Changes in fixed values by leaving out PO:')
  expect_output(print(outPO), 'Leave-one-out cross-validation score:')
  
  #ie no fish from dataset PO
  expect_setequal(outPO$Leaving_out_PO$names.fixed, c('bird_covariate', 'bird_intercept'))
  expect_setequal(names(outPO$Leaving_out_PO$bru_info$lhoods), 'PA_bird_PAresp')
  
  expect_setequal(outPO$Leaving_out_PO$source, 'PA')
  expect_false('PO' %in% names(outPO$Leaving_out_PO$species$speciesIn))
  expect_false('PO' %in% names(outPO$Leaving_out_PO$marks$marksIn))
  expect_length(outPO$Leaving_out_PO$optionsJoint$control.family, 1)
  

  })
  
  
