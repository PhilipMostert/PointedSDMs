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
  mesh <- fmesher::fm_mesh_2d_inla(boundary = fmesher::fm_as_segm(SpatialPoly), 
                                   max.edge = 2, crs = fmesher::fm_crs(projection))
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
  
  cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
  terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
  names(cov) <- 'covariate'
  
  obj <- startISDM(PO, PA, Projection = projection, Mesh = mesh, spatialCovariates = cov,
                 IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                 responsePA = responsePA, pointsIntercept = TRUE,pointsSpatial = NULL)
  
  spatMod <- fitISDM(data = obj,
                      options  = list(control.inla=list(int.strategy='eb')))
  
  ##Remove PO from the model
  outPO <- datasetOut(model = spatMod, dataset = 'PO', predictions = T)
  
  expect_setequal(class(outPO), c("datasetOut", "list"))
  expect_true(names(outPO) == 'Leaving_out_PO')
  expect_output(print(outPO),'Changes in fixed values by leaving out PO:')
  expect_output(print(outPO), 'Leave-one-out cross-validation score:')
  
  #ie no fish from dataset PO
  expect_setequal(outPO$Leaving_out_PO$names.fixed, c("covariate", "PA_intercept"))
  expect_setequal(names(outPO$Leaving_out_PO$bru_info$lhoods), 'PA_PAresp')
  
  expect_setequal(outPO$Leaving_out_PO$source, 'PA')
  expect_false('PO' %in% names(outPO$Leaving_out_PO$species$speciesIn))
  expect_false('PO' %in% names(outPO$Leaving_out_PO$marks$marksIn))
  expect_length(outPO$Leaving_out_PO$optionsJoint$control.family, 1)
  
  #Try copy spatial
  obj2 <- startISDM(PO, PA, Projection = projection, Mesh = mesh,
                   IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                   responsePA = responsePA, pointsIntercept = TRUE)
  
  spatMod2 <- fitISDM(data = obj2,
                     options  = list(control.inla=list(int.strategy='eb')))
  
  outPO2 <- datasetOut(model = spatMod2, dataset = 'PO', predictions = T)
  
  expect_true('PA_spatial' %in% names(outPO2$Leaving_out_PO$summary.random))
  expect_identical(deparse1(outPO2$Leaving_out_PO$componentsJoint), '~PA_intercept(1) + PA_spatial(main = geometry, model = PO_field) - 1')

  })
  
  
