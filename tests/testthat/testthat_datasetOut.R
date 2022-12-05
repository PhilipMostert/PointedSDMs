test_that('datasetOut is able to correctly remove the correct datasets and metadata from the integrated model.', {
  
  ##First generate a model:
   ##Just import a model??
  skip_on_cran()
  
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
  #iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  iPoints <- inlabru::ipoints(samplers = SpatialPoly)

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
                 IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
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
  
  
