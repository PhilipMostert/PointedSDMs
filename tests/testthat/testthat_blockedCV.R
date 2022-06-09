test_that('blockedCV completes spatial block cross-validation.', {
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
  
  obj <- intModel(PO, PA, Coordinates = coordnames, Projection = projection, Mesh = mesh,
                  IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                  responsePA = responsePA, markNames = NULL, markFamily = NULL,
                  speciesName = speciesName)
  
  obj$spatialBlock(k = 2, rows = 2, cols = 1)
  
  ##run model
  blocked <- blockedCV(data = obj, options  = list(control.inla=list(int.strategy='eb')))
  
  expect_setequal(class(blocked), c("blockedCV", "list"))
  expect_setequal(names(blocked), c( "DIC_fold_1", "DIC_fold_2",'Formula'))
  expect_equal(class(blocked$Formula), 'formula')
  expect_output(print(blocked), 'Spatial block cross-validation score:')
  expect_output(print(blocked), 'mean DIC score:')

})
