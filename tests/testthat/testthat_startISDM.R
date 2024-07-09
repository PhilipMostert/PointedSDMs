test_that('startISDM is able to initialize a specifyISDM object as well as correctly add the data to the integrated model', {
  
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
  PO$temp <- sample(x = 1:2, size = nrow(PO), replace = TRUE)
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PA) <- 'geometry'
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  PA$pointcov <- runif(n = nrow(PA))
  PA$temp <- sample(x = 1:2, size = nrow(PA), replace = TRUE)
  if (requireNamespace("INLA")) {
    mesh <<- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(SpatialPoly), 
                                max.edge = 2, crs = inlabru::fm_crs(projection))
    #iPoints <<- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
    iPoints <<- inlabru::fm_int(samplers = SpatialPoly, domain = mesh)
    
  }
  #iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  iPoints <- inlabru::fm_int(samplers = SpatialPoly, domain = mesh)
  ##Make PA a data.frame object

  responseCounts <- 'count'
  responsePA <- 'PAresp'
  trialName <- 'trial'
  pointCovs <- 'pointcov'
  
  cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
  terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
  names(cov) <- 'covariate'
  cov$cov2 <-  rgamma(n = terra::ncell(cov), shape = 2)
  
  
  obj <- startISDM(PO, PA, Projection = projection, Mesh = mesh,
                  IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                  responsePA = responsePA, spatialCovariates = NULL)
  
  #Test that object is created
  expect_true(all(class(obj) == c('specifyISDM', 'R6')))
  expect_setequal(names(obj$.__enclos_env__$private$modelData), c("PO", "PA"))
  expect_setequal(unlist(obj$.__enclos_env__$private$Family), c('cp', 'binomial'))
  expect_identical(iPoints, expected = obj$.__enclos_env__$private$IPS)
  
  ##Test object can be created from a list of data
  obj2 <- startISDM(list(PO, PA), Projection = projection, Mesh = mesh,
                          IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                          responsePA = responsePA, spatialCovariates = NULL)
  expect_setequal(names(obj2$.__enclos_env__$private$modelData), c("PO", "PA"))
  
  ##Test message if Temporal but pointsSpatial = 'copy'
  objMess <- expect_message(startISDM(list(PO, PA), Projection = projection, Mesh = mesh,
                            IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                            responsePA = responsePA, spatialCovariates = NULL, temporalName = 'temp',
                            pointsSpatial = 'correlate'),'Setting pointsSpatial to "shared" since it is required for the temporalModel.')
  
  ##Test error: data not sf
  PA2 <- as.data.frame(PA)
  expect_error(startISDM(PO, PA2, Projection = projection, Mesh = mesh,
                         IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                         responsePA = responsePA, spatialCovariates = NULL),'Datasets need to be sf objects.')
  
  ##Test error: pointsSpatial
  expect_error(startISDM(PO, PA2, Projection = projection, Mesh = mesh,
                         IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                         responsePA = responsePA, spatialCovariates = NULL, pointsSpatial = FALSE),'PointsSpatial needs to be one of: "shared", "copy", "individual", "correlate" or NULL.')
  expect_error(startISDM(PO, PA2, Projection = projection, Mesh = mesh,
                         IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                         responsePA = responsePA, spatialCovariates = NULL, pointsSpatial = c('shared', 'copy')),'PointsSpatial needs to be one of: "shared", "copy", "individual", "correlate" or NULL.')
  
  ##Test error: INLAmesh not an inla.mesh object
  expect_error(startISDM(PO, PA, Projection = projection, Mesh = list(),
                         IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                         responsePA = responsePA, spatialCovariates = NULL), 'Mesh needs to be a inla.mesh object.')
  #Try with cov
  objCov <- startISDM(PO, PA, Projection = projection, Mesh = mesh,
                        IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                        responsePA = responsePA, spatialCovariates = cov)
  expect_setequal(objCov$.__enclos_env__$private$spatcovsNames, c('covariate', 'cov2'))
  expect_setequal(objCov$.__enclos_env__$private$spatcovsClass, c('linear', 'linear'))
  
  #Try own formula
  objForm <- startISDM(PO, PA, Projection = projection, Mesh = mesh,
                      IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                      responsePA = responsePA, spatialCovariates = cov,
                      Formulas = list(covariateFormula = ~ covariate + I(covariate^2),
                                      biasFormula = ~ cov2))
  
  expect_equal(deparse1(objForm$.__enclos_env__$private$covariateFormula), "~covariate + I(covariate^2)")
  expect_equal(deparse1(objForm$.__enclos_env__$private$biasFormula), "~cov2")
  
  ##Expect output
  expect_output(obj$print(), 'Summary of startISDM data file:')
  expect_output(obj$print(), 'Summary of presence absence datasets')
  expect_output(obj$print(), 'Summary of presence only datasets')
  
  
})
