test_that('runModel runs a dataSDM object, and produces an INLA model with extra metadata.', {
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
                IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                responsePA = responsePA, markNames = NULL, markFamily = NULL,
                speciesName = speciesName, spatialCovariates = cov)
  
  ##run model
  spatMod <- runModel(data = obj, options  = list(control.inla=list(int.strategy='eb')))
  
  expect_setequal(class(spatMod), c('bruSDM', 'bru', 'iinla', 'inla'))
  expect_setequal(row.names(spatMod$summary.fixed), c( "fish_covariate", "bird_covariate",'fish_intercept', 'bird_intercept'))
  expect_equal(spatMod[['species']][['speciesVar']], 'species')
  expect_equal(spatMod[['species']][['speciesIn']][['PO']], 'fish')
  expect_equal(spatMod[['species']][['speciesIn']][['PA']], 'bird')
  
  expect_equal(deparse1(spatMod[['componentsJoint']]), "~-1 + shared_spatial(main = coordinates, model = shared_field) + fish_spatial(main = coordinates, model = fish_field) + bird_spatial(main = coordinates, model = bird_field) + fish_covariate(main = fish_covariate, model = \"linear\") + bird_covariate(main = bird_covariate, model = \"linear\") + fish_intercept(1) + bird_intercept(1)")
  
  
  expect_output(summary(spatMod), 'Summary for fish:')
  expect_output(summary(spatMod), 'Summary for bird:')
  
  expect_equal(spatMod[['spatCovs']][['name']], 'covariate')
  expect_equal(spatMod[['spatCovs']][['class']], c(covariate = 'numeric'))
  
  ##Run with marks
  obj2 <- intModel(PO, PA, Coordinates = coordnames, Projection = projection, Mesh = mesh,
                IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                responsePA = responsePA, markNames = c('factvar', 'binommark'), markFamily = c('multinomial', 'binomial'),
                speciesName = speciesName, spatialCovariates = cov, trialsMarks = 'marktrial')

  spatMod2 <- runModel(data = obj2,
                       options  = list(control.inla=list(int.strategy='eb')))
  
  expect_setequal(spatMod2[['source']], c('PO', 'PO', 'PA', 'PA'))
  expect_setequal(sapply(spatMod2$.args$control.family, function(x) x$link), c('log', 'log', 'cloglog', 'cloglog'))
  
  })
