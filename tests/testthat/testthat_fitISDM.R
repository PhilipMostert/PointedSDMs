test_that('fitISDM runs a dataSDM object, and produces an INLA model with extra metadata.', {
  skip_on_cran()
  
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
                IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                responsePA = responsePA, markNames = NULL, markFamily = NULL, speciesSpatial = 'individual',
                speciesName = speciesName, spatialCovariates = cov)
  
  ##run model
  spatMod <- fitISDM(data = obj, options  = list(control.inla=list(int.strategy='eb')))
  
  expect_setequal(class(spatMod), c('bruSDM', 'bru', 'iinla', 'inla'))
  expect_setequal(row.names(spatMod$summary.fixed), c( "fish_covariate", "bird_covariate",'fish_intercept', 'bird_intercept'))
  expect_equal(spatMod[['species']][['speciesVar']], 'species')
  expect_equal(spatMod[['species']][['speciesIn']][['PO']], 'fish')
  expect_equal(spatMod[['species']][['speciesIn']][['PA']], 'bird')
  
  expect_equal(deparse1(spatMod[['componentsJoint']]), "~-1 + shared_spatial(main = geometry, model = shared_field) + fish_PO_spatial(main = geometry, model = fish_PO_field) + bird_PA_spatial(main = geometry, model = bird_PA_field) + fish_covariate(main = fish_covariate, model = \"linear\") + bird_covariate(main = bird_covariate, model = \"linear\") + fish_intercept(1) + bird_intercept(1)")
  
  
  expect_output(summary(spatMod), 'Summary for fish:')
  expect_output(summary(spatMod), 'Summary for bird:')
  
  expect_equal(spatMod[['spatCovs']][['name']], 'covariate')
  expect_equal(spatMod[['spatCovs']][['class']], c(covariate = 'numeric'))
  
  ##Run with marks
  obj2 <- intModel(PO, PA, Coordinates = coordnames, Projection = projection, Mesh = mesh,
                IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts, 
                responsePA = responsePA, markNames = c('factvar', 'binommark'), markFamily = c('multinomial', 'binomial'),
                speciesName = speciesName, spatialCovariates = cov, trialsMarks = 'marktrial')

  spatMod2 <- fitISDM(data = obj2,
                       options  = list(control.inla=list(int.strategy='eb')))
  
  expect_setequal(spatMod2[['source']], c('PO', 'PO', 'PA', 'PA'))
  expect_setequal(sapply(spatMod2$.args$control.family, function(x) x$link), c('log', 'log', 'cloglog', 'cloglog'))
  
  })
