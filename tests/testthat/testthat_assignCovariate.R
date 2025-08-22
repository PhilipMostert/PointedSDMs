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
PO$species <- sample(x = c('fish1', 'fish2'), size = nrow(PO), replace = TRUE)
PO$temp <- sample(x = c(1,2), nrow(PO), replace = TRUE)
#Random presence absence dataset
PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
st_geometry(PA) <- 'geometry'
PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
#Add trial name
PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
PA$pointcov <- runif(n = nrow(PA))
PA$binommark <- sample(x = 0:1, size = nrow(PA), replace = TRUE)
PA$marktrial <- sample(x = 1:3, size = nrow(PA), replace = TRUE)
PA$species <- sample(x = c('bird1', 'bird2'), nrow(PA), replace = TRUE)
PA$temp <- sample(x = c(1,2), nrow(PA), replace = TRUE)
mesh <- fmesher::fm_mesh_2d_inla(boundary = fmesher::fm_as_segm(SpatialPoly), 
                                 max.edge = 2, crs = fmesher::fm_crs(projection))
#iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
iPoints <- fmesher::fm_int(samplers = SpatialPoly, domain = mesh)

spatialcovariates <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
terra::values(spatialcovariates) <- rgamma(n = terra::ncell(cov), shape = 2)
names(spatialcovariates) <- 'covariate'


spData <- list(PO, PA)
Check <- dataOrganize$new() 

Check$makeData(datapoints = spData, datanames = c('PO', 'PA'),
               coords = c('long', 'lat'), proj = projection, offsetname = NULL,
               pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
               speciesname = 'species', marks = NULL, temporalvar = 'temp',
               marktrialname = NULL, markfamily = NULL)


test_that('assignCovariate is able to correctly assign covariate information', {
  
  assign1 <- PointedSDMs:::assignCovariate(data = Check$Data, covariateEnv = environment(), 
                                covariateNames = 'covariate', projection = projection)
  
  expect_equal(class(assign1), 'list')
  expect_setequal(names(assign1), c('PO', 'PA'))
  
  expect_setequal(names(assign1$PO[[1]]), c(names(Check$Data$PO[[1]]), 'covariate'))
  expect_setequal(names(assign1$PA[[1]]), c(names(Check$Data$PA[[1]]), 'covariate'))
  
  expect_true(is.numeric(assign1$PO[[1]]$covariate))
  expect_true(is.numeric(assign1$PA[[1]]$covariate))
  
  assignIPS <- PointedSDMs:::assignCovariate(data = list(iPoints = iPoints), covariateEnv = environment(), 
                                           covariateNames = 'covariate', projection = projection, IPS = TRUE)
  
  
  expect_equal(class(assignIPS), c('sf', 'data.frame'))
  
  expect_setequal(names(assignIPS), c(names(iPoints), 'covariate'))

  expect_true(is.numeric(assignIPS$covariate))
  
  
  ##Check if multi-species works
  Check$makeSpecies(speciesname = 'species')
  
  assign2 <- PointedSDMs:::assignCovariate(data = Check$Data, covariateEnv = environment(), 
                                           covariateNames = 'covariate', projection = projection, 
                                           speciesName = 'species')
  
  expect_equal(class(assign2), 'list')
  expect_setequal(names(assign2), c('PO', 'PA'))
  
  expect_setequal(names(assign2$PO[[1]]), c(names(Check$Data$PO[[1]]), 'covariate'))
  expect_setequal(names(assign2$PA[[1]]), c(names(Check$Data$PA[[1]]), 'covariate'))
  
  expect_true(is.numeric(assign2$PO[[1]]$covariate))
  expect_true(is.numeric(assign2$PA[[1]]$covariate))
  
  assignIPSSpecies <- PointedSDMs:::assignCovariate(data = list(iPoints = iPoints), covariateEnv = environment(), 
                                             covariateNames = 'covariate', projection = projection, 
                                             IPS = TRUE, speciesName = 'species')
  
  
  expect_equal(class(assignIPSSpecies), c('sf', 'data.frame'))
  
  expect_setequal(names(assignIPSSpecies), c(names(iPoints), 'covariate'))
  
  expect_true(is.numeric(assignIPSSpecies$covariate))
  
  
})
  