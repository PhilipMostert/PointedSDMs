test_that('dataSet transforms the data to SpatialPointsDataFrames, and produces the relevant metadata.', {
  
  ##Generate arbitrary data:
  #Arbitrary projection
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
  PO$temp <- rpois(n = nrow(PO), lambda = 5)
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PA) <- 'geometry'
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  PA$pointcov <- runif(n = nrow(PA))
  PA$binommark <- sample(x = 0:1, size = nrow(PA), replace = TRUE)
  PA$marktrial <- sample(x = 2:5, size = nrow(PA), replace = TRUE)
  PA$species <- sample(x = c('bird1', 'bird2'), nrow(PA), replace = TRUE)
  PA$temp <- rpois(n = nrow(PA), lambda = 5)

  
  spData <- list(PO, PA)
  
  mod <- dataSet(datapoints = spData, datanames = c('PO', 'PA'),
                 coords = c('long', 'lat'), proj = projection,
                 pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                 speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'),
                 marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial', 'binomial'),
                 temporalvar = 'temp', offsetname = NULL)
  
  expect_setequal(names(mod), c("Data", "Family", "dataType", "varsIn", "Marks",
                                "marksType", "multinomVars", 'numObs'))
  expect_setequal(names(mod$Data), c('PO','PA'))
  
  expect_true(all(unlist(lapply(unlist(mod$Data, recursive = FALSE), function(x) inherits(x, 'sf')))))
  
  ##Should create a placeholder variable for the poresp + 
  #should keep marks +
  #should create new variables for the multinomial marks.
  expect_setequal(names(mod$Data$PO[[1]]), c("poresp", "numvar", "factvar", 'temp', 'geometry','._dataset_index_var_.',
                                                  "species", "factvar_phi", "factvar_response", 'speciesINDEX_VAR'))
  expect_true((all(mod$Data$PO[[1]]$factvar_phi == 1)))
  expect_true((all(mod$Data$PO[[1]]$factvar_response == 1)))
  expect_true(class(mod$Data$PO[[1]]$factvar) == 'character')
  
  expect_setequal(names(mod$Data$PA[[1]]), c("PAresp", "trial", "binommark", 'temp', 'geometry', '._dataset_index_var_.',
                                                  "marktrial", "species", "pointcov", 'speciesINDEX_VAR'))
  
  #Family for PO should be:
  # cp for the points;
  # uniform for the mark;
  # poisson for the multinomial mark.
  expect_setequal(mod$Family$PO, c('cp', 'uniform', 'poisson'))
  
  expect_setequal(mod$Family$PA, c('binomial', 'binomial'))
  
  expect_named(mod$dataType, c('PO', 'PA'))
  expect_setequal(mod$dataType, c('Present only', 'Present absence'))
  
  ##PO has no point covariates;
  ##PA has pointcov as a pointcovariate
  expect_true(is.null(unlist(mod$varsIn$PO)))
  expect_true(unlist(mod$varsIn$PA) == 'pointcov')
  
  expect_setequal(mod$Marks$PO, c('numvar', 'factvar'))
  expect_setequal(mod$Marks$PA, c('binommark'))
  
  expect_named(mod$marksType$PO, c('numvar', 'factvar'))
  expect_named(mod$marksType$PA, c('binommark'))
  
  expect_setequal(mod$marksType$PO, c('Uniform mark', 'Multinomial mark'))
  expect_setequal(mod$marksType$PA, c('Binomial mark'))
  
  expect_true('factvar' %in% mod$multinomVars)
  
  expect_true(mod$numObs[1] == nrow(PO))
  expect_true(mod$numObs[2] == nrow(PA))
  
  #Remove a dataset name
  expect_error(dataSet(datapoints = spData, datanames = c('PO'),
                       coords = colnames(PO@coords), proj = projection, offsetname = NULL,
                       pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                       speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'),
                       marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial', 'binomial')),
               'Number of dataset names needs to equal length of datasets.')
  
  #Remove a mark family
  expect_error(dataSet(datapoints = spData, datanames = c('PO','PA'),
                       coords = colnames(PO@coords), proj = projection, offsetname = NULL,
                       pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                       speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'),
                       marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial')),
               "Number of marks needs to equal the number of mark families.")
  
  })
