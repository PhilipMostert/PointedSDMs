test_that('dataSet transforms the data to SpatialPointsDataFrames, and produces the relevant metadata.', {
  
  ##Generate arbitrary data:
  #Arbitrary projection
  projection <- CRS('+proj=tmerc')
  
  #Make random shape to generate points on
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
  PO$species <- sample(x = c('fish1', 'fish2'), size = nrow(PO@coords), replace = TRUE)
  PO$temp <- rpois(n = nrow(PO@coords), lambda = 5)
  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
  PA$pointcov <- runif(n = nrow(PA@coords))
  PA$binommark <- sample(x = 2:5, size = nrow(PA@data), replace = TRUE)
  PA$marktrial <- sample(x = 0:1, size = nrow(PA@data), replace = TRUE)
  PA$species <- sample(x = c('bird1', 'bird2'), nrow(PA@data), replace = TRUE)
  PA$temp <- rpois(n = nrow(PA@coords), lambda = 5)
  ##Make PA a data.frame object
  PA <- data.frame(PA)
  
  spData <- list(PO, PA)
  
  mod <- dataSet(datapoints = spData, datanames = c('PO', 'PA'),
                 coords = colnames(PO@coords), proj = projection,
                 pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                 speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'),
                 marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial', 'binomial'),
                 temporalvar = 'temp', offsetname = NULL)
  
  expect_setequal(names(mod), c("Data", "Family", "dataType", "varsIn", "Marks",
                                "marksType", "multinomVars", 'numObs'))
  expect_setequal(names(mod$Data), c('PO','PA'))
  
  expect_true(all(lapply(unlist(mod$Data), function(x) inherits(x, 'Spatial'))))
  
  ##Should create a placeholder variable for the poresp + 
  #should keep marks +
  #should create new variables for the multinomial marks.
  expect_setequal(names(mod$Data$PO[[1]]@data), c("poresp", "numvar", "factvar", 'temp',
                                                  "species", "factvar_phi", "factvar_response"))
  expect_true((all(mod$Data$PO[[1]]@data$factvar_phi == 1)))
  expect_true((all(mod$Data$PO[[1]]@data$factvar_response == 1)))
  expect_true(class(mod$Data$PO[[1]]@data$factvar) == 'character')
  
  expect_setequal(names(mod$Data$PA[[1]]@data), c("PAresp", "trial", "binommark", 'temp',
                                                  "marktrial", "species", "pointcov"))
  
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
  
  expect_true(mod$numObs[1] == nrow(PO@data))
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
