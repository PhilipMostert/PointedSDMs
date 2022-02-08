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
PO$species <- sample(x = c('fish'), size = nrow(PO@coords), replace = TRUE)
#Random presence absence dataset
PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
#Add trial name
PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
PA$pointcov <- runif(n = nrow(PA@coords))
PA$binommark <- sample(x = 2:5, size = nrow(PA@data), replace = TRUE)
PA$marktrial <- sample(x = 0:1, size = nrow(PA@data), replace = TRUE)
PA$species <- sample(x = c('bird'), nrow(PA@data), replace = TRUE)
mesh <- PointedSDMs::MakeSpatialRegion(bdry = SpatialPoly, meshpars = list(max.edge = 2))
mesh <- mesh$mesh
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
spatCovs <- NULL #for now
speciesName <- 'species'
markSpatial <- TRUE
marksIntercept <- TRUE

test_that('dataSDMs initialize works as expected.', {
  
  check <<- dataSDM$new(coordinates = coordnames,
                       projection = projection,
                       Inlamesh = mesh,
                       responsepa = responsePA,
                       trialspa = trialName,
                       responsecounts = responseCounts,
                       marksnames = markNames,
                       trialsmarks = markTrial,
                       marksfamily = marksFamily,
                       pointcovariates = pointCovs,
                       marksspatial = markSpatial,
                       marksintercept = marksIntercept,
                       spatialcovariates = spatCovs,
                       speciesname = speciesName,
                       ips = iPoints,
                       spatial = TRUE,
                       intercepts = TRUE)
  ##Test private classes
  expect_true(all(check$.__enclos_env__$private$Coordinates%in%coordnames))
  expect_true(identical(check$.__enclos_env__$private$Projection,projection))
  expect_true(identical(check$.__enclos_env__$private$INLAmesh, mesh))
  
  expect_true(check$.__enclos_env__$private$responsePA == responsePA)
  expect_true(check$.__enclos_env__$private$responseCounts == responseCounts)
  
  expect_true(all(check$.__enclos_env__$private$markNames == markNames))
  expect_true(all(check$.__enclos_env__$private$markFamily == marksFamily))
  expect_true(!is.null(check$.__enclos_env__$private$trialsMarks))
  
  expect_true(check$.__enclos_env__$private$Spatial)
  expect_true(check$.__enclos_env__$private$Intercepts)
  
  #Remove one of the response variable names
  expect_error(dataSDM$new(coordinates = coordnames,
                           projection = projection,
                           Inlamesh = mesh,
                           responsepa = NULL,
                           responsecounts = responseCounts,
                           marksnames = markNames,
                           marksfamily = marksFamily,
                           pointcovariates = pointCovs,
                           spatialcovariates = spatCovs,
                           ips = iPoints,
                           spatial = FALSE,
                           intercepts = FALSE),
               'At least one of responseCounts and responsePA are NULL. Please provide both.')
  
  #Test mesh error
  expect_error(dataSDM$new(coordinates = coordnames,
                           projection = projection,
                           Inlamesh = list(),
                           responsepa = responsePA,
                           responsecounts = responseCounts,
                           marksnames = markNames,
                           marksfamily = marksFamily,
                           pointcovariates = pointCovs,
                           spatialcovariates = spatCovs,
                           ips = iPoints,
                           spatial = FALSE,
                           intercepts = FALSE),
               'INLAmesh needs to be an inla.mesh object.')
  
  ##Test coordinates error
  expect_error(dataSDM$new(coordinates = 'long',
                           projection = projection,
                           Inlamesh = mesh,
                           responsepa = responsePA,
                           responsecounts = responseCounts,
                           marksnames = markNames,
                           marksfamily = marksFamily,
                           pointcovariates = pointCovs,
                           spatialcovariates = spatCovs,
                           ips = iPoints,
                           spatial = FALSE,
                           intercepts = FALSE),
               'Coordinates needs to be a vector of length 2 containing the coordinate names.')
  
  #And unique error
  expect_error(dataSDM$new(coordinates = c('long','long'),
                           projection = projection,
                           Inlamesh = mesh,
                           responsepa = responsePA,
                           responsecounts = responseCounts,
                           marksnames = markNames,
                           marksfamily = marksFamily,
                           pointcovariates = pointCovs,
                           spatialcovariates = spatCovs,
                           ips = iPoints,
                           spatial = FALSE,
                           intercepts = FALSE),
               'Coordinates need to be unique values.')
  
  ##Check projection error
  expect_error(dataSDM$new(coordinates = coordnames,
                           projection = list(),
                           Inlamesh = mesh,
                           responsepa = responsePA,
                           responsecounts = responseCounts,
                           marksnames = markNames,
                           marksfamily = marksFamily,
                           pointcovariates = pointCovs,
                           spatialcovariates = spatCovs,
                           ips = iPoints,
                           spatial = FALSE,
                           intercepts = FALSE),
               'Projection needs to be a CRS object.')
  

})

test_that('addData can correctly add and store the relevent metadata properly.', {
  
  ##Add the datasets
  check$addData(PO,PA)
  
  ##ie also inclusive of marks + species
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO_fish_coordinates", "PO_fish_numvar", "PO_fish_factvar_response", "PA_bird_PAresp", "PA_bird_binommark"))
  expect_true(all(unlist(lapply(check$.__enclos_env__$private$modelData, function(x) class(x))) == c('bru_like', 'list')))
  #ie also inclusive of marks
  expect_equal(check$.__enclos_env__$private$dataSource, c(rep('PO',3), rep('PA',2)))
  
  ##Check link functions
  check$.__enclos_env__$private$optionsINLA
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'log', 'log', 'cloglog', 'cloglog'))
  ##Add a new counts dataset; first try without a species variable (compulsory since speciesName specified in the initialize...)
  Pcount <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  Pcount$count <- rpois(n = nrow(Pcount@coords), lambda = 2)
  expect_error(check$addData(Pcount),'The species variable name is required to be present in all the datasets.')
  
  #Generate species
  Pcount$species <- 'dog'
  check$addData(Pcount)
  ##Check that the data has been added into the model
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO_fish_coordinates", "PO_fish_numvar", "PO_fish_factvar_response", "PA_bird_PAresp", "PA_bird_binommark", "Pcount_dog_count"))
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'log', 'log', 'cloglog', 'cloglog', 'log'))
  
  #Add a new PA dataset with a different response variable name
  PA2 <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA2$newResponse <- sample(0:1, nrow(PA2@coords), replace = TRUE)
  PA2$species <- 'insect'
  check$addData(PA2, responsePA = 'newResponse')
  
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO_fish_coordinates", "PO_fish_numvar", "PO_fish_factvar_response", "PA_bird_PAresp", "PA_bird_binommark", "Pcount_dog_count", 'PA2_insect_PAresp'))
  expect_setequal(check$.__enclos_env__$private$dataSource, c("PO", "PO", "PO", "PA", "PA", "Pcount", "PA2"))
  ##The correct species_field is in the components, ie the object updated itself correctly
  expect_true('species_spatial(main = coordinates, model = speciesModel, group = species, ngroup = 4)' %in% check$.__enclos_env__$private$Components)
  expect_equal(unlist(check$.__enclos_env__$private$speciesIn), c(PO = 'fish', PA = 'bird', Pcount = 'dog', PA2 = 'insect'))
  
  })

#Make a random spatial covariate
cov <- sp::spsample(x = SpatialPoly, n = 100000, type = 'random')
cov$covariate <- rgamma(n = 100000, shape = 2)
cov <- sp::SpatialPixelsDataFrame(points = cov@coords,
                                  data = data.frame(covariate = cov$covariate),
                                  proj4string = projection,
                                  tolerance = 0.898631)

test_that('spatialCovariates can add spatial covariates to the model and succesfully update the formulas and components of the model.', {
  
  check$spatialCovariates(cov)
  
  #Is the new covariate in the formulas
  expect_true(all(unlist(lapply(check$.__enclos_env__$private$modelData, function(x) {
    
    'covariate'%in% x$include
    
  } ))))
  
  #covariate in components
  expect_true("fish_covariate(main = covariate, model = \"linear\") + bird_covariate(main = covariate, model = \"linear\") + dog_covariate(main = covariate, model = \"linear\") + insect_covariate(main = covariate, model = \"linear\")" %in% check$.__enclos_env__$private$Components)
  
})

test_that('addBias is able to add bias fields to the model as well as succesfully update the relevent formulas and components of the model.', {
  
  #Check adding bias to the present only dataset
  pcmatern <- inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(2, 0.01),
                                  prior.range = c(1, 0.05))
  check$addBias(datasetNames = 'PO', biasField = pcmatern)
  
  ##Check that bias field is in the formula of the PO dataset
  POdats <- grepl('^PO', names(check$.__enclos_env__$private$modelData))
  expect_true(all(unlist(lapply(check$.__enclos_env__$private$modelData[POdats], function(x) {
    
    'PO_bias_field' %in% x$include
    
  }))))
  expect_false(all(unlist(lapply(check$.__enclos_env__$private$modelData[!POdats], function(x) {
    
    'PO_bias_field' %in% x$include
    
  }))))
  expect_true('PO_bias_field(main = coordinates, model = biasField)' %in% check$.__enclos_env__$private$Components)
  
  expect_identical(check$.__enclos_env__$private$biasField, pcmatern)
  
  #Change bias field
  
})

test_that('changeFormula is able to change the formula of a dataset', {
  
  ##remove the covariate from the PO dataset
  check$changeFormula('PO', formula = ~ ., keepSpatial = TRUE, keepIntercepts = TRUE)
  
  expect_setequal(check$.__enclos_env__$private$modelData$PO_fish_coordinates$include_components, c("shared_spatial", "species_spatial", "fish_intercept" ))
  expect_setequal(check$.__enclos_env__$private$modelData$PA_bird_PAresp$include_components, c("species_spatial", "shared_spatial", "bird_intercept", "pointcov", "covariate"))
  
  #remove covariate for the numvar mark
  check$changeFormula('PO', markName = 'numvar', formula = ~ ., keepSpatial = TRUE, keepIntercepts = TRUE)
  expect_setequal(check$.__enclos_env__$private$modelData$PO_fish_numvar$include_components, c("shared_spatial", "species_spatial", "fish_intercept"))
  
})
