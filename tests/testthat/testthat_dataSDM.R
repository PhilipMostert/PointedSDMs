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
PO$temp <- sample(x = c(1,2), size = nrow(PO@coords), replace = TRUE)
#Random presence absence dataset
PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
#Add trial name
PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
PA$pointcov <- runif(n = nrow(PA@coords))
PA$binommark <- sample(x = 2:5, size = nrow(PA@data), replace = TRUE)
PA$marktrial <- sample(x = 0:1, size = nrow(PA@data), replace = TRUE)
PA$species <- sample(x = c('bird'), nrow(PA@data), replace = TRUE)
PA$temp <- sample(x = c(1,2), size = nrow(PA@coords), replace = TRUE)
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
speciesName <- 'species'
markSpatial <- TRUE
marksIntercept <- TRUE
speciesSpatial <- TRUE
temporalName <- 'temp'
temporalModel <- deparse(list(model = 'ar1'))

cov <- sp::spsample(x = SpatialPoly, n = 100000, type = 'random')
cov$covariate <- rgamma(n = 100000, shape = 2)
cov <- sp::SpatialPixelsDataFrame(points = cov@coords,
                                  data = data.frame(covariate = cov$covariate),
                                  proj4string = projection,
                                  tolerance = 0.898631)


test_that('dataSDMs initialize works as expected.', {
  
  check <<- dataSDM$new(coordinates = coordnames,
                       projection = projection,
                       Inlamesh = mesh, speciesspatial = speciesSpatial,
                       responsepa = responsePA,
                       trialspa = trialName,
                       responsecounts = responseCounts,
                       marksnames = markNames,
                       trialsmarks = markTrial,
                       marksfamily = marksFamily,
                       pointcovariates = pointCovs,
                       marksspatial = markSpatial,
                       marksintercept = marksIntercept,
                       spatialcovariates = cov,
                       speciesname = speciesName,
                       ips = iPoints,
                       spatial = TRUE, temporal = temporalName, 
                       intercepts = TRUE, temporalmodel = temporalModel)
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
  
  expect_equal(check$.__enclos_env__$private$temporalName, 'temp')
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
               'Mesh needs to be an inla.mesh object.')
  
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
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'log', 'log', 'cloglog', 'cloglog'))
  ##Add a new counts dataset; first try without a species variable (compulsory since speciesName specified in the initialize...)
  Pcount <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  Pcount$count <- rpois(n = nrow(Pcount@coords), lambda = 2)
  expect_error(check$addData(Pcount),'The species variable name is required to be present in all the datasets.')
  
  #Generate species
  Pcount$species <- 'dog'
  Pcount$temp <- sample(c(1,2), nrow(Pcount@data), TRUE)
  check$addData(Pcount)
  ##Check that the data has been added into the model
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO_fish_coordinates", "PO_fish_numvar", "PO_fish_factvar_response", "PA_bird_PAresp", "PA_bird_binommark", "Pcount_dog_count"))
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'log', 'log', 'cloglog', 'cloglog', 'log'))
  
  #Add a new PA dataset with a different response variable name
  PA2 <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA2$newResponse <- sample(0:1, nrow(PA2@coords), replace = TRUE)
  PA2$species <- 'insect'
  PA2$temp <- sample(c(1,2), nrow(PA2@data), TRUE)
  check$addData(PA2, responsePA = 'newResponse')
  
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO_fish_coordinates", "PO_fish_numvar", "PO_fish_factvar_response", "PA_bird_PAresp", "PA_bird_binommark", "Pcount_dog_count", 'PA2_insect_PAresp'))
  expect_setequal(check$.__enclos_env__$private$dataSource, c("PO", "PO", "PO", "PA", "PA", "Pcount", "PA2"))
  ##The correct species_field is in the components, ie the object updated itself correctly
  expect_true('shared_spatial(main = coordinates, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))' %in% check$.__enclos_env__$private$Components)
  expect_equal(unlist(check$.__enclos_env__$private$speciesIn), c(PO = 'fish', PA = 'bird', Pcount = 'dog', PA2 = 'insect'))
  
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
    
    'PO_biasField' %in% x$include
    
  }))))
  expect_false(all(unlist(lapply(check$.__enclos_env__$private$modelData[!POdats], function(x) {
    
    'PO_biasField' %in% x$include
    
  }))))
  expect_true('PO_biasField(main = coordinates, model = PO_bias_field)' %in% check$.__enclos_env__$private$Components)
  
  expect_identical(check$spatialFields$biasFields$PO, pcmatern)
  
  #Change bias field
  
})

test_that('updateFormula is able to change the formula of a dataset', {
  
  #Test error: species/mark/dataset all NULL
  expect_error(check$updateFormula(Formula = ~ covariate), 'At least one of: datasetName, speciesName, markName or allDatasets needs to be specified.')
  
  #Check the error regarding adding a new response variable.
  expect_error(check$updateFormula(datasetName = 'PO', Formula = notResponse ~ covariate), 'Please remove the response variable of the formula.')
  
  ##remove the covariate from the PO dataset
  check$updateFormula(datasetName = 'PO', Formula = ~ . - covariate)
  
  expect_setequal(check$.__enclos_env__$private$modelData$PO_fish_coordinates$include_components, c("fish_spatial", "shared_spatial", "fish_intercept", "PO_biasField"))
  expect_setequal(check$.__enclos_env__$private$modelData$PA_bird_PAresp$include_components, c("bird_covariate", "bird_spatial", "shared_spatial", "bird_intercept", "pointcov"))
  
  #remove covariate for the numvar mark
  check$updateFormula('PO', markName = 'numvar', Formula = ~ . - covariate, keepSpatial = TRUE, keepIntercepts = TRUE)
  expect_setequal(check$.__enclos_env__$private$modelData$PO_fish_numvar$include_components, c("numvar_spatial", "numvar_intercept", "PO_biasField"))
  
  #create completely new formula for PO
  check$updateFormula(datasetName = 'PO', newFormula = ~ . + covariate + I(covariate^2))
  expect_equal(deparse1(check$.__enclos_env__$private$modelData$PO_fish_coordinates$formula), 'coordinates ~ shared_spatial + fish_intercept + fish_spatial + PO_biasField + covariate + I(covariate^2)')

  })

test_that('changeComponents can change the components of the model', {
  
  #remove binmark_spatial from model
  check$changeComponents(removeComponent = 'binommark_spatial')
  expect_false('binommark_spatial(main = coordinates, model = binommark_field)'%in%check$.__enclos_env__$private$Components)
  
  ##Add it back into the components
  check$changeComponents(addComponent = 'binommark_spatial(main = coordinates, model = binommark_field)')
  expect_true('binommark_spatial(main = coordinates, model = binommark_field)'%in%check$.__enclos_env__$private$Components)
  
  #customize component for whatever reason
  check$changeComponents(addComponent = 'binommark_spatial(main = coordinates, model = different_field)')
  expect_false('binommark_spatial(main = coordinates, model = binommark_field)'%in%check$.__enclos_env__$private$Components)
  expect_true('binommark_spatial(main = coordinates, model = different_field)'%in%check$.__enclos_env__$private$Components)
  
  
})

test_that('priorsFixed can add the correct priors to the fixed effects', {
  
  #test errors
   #incorrect effect
  expect_error(check$priorsFixed(effect = 'notcovariate', mean.linear = 200, prec.linear = 20), 'Fixed effect provided not present in the model. Please add covariates using the "spatialCovariates" or "pointCovariates" argument in bruSDM.')
   #species not in model
  expect_error(check$priorsFixed(effect = 'covariate', species = 'monkey'), 'Species given is not available in the model.')
  
  #arbitrary mean and precision for insect_covariate
  check$priorsFixed(effect = 'covariate', species = 'insect', mean.linear = 200, prec.linear = 20)
  expect_true('insect_covariate(main = insect_covariate, model = "linear", mean.linear = 200, prec.linear = 20)' %in% check$.__enclos_env__$private$Components)
  
  
})

test_that('specifySpatial can correctly specify the spatial fields', {
  
  #Check errors:
   #give none of: sharedSpatial, species, mark, bias
  expect_error(check$specifySpatial(prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'At least one of sharedSpatial, dataset, species or mark needs to be provided.')
   #give wrong species name
  expect_error(check$specifySpatial(species = 'elephant', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'Species name provided is not currently in the model.')
   #give wrong mark name
  expect_error(check$specifySpatial(mark = 'height', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'Mark name provided is not currently in the model.')
   #give wrong bias field
  expect_error(check$specifySpatial(bias = 'PA', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)))
  
  #specify the sharedSpatial
  check$specifySpatial(sharedSpatial = TRUE, prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5))
  expect_equal(check$spatialFields$sharedField$sharedField$model, 'pcmatern')
  
  #remove the spatial field from all processes
  check$specifySpatial(sharedSpatial = TRUE, remove = TRUE)
  expect_false('shared_spatial(main = coordinates, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))'%in%check$.__enclos_env__$private$Components)
  expect_false(any(unlist(lapply(check$.__enclos_env__$private$modelData, function(x) 'shared_spatial' %in% x$include_components))))
    
})