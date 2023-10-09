#Make random shape to generate points on
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
PO$temp <- sample(x = c(1,2), size = nrow(PO), replace = TRUE)
#Random presence absence dataset

PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
#Add trial name
PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
PA$pointcov <- runif(n = nrow(PA))
PA$binommark <- sample(x = 0:1, size = nrow(PA), replace = TRUE)
PA$marktrial <- sample(x = 2:5, size = nrow(PA), replace = TRUE)
PA$species <- sample(x = c('bird'), nrow(PA), replace = TRUE)
PA$temp <- sample(x = c(1,2), size = nrow(PA), replace = TRUE)

if (requireNamespace("INLA")) {
mesh <<- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(SpatialPoly), 
                            max.edge = 2, crs = inlabru::fm_crs(projection))
#iPoints <<- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
iPoints <<- inlabru::fm_int(samplers = SpatialPoly, domain = mesh)

}

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
markTrial <- 'marktrial'
pointCovs <- 'pointcov'
speciesName <- 'species'
markSpatial <- TRUE
marksIntercept <- TRUE
speciesSpatial <- TRUE
temporalName <- 'temp'
temporalModel <- deparse(list(model = 'ar1'))
copyModel = deparse1(list(beta = list(fixed = FALSE)))

cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
names(cov) <- 'covariate'


test_that('dataSDMs initialize works as expected.', {
  
  skip_on_cran()

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
                       ips = iPoints, copymodel = copyModel,
                       spatial = 'shared', temporal = temporalName, 
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
  
  expect_equal(check$.__enclos_env__$private$Spatial, 'shared')
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
                           spatial = NULL,
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
                           spatial = NULL,
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
                           spatial = NULL,
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
                           spatial = NULL,
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
                           spatial = NULL,
                           intercepts = FALSE),
               'Projection needs to be a character object.')
  

})

test_that('addData can correctly add and store the relevant metadata properly.', {
  skip_on_cran()
  
  ##Add the datasets
  check$addData(PO,PA)
  
  ##ie also inclusive of marks + species
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO", "PA"))
  expect_true(all(unlist(lapply(check$.__enclos_env__$private$modelData, function(x) class(x))) == 'list'))
  #ie also inclusive of marks
  expect_equal(check$.__enclos_env__$private$dataSource, c(rep('PO',3), rep('PA',2)))
  
  ##Check link functions
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'log', 'log', 'cloglog', 'cloglog'))
  ##Add a new counts dataset; first try without a species variable (compulsory since speciesName specified in the initialize...)
  Pcount <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(Pcount) <- 'geometry'
  Pcount$count <- rpois(n = nrow(Pcount), lambda = 2)
  expect_error(check$addData(Pcount),'The species variable name is required to be present in all the datasets.')
  
  #Generate species
  Pcount$species <- 'dog'
  Pcount$temp <- sample(c(1,2), nrow(Pcount), TRUE)
  check$addData(Pcount)
  ##Check that the data has been added into the model
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO", "PA", "Pcount"))
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'log', 'log', 'cloglog', 'cloglog', 'log'))
  
  #Add a new PA dataset with a different response variable name
  PA2 <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  PA2$newResponse <- sample(0:1, nrow(PA2), replace = TRUE)
  PA2$species <- 'insect'
  PA2$temp <- sample(c(1,2), nrow(PA2), TRUE)
  check$addData(PA2, responsePA = 'newResponse')
  
  expect_setequal(names(check$.__enclos_env__$private$modelData), c("PO", "PA", "Pcount", "PA2"))
  expect_setequal(check$.__enclos_env__$private$dataSource, c("PO", "PO", "PO", "PA", "PA", "Pcount", "PA2"))
  ##The correct species_field is in the components, ie the object updated itself correctly
  expect_true('shared_spatial(main = geometry, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))' %in% check$.__enclos_env__$private$Components)
  expect_equal(unlist(check$.__enclos_env__$private$speciesIn), c(PO = 'fish', PA = 'bird', Pcount = 'dog', PA2 = 'insect'))
  
  })

test_that('addBias is able to add bias fields to the model as well as succesfully update the relevant formulas and components of the model.', {
  skip_on_cran()
  
  #Check adding bias to the present only dataset
  pcmatern <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(2, 0.01),
                                  prior.range = c(1, 0.05))
  check$addBias(datasetNames = 'PO', biasField = pcmatern)
  
  ##Check that bias field is in the formula of the PO dataset
  expect_true(unlist(lapply(check$.__enclos_env__$private$Formulas$PO$fish, function(x) {
    
    'PO_biasField' %in% x$RHS
    
  }))[1])

  expect_true("PO_biasField(main = geometry, model = PO_bias_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))"  %in% check$.__enclos_env__$private$Components)
  
})

test_that('updateFormula is able to change the formula of a dataset', {
  skip_on_cran()
  
  #Test error: species/mark/dataset all NULL
  expect_error(check$updateFormula(Formula = ~ covariate), 'At least one of: datasetName, speciesName, markName or allProcesses needs to be specified.')
  
  #Check the error regarding adding a new response variable.
  expect_error(check$updateFormula(datasetName = 'PO', Formula = notResponse ~ covariate), 'Please remove the response variable of the formula.')
  
  ##remove the covariate from the PO dataset
  check$updateFormula(datasetName = 'PO', Formula = ~ . - covariate)
  
  expect_setequal(check$.__enclos_env__$private$Formulas$PO$fish$geometry$RHS, c("fish_PO_spatial", "shared_spatial", "fish_intercept", "PO_biasField"))
  expect_setequal(check$.__enclos_env__$private$Formulas$PA$bird$PAresp$RHS, c("bird_covariate", "bird_PA_spatial", "shared_spatial", "bird_intercept", "pointcov"))
  
  })

test_that('changeComponents can change the components of the model', {
  skip_on_cran()
  
  #remove binmark_spatial from model
  check$changeComponents(removeComponent = 'binommark_spatial')
  expect_false('binommark_spatial(main = geometry, model = binommark_field)'%in%check$.__enclos_env__$private$Components)
  
  ##Add it back into the components
  check$changeComponents(addComponent = 'binommark_spatial(main = geometry, model = binommark_field)')
  expect_true('binommark_spatial(main = geometry, model = binommark_field)'%in%check$.__enclos_env__$private$Components)
  
  #customize component for whatever reason
  check$changeComponents(addComponent = 'binommark_spatial(main = geometry, model = different_field)')
  expect_false('binommark_spatial(main = geometry, model = binommark_field)'%in%check$.__enclos_env__$private$Components)
  expect_true('binommark_spatial(main = geometry, model = different_field)'%in%check$.__enclos_env__$private$Components)
  
  
})

test_that('priorsFixed can add the correct priors to the fixed effects', {
  skip_on_cran()
  
  #test errors
   #incorrect effect
  expect_error(check$priorsFixed(Effect = 'notcovariate', mean.linear = 200, prec.linear = 20), 'Fixed effect provided not present in the model. Please add covariates using the "spatialCovariates" or "pointCovariates" argument in intModel')
   #species not in model
  expect_error(check$priorsFixed(Effect = 'covariate', Species = 'monkey'), 'Species given is not available in the model.')
  
  #arbitrary mean and precision for insect_covariate
  check$priorsFixed(Effect = 'covariate', Species = 'insect', mean.linear = 200, prec.linear = 20)
  expect_true('insect_covariate(main = insect_covariate, model = "linear", mean.linear = 200, prec.linear = 20)' %in% check$.__enclos_env__$private$Components)
  
  
})

test_that('specifySpatial can correctly specify the spatial fields', {
  skip_on_cran()
  
  #Check errors:
   #give none of: sharedSpatial, species, mark, bias
  expect_error(check$specifySpatial(prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'At least one of sharedSpatial, datasetName, dataset, Species or Mark needs to be provided.')
   #give wrong species name
  expect_error(check$specifySpatial(Species = 'elephant', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'Species name provided is not currently in the model.')
   #give wrong mark name
  expect_error(check$specifySpatial(Mark = 'height', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'Mark name provided is not currently in the model.')
   #give wrong bias field
  expect_error(check$specifySpatial(bias = 'PA', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)))
  
  #specify the sharedSpatial
  check$specifySpatial(sharedSpatial = TRUE, prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5))
  expect_equal(check$spatialFields$sharedField$sharedField$model, 'pcmatern')
  
  #remove the spatial field from all processes
  check$specifySpatial(sharedSpatial = TRUE, Remove = TRUE)
  expect_false('shared_spatial(main = geometry, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))'%in%check$.__enclos_env__$private$Components)
  expect_false(any(unlist(lapply(check$.__enclos_env__$private$modelData, function(x) 'shared_spatial' %in% x$include_components))))
    
})

