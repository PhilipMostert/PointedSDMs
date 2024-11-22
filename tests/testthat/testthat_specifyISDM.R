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
PO$temp <- sample(x = c(1,2), size = nrow(PO), replace = TRUE)
PO$otherVar <- rgamma(n = nrow(PO), shape = 12)
#Random presence absence dataset

PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
#Add trial name
PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
PA$pointcov <- runif(n = nrow(PA))
PA$temp <- sample(x = c(1,2), size = nrow(PA), replace = TRUE)

Pcount <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
st_geometry(Pcount) <- 'geometry'
Pcount$count <- rpois(n = nrow(Pcount), lambda = 2)
Pcount$temp <- 2

if (requireNamespace("INLA")) {
  mesh <<- fmesher::fm_mesh_2d_inla(boundary = fmesher::fm_as_segm(SpatialPoly), 
                              max.edge = 2, crs = fmesher::fm_crs(projection))
  #iPoints <<- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  iPoints <<- fmesher::fm_int(samplers = SpatialPoly, domain = mesh)
  
}

responseCounts <- 'count'
responsePA <- 'PAresp'
trialName <- 'trial'
pointCovs <- 'pointcov'

temporalName <- 'temp'
temporalModel <- deparse(list(model = 'ar1'))
copyModel = deparse1(list(beta = list(fixed = FALSE)))

cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
names(cov) <- 'covariate'


test_that('specifyISDM initialize works as expected.', {
  
  skip_on_cran()
  
  check <<- specifyISDM$new(data = list(PO, PA, Pcount),
                        initialnames = c('PO', 'PA', 'Pcount'),
                        projection = projection,
                        Inlamesh = mesh,
                        responsepa = responsePA,
                        trialspa = trialName,
                        responsecounts = responseCounts,
                        pointcovariates = pointCovs,
                        spatialcovariates = cov,
                        formulas = NULL,
                        offset = NULL,
                        ips = iPoints, copymodel = copyModel,
                        spatial = 'shared', temporal = temporalName, 
                        intercepts = TRUE, temporalmodel = temporalModel)
  
  ##Test private classes
  expect_true(identical(check$.__enclos_env__$private$Projection,projection))
  expect_true(identical(check$.__enclos_env__$private$INLAmesh, mesh))
  
  expect_true(check$.__enclos_env__$private$responsePA == responsePA)
  expect_true(check$.__enclos_env__$private$responseCounts == responseCounts)
  expect_true(check$.__enclos_env__$private$trialsPA == trialName)
  
  expect_equal(check$.__enclos_env__$private$Spatial, 'shared')
  expect_true(check$.__enclos_env__$private$Intercepts)
  
  expect_equal(check$.__enclos_env__$private$temporalName, 'temp')
  expect_equal(check$.__enclos_env__$private$ptcovsClass, 'numeric', ignore_attr = TRUE)
  expect_equal(check$.__enclos_env__$private$spatcovsNames, 'covariate')
  expect_setequal(unlist(check$.__enclos_env__$private$optionsINLA), c('log', 'cloglog', 'log'))
  expect_true(inherits(check$spatialFields$sharedField$sharedField, 'inla.spde2'))
  
  expect_true('pointcov' %in% names(check$.__enclos_env__$private$IPS))
  expect_true(all(is.na(check$.__enclos_env__$private$IPS$pointcov)))
  expect_true('covariate' %in% names(check$.__enclos_env__$private$IPS))
  expect_true(inherits(check$.__enclos_env__$private$IPS$covariate, 'numeric'))
  
  checkNocovs <<- specifyISDM$new(data = list(PO, PA, Pcount),
                            initialnames = c('PO', 'PA', 'Pcount'),
                            projection = projection,
                            Inlamesh = mesh,
                            responsepa = responsePA,
                            trialspa = trialName,
                            responsecounts = responseCounts,
                            pointcovariates = NULL,
                            spatialcovariates = NULL,
                            formulas = NULL,
                            offset = NULL,
                            ips = iPoints, copymodel = copyModel,
                            spatial = 'shared', temporal = temporalName, 
                            intercepts = TRUE, temporalmodel = temporalModel)
  
  expect_false('pointcov' %in% names(checkNocovs$.__enclos_env__$private$IPS))
  expect_false('covariate' %in% names(checkNocovs$.__enclos_env__$private$IPS))
  
  #Remove one of the response variable names
  expect_error(specifyISDM$new(data = list(PO, PA, Pcount),
                               initialnames = c('PO', 'PA', 'Pcount'),
                               projection = projection,
                               Inlamesh = mesh,
                               responsepa = NULL,
                               trialspa = trialName,
                               responsecounts = responseCounts,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = NULL,
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'shared', temporal = temporalName, 
                               intercepts = TRUE, temporalmodel = temporalModel),
               'At least one of responseCounts and responsePA are NULL. Please provide both.')
  expect_error(specifyISDM$new(data = list(PO, PA, Pcount),
                               initialnames = c('PO', 'PA', 'Pcount'),
                               projection = projection,
                               Inlamesh = mesh,
                               responsepa = responsePA,
                               trialspa = trialName,
                               responsecounts = NULL,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = NULL,
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'shared', temporal = temporalName, 
                               intercepts = TRUE, temporalmodel = temporalModel),
               'At least one of responseCounts and responsePA are NULL. Please provide both.')
  
  #Test mesh error
  expect_error(specifyISDM$new(data = list(PO, PA, Pcount),
                               initialnames = c('PO', 'PA', 'Pcount'),
                               projection = projection,
                               Inlamesh = list(),
                               responsepa = NULL,
                               trialspa = trialName,
                               responsecounts = responseCounts,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = NULL,
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'shared', temporal = temporalName, 
                               intercepts = TRUE, temporalmodel = temporalModel),
               'Mesh needs to be an fm_mesh_2d object.')
  
  ##Check projection error
  expect_error(specifyISDM$new(data = list(PO, PA, Pcount),
                               initialnames = c('PO', 'PA', 'Pcount'),
                               projection = 84,
                               Inlamesh = mesh,
                               responsepa = NULL,
                               trialspa = trialName,
                               responsecounts = responseCounts,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = NULL,
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'shared', temporal = temporalName, 
                               intercepts = TRUE, temporalmodel = temporalModel),
               'Projection needs to be a character object.')
  
  ##Test changing PointSpatial
  checkInd <<- specifyISDM$new(data = list(PO, PA, Pcount),
                            initialnames = c('PO', 'PA', 'Pcount'),
                            projection = projection,
                            Inlamesh = mesh,
                            responsepa = responsePA,
                            trialspa = trialName,
                            responsecounts = responseCounts,
                            pointcovariates = pointCovs,
                            spatialcovariates = cov,
                            formulas = NULL,
                            offset = NULL,
                            ips = iPoints, copymodel = copyModel,
                            spatial = 'individual', temporal = temporalName, 
                            intercepts = TRUE, temporalmodel = temporalModel)
  
  expect_true(length(checkInd$spatialFields$sharedField) == 0)
  expect_setequal(names(checkInd$spatialFields$datasetFields), c('PO', 'PA', 'Pcount'))
  expect_true(all(unlist(lapply(checkInd$spatialFields$datasetFields, function(x) inherits(x, 'inla.spde2')))))
  
  checkCor <<- specifyISDM$new(data = list(PO, PA, Pcount),
                            initialnames = c('PO', 'PA', 'Pcount'),
                            projection = projection,
                            Inlamesh = mesh,
                            responsepa = responsePA,
                            trialspa = trialName,
                            responsecounts = responseCounts,
                            pointcovariates = pointCovs,
                            spatialcovariates = cov,
                            formulas = NULL,
                            offset = NULL,
                            ips = iPoints, copymodel = copyModel,
                            spatial = 'correlate', temporal = temporalName, 
                            intercepts = TRUE, temporalmodel = temporalModel)
  
  expect_true(length(checkCor$spatialFields$sharedField) == 1)
  expect_true(inherits(checkCor$spatialFields$sharedField$sharedField, 'inla.spde2'))
  expect_true('._dataset_index_var_.' %in% names(checkCor$.__enclos_env__$private$IPS))
  
  checkCop <<- specifyISDM$new(data = list(PO, PA, Pcount),
                            initialnames = c('PO', 'PA', 'Pcount'),
                            projection = projection,
                            Inlamesh = mesh,
                            responsepa = responsePA,
                            trialspa = trialName,
                            responsecounts = responseCounts,
                            pointcovariates = pointCovs,
                            spatialcovariates = cov,
                            formulas = NULL,
                            offset = NULL,
                            ips = iPoints, copymodel = copyModel,
                            spatial = 'copy', temporal = temporalName, 
                            intercepts = TRUE, temporalmodel = temporalModel)
   
  expect_true(length(checkCop$spatialFields$datasetFields) == 1)
  expect_true(inherits(checkCop$spatialFields$datasetFields$PO, 'inla.spde2'))
  
  ##Check offset
  PO$offset <- 0.5
  expect_error(specifyISDM$new(data = list(PO, PA, Pcount),
                  initialnames = c('PO', 'PA', 'Pcount'),
                  projection = projection,
                  Inlamesh = mesh,
                  responsepa = responsePA,
                  trialspa = trialName,
                  responsecounts = responseCounts,
                  pointcovariates = pointCovs,
                  spatialcovariates = cov,
                  formulas = NULL,
                  offset = 'offset',
                  ips = iPoints, copymodel = copyModel,
                  spatial = 'shared', temporal = temporalName, 
                  intercepts = TRUE, temporalmodel = temporalModel), 'The offset variable cannot be called "offset". Please choose a new name.')
  
  PO$offsetVar <- 0.5
  checkOff <<- specifyISDM$new(data = list(PO, PA, Pcount),
                               initialnames = c('PO', 'PA', 'Pcount'),
                               projection = projection,
                               Inlamesh = mesh,
                               responsepa = responsePA,
                               trialspa = trialName,
                               responsecounts = responseCounts,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = 'offsetVar',
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'shared', temporal = temporalName, 
                               intercepts = TRUE, temporalmodel = temporalModel)
  
  PO2 <- PO
  checknoSpat <<- specifyISDM$new(data = list(PO, PA, Pcount, PO2),
                               initialnames = c('PO', 'PA', 'Pcount', 'PO2'),
                               projection = projection,
                               Inlamesh = mesh,
                               boundary = SpatialPoly,
                               responsepa = responsePA,
                               trialspa = trialName,
                               responsecounts = responseCounts,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = 'offsetVar',
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'shared', temporal = NULL, 
                               intercepts = TRUE, temporalmodel = temporalModel)
  
})

test_that('addBias is able to add bias fields to the model as well as succesfully update the relevant formulas and components of the model.', {
  
  skip_on_cran()
  
  #Check adding bias to the present only dataset
  
  if (requireNamespace("INLA")) {
  
  pcmatern <- INLA::inla.spde2.pcmatern(mesh,
                                        prior.sigma = c(2, 0.01),
                                        prior.range = c(1, 0.05))
  
  }
  check$addBias(datasetNames = 'PO', biasField = pcmatern, copyModel = FALSE)
  
  expect_true("PO_biasField(main = geometry, model = PO_bias_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))"  %in% check$.__enclos_env__$private$Components)
  expect_equal(names(check$spatialFields$biasFields), 'PO')
  expect_true(inherits(check$spatialFields$biasFields$PO, 'inla.spde2'))
  
  ##Check that bias field is in the formula of the PO dataset
  expect_true(unlist(lapply(check$.__enclos_env__$private$Formulas$PO$PO, function(x) {
    
    'PO_biasField' %in% x$RHS
    
  }))[1])
  
  check$addBias(datasetNames = 'PO', biasField = pcmatern, temporalModel = list(model = 'rw1'), copyModel = FALSE)
  expect_true("PO_biasField(main = geometry, model = PO_bias_field, group = temp, ngroup = 2, control.group = list(model = \"rw1\"))" %in% check$.__enclos_env__$private$Components)
  
  #ADD bias for all PO
  checknoSpat$addBias(allPO = TRUE, copyModel = FALSE, shareModel = FALSE)
  expect_true(all(c("PO_biasField(main = geometry, model = PO_bias_field)", "PO2_biasField(main = geometry, model = PO2_bias_field)")%in% checknoSpat$.__enclos_env__$private$Components))
  expect_setequal(names(checknoSpat$spatialFields$biasFields), c('PO', 'PO2'))
  expect_true(all(unlist(lapply(checknoSpat$spatialFields$biasFields, function(x) inherits(x, 'inla.spde2')))))
  checknoSpat$spatialFields$biasFields <- list()
  #copy
  checknoSpat$addBias(allPO = TRUE, copyModel = TRUE)
  expect_true(all(c("PO_biasField(main = geometry, model = PO_bias_field)", "PO2_biasField(main = geometry, copy = \"PO_biasField\", hyper = list(beta = list(fixed = FALSE)))" )%in% checknoSpat$.__enclos_env__$private$Components))
  expect_setequal(names(checknoSpat$spatialFields$biasFields), c('PO', 'PO2'))
  expect_true(all(unlist(lapply(checknoSpat$spatialFields$biasFields, function(x) inherits(x, 'inla.spde2')))))
  checknoSpat$spatialFields$biasFields <- list()
  #shared
  expect_error(checknoSpat$addBias(allPO = TRUE, shareModel = TRUE), 'Only one of copyModel and shareModel may be TRUE.')
  checknoSpat$addBias(allPO = TRUE, shareModel = TRUE, copyModel = FALSE)
  expect_true("sharedBias_biasField(main = geometry, model = sharedBias_bias_field)" %in% checknoSpat$.__enclos_env__$private$Components)
  expect_true(names(checknoSpat$spatialFields$biasFields) == 'sharedBias')
  expect_true(inherits(checknoSpat$spatialFields$biasFields$sharedBias, 'inla.spde2'))
  

})

test_that('updateFormula is able to change the formula of a dataset', {
  
  skip_on_cran()
  
  #Test error: species/mark/dataset all NULL
  expect_error(check$updateFormula(Formula = ~ covariate), 'Please specify the name of the dataset for which you want to change the formula using datasetName.')
  
  #Check the error regarding adding a new response variable.
  expect_error(check$updateFormula(datasetName = 'PO', Formula = notResponse ~ covariate), 'Please remove the response variable of the formula.')
  
  #Check printing formula
  expect_true(inherits(check$updateFormula(datasetName = 'PO'), 'list'))
  
  ##remove the covariate from the PO dataset
  check$updateFormula(datasetName = 'PO', Formula = ~ . - covariate)
  
  expect_setequal(check$.__enclos_env__$private$Formulas$PO$PO$geometry$RHS, c('PO_intercept', "shared_spatial", "PO_biasField"))
  expect_setequal(check$.__enclos_env__$private$Formulas$PA$PA$PAresp$RHS, c('PA_intercept', "shared_spatial", "pointcov", 'covariate'))
  
  ##Completely change the formula 
  check$updateFormula(datasetName = 'PA', newFormula = ~ exp(PA_intercept + shared_spatial + pointcov + covariate + I(covariate^2)))
  expect_equal(check$.__enclos_env__$private$Formulas$PA$PA$PAresp$LHS, PAresp ~ exp(PA_intercept + shared_spatial + pointcov + covariate + 
                                                                                            I(covariate^2)), ignore_attr = TRUE)
  expect_true(is.null(check$.__enclos_env__$private$Formulas$PA$PA$PAresp$RHS))
  
  ##Check covariateFormula
  check <<- specifyISDM$new(data = list(PO, PA, Pcount),
                            initialnames = c('PO', 'PA', 'Pcount'),
                            projection = projection,
                            Inlamesh = mesh,
                            responsepa = responsePA,
                            trialspa = trialName,
                            responsecounts = responseCounts,
                            pointcovariates = pointCovs,
                            spatialcovariates = cov,
                            formulas = list(covariateFormula = ~ covariate + I(covariate^2)),
                            offset = NULL,
                            ips = iPoints, copymodel = copyModel,
                            spatial = 'shared', temporal = temporalName, 
                            intercepts = TRUE, temporalmodel = temporalModel)
  check$updateFormula(processLevel = TRUE, Formula = ~ . - I(covariate^2))
  expect_equal(deparse1(check$.__enclos_env__$private$covariateFormula), '~covariate')
  expect_true("Fixed__Effects__Comps(main = ~covariate - 1, model = \"fixed\")" %in% check$.__enclos_env__$private$Components)
  
  check$updateFormula(newFormula = ~exp(covariate), processLevel = TRUE)
  expect_equal(deparse1(check$.__enclos_env__$private$covariateFormula), '~exp(covariate)')
  expect_true("Fixed__Effects__Comps(main = ~exp(covariate) - 1, model = \"fixed\")" %in% check$.__enclos_env__$private$Components)
  
})

test_that('changeComponents can change the components of the model', {
  
  skip_on_cran()
  
  #remove binmark_spatial from model
  check$changeComponents(removeComponent = 'shared_spatial')
  expect_false('shared_spatial(main = geometry, model = shared_field, group = temp, 
    ngroup = 2, control.group = list(model = "ar1"))'%in%check$.__enclos_env__$private$Components)
  
  ##Add it back into the components
  check$changeComponents(addComponent = 'shared_spatial(main = geometry, model = shared_field, group = temp, 
    ngroup = 2, control.group = list(model = "ar1"))')
  expect_true('shared_spatial(main = geometry, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))'%in%check$.__enclos_env__$private$Components)
  
  #customize component for whatever reason
   #This won't work becuase different_field needs to be in the fitISDM environment?
  checknoSpat$changeComponents(addComponent = 'PO_biasField(main = geometry, model = different_field)')
  expect_false('PO_biasField(main = geometry, model = PO_bias_field)'%in%checknoSpat$.__enclos_env__$private$Components)
  expect_true('PO_biasField(main = geometry, model = different_field)'%in% checknoSpat$.__enclos_env__$private$Components)
  
  expect_output(check$changeComponents(), 'Model components:')
  
})

test_that('priorsFixed can add the correct priors to the fixed effects', {
  
  skip_on_cran()
  
  #incorrect effect
  expect_error(check$priorsFixed(Effect = 'notcovariate', mean.linear = 200, prec.linear = 20), 'Fixed effect provided not present in the model. Please add covariates using the "spatialCovariates" or "pointCovariates" argument in intModel')
  
  #arbitrary mean and precision for insect_covariate
  check$priorsFixed(Effect = 'covariate', mean.linear = 200, prec.linear = 20)
  expect_true("covariate(main = covariate, model = \"linear\", mean.linear = 200, prec.linear = 20)" %in% check$.__enclos_env__$private$Components)
  
  ##Incorrect datasetName
  expect_error(check$priorsFixed(Effect = 'Intercept', mean.linear = 1, prec.linear = 1, datasetName = 'PCounts'), 'datasetName is not the name of a dataset added to the model.')
  
  
  check$priorsFixed(Effect = 'Intercept', mean.linear = 1, prec.linear = 1, datasetName = 'PA')
  expect_true("PA_intercept(1, mean.linear = 1, prec.linear = 1)" %in% check$.__enclos_env__$private$Components)
  
})

test_that('specifySpatial can correctly specify the spatial fields', {
  
  skip_on_cran()
  
  #Check errors:
  #give none of: sharedSpatial, species, mark, bias
  expect_error(check$specifySpatial(prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)), 'At least one of sharedSpatial, datasetName or Bias needs to be provided.')
  #give wrong bias field
  expect_error(check$specifySpatial(Bias = 'PA', prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5)))
  
  #specify the sharedSpatial
  check$specifySpatial(sharedSpatial = TRUE, prior.range = c(1,0.1), prior.sigma = c(0.2, 0.5))
  expect_equal(check$spatialFields$sharedField$sharedField$model, 'pcmatern')
  
  #remove the spatial field from all processes
  check$specifySpatial(sharedSpatial = TRUE, Remove = TRUE)
  expect_false('shared_spatial(main = geometry, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))'%in%check$.__enclos_env__$private$Components)
  expect_false(any(unlist(lapply(check$.__enclos_env__$private$modelData, function(x) 'shared_spatial' %in% x$include_components))))
  
})

test_that('changeLink can correctly change the link function of a process', {
  
  skip_on_cran()
  
  expect_error(checknoSpat$changeLink(datasetName = 'PO2'))
  
  checknoSpat$changeLink(datasetName = 'PO2', Link = 'exp')

  expect_true(checknoSpat$.__enclos_env__$private$optionsINLA$control.family[[4]]$link == 'exp')
  
  checknoSpat$changeLink(datasetName = 'PA', Link = 'logit')
  
  expect_true(checknoSpat$.__enclos_env__$private$optionsINLA$control.family[[2]]$link == 'logit')
  
  })

test_that('spatialBlock can correctly block the spatial domain', {
  
  skip_on_cran()
  
  checknoSpat$spatialBlock(k = 2, rows_cols = c(1,2), plot = FALSE)
  
  expect_true(all(sapply(unlist(checknoSpat$.__enclos_env__$private$modelData, recursive = FALSE), function(x) '.__block_index__' %in% names(x))))
  
  expect_true(all(sapply(unlist(checknoSpat$.__enclos_env__$private$modelData, recursive = FALSE), function(x) x$ '.__block_index__' %in% 1:2)))
  
})

test_that('addSamplers can correctly add new samplers to a process', {
  skip_on_cran()
  
  wrongSamp <- 1:2
  
  expect_error(check$addSamplers(datasetName = 'PO', Samplers = wrongSamp), 'Samplers needs to be a sf object.')
  
  #Add random (wrong polygon)
  
  x <- c(16,  17,  24, 22, 16)
  y <- c(59, 55, 55.0341796875, 61, 59)
  xy <- cbind(x, y)
  SpatialPoly2 <- st_as_sf(st_sfc(st_polygon(list(xy)), crs = projection))
  
  check$addSamplers(datasetName = 'PO', Samplers = SpatialPoly2)
  expect_equal(check$.__enclos_env__$private$Samplers$PO, SpatialPoly2, ignore_attr = TRUE)
  
})

test_that('specifyRandom can correctly change the priors of the random effects', {
  
  skip_on_cran()
  #new again
  check <<- specifyISDM$new(data = list(PO, PA, Pcount),
                            initialnames = c('PO', 'PA', 'Pcount'),
                            projection = projection,
                            Inlamesh = mesh,
                            responsepa = responsePA,
                            trialspa = trialName,
                            responsecounts = responseCounts,
                            pointcovariates = pointCovs,
                            spatialcovariates = cov,
                            formulas = NULL,
                            offset = NULL,
                            ips = iPoints, copymodel = copyModel,
                            spatial = 'shared', temporal = temporalName, 
                            intercepts = TRUE, temporalmodel = temporalModel)
  
  
  check$specifyRandom(temporalModel = list(model = 'rw1'))
  expect_true("shared_spatial(main = geometry, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"rw1\"))" %in% 
                check$.__enclos_env__$private$Components)
  
  checkCop <<- specifyISDM$new(data = list(PO, PA, Pcount),
                               initialnames = c('PO', 'PA', 'Pcount'),
                               projection = projection,
                               Inlamesh = mesh,
                               responsepa = responsePA,
                               trialspa = trialName,
                               responsecounts = responseCounts,
                               pointcovariates = pointCovs,
                               spatialcovariates = cov,
                               formulas = NULL,
                               offset = NULL,
                               ips = iPoints, copymodel = copyModel,
                               spatial = 'copy', temporal = NULL, 
                               intercepts = TRUE, temporalmodel = temporalModel)
  
  checkCop$addBias(datasetNames = c('PA', 'PO', 'Pcount'))
  
  
  checkCop$specifyRandom(copyModel = list(beta = list(fixed = TRUE)))
  expect_true(all(c("PA_spatial(main = geometry, copy = \"PO_spatial\", hyper = list(beta = list(fixed = TRUE)))", "Pcount_spatial(main = geometry, copy = \"PO_spatial\", hyper = list(beta = list(fixed = TRUE)))") %in% checkCop$.__enclos_env__$private$Components))
  
  ##Add bias to all
  checkCop$specifyRandom(copyBias = list(hyper = list(beta = list(initial = 2))))
  expect_true(all(c("PO_biasField(main = geometry, copy = \"PA_biasField\", hyper = list(hyper = list(beta = list(initial = 2))))",
                    "Pcount_biasField(main = geometry, copy = \"PA_biasField\", hyper = list(hyper = list(beta = list(initial = 2))))") %in% checkCop$.__enclos_env__$private$Components))
  
})

