test_that('fitISDM runs a dataSDM object, and produces an INLA model with extra metadata.', {
  skip_on_cran()
  ##Re do this for startISDM + startSpecies + startMarks
  
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
  PA$species <- sample(x = c('bird'), nrow(PA), replace = TRUE)
  mesh <- fmesher::fm_mesh_2d_inla(boundary = fmesher::fm_as_segm(SpatialPoly), 
                             max.edge = 2, crs = fmesher::fm_crs(projection))
  #iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
  iPoints <- fmesher::fm_int(samplers = SpatialPoly, domain = mesh)
  
  responseCounts <- 'count'
  responsePA <- 'PAresp'
  trialName <- 'trial'
  pointCovs <- 'pointcov'
  speciesName <- 'species'
  
  cov <- terra::rast(st_as_sf(SpatialPoly), crs = projection)
  terra::values(cov) <- rgamma(n = terra::ncell(cov), shape = 2)
  names(cov) <- 'covariate'
  cov$cov2 <-  rgamma(n = terra::ncell(cov), shape = 2)
  
  #Start normal ISDM
  obj <- startISDM(PO, PA, Projection = projection, Mesh = mesh,responsePA = responsePA,
                   IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                   spatialCovariates = cov)
  
  ##run model
  spatMod <- fitISDM(data = obj, options  = list(control.inla=list(int.strategy='eb', diagonal = 100)))
  
  expect_setequal(class(spatMod), c('modISDM', 'bru', 'iinla', 'inla'))
  expect_contains(c(row.names(spatMod$summary.fixed)), c("covariate", "cov2", "PO_intercept", "PA_intercept"))
  expect_equal(spatMod[['species']][['speciesVar']], NULL)
  expect_equal(deparse1(spatMod[['componentsJoint']]),"~-1 + PO_spatial(main = geometry, model = PO_field) + PA_spatial(main = geometry, copy = \"PO_spatial\", hyper = list(beta = list(fixed = FALSE))) + covariate(main = covariate, model = \"linear\") + cov2(main = cov2, model = \"linear\") + PO_intercept(1) + PA_intercept(1)")
  
  expect_equal(spatMod[['spatCovs']][['name']], c('covariate', 'cov2'))
  expect_equal(spatMod[['spatCovs']][['class']], c(covariate = 'numeric', cov2 = 'numeric'))
  
  ##Try with custom formula
  obj2 <- startISDM(PO, PA, Projection = projection, Mesh = mesh, responsePA = responsePA,
                   IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                   spatialCovariates = cov, Formulas = list(covariateFormula = ~ covariate + I(covariate^2),
                                                            biasFormula = ~ cov2))
  spatMod2 <- fitISDM(data = obj2, options  = list(control.inla=list(int.strategy='eb', diagonal = 100)))

  expect_true(all(c('Fixed__Effects__Comps', 'Bias__Effects__Comps', 'PO_spatial', 'PA_spatial') %in% names(spatMod2$summary.random)))
  expect_setequal(spatMod2$summary.random$Fixed__Effects__Comps$ID, c('covariate', 'I(covariate^2)'))
  expect_equal(spatMod2$summary.random$Bias__Effects__Comps$ID, 'cov2')

  ##Run with species
  obj3 <- startSpecies(PO, PA, Projection = projection, Mesh = mesh, responsePA = responsePA,
                    IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                    spatialCovariates = cov, speciesName = 'species')

  spatMod3 <- fitISDM(data = obj3,
                       options  = list(control.inla=list(int.strategy='eb', diagonal = 100)))
  
  expect_setequal(spatMod3$species$speciesIn$PO, c('fish'))
  expect_setequal(spatMod3$species$speciesIn$PA, c('bird'))
  expect_equal(spatMod3$species$speciesVar, 'species')
  expect_true(spatMod3$species$speciesEffects$Intercepts)
  expect_true(spatMod3$species$speciesEffects$Environmental)
  
  expect_setequal(names(spatMod3$summary.random), c('PO_spatial', 'PA_spatial', 'speciesShared', 'species_intercepts'))
  expect_setequal(sapply(spatMod3$.args$control.family, function(x) x$link), c('log', 'log', 'cloglog', 'cloglog'))
  
  
  #Try with different arguments
  obj4 <- startSpecies(PO, PA, Projection = projection, Mesh = mesh,responsePA = responsePA,
                       IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                       spatialCovariates = cov, speciesName = 'species', speciesIntercept = FALSE,
                       speciesEnvironment = FALSE)
  
  spatMod4 <- fitISDM(data = obj4,
                      options = list(control.inla=list(int.strategy='eb', diagonal = 100)))
  expect_false(spatMod4$species$speciesEffects$Intercepts)
  expect_false(spatMod4$species$speciesEffects$Environmental)
  expect_setequal(row.names(spatMod4$summary.fixed), c('covariate', 'cov2', 'PO_intercept', 'bird_intercept', 'fish_intercept', 'PA_intercept'))
  
  #Try with formulas
  obj5 <- startSpecies(PO, PA, Projection = projection, Mesh = mesh, responsePA = responsePA,
                       IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                       spatialCovariates = cov, speciesName = 'species',  
                       Formulas = list(covariateFormula = ~ covariate + I(covariate^2),
                                       biasFormula = ~ cov2))
  spatMod5 <- fitISDM(data = obj5,
                      options = list(control.inla=list(int.strategy='eb', diagonal = 100)))
  expect_setequal(names(spatMod5$summary.random), c("PO_spatial", "speciesShared", "Bias__Effects__Comps",
                              "bird_Fixed__Effects__Comps","species_intercepts", "fish_Fixed__Effects__Comps", "PA_spatial"))
  
  ##Try with Marks
  obj6 <- startMarks(PO, Projection = projection, Mesh = mesh, responsePA = responsePA,
                     IPS = iPoints, trialsPA = trialName, responseCounts = responseCounts,
                     markNames = 'numvar', markFamily = 'gaussian')
  
  spatMod6 <- fitISDM(obj6)
  expect_setequal(names(spatMod6$summary.random), c('shared_spatial', 'numvar_spatial', 'PO_numvar_spatial'))
  expect_setequal(row.names(spatMod6$summary.fixed), c('PO_intercept', 'numvar_intercept'))
  
  })
