##Load in the necessary data
Check <- dataOrganize$new() 

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
PA$marktrial <- sample(x = 2:5, size = nrow(PA), replace = TRUE)
PA$species <- sample(x = c('bird1', 'bird2'), nrow(PA), replace = TRUE)
PA$temp <- sample(x = c(1,2), nrow(PA), replace = TRUE)
mesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(SpatialPoly), 
                           max.edge = 2, crs = inlabru::fm_crs(projection))
#iPoints <- inlabru::ipoints(samplers = SpatialPoly, domain = mesh)
iPoints <- inlabru::fm_int(samplers = SpatialPoly, domain = mesh)
##Make PA a data.frame object
PA$long <- st_coordinates(PA)[,1]
PA$lat <- st_coordinates(PA)[,2]
st_geometry(PA) <- NULL
PA <- data.frame(PA)

spData <- list(PO, PA)

test_that('The internal function makeData returns a list of SpatialPointDataFrame objects as
          well as the relevant metadata to be used in the integrated model.', {
    
    Check$makeData(datapoints = spData, datanames = c('PO', 'PA'),
                   coords = c('long', 'lat'), proj = projection, offsetname = NULL,
                   pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                   speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'), temporalvar = 'temp',
                   marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial', 'binomial'))
            
    expect_setequal(names(Check$Data), c('PO','PA'))
            
    expect_true(all(unlist(sapply(unlist(Check$Data, recursive = FALSE), function(x) inherits(x, 'sf')))))
            
    ##Should create a placeholder variable for the poresp + 
    #should keep marks +
    #should create new variables for the multinomial marks.
    expect_setequal(names(Check$Data$PO[[1]]), c("poresp", "numvar", "factvar", 'temp', 'geometry',
                                                    "species", "factvar_phi", "factvar_response"))
    expect_true((all(Check$Data$PO[[1]]$factvar_phi == 1)))
    expect_true((all(Check$Data$PO[[1]]$factvar_response == 1)))
    expect_true(class(Check$Data$PO[[1]]$factvar) == 'character')
    
    expect_setequal(names(Check$Data$PA[[1]]), c("PAresp", "trial", "binommark", 'temp', 'geometry',
                                                    "marktrial", "species", "pointcov"))
    #Family for PO should be:
    # cp for the points;
    # uniform for the mark;
    # poisson for the multinomial mark.
    expect_setequal(Check$Family$PO, c('cp', 'uniform', 'poisson'))
    
    expect_setequal(Check$Family$PA, c('binomial', 'binomial'))
    
    expect_named(Check$dataType, c('PO', 'PA'))
    expect_setequal(Check$dataType, c('Present only', 'Present absence'))
    
    ##PO has no point covariates;
    ##PA has pointcov as a pointcovariate
    expect_true(is.null(unlist(Check$varsIn$PO)))
    expect_true(unlist(Check$varsIn$PA) == 'pointcov')
    
    expect_setequal(Check$Marks$PO, c('numvar', 'factvar'))
    expect_setequal(Check$Marks$PA, c('binommark'))
    
    expect_named(Check$marksType$PO, c('numvar', 'factvar'))
    expect_named(Check$marksType$PA, c('binommark'))
    
    expect_setequal(Check$marksType$PO, c('Uniform mark', 'Multinomial mark'))
    expect_setequal(Check$marksType$PA, c('Binomial mark'))
    
    expect_true('factvar' %in% Check$multinomVars)
    
    expect_true(Check$numObs[1] == nrow(PO))
    expect_true(Check$numObs[2] == nrow(PA))
    
    #Ie there are three processes in PO: the points + 2 marks
    expect_length(Check$dataSource[[1]], 3)
    
    #Ie there are two processes in PO: the points + 21marks
    expect_length(Check$dataSource[[2]],2)
    
    #Remove a dataset name
    expect_error(Check$makeData(datapoints = spData, datanames = c('PO'),
                                coords = colnames(PO@coords), proj = projection,
                                pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                                speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'),
                                marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial', 'binomial')),
                                'Number of dataset names needs to equal length of datasets.')
    
    #Remove a mark family
    expect_error(Check$makeData(datapoints = spData, datanames = c('PO','PA'),
                                coords = colnames(PO@coords), proj = projection,
                                pointcovnames = 'pointcov', paresp = 'PAresp', countsresp = 'counts', trialname = 'trial',
                                speciesname = 'species', marks = c('numvar', 'factvar', 'binommark'),
                                marktrialname = 'marktrial', markfamily = c('uniform', 'multinomial')),
                                "Number of marks needs to equal the number of mark families.")
          })

test_that('makeSpecies is able to split the data up by species.', {
  
  Check$makeSpecies(speciesname = 'species')
  
  expect_setequal(names(Check$Data$PO), c('PO_fish1', 'PO_fish2'))
  expect_setequal(names(Check$Data$PA), c('PA_bird1', 'PA_bird2'))
  
  expect_setequal(Check$SpeciesInData$PO, c('fish1', 'fish2'))
  expect_setequal(Check$SpeciesInData$PA, c('bird1', 'bird2'))
  
  #When converting factors to numeric, we expect them to be ordered alphabetically
  ## So bird1 and bird2 should get 1 and 2; fish1 and fish2 should get 3 and 4
  expect_true(all(unlist(Check$speciesNumeric$species$PO)%in%c(3,4)))
  expect_true(all(unlist(Check$speciesNumeric$species$PA)%in%c(1,2)))
  
  #ie multiply the length of process by #species = 2
  expect_length(Check$dataSource[[1]], 3 * 2)
  expect_length(Check$dataSource[[2]], 2 * 2)
  
  
})

test_that('makeMultinom is able to organize and create usable multinomial data for INLA', {
  
  #This function was implicitly checked with makeSpecies?
  
  Check$makeMultinom(multinomVars = Check$multinomVars, return = 'marks', oldVars = NULL)
  
  expect_setequal(unlist(Check$multinomIndex$factvar$PO),c('a','b'))
  #No factor var present in PA, so should expect NA
  expect_true(all(is.na(unlist(Check$multinomIndex$factvar$PA))))
  
  expect_true(all(as.numeric(factor(unlist(Check$multinomIndex$factvar$PO))) == unlist(Check$multinomNumeric$factvar$PO)))
  
  #The factor variable should now be numeric in the data; the index is stored in multinoIndex
  expect_true(all(sapply(Check$Data$PO, function(x) class(x$factvar) == 'numeric')))
  
  })

test_that('makeFormulas is able to make the correct formulas for the different processes
          based on their response variables, and the available covariates.', {
            
            
            #Spatcovs is the names of the spatial covariates
            #specnesname is the name of the species variable
            Check$makeFormulas(spatcovs = 'spatcovs', speciesname = 'species', markintercept = TRUE,
                               paresp = 'PAresp', countresp = 'counts', marksspatial = TRUE, speciesspatial = 'individual',
                               marks = c('numvar', 'factvar', 'binommark'), temporalname = 'temp',
                               spatial = 'shared', intercept = TRUE, pointcovs = 'pointcov')
            
            expect_setequal(names(Check$Formulas), c('PO', 'PA'))
            expect_setequal(names(Check$Formulas$PO), c('fish1','fish2'))
            expect_setequal(names(Check$Formulas$PA), c('bird1','bird2'))
            expect_setequal(names(Check$Formulas$PO$fish1), c('geometry', 'numvar', 'factvar_response'))
            expect_setequal(names(Check$Formulas$PO$fish2), c('geometry', 'numvar', 'factvar_response'))
            expect_setequal(names(Check$Formulas$PA$bird1), c('PAresp', 'binommark'))
            expect_setequal(names(Check$Formulas$PA$bird2), c('PAresp', 'binommark'))
            
            
            
            expect_equal(deparse1(Check$Formulas$PO$fish1$geometry$LHS), 
                         'geometry ~ .')
            expect_equal(deparse1(Check$Formulas$PO$fish1$numvar$LHS), 
                        'numvar ~ .')
            expect_equal(deparse1(Check$Formulas$PO$fish1$factvar_response$LHS), 
                         'factvar_response ~ .')
            
            expect_equal(deparse1(Check$Formulas$PA$bird1$PAresp$LHS), 
                        'PAresp ~ .')
            expect_equal(deparse1(Check$Formulas$PA$bird1$binommark$LHS), 
                        'binommark ~ .')
            
            expect_setequal(Check$Formulas$PO$fish1$geometry$RHS,
                          c("fish1_spatcovs", "fish1_PO_spatial", "shared_spatial", "fish1_intercept"))
            expect_setequal(Check$Formulas$PO$fish1$numvar$RHS,
                            c("fish1_spatcovs", "numvar_intercept", "numvar_spatial"))
            expect_setequal(Check$Formulas$PO$fish1$factvar_response$RHS,
                            c("fish1_spatcovs",
                              "factvar_spatial", "factvar", "factvar_phi"))
            
            expect_setequal(Check$Formulas$PA$bird1$PAresp$RHS,
                           c("bird1_spatcovs", "bird1_PA_spatial", "shared_spatial", "bird1_intercept", "pointcov"))
            expect_setequal(Check$Formulas$PA$bird2$binommark$RHS,
                            c("bird2_spatcovs", "binommark_spatial", "binommark_intercept"))
            
            ##Change terms
             #Set spatial and intercept to FALSE
            Check$makeFormulas(spatcovs = 'spatcovs', speciesname = 'species', temporalname = 'temp',
                               paresp = 'PAresp', countresp = 'counts', marksspatial = FALSE, speciesspatial = NULL,
                               marks = c('numvar', 'factvar', 'binommark'), markintercept = FALSE,
                               spatial = NULL, intercept = FALSE, pointcovs = 'pointcov')
          
            
            expect_setequal(Check$Formulas$PO$fish1$geometry$RHS,
                            c("fish1_spatcovs"))
            expect_setequal(Check$Formulas$PO$fish1$numvar$RHS,
                            c("fish1_spatcovs"))
            expect_setequal(Check$Formulas$PO$fish1$factvar_response$RHS,
                            c("fish1_spatcovs", "factvar", "factvar_phi"))
            
            expect_setequal(Check$Formulas$PA$bird1$PAresp$RHS,
                            c("bird1_spatcovs", "pointcov"))
            expect_setequal(Check$Formulas$PA$bird2$binommark$RHS,
                            c("bird2_spatcovs"))
            
            ##Change terms
            #Set spatcovs to NULL
            Check$makeFormulas(spatcovs = NULL, speciesname = 'species', marksspatial = TRUE, speciesspatial = 'individual',
                               paresp = 'PAresp', countresp = 'counts', markintercept = FALSE,
                               marks = c('numvar', 'factvar', 'binommark'), temporalname = 'temp',
                               spatial = 'shared', intercept = TRUE, pointcovs = 'pointcov')
            
            expect_setequal(Check$Formulas$PO$fish1$geometry$RHS,
                            c("fish1_PO_spatial", "shared_spatial", "fish1_intercept"))
            expect_setequal(Check$Formulas$PO$fish1$numvar$RHS,
                            c("numvar_spatial"))
            expect_setequal(Check$Formulas$PO$fish1$factvar_response$RHS,
                            c("factvar_spatial", "factvar", "factvar_phi"))
            
            expect_setequal(Check$Formulas$PA$bird1$PAresp$RHS,
                            c("bird1_PA_spatial", "shared_spatial", "bird1_intercept", "pointcov"))
            expect_setequal(Check$Formulas$PA$bird2$binommark$RHS,
                            c("binommark_spatial"))
            
            ##Try copy model
            Check$makeFormulas(spatcovs = NULL, speciesname = 'species', marksspatial = TRUE, speciesspatial = 'individual',
                                     paresp = 'PAresp', countresp = 'counts', markintercept = FALSE,
                                     marks = c('numvar', 'factvar', 'binommark'), temporalname = 'temp',
                                     spatial = 'copy', intercept = TRUE, pointcovs = 'pointcov')
            
            expect_setequal(Check$Formulas$PO$fish2$geometry$RHS,
                            c('PO_spatial',"fish2_PO_spatial", "fish2_intercept"))
            expect_setequal(Check$Formulas$PO$fish1$geometry$RHS,
                            c("PO_spatial", "fish1_intercept", "fish1_PO_spatial"))
            expect_setequal(Check$Formulas$PA$bird2$PAresp$RHS,
                            c("PA_spatial", "bird2_intercept", 'pointcov', "bird2_PA_spatial"))
            expect_setequal(Check$Formulas$PA$bird1$PAresp$RHS,
                            c("PA_spatial", "bird1_intercept", "pointcov", "bird1_PA_spatial"))
            
            })

#Change back to original
Check$makeFormulas(spatcovs = 'spatcovs', speciesname = 'species', marksspatial = TRUE,
                   paresp = 'PAresp', countresp = 'counts', markintercept = TRUE, speciesspatial = 'individual',
                   marks = c('numvar', 'factvar', 'binommark'), temporalname = 'temp',
                   spatial = 'shared', intercept = TRUE, pointcovs = 'pointcov')

test_that('makeComponents is able to make the correct components for all the processes
          based on the predictors and spatial effects available.', {
            
    comps <- Check$makeComponents(spatial = 'shared', intercepts = TRUE, datanames = c('PO','PA'), marksintercept = TRUE,
                         marks = c('numvar', 'factvar', 'binommark'), temporalmodel = deparse(list(model = "ar1")),
                         multinomnames = 'factvar', pointcovariates = 'pointcov', marksspatial = TRUE, offsetname = NULL,
                         speciesname = 'species', covariatenames = 'spatcovs', temporalname = 'temp', speciesspatial = 'individual',
                         covariateclass = 'numeric', numtime = 2, copymodel = Check$.__enclos_env__$private$copyModel)
    
    expect_setequal(comps,c("shared_spatial(main = geometry, model = shared_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))",
                            "fish2_PO_spatial(main = geometry, model = fish2_PO_field)",                                                                   
                            "fish1_PO_spatial(main = geometry, model = fish1_PO_field)",                                                                   
                            "bird2_PA_spatial(main = geometry, model = bird2_PA_field)",                                                                   
                            "bird1_PA_spatial(main = geometry, model = bird1_PA_field)",                                                                   
                            "numvar_spatial(main = geometry, model = numvar_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))"  ,                                                                 
                            "factvar_spatial(main = geometry, model = factvar_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))",                                                               
                            "binommark_spatial(main = geometry, model = binommark_field, group = temp, ngroup = 2, control.group = list(model = \"ar1\"))",                                                           
                            "fish2_spatcovs(main = fish2_spatcovs, model = \"numeric\")",                                                               
                            "fish1_spatcovs(main = fish1_spatcovs, model = \"numeric\")",                                                               
                            "bird2_spatcovs(main = bird2_spatcovs, model = \"numeric\")",                                                               
                            "bird1_spatcovs(main = bird1_spatcovs, model = \"numeric\")",                                                               
                            "pointcov",                                                                                                                 
                            "fish2_intercept(1)",                                                                                                       
                            "fish1_intercept(1)",                                                                                                       
                            "bird2_intercept(1)",                                                                                                      
                            "bird1_intercept(1)",                                                                                                       
                            "factvar(main = factvar, model = \"iid\",constr = FALSE, fixed=TRUE)",                                                      
                            "factvar_phi(main = factvar_phi, model = \"iid\", initial = -10, fixed = TRUE)",                                            
                            "numvar_intercept(1)",                                                                                                      
                            "binommark_intercept(1)"))
    
    ## Change arguments
    #spatial and intercepts == FALSE
    comps2 <- Check$makeComponents(spatial = NULL, intercepts = FALSE, datanames = c('PO','PA'),
                                  marks = c('numvar', 'factvar', 'binommark'), marksspatial = FALSE, offsetname = NULL,
                                  multinomnames = 'factvar', pointcovariates = 'pointcov', marksintercept = FALSE,
                                  speciesname = 'species', covariatenames = 'spatcovs', speciesspatial = 'individual',
                                  covariateclass = 'numeric', numtime =  2,  copymodel = Check$.__enclos_env__$private$copyModel)
    
    expect_setequal(comps2,c("fish2_PO_spatial(main = geometry, model = fish2_PO_field)",                       
                             "fish1_PO_spatial(main = geometry, model = fish1_PO_field)",                       
                             "bird2_PA_spatial(main = geometry, model = bird2_PA_field)",                       
                             "bird1_PA_spatial(main = geometry, model = bird1_PA_field)",                       
                             "fish2_spatcovs(main = fish2_spatcovs, model = \"numeric\")",                   
                             "fish1_spatcovs(main = fish1_spatcovs, model = \"numeric\")",                   
                             "bird2_spatcovs(main = bird2_spatcovs, model = \"numeric\")",                   
                             "bird1_spatcovs(main = bird1_spatcovs, model = \"numeric\")",                   
                             "pointcov",                                                                     
                             "factvar(main = factvar, model = \"iid\",constr = FALSE, fixed=TRUE)",          
                             "factvar_phi(main = factvar_phi, model = \"iid\", initial = -10, fixed = TRUE)"))
    
    #checkComponents with a copy model
    compsCopy <- Check$makeComponents(spatial = 'copy', intercepts = FALSE, datanames = c('PO','PA'),
                         marks = c('numvar', 'factvar', 'binommark'), marksspatial = FALSE, offsetname = NULL,
                         multinomnames = 'factvar', pointcovariates = 'pointcov', marksintercept = FALSE,
                         speciesname = 'species', covariatenames = 'spatcovs', speciesspatial = 'individual', temporalname = NULL,
                         covariateclass = 'numeric', numtime =  NULL,  copymodel = "list(beta = list(fixed = FALSE))")
    
    expect_setequal(compsCopy,c("PO_spatial(main = geometry, model = PO_field)",                                          
                             "PA_spatial(main = geometry, copy = \"PO_spatial\", hyper = list(beta = list(fixed = FALSE)))",
                             "fish1_PO_spatial(main = geometry, model = fish1_PO_field)",                                  
                             "fish2_PO_spatial(main = geometry, model = fish2_PO_field)",                                       
                             "bird1_PA_spatial(main = geometry, model = bird1_PA_field)",                                       
                             "bird2_PA_spatial(main = geometry, model = bird2_PA_field)",                                       
                             "fish1_spatcovs(main = fish1_spatcovs, model = \"numeric\")",                                   
                             "fish2_spatcovs(main = fish2_spatcovs, model = \"numeric\")",                                   
                             "bird1_spatcovs(main = bird1_spatcovs, model = \"numeric\")",                                   
                             "bird2_spatcovs(main = bird2_spatcovs, model = \"numeric\")",                                   
                             "pointcov",                                                                                     
                             "factvar(main = factvar, model = \"iid\",constr = FALSE, fixed=TRUE)",                          
                             "factvar_phi(main = factvar_phi, model = \"iid\", initial = -10, fixed = TRUE)"))
    
            
          })

