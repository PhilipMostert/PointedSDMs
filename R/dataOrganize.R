#' @title R6 class to assist in reformatting the data to be used in dataSDM.
#' @description Internal functions used to temporarily store data and other information before adding to dataSDM.
#' @export

#' @importFrom R6 R6Class

dataOrganize <- R6::R6Class(classname = 'dataOrganize',
                            lock_objects = FALSE)

dataOrganize$set('public', 'Data', list())

dataOrganize$set('public', 'Family', list())

dataOrganize$set('public', 'dataType', c())

dataOrganize$set('public', 'Marks', list())

dataOrganize$set('public', 'marksType', list())

dataOrganize$set('public', 'numObs', c())

dataOrganize$set('public', 'varsIn', list())

dataOrganize$set('public', 'multinomVars', c())

dataOrganize$set('public', 'multinomIndex', list())

dataOrganize$set('public', 'multinomNumeric', list())

dataOrganize$set('public', 'SpeciesInData', list())

dataOrganize$set('public', 'Formulas', list())

dataOrganize$set('public', 'dataSource', list())

#' @description Function to form the data into sf objects, as well as extract other information about the model required.
#' @param datapoints A list of datasets as either sf, data.frame or SpatialPoints objects.
#' @param datanames A vector of the names of the datasets.
#' @param coords Names of the coordinates used in the model.
#' @param proj The projection reference system used in the model.
#' @param marktrialname Name of the trial variable used by the binomial marks.
#' @param paresp Name of the response variable used by the presence absence datasets.
#' @param countsresp Name of the response variable used by the counts data.
#' @param trialname Name of the trial variable used by the presence absence datasets.
#' @param speciesname Name of the species name variable.
#' @param marks Name of the marks considered in the model.
#' @param pointcovnames Name of the point covariates used in the model.
#' @param markfamily A vector describing the distribution of the marks.
#' @param temporalvar The name of the temporal variable.
#' @param offsetname The name of the offset 

dataOrganize$set('public', 'makeData', function(datapoints, datanames, coords, proj, marktrialname,
                                                paresp, countsresp, trialname, speciesname, marks,
                                                pointcovnames, markfamily, temporalvar, offsetname) {
  
  ##REMOVE ALL sp AND DATAFRAME METHODS
  
  dataMade <- dataSet(datapoints = datapoints, datanames = datanames,
                      coords = c('CoordLoc1', 'CoordLoc2'), proj = proj, marks = marks,
                      paresp = paresp, countsresp = countsresp, pointcovnames = pointcovnames,
                      trialname = trialname, speciesname = speciesname, temporalvar = temporalvar,
                      marktrialname = marktrialname, markfamily = markfamily, offsetname = offsetname)
  
  self$Data <- dataMade[['Data']]
  self$Family <- dataMade[['Family']]
  self$dataType <- dataMade[['dataType']]
  self$varsIn <- dataMade[['varsIn']]
  self$Marks <- dataMade[['Marks']]
  self$marksType <- dataMade[['marksType']]
  self$numObs <- dataMade[['numObs']]
  self$multinomVars <- dataMade[['multinomVars']]
  
  if (is.null(speciesname)) {
    
    for (data in 1:length(self$Data))
      
      names(self$Data[[data]]) <- names(self$Data)[data]
  
  }
  
  self$dataSource <- lapply(seq_along(dataMade[['Family']]), function(x, name, i) rep(name[i], each = length(x[[i]])),
                            x = dataMade[['Family']], name = names(dataMade[['Family']]))
  
  
})

#' @description Function used to separate the datasets by species.
#' @param speciesname The name of the species variable.
#' @param repl Are species effects replicate. Defaults to \code{FALSE}.

dataOrganize$set('public', 'makeSpecies', function(speciesname, repl = FALSE) {
  
  self$SpeciesInData <- vector(mode = 'list', length = length(self$Data))
  
  ##Need #unique species in each
  numUnique <- c()
  
  for (dataset in 1:length(self$Data)) {
    
    allSpecies <- data.frame(self$Data[[dataset]][[1]])[, speciesname]
    self$SpeciesInData[[dataset]] <- uniqueSpecies <- unique(allSpecies)
    numUnique[dataset] <- length(uniqueSpecies)
    names(self$SpeciesInData)[[dataset]] <- datasetName <- names(self$Data)[[dataset]]
    speciesList <- list()
    
    for (species in uniqueSpecies) {
      
      speciesList[[species]] <- self$Data[[dataset]][[1]][data.frame(self$Data[[dataset]][[1]])[,speciesname] == species,]
      
    }
    
    self$Data[[dataset]] <- speciesList
    
    names(self$Data[[dataset]]) <- paste0(datasetName,'_',uniqueSpecies)
    
  }
  
  self$dataSource <- mapply(rep, self$dataSource, each = numUnique)
  
  ##Make species index here ...
  
  self$makeMultinom(multinomVars = speciesname,
                    return = 'species', oldVars = NULL, repl = repl)
  
})

#' @description Function to make multinomial marks, if any are included in the model.
#' @param multinomVars Name of the multinomial marks.
#' @param return What to return: species or marks.
#' @param oldVars If any multinomial marks were included in a previous iteration, what where their names.
#' @param repl Species replicate model included. Defaults to \code{FALSE}.

dataOrganize$set('public', 'makeMultinom', function(multinomVars, return, oldVars, repl = FALSE) {
  
  #if (return %in% c('marks','species')) stop('Something went wrong.')
  #Need to change the multinomial vars into numeric
  # so create a generic index across datasets.
  multinomIndex <- vector(mode = 'list', length = length(multinomVars))
  names(multinomIndex) <- multinomVars
  multinomNumeric <- vector(mode = 'list', length(multinomVars))
  names(multinomNumeric) <- multinomVars
  
  for (var in multinomVars) {
    
    multinomIndex[[var]] <- vector(mode = 'list', length(self$Data))
    names(multinomIndex[[var]]) <- names(self$Data)
    multinomNumeric[[var]] <- vector(mode = 'list', length(self$Data))
    names(multinomNumeric[[var]]) <- names(self$Data)
    
    for (dataset in names(self$Data)) {
      
      #multinomIndex[[var]][[dataset]] <- vector(mode = 'list', length(self$Data[[dataset]]))
      #multinomNumeric[[var]][[dataset]] <- vector(mode = 'list', length(self$Data[[dataset]]))
      if (is.null(names(self$Data[[dataset]]))) spIndex <- 1
      else spIndex <- length(names(self$Data[[dataset]]))
      
      for (species in 1:spIndex) {
        
        if (var %in% names(self$Data[[dataset]][[species]])) {
          
          multinomIndex[[var]][[dataset]][[species]] <- as.character(data.frame(self$Data[[dataset]][[species]])[, var])
          
        }
        else multinomIndex[[var]][[dataset]][[species]] <- NA 
        
      }
     
      
    }
    
    if (!is.null(oldVars)) {
      
      if (var %in% names(oldVars)) {
        
        lenOld <- length(oldVars[[var]])
        
        facts <- as.numeric(as.factor(c(oldVars[[var]],unlist(multinomIndex[[var]]))))
        facts <- facts[-(1:lenOld)]
        
      } else facts <- as.numeric(as.factor(unlist(multinomIndex[[var]])))
      
    }
    else facts <- as.numeric(as.factor(unlist(multinomIndex[[var]])))
    
    multinomNumeric[[var]] <- relist(facts,  
                                     multinomIndex[[var]])
    
    ##Now need to create index back in
    for (x in names(multinomNumeric[[var]])) {
      
      for (y in 1:length((multinomNumeric[[var]][[x]]))) {
        
        if (!all(is.na((multinomNumeric[[var]][[x]][[y]])))) {
          
          self$Data[[x]][[y]][,var] <- multinomNumeric[[var]][[x]][[y]]
          if (repl) self$Data[[x]][[y]][,'speciesSpatialGroup'] <- multinomNumeric[[var]][[x]][[y]]
          
        }
      }
    }
  }
  
  if (return == 'species') {
    
    self$speciesNumeric <- multinomNumeric
    self$speciesIndex <- multinomIndex
    
  }
  else 
    if (return == 'time') {
      
      self$timeNumeric <- multinomNumeric
      self$timeIndex <- multinomIndex
      
    }
  else {
    
    self$multinomNumeric <- multinomNumeric
    self$multinomIndex <- multinomIndex
    
  }
  
})

#' @description Function used to create formulas for the processes.
#' @param spatcovs Names of the spatial covariates used in the model.
#' @param speciesname Name of the species variable.
#' @param paresp Name of the presence absence response variable.
#' @param countresp Name of the count data response variable.
#' @param marks Name of the marks used in the model.
#' @param temporalname Name of the temporal variable used in the model.
#' @param spatial Logical: are spatial effects run in the model.
#' @param marksspatial Logical: should spatial fields be included for the marks.
#' @param speciesintercept Logical: should specific species intercept terms be created for the species.
#' @param speciesenvironment Logical: should species specific environmental terms be created.
#' @param intercept Logical: are intercepts run in the model.
#' @param markintercept Logical: are intercepts run for the marks in the model.
#' @param pointcovs Name of the point covariates.
#' @param speciesspatial Logical: should the species have spatial fields.
#' @param speciesindependent Logical: make independent species effects.
#' @param biasformula Terms to include for PO data.
#' @param covariateformula Terms to include for the covariate formula.

dataOrganize$set('public', 'makeFormulas', function(spatcovs, speciesname,
                                                    paresp, countresp, marks, marksspatial, 
                                                    speciesintercept, speciesenvironment,
                                                    spatial, intercept, temporalname, speciesindependent,
                                                    markintercept, pointcovs, speciesspatial,
                                                    biasformula, covariateformula) {
  
  #if (length(self$multinomVars) != 0) marks[marks %in% self$multinomVars] <- paste0(marks[marks %in% self$multinomVars],'_response')
  
  if (!is.null(marks)) {
    
    if (length(self$multinomVars) != 0) marksResps <- c(paste0(marks[marks %in% self$multinomVars],'_response'), marks[!marks %in% self$multinomVars])
    else marksResps <- marks
    
  }
  else marksResps <- NULL
  
  if (!is.null(biasformula)) {
    
    biasTerms <- labels(terms(biasformula))
    spatcovs <- spatcovs[!spatcovs %in% biasTerms]
    if (identical(spatcovs, character(0))) spatcovs <- NULL
    
  }
  
  formulas <- vector(mode = 'list', length = length(self$Data))
  names(formulas) <- names(self$Data)
  
  for (dataset in 1:length(self$Data)) {
    
    pointsResponse <- lapply(self$Data[[dataset]], function(x) {
      
      if (paresp %in% names(x)) names(x)[names(x) %in% c(paresp, marksResps)]
      else 
        if (countresp %in% names(x)) names(x)[names(x) %in% c(countresp, marksResps)]
      else {
        
        marksIn <- names(x)[names(x) %in% marksResps]
        responses <- c('geometry', marksIn)
        responses      
        
      }
      
    })
    
    datasetCovs <- self$varsIn[[dataset]] 
    formulas[[dataset]] <- vector(mode = 'list', length = length(self$Data[[dataset]]))
    
    for (species in 1:length(self$Data[[dataset]])) {
      
      if (!is.null(speciesname)) {
        
        ##Fix this later depending on weather we are running INLA grouped model for the species or not...
        if (length(self$speciesIndex[[speciesname]][[dataset]][[species]]) != 0) speciesIn <- unique(self$speciesIndex[[speciesname]][[dataset]][[species]])   
        else speciesIn <- unique(self$Data[[dataset]][[species]][,speciesname])
        names(formulas[[dataset]])[species] <- speciesIn
        #length_index <- length(pointsResponse[[species]])
        length_index <- length(speciesIn)
        
      } 
      else {
        
        speciesIn <- NULL
        names(formulas[[dataset]]) <- names(formulas)[[dataset]]
        length_index <- 1
      }
      
      formulas[[dataset]][[species]] <- vector(mode = 'list', length = length(pointsResponse[[species]]))
      
      for (response in 1:length_index) {
        
        for (j in 1:length(pointsResponse[[response]])) {
          
          #if (!is.null(temporalname)) temp <- paste0(temporalname,'_effect')
          #else temp <- NULL
          
          if (!is.null(speciesIn)) {
             if (pointsResponse[[response]][j] %in% c('geometry', paresp, countresp)) {
              ##Change this part ot the speciesIn: not sure what the one below does...
              if (!is.null(speciesspatial)) {
                
                if (speciesspatial == 'shared' || speciesspatial == 'replicate') speciesspat <- 'speciesShared'
                else {
                
                  
                if (speciesindependent) speciesspat <- paste0(speciesIn,'_spatial')
                else speciesspat <- paste0(do.call(paste0, expand.grid(paste0(speciesIn,'_'),
                                                                                              names(self$Data)[[dataset]])),
                                                                  '_spatial') ## new argument called speciesSpatial??
                
                }
                
              }
              else speciesspat <- NULL

            
            #else {  
            #  
            #  if (spatial) spat <- c(paste0(speciesIn,'_spatial'), 'shared_spatial')
            #  else spat <- NULL
              
            #}
            if (!is.null(speciesintercept)) {
              
              if (speciesintercept) spint <- paste0(speciesname, '_intercepts')
              else spint <- paste0(speciesIn,'_intercept') 
              
            } else spint <- NULL
              
            
            if (intercept) {
              
              int <- paste0(names(self$Data)[[dataset]],'_intercept') 
              
            } else int <- NULL  
            
             }
            else {
              
              speciesspat <- NULL
              spint <- NULL

            }
          } 
          else {
            
            speciesspat <- NULL
            spint <- NULL
            if (intercept) int <- paste0(names(self$Data)[[dataset]], '_intercept')
            else int <- NULL
            
          }
        
          if (!is.null(spatial)) {
            
            
           if (spatial %in% c('shared', 'correlate')) spat <- 'shared_spatial'
            else 
              if (spatial %in% c('individual', 'copy')) spat <- paste0(names(self$Data)[[dataset]], '_spatial')
              
              } else spat <- NULL
            
          

          
          if (!is.na(datasetCovs)) {
            #Species specific? dataset specific?
            if (pointsResponse[[response]][j] %in% c(paresp, countresp, 'geometry')) addcovs <- unlist(datasetCovs)
            else addcovs <- NULL
            
          }
          else addcovs <- NULL
          
          if (!is.null(spatcovs)) {
            
            if (!is.null(biasformula)) {
              
              if (!pointsResponse[[response]][j] %in% c(paresp, countresp, marks)) {
                
                if (speciesenvironment) biascov <- 'Bias__Effects__Comps'#paste0(speciesIn, '_Bias__Effects__Comps')
                else biascov <- 'Bias__Effects__Comps'
              
                } 
              else biascov <- NULL
              
            } else biascov <- NULL
            
            if (!is.null(covariateformula)) {
              
              if (speciesenvironment) covs <- paste0(speciesIn, '_Fixed__Effects__Comps')
              else covs <- 'Fixed__Effects__Comps'
              
            }
            
            else {
              
              if (!is.null(speciesname)) {
              
              if (speciesenvironment) covs <- paste0(speciesIn, '_', spatcovs)
              else covs <- spatcovs
              
            }
            else covs <- spatcovs
            
          } 
          
          } else {
            
            covs <- NULL
            biascov <- NULL
            
          }
          
          if (!is.null(marks)) {
            
            if (pointsResponse[[response]][j] %in% c('geometry', paresp, countresp)) {
              
              markspat <- NULL
              marksint <- NULL
              
            }
            else {
              
              #if (!is.null(spatial)) spat <- NULL
              spat <- paste0(names(self$Data)[[dataset]], '_', pointsResponse[[response]][j], '_spatial')
              
              if (marksspatial) {
                
                if (pointsResponse[[response]][j] %in% paste0(marks,'_response')) markspat <- paste0(marks[pointsResponse[[response]][j] == paste0(marks,'_response')],'_spatial')
                else markspat <- paste0(pointsResponse[[response]][j], '_spatial')
                
                
              } else markspat <- NULL
              
              if (intercept) int <- NULL
              
              if (markintercept) {
                
                if (length(self$multinomVars) != 0) marks_intercepts <- marks[!marks %in% self$multinomVars]
                else marks_intercepts <- marks
                
                #if (pointsResponse[[response]][j] %in% paste0(marks_intercepts,'_response')) marksint <- paste0(marks[pointsResponse[[response]][j] == paste0(marks,'_response')],'_intercept')
                #else marksint <- paste0(pointsResponse[[response]][j], '_intercept')
                #else marksint <- NULL
                if (pointsResponse[[response]][j] %in% marks_intercepts) marksint <- paste0(pointsResponse[[response]][j], '_intercept')
                else marksint <- NULL
                
                
              } else marksint <- NULL
              
            }
            
          }
          else {
            
            markspat <- NULL
            marksint <- NULL
            
          }
          
          RHS <- c(covs, spat, int, addcovs, markspat, marksint, speciesspat, biascov, spint) # temp
          
          if (pointsResponse[[response]][j] %in% paste0(self$multinomVars,'_response')) { #paste multinomvar and phi # Need to convert multinomvar to numeric
            
            mnVar <- self$multinomVars[pointsResponse[[response]][j] %in% paste0(self$multinomVars,'_response')]
            RHS <- c(RHS, mnVar, paste0(mnVar, '_phi'))
            
          }
          
          #formulas[[dataset]][[species]][[j]] <- formula(paste0(pointsResponse[[response]][j], ' ~ ', paste0(RHS, collapse = ' + ')))
          ##Test
          formulas[[dataset]][[species]][[j]] <- vector(mode = 'list', length = 2)
          names(formulas[[dataset]][[species]][[j]]) <- c('LHS', 'RHS')
          formulas[[dataset]][[species]][[j]][[1]] <- formula(paste0(pointsResponse[[response]][j], ' ~ .'))
          formulas[[dataset]][[species]][[j]][[2]] <- RHS
          names(formulas[[dataset]][[species]])[j] <- pointsResponse[[response]][j]
          
        }
      }
    }
  }
  self$Formulas <- formulas
  
})

#' @description Function used to make components for the model.
#' @param spatial Logical: are spatial effects run in the model.
#' @param intercepts Logical: are intercepts run in the model.
#' @param datanames Names of the datasets used in the model.
#' @param marks Names of the marks used in the model.
#' @param speciesname Name of the species variable.
#' @param multinomnames Names of the multinomial marks.
#' @param pointcovariates Names of the point covariates.
#' @param covariatenames Names of the spatially varying covariates.
#' @param covariateclass The classes of the spatially varying covariates.
#' @param temporalname Name of the temporal variable used in the model.
#' @param marksspatial Logical: should spatial fields be included for the marks.
#' @param marksintercept Logical: should intercepts be included for the marks.
#' @param speciesintercept: Logical: should intercept terms be created for each species.
#' @param numtime Number of time increments included in the model.
#' @param speciesspatial Logical: Should the species be run with spatial fields.
#' @param speciesenvironment Logical: Should the species have their own environmental effects.
#' @param offsetname Name of the offset column in the datasets.
#' @param copymodel List of the hyper parameters for the \code{copy} model.
#' @param speciesindependent Logical: should species effects be made independent.
#' @param biasformula Terms to include for PO data.
#' @param covariateformula Terms to include for the covariate formula.
#' @param marksCopy Names of where each mark occurs in the dataset.

dataOrganize$set('public', 'makeComponents', function(spatial, intercepts, 
                                                      datanames, marks, speciesname,
                                                      multinomnames, pointcovariates,
                                                      covariatenames,covariateclass,
                                                      marksspatial,
                                                      marksintercept,
                                                      temporalname,
                                                      speciesspatial,
                                                      speciesenvironment,
                                                      speciesintercept,
                                                      numtime,temporalmodel,
                                                      offsetname,
                                                      copymodel,
                                                      speciesindependent,
                                                      biasformula,
                                                      covariateformula,
                                                      marksCopy) {
  
  if (!is.null(biasformula)) {
    
    biasTerms <- labels(terms(biasformula))
    removeIndex <- !covariatenames %in% biasTerms
    covariatenames <- covariatenames[removeIndex]
    covariateclass <- covariateclass[removeIndex]
    
    if (identical(covariatenames, character(0))) {
      
      covariatenames <- NULL
      covariateclass <- NULL
      
    }
    
    
  }
  ##Copy for marks fields???
  if (length(self$SpeciesInData) != 0) species <- unique(unlist(self$SpeciesInData))
  else species = NULL
  
  if (!is.null(spatial)) {
   
    if (spatial == 'shared') {
    
    if (!is.null(temporalname)) spat <- paste0('shared_spatial(main = geometry, model = shared_field, group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
    else spat <- paste0('shared_spatial(main = geometry, model = shared_field)')
    
    }
    else 
      if (spatial == 'copy') {
      
     mainName <- datanames[[1]]
     
     if (!is.null(temporalname)) {
       
       spatMain <- paste0(mainName, '_spatial(main = geometry, model = ', paste0(mainName,'_field'), ', group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
       spatCopy <-  paste0(datanames[datanames != mainName], '_spatial(main = geometry, copy = \"', paste0(mainName,'_spatial'),'\", group = ', temporalname, ', control.group = ', temporalmodel,', hyper = ', copymodel,')')
       
     }
     else {
       
       spatMain <-  paste0(mainName, '_spatial(main = geometry, model = ', paste0(mainName,'_field'),')')
       spatCopy <-  paste0(datanames[datanames != mainName], '_spatial(main = geometry, copy = \"', paste0(mainName,'_spatial'),'\", hyper = ', copymodel,')')
       
     }
     
     spat <- c(spatMain, spatCopy)    
        
      }
     else
       if (spatial == 'correlate') {
         
         spat <- paste0('shared_spatial(main = geometry, model = shared_field, group = ._dataset_index_var_., control.group = list(model = "exchangeable"))')
         
       }
    else {
      
    if (!is.null(temporalname)) spat <- paste0(datanames, '_spatial(main = geometry, model = ', paste0(datanames,'_field'), ', group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
    else spat <-  paste0(datanames, '_spatial(main = geometry, model =', paste0(datanames,'_field'),')')
      
    }
    
  } else spat <- NULL
  
  if (!is.null(species)) {
      
      if (!is.null(speciesintercept)) {
        
        if (speciesintercept) {
          
          if (intercepts) spint <- paste0(speciesname, '_intercepts(main = ', speciesname, ', model = "iid", constr = TRUE, hyper = list(prec = list(prior = "loggamma", param = c(1, 5e-05))))')
          else spint <- paste0(speciesname, '_intercepts(main = ', speciesname, ', model = "iid", constr = FALSE, hyper = list(prec = list(prior = "loggamma", param = c(1, 5e-05))))')
        }
        else spint <- paste0(species, '_intercept(1)')
        
      } else spint <- NULL
    
    if (!is.null(speciesspatial)) { ## Then if copy or individual ...

      #if (length(speciesspatial) > 0) {
      if (speciesspatial == 'shared') speciesSpat <- 'speciesShared(main = geometry, model = speciesField)'
      
      else
        if (speciesspatial == 'replicate') {
          
          speciesSpat <- 'speciesShared(main = geometry, model = speciesField, group = speciesSpatialGroup, control.group = list(model = "iid"))'
          
          
        } else {
      #speciesSpat <- paste0(species,'_spatial(main = coordinates, model = ',paste0(species,'_spdeModel'),')', collapse = ' + ')
 
      #} 
      #else {
      #ie we are assuming and INLA grouped model
      if (length(self$speciesIndex) != 0) {
        ##Change the species part to model = paste0(speciesname) ##where speciesname = species
         # but keep the speciesSpat framework for the temporal part of the model
        #speciesSpat <- paste0(speciesname, '_spatial(main = coordinates, model = speciesModel, group = ',speciesname,', ngroup = ', numspecies,')')
        #if (speciesspatial == 'individual' || length(species) == 1) {
        if (speciesspatial == 'individual') {  
          if (speciesindependent) speciesSpat <- paste0(species,'_spatial(main = geometry, model = ',paste0(species,'_field)'))
          else {
          
          speciesComps <- list()
          
          for (spec in species) {
            
          specInData <- sapply(self$SpeciesInData, FUN = function(x) spec %in% x)
            
          speciesComps[[spec]] <- paste0(do.call(paste0, expand.grid(paste0(spec,'_'), names(specInData[specInData]))),'_spatial(main = geometry, model = ',
                                                                  paste0(do.call(paste0, expand.grid(paste0(spec,'_'),
                                                                                                     names(specInData[specInData]))),'_field)'))
          
          }
          
          speciesSpat <- unlist(speciesComps)
          
          }
          
        }
        else {
        
          if (speciesindependent) {
            
            speciesOne <- paste0(species[1],'_spatial(main = geometry, model = ',paste0(species[1],'_field)'))
            if (length(species) > 1) speciesOther <- paste0(species[-1],'_spatial(main = geometry, copy = \"', species[1],'_spatial\",  hyper = list(beta = list(fixed = FALSE)))')
            else speciesOther <- NULL
          
          }
          else {
          
          speciesOne <- list()
          speciesOther <- list()
          
          for (spec in species) {
   
            specInData <- sapply(self$SpeciesInData, FUN = function(x) spec %in% x)
        
            if (sum(specInData) == 1) {
              
              speciesOne[[spec]] <- paste0(spec,'_', names(specInData[specInData]),'_spatial(main = geometry, model = ',
                                  paste0(spec,'_', names(specInData[specInData]),'_field)'))
              speciesOther[[spec]] <- NULL
              
            }
            
            else {
              
              dataForSpec <- names(specInData[specInData])
              
              speciesOne[[spec]] <- paste0(paste0(spec, '_', dataForSpec[1]),'_spatial(main = geometry, model = ',
                                           paste0(spec, '_',dataForSpec[1]),'_field)')
              
              speciesOther[[spec]] <- paste0(do.call(paste0, expand.grid(paste0(spec, '_'), dataForSpec[-1])),
                                             '_spatial(main = geometry, copy = \"', do.call(paste0, expand.grid(paste0(spec, '_'), dataForSpec[1]))
                                             ,'_spatial\",  hyper = list(beta = list(fixed = FALSE)))')
              
              
            }
            
            
            
          }
          
          }
          speciesSpat <- c(unlist(speciesOne), unlist(speciesOther))
          
          
        }
      }
      }
      #else speciesSpat <- paste0(species,'_spatial(main = coordinates, model = speciesModel)', collapse = ' + ') #change this to speciesModel
      
      #}
      
    }
    else speciesSpat <- NULL
    
  } 
  else {
    
    speciesSpat <- NULL
    spint <- NULL
    
  }
  
  #if (!is.null(temporalname)) {
  #  
  #  tempSpat <- paste0(temporalname,'_effect(main = coordinates, model = temporal_field, group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
    
  # } else tempSpat <- NULL
  
  #Covariates: Dataset covariates #Done
  #            Species covariates #Done
  #            Marks covariates? #Surely if we are grouping by species, then mark covariates will be the same as dataset/species??
  #            Point covariates #Done
  
  if (!is.null(marks)) {
    
    if (marksspatial) {
      
      if (!is.null(spatial)) {
        
        marksCopy <- mapply(FUN = function(x, names) paste0(names,'_', x), x = marksCopy, names = names(marksCopy))
        marksCopySpat <- list()
        
        if (spatial  %in% c('shared', 'correlate')) marksCopySpat <- paste0(unlist(marksCopy), '_spatial(main = geometry, copy = \"shared_spatial\")')#if shared or correlate
        else
          if (spatial == 'individual') {
            for (dat in names(marksCopy)) {
            
            marksCopySpat[[dat]] <- paste0(unlist(marksCopy[dat]), '_spatial(main = geometry, copy = \"', dat, '_spatial\")') #xx are datasetnames
            
            }
          }
        else {
          
          marksCopySpat <- paste0(unlist(marksCopy), '_spatial(main = geometry, copy = \"', mainName, '_spatial\")')
          
          
        }

        if (inherits(marksCopySpat, 'list')) marksCopySpat <- unlist(marksCopySpat)
        
      } else marksCopySpat <- NULL
      
      if (!is.null(temporalname)) marksSpat <- paste0(marks, '_spatial(main = geometry, model = ', paste0(marks,'_field'), ', group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
      else marksSpat <- paste0(marks, '_spatial(main = geometry, model = ', paste0(marks,'_field)'))
      
    }
    else {
      
      marksSpat <- NULL
      marksCopySpat <- NULL
      
    }
    
    if (marksintercept) {
      
      if (length(self$multinomVars) != 0)  marks_intercepts <- marks[!marks %in% self$multinomVars]
      else marks_intercepts <- marks
      
      if (identical(marks_intercepts, 'character(0)')) marksInt <- NULL
      else marksInt <- paste0(marks_intercepts, '_intercept(1)')
      
    } else marksInt <- NULL
    
  }
  else {
    
    marksSpat <- NULL
    marksInt <- NULL
    marksCopySpat <- NULL
  }
  
  if (!is.null(covariatenames)) {
    
    ##IF bias covariates
    
    if (!is.null(biasformula)) bias <- makeFormulaComps(form = biasformula, species = FALSE, speciesnames = species, type = 'Bias')
    else bias <- NULL
      
    
     #IF covariateFormula
    
    if (!is.null(covariateformula)) {
      
      covs <-  makeFormulaComps(form = covariateformula, species = speciesenvironment, speciesnames = species, type = 'Covariate')
      
    }
    
    else {
      
      if (!is.null(species) && speciesenvironment) {
      
      speciesCovs <- apply(expand.grid(paste0(species,'_'), covariatenames), MARGIN = 1, FUN = paste0,collapse = '')
      speciesCovClass <- rep(covariateclass, each = length(species))
      covs <- paste0(speciesCovs, '(main = ', speciesCovs, ', model = \"',speciesCovClass,'\")') #, collapse = ' + '
      
    }
    else covs <- paste0(covariatenames, '(main = ', covariatenames, ', model = \"',covariateclass,'\")') # , collapse = ' + '
    
    }
  } 
  else {
    
    covs <- NULL
    bias <- NULL
    
  }
  
  if (!is.null(pointcovariates)) {
    
    covsPoints <- pointcovariates#paste0(pointcovariates, collapse = ' + ')
    
  } else covsPoints <- NULL
  
  if (!is.null(offsetname)) {
    
    if (any(offsetname %in% unlist(self$varsIn))) {
     
      offsetIn <- offsetname[offsetname %in% unlist(self$varsIn)]
      
      offsetTerm <- paste0(offsetIn,'(log(',offsetIn,'), model = "offset")')
      
    } else offsetTerm <- NULL
    
  } else offsetTerm <- NULL
  
  #Intercepts: Dataset intercepts
  #            species intercepts
  #            marks intercepts # Am I doing this? Or are we doing dataset/species mark intercepts
  
  if (intercepts) {
    
  int <- paste0(datanames, '_intercept(1)')
  
  if (!is.null(marks)) intMarks <- paste0(marks, '_intercept(1)')
  else intMarks <- NULL
    
  } 
  else {
    
    int <- NULL
    intMarks <- NULL
    
  } 
  
  if (!is.null(multinomnames)) {
    
    multinomVars <- paste0(unique(multinomnames),'(main = ', unique(multinomnames), ', model = "iid",constr = FALSE, fixed=TRUE)', collapse = ' + ')
    multinomPhi <- paste0(unique(multinomnames),'_phi(main = ', unique(multinomnames), '_phi, model = "iid", initial = -10, fixed = TRUE)', collapse = ' + ')
    
  }
  else {
    
    multinomVars <- NULL
    multinomPhi <- NULL
    
  }
  
  RHS <- c(spat, speciesSpat, marksSpat, covs, covsPoints, int, multinomVars, multinomPhi, marksInt, offsetTerm, bias, spint,marksCopySpat)
  
  RHS
  
})

#' @description Function to make the datasets into likelihoods.
#' @param mesh An \code{fm_mesh_2d} object.
#' @param ips Integration points used.
#' @param paresp The response variable name for the presence absence datasets.
#' @param ntrialsvar The trials variable name for the presence absence datasets.
#' @param markstrialsvar The trial variable name for the binomial marks.
#' @param speciesname The name of the species variable used.

dataOrganize$set('public', 'makeLhoods', function(mesh, ips,
                                                  paresp, ntrialsvar,
                                                  markstrialsvar, speciesname) {
  
  
  Likelihoods <- list()
  #formulas <- unlist(self$Formulas)
  
  for (dataset in 1:length(self$Data)) {
    
    for (species in 1:length(self$Data[[dataset]])) {
      
      Ntrialsvar <- list()
      
      if (!is.null(ntrialsvar)) {
        
        if (ntrialsvar %in% names(self$Data[[dataset]][[species]]@data)) {
          
          Ntrialsvar[[1]] <- self$Data[[dataset]][[species]]@data[,ntrialsvar]
          
        } else Ntrialsvar[[1]] <- NA
        
      } else Ntrialsvar[[1]] <- NA
      
      if (!is.null(markstrialsvar)) {
        
        if (markstrialsvar %in% names(self$Data[[dataset]][[species]]@data)) {
          
          Ntrialsvar[[2]] <- self$Data[[dataset]][[species]]@data[,markstrialsvar]
          
        } else Ntrialsvar[[2]] <- NA
        
      } else Ntrialsvar[[2]] <- NA
      
      for (process in 1:length(self$Family[[dataset]])) {
        
        Likindex <- length(Likelihoods) + 1
        
        if (self$Family[[dataset]][process] == 'binomial') {
          
          if (!is.na(Ntrialsvar[[1]]) || !is.na(Ntrialsvar[[2]])) {
            
            if (as.character(self$Formulas[[dataset]][[species]][[process]][['LHS']])[2] == paresp) Ntrials <- Ntrialsvar[[1]]
            
            else Ntrials <- Ntrialsvar[[2]]
            
          } else Ntrials <- 1
        } 
        else Ntrials <- 1
        
        if (!is.null(speciesname)) {
          
          if (length(self$speciesIndex) != 0) {
            
            speciesRep <- data.frame(rep(unique(self$Data[[dataset]][[species]]@data[,speciesname]), nrow(ips@coords)))
            names(speciesRep) <- speciesname
            
            IPS <- ips
            IPS@data <- cbind(ips@data, speciesRep)
            
          }
          
        }
        else IPS <- ips
        ##rm formulas for now

        # bru_like_list will use the like-tag, not the list names.
        if (is.null(names(self$Data[[dataset]])[species])) nameGive <- names(self$Data)[[dataset]]
        else nameGive <- names(self$Data[[dataset]])[species]
        
        like_name <- paste0(nameGive, '_', as.character(self$Formulas[[dataset]][[species]][[process]][['LHS']])[2])
        
        Likelihoods[[Likindex]] <- inlabru::like(formula = self$Formulas[[dataset]][[species]][[process]][['LHS']],
                                                 include = self$Formulas[[dataset]][[species]][[process]][['RHS']],
                                                 data = self$Data[[dataset]][[species]], 
                                                 Ntrials = Ntrials,
                                                 ips = IPS,
                                                 family = self$Family[[dataset]][process],
                                                 tag = like_name)
        
        names(Likelihoods)[[Likindex]] <- like_name
        
      }
      
    }
    
  }
  
  self$Data <- Likelihoods
  
})

