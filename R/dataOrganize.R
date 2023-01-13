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

#' @description Function to form the data into SpatialPointsDataFrames, as well as extract other information about the model required.
#' @param datapoints A list of datasets as either data.frame or SpatialPoints objects.
#' @param datanames A vector of the names of the datasets.
#' @param coords Names of the coordinates used in the model.
#' @param proj The projection reference system used in the model.
#' @param marktrialname Name of the trial variable used by the binomial marks.
#' @param paresp Name of the response variable used by the presence absence datasets.
#' @param countsresp Name of the respons variable used by the counts data.
#' @param trialname Name of the trial variable used by the presence absence datasets.
#' @param speciesname Name of the species name variable.
#' @param marks Name of the marks considered in the model.
#' @param pointcovnames Name of the point covariates used in the model.
#' @param markfamily A vector describing the distribution of the marks.

dataOrganize$set('public', 'makeData', function(datapoints, datanames, coords, proj, marktrialname,
                                                paresp, countsresp, trialname, speciesname, marks,
                                                pointcovnames, markfamily, temporalvar, offsetname) {
  
  dataMade <- dataSet(datapoints = datapoints, datanames = datanames,
                      coords = coords, proj = proj, marks = marks,
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

dataOrganize$set('public', 'makeSpecies', function(speciesname) {
  
  self$SpeciesInData <- vector(mode = 'list', length = length(self$Data))
  
  ##Need #unique species in each
  numUnique <- c()
  
  for (dataset in 1:length(self$Data)) {
    
    allSpecies <- self$Data[[dataset]][[1]]@data[, speciesname]
    self$SpeciesInData[[dataset]] <- uniqueSpecies <- unique(allSpecies)
    numUnique[dataset] <- length(uniqueSpecies)
    names(self$SpeciesInData)[[dataset]] <- datasetName <- names(self$Data)[[dataset]]
    speciesList <- list()
    
    for (species in uniqueSpecies) {
      
      speciesList[[species]] <- self$Data[[dataset]][[1]][self$Data[[dataset]][[1]]@data[,speciesname] == species,]
      
    }
    
    self$Data[[dataset]] <- speciesList
    
    names(self$Data[[dataset]]) <- paste0(datasetName,'_',uniqueSpecies)
    
  }
  
  self$dataSource <- mapply(rep, self$dataSource, each = numUnique)
  
  ##Make species index here ...
  
  self$makeMultinom(multinomVars = speciesname,
                    return = 'species', oldVars = NULL)
  
})

#' @description Function to make multinomial marks, if any are included in the model.
#' @param multinomVars Name of the multinomial marks.
#' @param return What to return: species or marks.
#' @param oldVars If any multinomial marks were included in a previous iteration, what where their names.

dataOrganize$set('public', 'makeMultinom', function(multinomVars, return, oldVars) {
  
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
        
        if (var %in% names(self$Data[[dataset]][[species]]@data)) {
          
          multinomIndex[[var]][[dataset]][[species]] <- as.character(self$Data[[dataset]][[species]]@data[, var])
          
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
          
          self$Data[[x]][[y]]@data[,var] <- multinomNumeric[[var]][[x]][[y]]
          
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
#' @param intercept Logical: are intercepts run in the model.
#' @param markintercept Logical: are intercepts run for the marks in the model.
#' @param pointcovs Name of the point covariates.
#' @param speciesspatial Logical: should the species have spatial fields.

dataOrganize$set('public', 'makeFormulas', function(spatcovs, speciesname,
                                                    paresp, countresp, marks, marksspatial,
                                                    spatial, intercept, temporalname,
                                                    markintercept, pointcovs, speciesspatial) {
  
  #if (length(self$multinomVars) != 0) marks[marks %in% self$multinomVars] <- paste0(marks[marks %in% self$multinomVars],'_response')
  
  if (!is.null(marks)) {
    
    if (length(self$multinomVars) != 0) marksResps <- c(paste0(marks[marks %in% self$multinomVars],'_response'), marks[!marks %in% self$multinomVars])
    else marksResps <- marks
    
  }
  else marksResps <- NULL
  
  formulas <- vector(mode = 'list', length = length(self$Data))
  names(formulas) <- names(self$Data)
  
  for (dataset in 1:length(self$Data)) {
    
    pointsResponse <- lapply(unlist(self$Data[[dataset]]), function(x) {
      
      if (paresp %in% names(x)) names(x)[names(x) %in% c(paresp, marksResps)]
      else 
        if (countresp %in% names(x)) names(x)[names(x) %in% c(countresp, marksResps)]
      else {
        
        marksIn <- names(x)[names(x) %in% marksResps]
        responses <- c('coordinates', marksIn)
        responses      
        
      }
      
    })
    
    datasetCovs <- self$varsIn[[dataset]] 
    formulas[[dataset]] <- vector(mode = 'list', length = length(self$Data[[dataset]]))
    
    for (species in 1:length(self$Data[[dataset]])) {
      
      if (!is.null(speciesname)) {
        
        ##Fix this later depending on weather we are running INLA grouped model for the species or not...
        if (length(self$speciesIndex[[speciesname]][[dataset]][[species]]) != 0) speciesIn <- unique(self$speciesIndex[[speciesname]][[dataset]][[species]])   
        else speciesIn <- unique(self$Data[[dataset]][[species]]@data[,speciesname])
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
             if (pointsResponse[[response]][j] %in% c('coordinates', paresp, countresp)) {
              ##Change this part ot the speciesIn: not sure what the one below does...
              if (speciesspatial) speciesspat <- paste0(speciesIn,'_spatial') ## new argument called speciesSpatial??
              else speciesspat <- NULL

            
            #else {  
            #  
            #  if (spatial) spat <- c(paste0(speciesIn,'_spatial'), 'shared_spatial')
            #  else spat <- NULL
              
            #}
            
            if (intercept) int <- paste0(speciesIn,'_intercept') 
            else int <- NULL
             }
            else {
              
              speciesspat <- NULL

            }
          } 
          else {
            
            speciesspat <- NULL
            
            if (intercept) int <- paste0(names(self$Data)[[dataset]], '_intercept')
            else int <- NULL
            
          }
        
          if (!is.null(spatial)) {
            
            
           if (spatial == 'shared') spat <- 'shared_spatial'
            else 
              if (spatial %in% c('individual', 'copy')) spat <- paste0(names(self$Data)[[dataset]], '_spatial')
              
              } else spat <- NULL
            
          

          
          if (!is.na(datasetCovs)) {
            #Species specific? dataset specific?
            if (pointsResponse[[response]][j] %in% c(paresp, countresp, 'coordinates')) addcovs <- unlist(datasetCovs)
            else addcovs <- NULL
            
          }
          else addcovs <- NULL
          
          if (!is.null(spatcovs)) {
            
            if (!is.null(speciesname)) covs <- paste0(speciesIn, '_', spatcovs)
            else covs <- spatcovs
            
          }
          else covs <- NULL
          
          if (!is.null(marks)) {
            
            if (pointsResponse[[response]][j] %in% c('coordinates', paresp, countresp)) {
              
              markspat <- NULL
              marksint <- NULL
              
            }
            else {
              
              if (!is.null(spatial)) spat <- NULL
              
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
          
          RHS <- c(covs, spat, int, addcovs, markspat, marksint, speciesspat) # temp
          
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
#' @param numtime Number of time increments included in the model.
#' @param speciesspatial Logical: Should the species be run with spatial fields.
#' @param offsetname Name of the offset column in the datasets.
#' @param copymodel List of the hyper parameters for the \code{copy} model.

dataOrganize$set('public', 'makeComponents', function(spatial, intercepts, 
                                                      datanames, marks, speciesname,
                                                      multinomnames, pointcovariates,
                                                      covariatenames,covariateclass,
                                                      marksspatial,
                                                      marksintercept,
                                                      temporalname,
                                                      speciesspatial,
                                                      numtime,temporalmodel,
                                                      offsetname,
                                                      copymodel) {
  ##Copy for marks fields???
  if (length(self$SpeciesInData) != 0) species <- unique(unlist(self$SpeciesInData))
  else species = NULL
  
  if (!is.null(spatial)) {
    
    if (spatial == 'shared') {
    
    if (!is.null(temporalname)) spat <- paste0('shared_spatial(main = coordinates, model = shared_field, group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
    else spat <- paste0('shared_spatial(main = coordinates, model = shared_field)')
    
    }
    else 
      if (spatial == 'copy') {
      
     mainName <- datanames[[1]]
     
     if (!is.null(temporalname)) spatMain <- paste0(mainName, '_spatial(main = coordinates, model = ', paste0(mainName,'_field'), ', group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
     else spatMain <-  paste0(mainName, '_spatial(main = coordinates, model = ', paste0(mainName,'_field'),')')
     spatCopy <-  paste0(datanames[datanames != mainName], '_spatial(main = coordinates, copy = \"', paste0(mainName,'_spatial'),'\", hyper = ', copymodel,')')
     spat <- c(spatMain, spatCopy)    
        
      }
    else {
      
    if (!is.null(temporalname)) spat <- paste0(datanames, '_spatial(main = coordinates, model = ', paste0(datanames,'_field'), ', group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
    else spat <-  paste0(datanames, '_spatial(main = coordinates, model =', paste0(datanames,'_field'),')')
      
    }
    
  } else spat <- NULL
  
  if (!is.null(species)) {
    
    if (speciesspatial) {
      
      #if (length(speciesspatial) > 0) {
      
      #speciesSpat <- paste0(species,'_spatial(main = coordinates, model = ',paste0(species,'_spdeModel'),')', collapse = ' + ')
      
      #} 
      #else {
      
      #ie we are assuming and INLA grouped model
      if (length(self$speciesIndex) != 0) {
        ##Change the species part to model = paste0(speciesname) ##where speciesname = species
         # but keep the speciesSpat framework for the temporal part of the model
        #speciesSpat <- paste0(speciesname, '_spatial(main = coordinates, model = speciesModel, group = ',speciesname,', ngroup = ', numspecies,')')
        speciesSpat <- paste0(species,'_spatial(main = coordinates, model = ',paste0(species,'_field)'))
        
      }
      #else speciesSpat <- paste0(species,'_spatial(main = coordinates, model = speciesModel)', collapse = ' + ') #change this to speciesModel
      
      #}
      
    }
    else speciesSpat <- NULL
    
  } 
  else speciesSpat <- NULL
  
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
      
      if (!is.null(temporalname)) marksSpat <- paste0(marks, '_spatial(main = coordinates, model = ', paste0(marks,'_field'), ', group = ', temporalname, ', ngroup = ', numtime,', control.group = ', temporalmodel,')')
      else marksSpat <- paste0(marks, '_spatial(main = coordinates, model = ', paste0(marks,'_field)'))
      
    }
    else marksSpat <- NULL
    
    if (marksintercept) {
      
      if (length(self$multinomVars) != 0)  marks_intercepts <- marks[!marks %in% self$multinomVars]
      else marks_intercepts <- marks
      
      marksInt <- paste0(marks_intercepts, '_intercept(1)')
      
    } else marksInt <- NULL
    
  }
  else {
    
    marksSpat <- NULL
    marksInt <- NULL
    
  }
  
  if (!is.null(covariatenames)) {
    
    if (!is.null(species)) {
      
      speciesCovs <- apply(expand.grid(paste0(species,'_'), covariatenames), MARGIN = 1, FUN = paste0,collapse='')
      speciesCovClass <- rep(covariateclass, each = length(species))
      covs <- paste0(speciesCovs, '(main = ', speciesCovs, ', model = \"',speciesCovClass,'\")') #, collapse = ' + '
      
    }
    else covs <- paste0(covariatenames, '(main = ', covariatenames, ', model = \"',covariateclass,'\")') # , collapse = ' + '
    
  } 
  else covs <- NULL
  
  if (!is.null(pointcovariates)) {
    
    covsPoints <- pointcovariates#paste0(pointcovariates, collapse = ' + ')
    
  } else covsPoints <- NULL
  
  if (!is.null(offsetname)) {
    
    if (offsetname %in% unlist(self$varsIn)) {
      
      offsetTerm <- paste0(offsetname,'(log(',offsetname,'), model = "offset")')
      
    } else offsetTerm <- NULL
    
  } else offsetTerm <- NULL
  
  #Intercepts: Dataset intercepts
  #            species intercepts
  #            marks intercepts # Am I doing this? Or are we doing dataset/species mark intercepts
  
  if (intercepts) {
    
    if (!is.null(species)) int <- paste0(species, '_intercept(1)')
    else int <- paste0(datanames, '_intercept(1)')
    
    if (!is.null(marks)) intMarks <- paste0(marks, '_intercept(1)')
    else intMarks <- NULL
    
  } 
  else {
    
    int <- NULL
    intMarks <- NULL
    
  } 
  
  if (!is.null(multinomnames)) {
    
    multinomVars <- paste0(multinomnames,'(main = ', multinomnames, ', model = "iid",constr = FALSE, fixed=TRUE)', collapse = ' + ')
    multinomPhi <- paste0(multinomnames,'_phi(main = ', multinomnames, '_phi, model = "iid", initial = -10, fixed = TRUE)', collapse = ' + ')
    
  }
  else {
    
    multinomVars <- NULL
    multinomPhi <- NULL
    
  }
  
  RHS <- c(spat, speciesSpat, marksSpat, covs, covsPoints, int, multinomVars, multinomPhi, marksInt, offsetTerm)
  
  RHS
  
})

#' @description Function to make the datasets into likelihoods.
#' @param mesh An inla.mesh object.
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
        
        Likelihoods[[Likindex]] <- inlabru::like(formula = self$Formulas[[dataset]][[species]][[process]][['LHS']],
                                                 include = self$Formulas[[dataset]][[species]][[process]][['RHS']],
                                                 data = self$Data[[dataset]][[species]], 
                                                 Ntrials = Ntrials,
                                                 mesh = mesh,
                                                 ips = IPS,
                                                 family = self$Family[[dataset]][process])
        
        if (is.null(names(self$Data[[dataset]])[species])) nameGive <- names(self$Data)[[dataset]]
        else nameGive <- names(self$Data[[dataset]])[species]
        
        names(Likelihoods)[[Likindex]] <- paste0(nameGive, '_', as.character(self$Formulas[[dataset]][[species]][[process]][['LHS']])[2])
        
      }
      
    }
    
  }
  
  self$Data <- Likelihoods
  
})

