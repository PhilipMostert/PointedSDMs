#' @title runModel: function used to run the integrated model.
#' 
#' @description This function takes a \code{bruSDM} object and produces an \code{inlabru} model object with additional lists and meta-data added.
#' 
#' @param data A bruSDM data file to be used in the integrated model.
#' @param options A list of INLA options used in the model. Defaults to \code{list()}.
#' 
#' @export

runModel <- function(data, options = list()) {
  
  if (!inherits(data, 'dataSDM')) stop('data needs to be a dataSDM object.')

  ## Need to assign all relevant variables into this environment
  #$.__enclos_env__$private$
  
  #Things needed:
   #proj ## why do I need proj ...
   #spatialmodel (for points and for species) #Do we need to get the environment for these??
   #spatialcovariates
   #components
  
   ## Things to do here:
   # Assign all relevant variables to this environment
    # Most notably the spatial covariates: get each covariate as its own SPpixeldataframe

  data2ENV(data = data, env = environment())
  #stop(return(ls()))
  ## Get all components in formula; get all components but without the ()
   # if not in formulas then remove from components
  
  formula_terms <- unique(unlist(lapply(data$.__enclos_env__$private$modelData, function(x) {
    
    if (is.null(x$include_components))  attributes(terms(x$formula))[['term.labels']]
    else x$include_components
    
  })))
  
  comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
  
  comp_keep <- comp_terms %in% formula_terms
  
  componentsJoint <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))

  ##in case there are duplicates, will it cause an error??
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  allLiks <- do.call(like_list, data$.__enclos_env__$private$modelData)

  ##For now just set bru_max_iter = 1
  
  optionsJoint <- append(data$.__enclos_env__$private$optionsINLA, options)
  
  ##For now set bru_max_iter to 1
  optionsJoint$bru_max_iter <- 1

  inlaModel <- inlabru::bru(components = componentsJoint,
                                 allLiks, options = optionsJoint)

  if (length(data$.__enclos_env__$private$multinomVars) != 0) {
    
    for (name in names(data$.__enclos_env__$private$multinomIndex)) {
   
      inlaModel$summary.random[[name]]['ID'] <- unique(data$.__enclos_env__$private$multinomIndex[[name]])

    }
    
  }

  inlaModel[['componentsJoint']] <- componentsJoint
  inlaModel[['optionsJoint']] <- optionsJoint
  inlaModel[['source']] <- as.vector(unlist(data$.__enclos_env__$private$dataSource))
  inlaModel[['spatCovs']] <- list(name = data$.__enclos_env__$private$spatcovsNames,
                                  class = data$.__enclos_env__$private$ptcovsClass)
  inlaModel[['species']] <- list(speciesIn = data$.__enclos_env__$private$speciesIn,
                                 speciesVar = data$.__enclos_env__$private$speciesName)
  inlaModel[['dataType']] <- c(na.omit(data$.__enclos_env__$private$printSummary$Type),
                               na.omit(unlist(unname(data$.__enclos_env__$private$printSummary$marksType))))
  inlaModel[['marks']] <- list(marksIn = data$.__enclos_env__$private$printSummary$Marks,
                                      multinomVars = data$.__enclos_env__$private$multinomVars)
  inlaModel[['biasData']] <- names(data$spatialFields$biasField)
  inlaModel[['spatial']] <- list(points = pointsSpatial,
                                 species = speciesSpatial,
                                 marks = marksSpatial)
  inlaModel[['temporal']] <- list(temporalIn = data$.__enclos_env__$private$temporalVars, # I think ... do we need all these vars??? Can we just use unique unlist( ... )
                                  temporalVar = data$.__enclos_env__$private$temporalName)
  
  class(inlaModel) <- c('bruSDM', class(inlaModel))
  
  inlaModel
  
  }

