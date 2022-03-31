#' @title runModel: function used to run the integrated model.
#' 
#' @description This function takes a \code{bruSDM} object and produces an \code{inlabru} model object with additional lists and meta-data added.
#' 
#' @param data A bruSDM data file to be used in the integrated model.
#' @param options A list of INLA options used in the model. Defaults to \code{list()}.
#' 
#' @return An inlabru model with additional information attached.
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' #Create dataSDM object
#' 
#' dataObject <- bruSDM(...)
#' 
#' #Run the joint model
#' 
#' joint_model <- runModel(dataObject, 
#'                         options = list(control.inla = list(int.strategy = 'eb')))
#' 
#' #Print summary of model
#' 
#' joint_model
#' 
#' }
#' 
#' @export

runModel <- function(data, options = list()) {
  
  if (!inherits(data, 'dataSDM')) stop('data needs to be a dataSDM object.')

  data2ENV(data = data, env = environment())
  
  if (data$.__enclos_env__$private$Spatial) pointsSpatial <- TRUE
  else pointsSpatial <- FALSE
  
  if (!is.null(names(data$spatialFields$markFields))) marksSpatial <- TRUE
  else marksSpatial <- FALSE
  
  if (!is.null(names(data$spatialFields$speciesFields))) speciesSpatial <- TRUE
  else speciesSpatial <- FALSE

  formula_terms <- unique(unlist(lapply(unlist(unlist(dataObj$.__enclos_env__$private$Formulas, recursive = F), recursive = F), function(x) {
    
    if (is.null(x$RHS))  attributes(terms(x$LHS))[['term.labels']]
    else x$RHS
    
  }))
  )
  
  comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
  
  comp_keep <- comp_terms %in% formula_terms
  
  componentsJoint <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))

  ##in case there are duplicates, will it cause an error??
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  allLiks <- do.call(inlabru::like_list,
                     makeLhoods(data = data$.__enclos_env__$private$modelData,
                     formula = data$.__enclos_env__$private$Formulas,
                     family = data$.__enclos_env__$private$Family,
                     mesh = data$.__enclos_env__$private$INLAmesh,
                     ips = data$.__enclos_env__$private$IPS,
                     paresp = data$.__enclos_env__$private$responsePA,
                     ntrialsvar = data$.__enclos_env__$private$trialsPA,
                     markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                     speciesname = data$.__enclos_env__$private$speciesName,
                     speciesindex = data$.__enclos_env__$private$speciesIndex))
  
  optionsJoint <- append(data$.__enclos_env__$private$optionsINLA, options)
  
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

