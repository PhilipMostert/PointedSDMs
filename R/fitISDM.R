#' @title \emph{fitISDM}: function used to run the integrated model.
#' 
#' @description This function takes a \code{intModel} object and produces an \code{inlabru} model object with additional lists and meta-data added.
#' 
#' @param data A intModel object to be used in the integrated model.
#' @param options A list of INLA options used in the model. Defaults to \code{list()}.
#' 
#' @import inlabru
#' @import stats
#' 
#' @return An inlabru model with additional lists containing some more metadata attached.
#' 
#' @examples 
#' 
#' \dontrun{
#'  
#'  if (requireNamespace('INLA')) {
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  proj <- CRS("+proj=longlat +ellps=WGS84")
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- fitISDM(organizedData, 
#'   options = list(control.inla = list(int.strategy = 'eb')))
#'    
#'   #Print summary of model
#'   modelRun
#'    
#'  }
#'}
#' 
#' @export

fitISDM <- function(data, options = list()) {
  
  if (!inherits(data, 'dataSDM')) stop('data needs to be a dataSDM object.')
  
  if (is.null(data$.__enclos_env__$private$INLAmesh)) stop('An inla.mesh object is required before any model is run.')
  
  data2ENV(data = data, env = environment())
  
  if (!is.null(data$.__enclos_env__$private$Spatial)) pointsSpatial <- data$.__enclos_env__$private$Spatial
  else pointsSpatial <- FALSE
  
  if (!is.null(names(data$spatialFields$markFields))) marksSpatial <- TRUE
  else marksSpatial <- FALSE
  
  if (!is.null(names(data$spatialFields$speciesFields))) speciesSpatial <- TRUE
  else speciesSpatial <- FALSE
  
  formula_terms <- unique(unlist(lapply(unlist(unlist(data$.__enclos_env__$private$Formulas, recursive = F), recursive = F), function(x) {
    
    if (is.null(x$RHS)) {
      
      if (length(attributes(terms(x$LHS))[['term.labels']]) != 1) x$RHS
      else {
        
        getTerms <- gsub('^.|[()]','',as.character(x$LHS)[3], perl = F)
        getTerms <- unlist(strsplit(getTerms, split = ' '))
        getTerms[!getTerms %in% c('+', '-', '/', '*', ':')]
        
      }
      
    }
    else x$RHS
    
  }))
  )
  
  comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
  
  #Will need to change this to say comp_terms %in% c(formula_terms, bias_terms)
  comp_keep <- comp_terms %in% c(formula_terms, 
                                 paste0(data$.__enclos_env__$private$dataSource, '_samplers_field'),
                                 paste0(data$.__enclos_env__$private$dataSource, 'samplers'))
  
  componentsJoint <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
  
  ##in case there are duplicates, will it cause an error??
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  if (!is.null(data$.__enclos_env__$private$temporalName)) {
    
    numTime <- length(unique(unlist(data$.__enclos_env__$private$temporalVars)))
    
    newIPS <- rep(list(data$.__enclos_env__$private$IPS), numTime)
    
    newIPS <- do.call(sp::rbind.SpatialPointsDataFrame, newIPS)
    
    newIPS@data[, data$.__enclos_env__$private$temporalName] <- rep((1:length(numTime)), each = nrow(data$.__enclos_env__$private$IPS@data))
    
    newIPS@proj4string <- data$.__enclos_env__$private$Projection
    
    data$.__enclos_env__$private$IPS <- newIPS
    
  }
  
  allLiks <- do.call(inlabru::like_list,
                     makeLhoods(data = data$.__enclos_env__$private$modelData,
                                formula = data$.__enclos_env__$private$Formulas,
                                family = data$.__enclos_env__$private$Family,
                                mesh = data$.__enclos_env__$private$INLAmesh,
                                ips = data$.__enclos_env__$private$IPS,
                                samplers = data$.__enclos_env__$private$Samplers,
                                paresp = data$.__enclos_env__$private$responsePA,
                                ntrialsvar = data$.__enclos_env__$private$trialsPA,
                                markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                                speciesname = data$.__enclos_env__$private$speciesName,
                                speciesindex = data$.__enclos_env__$private$speciesIndex))
  
  if (length(data$.__enclos_env__$private$biasData) > 0) {
    
    biasLikes <- list()
    
    for (bias in names(data$.__enclos_env__$private$biasData)) {
      
      biasLikes[[paste0(bias, '_samplers')]] <- inlabru::like(formula = coordinates ~ .,
                                                              family = 'cp',
                                                              data = do.call(sp::rbind.SpatialPointsDataFrame, data$.__enclos_env__$private$modelData[[bias]]),
                                                              samplers = data$.__enclos_env__$private$biasData[[bias]],
                                                              ips = data$.__enclos_env__$private$IPS,
                                                              domain = list(coordinates = data$.__enclos_env__$private$INLAmesh),
                                                              include = c(paste0(bias, '_samplers_field'), paste0(bias,'_samplers'), data$.__enclos_env__$private$spatcovsNames))
      
      
    }
    
    biasLikes <- do.call(inlabru::like_list, biasLikes)
    
    allLiks <- inlabru::like_list(do.call(append, list(allLiks, biasLikes)))
    
    data$.__enclos_env__$private$optionsINLA$control.family <- append(data$.__enclos_env__$private$optionsINLA$control.family,
                                                                      lapply(1:length(biasLikes), function(x) list(link = 'log')))
    
  }
  
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

