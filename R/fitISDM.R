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
#'  proj <- "+proj=longlat +ellps=WGS84"
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
  
  if (!inherits(data, 'dataSDM') && !inherits(data, 'specifySpecies') && !inherits(data, 'specifyISDM') &&
      !inherits(data, 'specifyMarks')) stop('data needs to be either a specifySpecies, specifyISDM or specifyMarks object.')
  
  if (is.null(data$.__enclos_env__$private$INLAmesh)) stop('An fm_mesh_2d object is required before any model is run.')
  
  data2ENV(data = data, env = environment())
  
  if (!is.null(data$.__enclos_env__$private$Spatial)) pointsSpatial <- data$.__enclos_env__$private$Spatial
  else pointsSpatial <- FALSE
  
  if (!is.null(names(data$spatialFields$markFields))) marksSpatial <- TRUE
  else marksSpatial <- FALSE
  
  if(!is.null(data$.__enclos_env__$private$speciesSpatial)) speciesSpatial <- data$.__enclos_env__$private$speciesSpatial
  else speciesSpatial <- FALSE
  
  formula_terms <- unique(unlist(lapply(unlist(unlist(data$.__enclos_env__$private$Formulas, recursive = F), recursive = F), function(x) {
    
    if (is.null(x$RHS)) {
      
      if (length(attributes(terms(x$LHS))[['term.labels']]) != 1) attributes(terms(x$LHS))[['term.labels']]
      else {
        
        getTerms <- unlist(strsplit(as.character(x$LHS)[3], split = ' '))
        getTerms <- gsub('.*\\(', '',getTerms)
        getTerms <- gsub('\\)', '', getTerms)
        getTerms[!getTerms %in% c('+', '-', '/', '*', ':')]
        
      }
      
    }
    else x$RHS
    
  }))
  )
  
  if (!is.null(data$.__enclos_env__$private$temporalName)) {
    
    if (!data$.__enclos_env__$private$temporalName %in% names(data$.__enclos_env__$private$IPS)) {
    
    numTime <- length(unique(unlist(data$.__enclos_env__$private$temporalVars)))
    
    newIPS <- rep(list(data$.__enclos_env__$private$IPS), numTime)
    
    newIPS <- do.call(rbind, newIPS)
    
    newIPS[, data$.__enclos_env__$private$temporalName] <- rep(1:numTime, each = nrow(data$.__enclos_env__$private$IPS))
    
    newIPS <- st_transform(newIPS, data$.__enclos_env__$private$Projection)
    
    data$.__enclos_env__$private$IPS <- newIPS
    
    }
  }
  
  comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
  
  #Will need to change this to say comp_terms %in% c(formula_terms, bias_terms)
  comp_keep <- comp_terms %in% c(formula_terms, 
                                 paste0(data$.__enclos_env__$private$dataSource, '_samplers_field'),
                                 paste0(data$.__enclos_env__$private$dataSource, 'samplers'))
  
  componentsJoint <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
  
  ##in case there are duplicates, will it cause an error??
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  
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
                                speciesindex = data$.__enclos_env__$private$speciesIndex,
                                pointcovs = c(data$.__enclos_env__$private$pointCovariates, data$.__enclos_env__$private$Offset)))
  
  if (length(data$.__enclos_env__$private$biasData) > 0) {
    
    biasLikes <- list()
    
    for (bias in names(data$.__enclos_env__$private$biasData)) {
      
      biasLikes[[paste0(bias, '_samplers')]] <- inlabru::like(formula = coordinates ~ .,
                                                              family = 'cp',
                                                              data = do.call(sp::rbind.SpatialPointsDataFrame, data$.__enclos_env__$private$modelData[[bias]]),
                                                              samplers = data$.__enclos_env__$private$biasData[[bias]],
                                                              ips = data$.__enclos_env__$private$IPS,
                                                              domain = list(coordinates = data$.__enclos_env__$private$INLAmesh),
                                                              include = c(paste0(bias, '_samplers_field'), paste0(bias,'_samplers'), data$.__enclos_env__$private$spatcovsNames),
                                                              tag = paste0(bias, '_samplers'))
      
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
  
  if (is.null(data$.__enclos_env__$private$speciesIntercepts)) data$.__enclos_env__$private$speciesIntercepts <- FALSE
  if (data$.__enclos_env__$private$speciesIntercepts) row.names(inlaModel$summary.random[[paste0(data$.__enclos_env__$private$speciesName, '_intercepts')]]) <- data$.__enclos_env__$private$speciesTable[['species']]
  
  
  inlaModel[['componentsJoint']] <- componentsJoint
  inlaModel[['optionsJoint']] <- optionsJoint
  inlaModel[['source']] <- as.vector(unlist(data$.__enclos_env__$private$dataSource))
  inlaModel[['spatCovs']] <- list(name = data$.__enclos_env__$private$spatcovsNames,
                                  class = data$.__enclos_env__$private$ptcovsClass,
                                  env = data$.__enclos_env__$private$spatcovsEnv,
                                  covariateFormula = data$.__enclos_env__$private$covariateFormula,
                                  biasFormula = data$.__enclos_env__$private$biasFormula)
  inlaModel[['species']] <- list(speciesIn = data$.__enclos_env__$private$speciesIn,
                                 speciesVar = data$.__enclos_env__$private$speciesName,
                                 speciesEffects = list(Intercepts = data$.__enclos_env__$private$speciesIntercepts,
                                                       Environmental = data$.__enclos_env__$private$speciesEnvironment),
                                 speciesTable = data$.__enclos_env__$private$speciesTable)
  inlaModel[['dataType']] <- c(na.omit(data$.__enclos_env__$private$printSummary$Type),
                               na.omit(unlist(unname(data$.__enclos_env__$private$printSummary$marksType))))
  inlaModel[['marks']] <- list(marksIn = data$.__enclos_env__$private$printSummary$Marks,
                               multinomVars = data$.__enclos_env__$private$multinomVars)
  inlaModel[['biasData']] <- list(Fields = names(data$spatialFields$biasField), Comps = data$.__enclos_env__$private$biasFormula,
                                  Copy = data$.__enclos_env__$private$biasCopy)
  inlaModel[['spatial']] <- list(points = pointsSpatial,
                                 species = speciesSpatial,
                                 marks = marksSpatial)
  inlaModel[['temporal']] <- list(temporalIn = data$.__enclos_env__$private$temporalVars, # I think ... do we need all these vars??? Can we just use unique unlist( ... )
                                  temporalVar = data$.__enclos_env__$private$temporalName)
  
  #Do a switch here to assign the correct class
  
  
  class(inlaModel) <- switch(class(data)[1],
                             specifyISDM = c('modISDM', class(inlaModel)),
                             specifySpecies = c('modSpecies', class(inlaModel)),
                             specifyMarks = c('modMarks', class(inlaModel)),
                             dataSDM = c('bruSDM', class(inlaModel)))
  
  inlaModel
  
}

