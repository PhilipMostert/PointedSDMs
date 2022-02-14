#' @title runModel: function used to run the integrated model.
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

  
  if (!is.null(data$.__enclos_env__$private$spatcovsNames)) {
    
    spatCovs <- get(data$.__enclos_env__$private$spatcovsObj, 
                    envir = data$.__enclos_env__$private$spatcovsEnv)
    
    if (!inherits(spatCovs, 'Spatial')) spatCovs <- as(spatCovs, 'SpatialPixelsDataFrame')
    
    if (is.null(data$.__enclos_env__$private$speciesIn)) {
    
    for (name in data$.__enclos_env__$private$spatcovsNames) {

    pixelsDF <- sp::SpatialPixelsDataFrame(points = spatCovs@coords,
                                           data = data.frame(spatCovs@data[,name]),
                                           proj4string = data$.__enclos_env__$private$Projection)
    names(pixelsDF@data) <- name
    #Note that if species is non-null we would also need to paste the species name to this.
    assign(name, pixelsDF)
    
    }
  }
  else {
    
    for (species in unique(unlist(data$.__enclos_env__$private$speciesIn))) {
      
      for (name in data$.__enclos_env__$private$spatcovsNames) {

          pixelsDF <- sp::SpatialPixelsDataFrame(points = spatCovs@coords,
                                                 data = data.frame(spatCovs@data[,name]),
                                                 proj4string = data$.__enclos_env__$private$Projection)
          names(pixelsDF@data) <- paste0(species,'_',name)
          #Note that if species is non-null we would also need to paste the species name to this.
          assign(paste0(species,'_',name), pixelsDF)

        }
      }
    }
  }

  if (!is.null(data$.__enclos_env__$private$speciesIn)) {
   ## if speciesSPatial TRUE
   for (species in names(data$spatialFields$speciesFields)) {
 
     assign(paste0(species,'_field'), data$spatialFields$speciesFields[[species]])
     
   }
    
    speciesSpatial <- TRUE
    
  } else speciesSpatial <- FALSE
  
  if (data$.__enclos_env__$private$Spatial) {
    
  assign('spdeModel', data$spatialFields$sharedField)
  pointsSpatial <- TRUE
  
  } else pointsSpatial <- FALSE
  
  if (!is.null(data$.__enclos_env__$private$markNames)) {
    
    if (data$.__enclos_env__$private$marksSpatial) {
      
      for (mark in names(data$spatialFields$markFields)) {
        
        assign(paste0(mark,'_field'), data$spatialFields$markFields[[mark]])
      
      }
      
    marksSpatial <- TRUE  
 
    } else marksSpatial <- FALSE  
    
  } else marksSpatial <- FALSE  
  
  if (length(data$spatialFields$biasFields) != 0) {
    
    for (name in names(data$spatialFields$biasFields)) {
    
      assign(paste0(name,'_biasField'), data$spatialFields$biasFields[[name]])
      
    }
  }


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
  inlaModel[['biasData']] <- names(data$.__enclos_env__$private$biasField)
  inlaModel[['spatial']] <- list(points = pointsSpatial,
                                 species = speciesSpatial,
                                 marks = marksSpatial)
  
  class(inlaModel) <- c('bruSDM', class(inlaModel))
  
  inlaModel
  
  }

