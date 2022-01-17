#' @title runModel: function used to run the integrated model.
#' @param data A bruSDM data file to be used in the integrated model.
#' @param options A list of INLA options used in the model. ADD
#' 
#' @export

runModel <- function(data) {
  
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

  if (!is.null(data$.__enclose_env__$private$speciesIn)) {

    #if species model is non null
    #for (species in unique(unlist(data$.__enclose_env__$private$speciesIn))) {
      
    #  assign(paste0(species,'_spatial'), the spde model)
     
    #}
    
  }
  
  if (data$.__enclos_env__$private$Spatial) {
    
    assign('spdeModel', data$.__enclos_env__$private$pointsField)
   
    if (!is.null(data$.__enclos_env__$private$markNames)) assign('markModel', data$.__enclos_env__$private$markField)
     
    if (!is.null(data$.__enclos_env__$private$speciesName)) assign('speciesModel', data$.__enclos_env__$private$speciesField)
    
  }
  
  if (!is.null(data$.__enclos_env__$private$biasField)) assign('biasField', data$.__enclos_env__$private$biasField)

  ##Do the same for the species but run through a for loop and use assign.

  componentsJoint <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components, collapse = ' + ')))

  ##in case there are duplicates, will it cause an error??
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  allLiks <- do.call(like_list, data$.__enclos_env__$private$modelData)

  ##For now just set bru_max_iter = 1
  
  optionsJoint <- data$.__enclos_env__$private$optionsINLA
  #optionsJoint$bru_max_iter <- 1

  inlaModel <- inlabru::bru(components = componentsJoint,
                            allLiks, options = optionsJoint)

  if (length(data$.__enclos_env__$private$multinomVars) != 0) {
    
    for (name in names(data$.__enclos_env__$private$multinomIndex)) {
   
      inlaModel$summary.random[[name]]['ID'] <- unique(data$.__enclos_env__$private$multinomIndex[[name]])

    }
    
  }

  inlaModel[['componentsJoint']] <- componentsJoint
  inlaModel[['optionsJoint']] <- optionsJoint
  inlaModel[['source']] <- data$.__enclos_env__$private$dataSource
  inlaModel[['spatCovs']] <- list(name = data$.__enclos_env__$private$spatcovsNames,
                                  class = data$.__enclos_env__$private$ptcovsClass)
  inlaModel[['species']] <- list(speciesIn = data$.__enclos_env__$private$speciesIn,
                                 speciesVar = data$.__enclos_env__$private$speciesName)
  inlaModel[['dataType']] <- na.omit(c(data$.__enclos_env__$private$printSummary$Type,
                               unlist(unname(data$.__enclos_env__$private$printSummary$marksType))))
  inlaModel[['multinomVars']] <- data$.__enclos_env__$private$multinomVars
  
  class(inlaModel) <- c('bruSDM', class(inlaModel)) # ie I'm writing an s3 class for this?
  
  inlaModel
  
  }

