#' @description Internal function: used to assign objects specified in bruSDM to the dataSDM/blockedCV function environments.
#' 
#' @param data bruSDM data file to be used in the integrated model.
#' @param env Environment where the objects should be assigned.
#' 
#' @export

data2ENV <- function(data, env) {
  
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
        assign(name, pixelsDF, envir = env)
        
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
          assign(paste0(species,'_',name), pixelsDF, envir = env)
          
        }
      }
    }
  }
  
  if (!is.null(data$.__enclos_env__$private$speciesIn)) {
    
    if (data$.__enclos_env__$private$speciesSpatial) {
      
      for (species in names(data$spatialFields$speciesFields)) {
        
        assign(paste0(species,'_field'), data$spatialFields$speciesFields[[species]], envir = env)
        
      }
      
      speciesSpatial <- TRUE
      
    }  else speciesSpatial <- FALSE
    
  } else speciesSpatial <- FALSE
  
  if (data$.__enclos_env__$private$Spatial) {
    
    assign('shared_field', data$spatialFields$sharedField, envir = env)
    pointsSpatial <- TRUE
    
  } else pointsSpatial <- FALSE
  
  if (!is.null(data$.__enclos_env__$private$markNames)) {
    
    if (data$.__enclos_env__$private$marksSpatial) {
      
      for (mark in names(data$spatialFields$markFields)) {
        
        assign(paste0(mark,'_field'), data$spatialFields$markFields[[mark]], envir = env)
        
      }
      
      marksSpatial <- TRUE  
      
    } else marksSpatial <- FALSE  
    
  } else marksSpatial <- FALSE  
  
  if (length(data$spatialFields$biasFields) != 0) {
    
    for (name in names(data$spatialFields$biasFields)) {
      
      assign(paste0(name,'_bias_field'), data$spatialFields$biasFields[[name]], envir = env)
      
    }
  }
  
  
}