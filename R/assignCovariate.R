#' @description Function to assign environmental data to the observed data and integration points.
#' @param data Observation data or integration points.
#' @param covariateEnv Environment for the covariate layers.
#' @param covariateNames Names of the covariate layers.
#' @param timeVariable Name of the temporal variable.
#' @param timeData The temporal data associated with the observations.
#' @param speciesName Name of the species variable name.
#' @param IPS Logical: are covariates being assigned for integration points.
#' @param projection The required projection.
#' 
#' @import terra
#' @importFrom dplyr bind_rows
#' 

assignCovariate <- function(data, covariateEnv, covariateNames,
                            timeVariable = NULL, timeData = NULL, 
                            speciesName = NULL, IPS = FALSE, projection) {
  
  dataColNames <- lapply(data, function(x) lapply(x, names))
  datasetNames <- names(data)
  #will this work?
  #allCovs <- covariateNames
  ##FIX THIS IN INITIAL
  
  #biasCovs <- if(is.null(private$biasFormula)) NULL else labels(terms(private$biasFormula))
  #modelCovs <- allCovs[! allCovs %in% biasCovs]
  
  fullGeom <- bind_rows(lapply(data, bind_rows))
  ## do this per layer if list
  #Maybe easiest to do if temporal
  
  CovsGeom <- get('spatialcovariates',envir = covariateEnv)
  
  if (!is.null(timeVariable) && inherits(CovsGeom, 'list')) {
    
    if (is.null(names(CovsGeom))) names(CovsGeom) <- covariateNames
    
    uniTempVars <- unique(unlist(timeData))
    
    fullGeomCovs <- data.frame(matrix(NA, nrow = nrow(fullGeom),
                                      ncol = length(names(CovsGeom))))
    
    names(fullGeomCovs) <- names(CovsGeom)
    
    for (cov in covariateNames) {
      
      if (terra::nlyr(CovsGeom[[cov]]) == 1) {
        
        fullGeomCovs[[cov]] <- terra::extract(terra::project(CovsGeom[[cov]],
                                                             projection),
                                              fullGeom, ID = FALSE)[,1]
        
        if (any(is.na(fullGeomCovs[[cov]]))) fullGeomCovs[[cov]][is.na(fullGeomCovs[[cov]])] <- nearestValue(matrix(st_coordinates(fullGeom[is.na(fullGeomCovs[[cov]]),])[,c("X","Y")], ncol = 2), 
                                                                                                             terra::project(CovsGeom[[cov]], 
                                                                                                                            projection)[cov])
      }
      
      else {
        
        if (length(names(CovsGeom[[cov]])) == length(uniTempVars)) {
          
          if (!all(names(CovsGeom[[cov]]) %in% uniTempVars)) {
            
            warning(paste('Not all names for covariate', cov, 'are the same as the temporal variables for the observations. Will assume the covariate is ordered chronologically.'))
            names(CovsGeom[[cov]]) <- uniTempVars
            
          }
          
          for (timeP in uniTempVars) {
            
            ##Fix with using the layer argument in terra::extract ie layer = fullGeom[[timeVariable]]
             #But then figure out how to deal with NA values
            
            timeIn <- unlist(timeData) == timeP
            
            fullGeomCovs[timeIn, cov] <- terra::extract(terra::project(terra::subset(CovsGeom[[cov]], timeP),
                                                                       projection), 
                                                        fullGeom[timeIn,], 
                                                        ID = FALSE)[,1]
            
            if (any(is.na(fullGeomCovs[timeIn, cov]))) fullGeomCovs[timeIn, cov][is.na(fullGeomCovs[timeIn, cov])] <- nearestValue(matrix(st_coordinates(fullGeom[timeIn & is.na(fullGeomCovs[[cov]]),])[,c("X","Y")], 
                                                                                                                                          ncol = 2), 
                                                                                                                                   terra::project(terra::subset(CovsGeom[[cov]], timeP), 
                                                                                                                                                  projection))
            
          }
          
        }
        else {
          
          warning(paste0('Covariate ', cov, ' is a temporal covariate which does not contain the same temporal variables as the observation data, or has a length identical to the number of temporal locations. Will extract data from the first layer only.'))
          
          fullGeomCovs[[cov]] <- terra::extract(terra::project(CovsGeom[[cov]][[1]], projection), fullGeom, ID = FALSE)[,1]
          
          if (any(is.na(fullGeomCovs[[cov]]))) fullGeomCovs[[cov]][is.na(fullGeomCovs[[cov]])] <- nearestValue(matrix(st_coordinates(fullGeom[is.na(fullGeomCovs[[cov]]),])[,c("X","Y")], ncol = 2), 
                                                                                                               terra::project(CovsGeom[[cov]][[1]], 
                                                                                                                              projection)[1])
          
        }
        
      }
      
      
    }
    
  }
  else {
    
    if (inherits(CovsGeom, 'list')) {
      
      names(CovsGeom) <- NULL
      CovsGeom <- do.call(c, CovsGeom)
      if (is.null(names(CovsGeom))) names(CovsGeom) <- covariateNames
      
    }
    
    fullGeomCovs <- terra::extract(terra::project(CovsGeom, 
                                                  projection), 
                                   fullGeom, ID = FALSE)
    
  }
  
  ##Redo this for temporalcovs
  
  if(any(is.na(fullGeomCovs))){
    naRows <- lapply(fullGeomCovs, function(x) which(is.na(x)))  # identify missing rows
    naCovs <- names(naRows)[sapply(naRows, length) > 0]  # identify covs with missing data
    for(cov in naCovs){  # fill missing values for rows/covs using nearest neighbour 
      fullGeomCovs[naRows[[cov]], cov] <- 
        nearestValue(matrix(st_coordinates(fullGeom[naRows[[cov]],])[,c("X","Y")], ncol = 2), 
                     terra::project(CovsGeom, 
                                    projection)[cov])
      
    }
  }
  
  fullGeom <- cbind(fullGeom, fullGeomCovs) 

  if (!IPS) {
  # split by dataset 
  splitVar <- split(fullGeom, factor(fullGeom$._dataset_index_var_., levels = unique(fullGeom$._dataset_index_var_.)))
  names(splitVar) <- datasetNames
  # 
  
  ##if data colnames is a list in a list
  
  if (is.null(speciesName)) fullGeom <- lapply(splitVar, FUN = function(x) {
    
    ds <- names(dataColNames)[x$._dataset_index_var_.[1]]
    colsKeep <- dataColNames[[ds]][[1]]
    x <- list(x[,c(colsKeep, covariateNames)])
    names(x) <- names(dataColNames[[ds]])
    x
    
  })
  else fullGeom <- lapply(splitVar, function(x){
    
    xSplit <- split(x, factor(x[[paste0(speciesName,'INDEX_VAR')]], levels = unique(x[[paste0(speciesName,'INDEX_VAR')]])))
    ds <- names(dataColNames)[x$._dataset_index_var_.[1]]
    
    xSplit <- lapply(xSplit, function(x2){
      
      colsKeep <- dataColNames[[ds]][[paste0(ds, "_", x2[[paste0(speciesName,'INDEX_VAR')]][1])]]
      x2 <- x2[,c(colsKeep, covariateNames)]
      
      return(x2)
    })
    
    names(xSplit) <- names(dataColNames[[ds]])
    xSplit
  })
  }

  
  fullGeom
  
  
}