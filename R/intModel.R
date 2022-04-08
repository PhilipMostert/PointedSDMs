#' @title intModel: integrated model data object.
#' 
#' @description This function is used to create a \code{dataSDM} object containing all the relevant data and meta-data to be used in the integrated distribution model.
#' 
#' @param ... The datasets to be used in the model. May come as either a \code{data.frame} or \code{SpatialPoints*}, or as a list of datasets inheriting the mentioned classes.
#' @param spatialCovariates The spatial covariates used in the model. These covariates must be measured at every location (pixel) in the study area, and must be a \code{Raster*} or \code{SpatialPixelsDataFrame} object.
#' @param Coordinates A vector of length 2 containing the names of the coordinate variables used in the model.
#' @param Projection Projection used for the data and spatial covariates. Must be of class \code{CRS}.
#' @param Boundary ##NOT USED YET.
#' @param Mesh An \code{inla.mesh} object.
#' @param IPS Integration points to be used in the model. Defaults to \code{NULL} which will create integration points from the \code{inla.mesh} object.
#' @param speciesSpatial Logical argument: should the species have their own spatial fields. Defaults to \code{TRUE}.
#' @param markNames A vector with the mark names to be included in the integrated model. These marks must be in the same data objects as the points.
#' @param markFamily A vector with the statistical families of the marks. Must be the same length as markNames, and the position of the mark in the vector \code{markName} is associated with the position of the family in \code{markFamily}. Defaults to \code{NULL} which assigns each mark as "Gaussian".
#' @param pointCovariates The non-spatial covariates to be included in the integrated model. These covariates must be included in the same data object as the points.
#' @param pointsIntercept Logical argument: should the points be modeled with intercepts. Defaults to \code{TRUE}.
#' @param marksIntercept Logical argument: should the marks be modeled with intercepts. Defaults to \code{TRUE}.
#' @param pointsSpatial Logical argument: should the points have a shared spatial field. Defaults to \code{TRUE}.
#' @param marksSpatial Logical argument: should the marks have their own spatial field. Defaults to \code{TRUE}.
#' @param responseCounts Name of the response variable in the counts datasets. This variable name needs to be standardized across all counts datasets used in the integrated model. Defaults to \code{'counts'}.
#' @param responsePA Name of the response variable in the presence absence datasets. This variable name needs to be standardized across all present absence datasets. Defaults to \code{'present'}.
#' @param trialsPA Name of the trials response variable for the presence absence datasets. Defaults to \code{NULL}.
#' @param trialsMarks Name of the trials response variable for the binomial marks. Defaults to \code{NULL}.
#' @param speciesName Name of the species variable name. Specifying this argument calculates covariate values for the individual species, as well as including a spatial group model for the species. Defaults to \code{NULL}
#' @param temporalName Name of the temporal variable in the model. This variable is required to be in all the datasets. Defaults to \code{NULL}.
#' @param temporalModel List of model specifications given to the control.group argument in the time effect component.
#' 
#' @return A dataSDM object. Use ?dataSDM to see all the slot functions related to this object.
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' #Specify whats needed
#' 
#' dataset1 <- data.frame(...)
#' dataset2 <- sp::SpatialPointsDataFrame(...)
#' 
#' mesh <- inla.mesh.2d(...)
#' spatiialCovariates <- raster::raster(...)
#' Projection <- CRS(...)
#' Coords <- c('Long', 'Lat')
#' 
#' #Run function
#' 
#' object <- intModel(dataset1, dataset2, Mesh = Mesh,
#'                   spatialCovariates = spatialCovariates, Coords = Coords,
#'                   Projection = Projection)
#'                  
#' #Return summary of data
#' 
#' object
#' 
#' }
#' 
#' @export


##Need to remove the field things here...
 # + everything in dataSDM and dataOrganize
 #Then go to unit tests
intModel <- function(..., spatialCovariates = NULL, Coordinates,
                   Projection, Boundary = NULL, Mesh, IPS = NULL,
                   speciesSpatial = TRUE,
                   markNames = NULL, markFamily = NULL,
                   pointCovariates = NULL, pointsIntercept = TRUE, marksIntercept = TRUE,
                   pointsSpatial = TRUE, marksSpatial = TRUE,
                   responseCounts = 'counts', responsePA = 'present', trialsPA = NULL,
                   trialsMarks = NULL, speciesName = NULL, temporalName = NULL, temporalModel = list(model = 'ar1')) {
  
  if (length(Coordinates) != 2) stop('Coordinates needs to be a vector of length 2 containing the coordinate names.')
  
  if (Coordinates[1] == Coordinates[2]) stop('Coordinates need to be unique values.')
  
  if (class(Projection) != 'CRS') stop('Projection needs to be a CRS object.')
  
  if (class(Mesh) != 'inla.mesh') stop('Mesh needs to be a inla.mesh object.')
  
  if (!pointsSpatial && !is.null(temporalName)) {
    
    pointsSpatial <- TRUE
    message('Setting pointsSpatial to TRUE since it is required for the temporalModel.')
    
  }
  
  #if (length(list(...)) > 0)  
  dataPoints <- list(...)
  
  #if (length(dataPoints) == 0) stop('Please provide data in the ... argument.')
  if (length(dataPoints) > 0) {
    
  datasetClass <- unlist(lapply(dataPoints, class))
##Need something more here if object is list: get names from list
  if (length(datasetClass) == 1 && datasetClass == "list") {

    dataNames <- NULL
    dataPoints <- unlist(dataPoints, recursive = FALSE)
    datasetClass <- lapply(dataPoints, class)
    dataList <- TRUE
    
  }
  else dataList <- FALSE

  if (!all(unlist(datasetClass) %in% c("SpatialPointsDataFrame", "SpatialPoints", "data.frame"))) stop("Datasets need to be either a SpatialPoints* object or a data frame.")
  
  if (dataList) {
    
    if (is.null(dataNames)) {
      
      initialnames <- setdiff(gsub("list[(]|[)]", "", as.character(match.call(expand.dots = TRUE))), 
                           gsub("list[(]|[)]", "", as.character(match.call(expand.dots = FALSE))))
      initialnames <- unlist(strsplit(x = initialnames, split = ", "))
      
      if (length(initialnames) != length(dataPoints)) {
        
        initialnames <- names(list(...)[[1]])
        
        if (any(is.null(initialnames))) {
        warning("Issues with naming the datasets from a list. Will create generic dataset names.", 
                immediate. = FALSE)
        
        cat("\n")
        initialnames <- paste0("dataset_", seq_len(length(datasetClass)))
        
      }
      
    }
    
    }
  }
  else {
    
    initialnames <- setdiff(as.character(match.call(expand.dots = TRUE)), 
                            as.character(match.call(expand.dots = FALSE)))
  }
  
  } else initialnames <- NULL
  
  if (!is.null(temporalName)) temporalModel <- deparse(temporalModel)
  
  bruData <- dataSDM$new(coordinates = Coordinates, projection = Projection,
                         Inlamesh = Mesh, initialnames = initialnames,
                         responsecounts = responseCounts,
                         responsepa = responsePA,
                         marksnames = markNames,
                         marksfamily = markFamily,
                         pointcovariates = pointCovariates,
                         trialspa = trialsPA,
                         trialsmarks = trialsMarks,
                         spatial = pointsSpatial,
                         marksspatial = marksSpatial,
                         speciesname = speciesName,
                         intercepts = pointsIntercept,
                         marksintercepts = marksIntercept,
                         spatialcovariates = spatialCovariates,
                         boundary = Boundary,
                         ips = IPS,
                         temporal = temporalName,
                         temporalmodel = temporalModel,
                         speciesspatial = speciesSpatial)
  
  if (length(list(...)) == 0) warning('No point data added. You can add data to this object with `$.addData()`.')
 
  else {
    
    if (is.null(responseCounts) || is.null(responsePA)) stop('One of responseCounts and responsePA are NULL. Both are required arguments.')
    
    #if (!is.null(responsePA) && is.null(trialsPA)) warning('Present absence response name given but trials name is NULL.')

   bruData$addData(..., responseCounts = responseCounts,
                    responsePA = responsePA, trialsPA = trialsPA,
                    markNames = markNames, pointCovariates = pointCovariates,
                    trialsMarks = trialsMarks, speciesName = speciesName,
                    temporalName = temporalName)
    
  }
  
  bruData
  
}