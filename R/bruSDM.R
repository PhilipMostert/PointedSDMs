#' bruSDM function to ...
#' 
#' @param ...
#' @param spatialCovariates
#' @param Coordinates
#' @param Projection
#' @param Boundary
#' @param INLAmesh
#' @param pointsField
#' @param speciesField
#' @param markNames
#' @param pointCovariates
#' @param responseCounts
#' @param responsePA
#' @param trialsPA
#' @param trialsMarks
#' @param speciesName
#' 
#' @export


##Is markFamily even being used??
bruSDM <- function(..., spatialCovariates = NULL, Coordinates,
                   Projection, Boundary = NULL, INLAmesh, IPS = NULL,
                   pointsField = NULL, speciesField = NULL,
                   markNames = NULL, markFamily = NULL, marksField = NULL,
                   pointCovariates = NULL, Intercepts = TRUE, Spatial = TRUE, 
                   responseCounts = NULL, responsePA = NULL, trialsPA = NULL,
                   trialsMarks = NULL, speciesName = NULL) {
  
  if (length(Coordinates) != 2) stop('Coordinates needs to be a vector of length 2 containing the coordinate names.')
  
  if (Coordinates[1] == Coordinates[2]) stop('Coordinates need to be unique values.')
  
  if (class(Projection) != 'CRS') stop('Projection needs to be a CRS object.')
  
  if (class(INLAmesh) != 'inla.mesh') stop('INLAmesh needs to be a inla.mesh object.')
  
  #if (length(list(...)) > 0)  
  dataPoints <- list(...)
  
  #if (length(dataPoints) == 0) stop('Please provide data in the ... argument.')
  if (length(dataPoints) > 0) {
    
  datasetClass <- lapply(dataPoints, class)
  
  if (length(datasetClass) == 1 && class(datasetClass) == "list") {
    
    dataNames <- NULL
    dataPoints <- unlist(dataPoints)
    datasetClass <- lapply(dataPoints, class)
    dataList <- TRUE
    
  }
  else dataList <- FALSE
  
  if (any(!unlist(datasetClass) %in% c("SpatialPointsDataFrame", "SpatialPoints", "data.frame"))) stop("Datasets need to be either a SpatialPoints* object or a data frame.")
  
  if (dataList) {
    
    if (is.null(dataNames)) {
      
      initialnames <- setdiff(gsub("list[(]|[)]", "", as.character(match.call(expand.dots = TRUE))), 
                           gsub("list[(]|[)]", "", as.character(match.call(expand.dots = FALSE))))
      initialnames <- unlist(strsplit(x = initialnames, split = ", "))
      
      if (length(initialnames) != length(dataPoints)) {
        
        warning("Issues with naming the datasets from a list. Will create generic dataset names.", 
                immediate. = FALSE)
        
        cat("\n")
        initialnames <- paste0("dataset_", seq_len(length(datasets)))
        
      }
      
    }
    
  }
  else {
    
    initialnames <- setdiff(as.character(match.call(expand.dots = TRUE)), 
                            as.character(match.call(expand.dots = FALSE)))
  }
  
  } else initialnames <- NULL
  
  bruData <- dataSDM$new(coordinates = Coordinates, projection = Projection,
                         Inlamesh = INLAmesh, initialnames = initialnames,
                         responsecounts = responseCounts,
                         responsepa = responsePA,
                         marksnames = markNames,
                         marksfamily = markFamily,
                         pointcovariates = pointCovariates,
                         trialspa = trialsPA,
                         trialsmarks = trialsMarks,
                         spatial = Spatial,
                         speciesname = speciesName,
                         intercepts = Intercepts,
                         spatialcovariates = spatialCovariates,
                         boundary = Boundary,
                         ips = IPS)
  
  if (length(list(...)) == 0) warning('No point data added. You can add data to this object with `$.addData()`.')
 
  else {
    
    if (is.null(responseCounts) || is.null(responsePA)) stop('One of responseCounts and responsePA are NULL. Both are required arguments.')
    
    #if (!is.null(responsePA) && is.null(trialsPA)) warning('Present absence response name given but trials name is NULL.')

   bruData$addData(..., responseCounts = responseCounts,
                    responsePA = responsePA, trialsPA = trialsPA,
                    markNames = markNames, pointCovariates = pointCovariates,
                    trialsMarks = trialsMarks, speciesName = speciesName,
                    pointsField = pointsField, speciesField = speciesField,
                    marksField = marksField)
    
  }
  
  bruData
  
}