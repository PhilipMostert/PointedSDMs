#' @title \emph{intModel}: Function used to initialize the integrated species distribution model.
#' 
#' @description This function is used to create an object containing all the data, metadata and relevant components required for the integrated species distribution model and \pkg{INLA} to work.
#' As a result, the arguments associated with this function are predominantly related to describing variable names within the datasets that are relevant, and arguments related to what terms should be included in the formula for the integrated model. The output of this function is an \code{R6} object, and so there are a variety of public methods within the output of this function which can be used to further specify the model (see \code{?dataSDM} for a comprehensive description of these public methods).
#' 
#' @param ... The datasets to be used in the model. May come as either \code{data.frame} or \code{SpatialPoints*} objects, or as a list of objects with these classes. The classes of the datasets do not necessarily need to be standardized, however the variable names within them often have to be.
#' @param spatialCovariates The spatial covariates used in the model. These covariates must be measured at every location (pixel) in the study area, and must be a \code{Raster*} or \code{SpatialPixelsDataFrame} object. Can be either \code{numeric} or \code{factor}\\code{character} data.
#' @param Coordinates A vector of length 2 containing the names (class \code{character}) of the coordinate variables used in the model.
#' @param Projection The coordinate reference system used by both the spatial points and spatial covariates. Must be of class \code{\link[sp]{CRS}} (see \code{CRS} from the \pkg{sp} package for more details).
#' @param Boundary ##NOT USED YET.
#' @param Mesh An \code{inla.mesh} object required for the spatial random fields and the integration points in the model (see \code{\link[INLA]{inla.mesh.2d}} from the \pkg{INLA} package for more details). 
#' @param IPS The integration points to be used in the model (that is, the points on the map where the intensity of the model is calculated). See \code{\link[inlabru]{ipoints}} from the \pkg{inlabru} package for more details regarding these points; however defaults to \code{NULL} which will create integration points from the \code{Mesh} object.
#' @param speciesSpatial Logical argument: should the species have their own individual spatial fields in the integrated model. Defaults to \code{TRUE}.
#' @param markNames A vector with the mark names (class \code{character}) to be included in the integrated model. Marks are variables which are used to describe the individual points in the model (for example, in the field of ecology the size of the species or its feeding type could be considered). Defaults to \code{NULL}, however if this argument is non-\code{NULL}, the model run will become a marked point process. The marks must be included in the same data object as the points.
#' @param markFamily A vector with the statistical families (class \code{character}) assumed for the marks. Must be the same length as markNames, and the position of the mark in the vector \code{markName} is associated with the position of the family in \code{markFamily}. Defaults to \code{NULL} which assigns each mark as "Gaussian".
#' @param pointCovariates The non-spatial covariates to be included in the integrated model (for example, in the field of ecology the distance to the nearest road or time spent sampling could be considered). These covariates must be included in the same data object as the points.
#' @param Offset Name of the offset variable (class \code{character}) in the datasets. Defaults to \code{NULL}; if the argument is non-\code{NULL}, the variable name needs to be standardized across datasets (but does not need to be included in all datasets). The offset variable will be transformed onto the log-scale in the integrated model.
#' @param pointsIntercept Logical argument: should the points be modeled with intercepts. Defaults to \code{TRUE}.
#' @param marksIntercept Logical argument: should the marks be modeled with intercepts. Defaults to \code{TRUE}.
#' @param pointsSpatial Argument to determine whether the spatial field is shared between the datasets, or if each dataset has its own unique spatial field. May take on the values: \code{"shared"}, \code{"individual"} or \code{NULL} if no spatial field is required for the model. Defaults to \code{"shared"}.
#' @param marksSpatial Logical argument: should the marks have their own spatial field. Defaults to \code{TRUE}.
#' @param responseCounts Name of the response variable in the counts/abundance datasets. This variable name needs to be standardized across all counts datasets used in the integrated model. Defaults to \code{'counts'}.
#' @param responsePA Name of the response variable (class \code{character}) in the presence absence/detection non-detection datasets. This variable name needs to be standardized across all present absence datasets. Defaults to \code{'present'}.
#' @param trialsPA Name of the trials response variable (class \code{character}) for the presence absence datasets. Defaults to \code{NULL}.
#' @param trialsMarks Name of the trials response variable (class \code{character}) for the binomial marks (if included). Defaults to \code{NULL}.
#' @param speciesName Name of the species variable name (class \code{character}). Specifying this argument turns the model into a stacked species distribution model, and calculates covariate values for the individual species, as well as a species group model in the shared spatial field. Defaults to \code{NULL}
#' @param temporalName Name of the temporal variable (class \code{character}) in the model. This variable is required to be in all the datasets. Defaults to \code{NULL}.
#' @param temporalModel List of model specifications given to the control.group argument in the time effect component. Defaults to \code{list(model = 'ar1')}; see \code{\link[INLA]{control.group}} from the \pkg{INLA} package for more details.
#' 
#' @return A \code{\link{dataSDM}} object (class \code{R6}). Use \code{?dataSDM} to get a comprehensive description of the slot functions associated with this object.
#' 
#' @note The idea with this function is to describe the full model: that is, all the covariates and spatial effects will appear in all the formulas for the datasets and species.
#' If some of these terms should not be included in certain observation models in the integrated model, they can be thinned out using the \code{.$updateFormula} function.
#' Note: the point covariate and mark terms will only be included in the formulas for where they are present in a given dataset, and so these terms do not need to be thinned out if they are not required by certain observation models.
#' 
#' @import methods
#' 
#' @examples 
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
#'  #Set base model up
#'  baseModel <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj, responsePA = 'Present')
#'  
#'  #Print summary
#'  baseModel
#'  
#'  #Set up model with dataset specific spatial fields
#'  
#'  indSpat <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                      Projection = proj, pointsSpatial = 'individual', responsePA = 'Present')
#'                      
#'  #Model with offset variable
#'  offSet <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                      Projection = proj, Offset = 'area', responsePA = 'Present')
#'                      
#'  #Assume area as a mark
#'  markModel <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                      Projection = proj, markNames = 'area', markFamily = 'gamma',
#'                      responsePA = 'Present')
#'                       
#'  }
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
                     Offset = NULL, pointsSpatial = 'shared', marksSpatial = TRUE,
                     responseCounts = 'counts', responsePA = 'present', trialsPA = NULL,
                     trialsMarks = NULL, speciesName = NULL, temporalName = NULL, temporalModel = list(model = 'ar1')) {
  
  ##Things to do: make pointsSpatial one of: 'shared', 'individual', or NULL
   #Need to make a whole new list for all the individual spatial field effects...
  
  if (length(Coordinates) != 2) stop('Coordinates needs to be a vector of length 2 containing the coordinate names.')
  
  if (Coordinates[1] == Coordinates[2]) stop('Coordinates need to be unique values.')
  
  if (!inherits(Projection, 'CRS')) stop('Projection needs to be a CRS object.')
  
  if (!inherits(Mesh, 'inla.mesh')) stop('Mesh needs to be a inla.mesh object.')
  
  if (!is.null(pointsSpatial)) {
    
  if (length(pointsSpatial) != 1) stop('PointsSpatial needs to be one of: "shared", "individual" or NULL.')
    
  if (!pointsSpatial %in% c('shared', 'individual')) stop('PointsSpatial needs to be one of: "shared", "individual" or NULL.')
    
  } 
  
  if (pointsSpatial != 'shared' && !is.null(temporalName)) {
    
    pointsSpatial <- 'shared'
    message('Setting pointsSpatial to "shared" since it is required for the temporalModel.')
    
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
  
  if (!is.null(temporalName)) temporalModel <- deparse1(temporalModel)
  
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
                         speciesspatial = speciesSpatial,
                         offset = Offset)
  
  if (length(list(...)) == 0) warning('No point data added. You can add data to this object with `$.addData()`.')
  
  else {
    
    if (is.null(responseCounts) || is.null(responsePA)) stop('One of responseCounts and responsePA are NULL. Both are required arguments.')
    
    #if (!is.null(responsePA) && is.null(trialsPA)) warning('Present absence response name given but trials name is NULL.')
    
    bruData$addData(..., responseCounts = responseCounts,
                    responsePA = responsePA, trialsPA = trialsPA,
                    markNames = markNames, pointCovariates = pointCovariates,
                    trialsMarks = trialsMarks, speciesName = speciesName,
                    temporalName = temporalName, Offset = Offset)
    
  }
  
  bruData
  
}