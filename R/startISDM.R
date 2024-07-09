#' @title \emph{startISDM}: Function used to initialize the integrated species distribution model.
#' 
#' @description This function is used to create an object containing all the data, metadata and relevant components required for the integrated species distribution model and \pkg{INLA} to work.
#' As a result, the arguments associated with this function are predominantly related to describing variable names within the datasets that are relevant, and arguments related to what terms should be included in the formula for the integrated model. The output of this function is an \code{R6} object, and so there are a variety of public methods within the output of this function which can be used to further specify the model (see \code{?specifyISDM} or \code{.$help()} for a comprehensive description of these public methods).
#' 
#' @param ... The datasets to be used in the model. Must come as either \code{sf} objects, or as a list of named \code{sf} objects. 
#' @param spatialCovariates The spatial covariates used in the model. These covariates must be measured at every location (pixel) in the study area, and must be a \code{SpatialRaster} object. Can be either \code{numeric}, \code{factor} or \code{character} data. Defaults to \code{NULL} which includes no spatial effects in the model.
#' @param Projection The coordinate reference system used by both the spatial points and spatial covariates. Must be of class \code{character}.
#' @param Mesh An \code{inla.mesh} object required for the spatial random fields and the integration points in the model (see \code{\link[INLA]{inla.mesh.2d}} from the \pkg{INLA} package for more details). 
#' @param IPS The integration points to be used in the model (that is, the points on the map where the intensity of the model is calculated). See \code{\link[inlabru]{fm_int}} from the \pkg{inlabru} package for more details regarding these points; however defaults to \code{NULL} which will create integration points from the \code{Mesh} and \code{Boundary }objects.
#' @param Boundary A \code{sf} object of the study area. If not missing, this object is used to help create the integration points.
#' @param pointCovariates The non-spatial covariates to be included in the integrated model (for example, in the field of ecology the distance to the nearest road or time spent sampling could be considered). These covariates must be included in the same data object as the points, and do not necessarily need to be present in all datasets.
#' @param Offset Name of the offset variable (class \code{character}) in the datasets. Defaults to \code{NULL}; if the argument is non-\code{NULL}, the variable name needs to be standardized across datasets (but does not need to be included in all datasets). The offset variable will be transformed onto the log-scale in the integrated model.
#' @param pointsIntercept Logical argument: should the points be modeled with intercepts. Defaults to \code{TRUE}.
#' @param pointsSpatial Argument to determine whether the spatial field is shared between the datasets, or if each dataset has its own unique spatial field. The datasets may share a spatial field with \pkg{INLA}'s "copy" feature if the argument is set to \code{copy}. May take on the values: \code{"shared"}, \code{"individual"}, \code{"copy"} or \code{NULL} if no spatial field is required for the model. Defaults to \code{"copy"}.
#' @param responseCounts Name of the response variable in the counts/abundance datasets. This variable name needs to be standardized across all counts datasets used in the integrated model. Defaults to \code{'counts'}.
#' @param responsePA Name of the response variable (class \code{character}) in the presence absence/detection non-detection datasets. This variable name needs to be standardized across all present absence datasets. Defaults to \code{'present'}.
#' @param trialsPA Name of the trials response variable (class \code{character}) for the presence absence datasets. Defaults to \code{NULL}.
#' @param temporalName Name of the temporal variable (class \code{character}) in the model. This variable is required to be in all the datasets. Defaults to \code{NULL}.
#' @param Formulas A named list with two objects. The first one, \code{covariateFormula}, is a formula for the covariates and their transformations for the distribution part of the model. Defaults to \code{NULL} which includes all covariates specified in \code{spatialCovariates} into the model. The second, \code{biasFormula}, specifies which covariates are used for the PO datasets. Defaults to \code{NULL} which includes no covariates for the PO datasets.
#' 
#' @return A \code{\link{specifyISDM}} object (class \code{R6}). Use \code{?specifyISDM} to get a comprehensive description of the slot functions associated with this object.
#' 
#' @note The idea with this function is to describe the full model: that is, all the covariates and spatial effects will appear in all the formulas for the datasets and species.
#' If some of these terms should not be included in certain observation models in the integrated model, they can be thinned out using the \code{.$updateFormula} function.
#' Note: the point covariate will only be included in the formulas for where they are present in a given dataset, and so these terms do not need to be thinned out if they are not required by certain observation models.
#' 
#' @import methods
#' 
#' @examples 
#' 
#' \dontrun{
#'  if (requireNamespace('INLA')) {
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  proj <- "+proj=longlat +ellps=WGS84"
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set base model up
#'  baseModel <- startISDM(data, Mesh = mesh,
#'                              Projection = proj, responsePA = 'Present')
#'  
#'  #Print summary
#'  baseModel
#'  
#'  #Set up model with dataset specific spatial fields
#'  indSpat <- startISDM(data, Mesh = mesh,
#'                      Projection = proj, pointsSpatial = 'individual', responsePA = 'Present')
#'                      
#'  #Model with offset variable
#'  offSet <- startISDM(data, Mesh = mesh,
#'                      Projection = proj, Offset = 'area', responsePA = 'Present')
#'  
#'  #Add own formula
#'  quadModel <- startISDM(data, Mesh = mesh, 
#'                         Projection = proj, pointsSpatial = 'copy', responsePA = 'Present',
#'                         Formulas = list(covariateFormula = ~ NPP + I(NPP^2)))                    
#'                       
#'  }
#'  }
#' 
#' @export

startISDM <- function(..., spatialCovariates = NULL, 
                      Projection, Mesh, 
                      IPS = NULL, Boundary = NULL, pointCovariates = NULL,
                      Offset = NULL, pointsIntercept = TRUE, pointsSpatial = 'copy', 
                      responseCounts = 'counts', 
                      responsePA = 'present', 
                      trialsPA = NULL, temporalName = NULL,
                      Formulas = list(covariateFormula = NULL,
                                      biasFormula = NULL)) {

  #Are we keeping the argument names the same?
  
  if (!is.null(Boundary) && !inherits(Boundary, c('sf', 'sfc'))) stop('Boundary needs to be an sf object.')
  
  if (!inherits(Projection, 'character')) stop('Projection needs to be a character object.')
  
  if (!inherits(Mesh, 'inla.mesh')) stop('Mesh needs to be a inla.mesh object.')
  
  if (!is.null(pointsSpatial)) {
    
    if (length(pointsSpatial) != 1 || !pointsSpatial %in% c('shared', 'individual', 'copy', 'correlate')) stop('PointsSpatial needs to be one of: "shared", "copy", "individual", "correlate" or NULL.')
    
  } 
  
  #copy here?
   #Is this TRUE
  if (pointsSpatial %in% c('correlate') && !is.null(temporalName)) {
    
    pointsSpatial <- 'shared'
    message('Setting pointsSpatial to "shared" since it is required for the temporalModel.')
    
  }
  
  dataPoints <- list(...)
  
  #if (length(dataPoints) == 0) stop('Please provide data in the ... argument.')
  if (length(dataPoints) > 0) {
    
    datasetClass <- unlist(lapply(dataPoints, class))
    ##Need something more here if object is list: get names from list
    if (length(datasetClass) == 1 && datasetClass == "list") {
      
      if (!is.null(names(dataPoints[[1]]))) {
        
        dataList <- TRUE
        dataNames <- names(dataPoints[[1]])
        dataPoints <- unlist(dataPoints, recursive = FALSE)
        datasetClass <- lapply(dataPoints, class)

      }
      
      else {
        
        dataNames <- NULL
        dataPoints <- unlist(dataPoints, recursive = FALSE)
        datasetClass <- lapply(dataPoints, class)
        dataList <- TRUE
        
      }
      
    }
    else {
      
      dataList <- FALSE
      datasetNames <- NULL
      
    }
    
    isSF <- unlist(lapply(dataPoints, function(x) inherits(x, 'sf')))
    if (!all(isSF)) stop("Datasets need to be sf objects.")
    
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
      else initialnames <- dataNames
    }
    else {
      
      initialnames <- setdiff(as.character(match.call(expand.dots = TRUE)), 
                              as.character(match.call(expand.dots = FALSE)))
    }
    
  } else {
    
    initialnames <- NULL
    
  }
  
  if (length(initialnames) == 1 & !is.null(pointsSpatial)) pointsSpatial <- 'shared'
  
  if (!is.null(temporalName)) temporalModel <- deparse1(list(model = 'ar1'))
  else temporalModel <- NULL
  
  if (is.null(pointsSpatial) || pointsSpatial == 'shared') copyModel <- NULL
  else copyModel <- deparse1(list(beta = list(fixed = FALSE)))
  
  if (!is.null(Formulas$covariateFormula) &&
      !is.null(Formulas$biasFormula)) {
    
    ##Check that there is nothing overlapping here -- remove anything that is
    
  }
  
  bruData <- specifyISDM$new(data = dataPoints, projection = Projection,
                             Inlamesh = Mesh, initialnames = initialnames,
                             responsecounts = responseCounts,
                             responsepa = responsePA,
                             pointcovariates = pointCovariates,
                             trialspa = trialsPA,
                             spatial = pointsSpatial,
                             intercepts = pointsIntercept,
                             spatialcovariates = spatialCovariates,
                             boundary = Boundary,
                             ips = IPS,
                             temporal = temporalName,
                             temporalmodel = temporalModel,
                             offset = Offset,
                             copymodel = copyModel,
                             formulas = Formulas)
  
  bruData
  
  
  
}
