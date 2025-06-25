#' @title \emph{startSpecies}: Function used to initialize a multi-species integrated species distribution model.
#' 
#' @description This function is used to create an object containing all the data, metadata and relevant components required for the multi-species integrated species distribution model and \pkg{INLA} to work.
#' As a result, the arguments associated with this function are predominantly related to describing variable names within the datasets that are relevant, and arguments related to what terms should be included in the formula for the integrated model. The output of this function is an \code{R6} object, and so there are a variety of public methods within the output of this function which can be used to further specify the model (see \code{?specifySpecies} or \code{.$help()} for a comprehensive description of these public methods).
#' 
#' @param ... The datasets to be used in the model. Must come as either \code{sf} objects, or as a list of named \code{sf} objects. 
#' @param spatialCovariates The spatial covariates used in the model. These covariates must be measured at every location (pixel) in the study area, and must be a \code{SpatialRaster} object. Can be either \code{numeric}, \code{factor} or \code{character} data. Defaults to \code{NULL} which includes no spatial effects in the model.
#' @param Projection The coordinate reference system used by both the spatial points and spatial covariates. Must be of class \code{character}.
#' @param Mesh An \code{fm_mesh_2d} object required for the spatial random fields and the integration points in the model (see \code{\link[fmesher]{fm_mesh_2d_inla}} from the \pkg{fmesher} package for more details). 
#' @param IPS The integration points to be used in the model (that is, the points on the map where the intensity of the model is calculated). See \code{\link[fmesher]{fm_int}} from the \pkg{fmesher} package for more details regarding these points; however defaults to \code{NULL} which will create integration points from the \code{Mesh} object.
#' @param Boundary A \code{sf} object of the study area. If not missing, this object is used to help create the integration points.
#' @param speciesSpatial Argument to specify if each species should have their own spatial effect with different hyperparameters to be estimated using \pkg{INLA}'s "replicate" feature, of if a the field's should be estimated per species copied across datasets using \pkg{INLA}'s "copy" feature. Possible values include: \code{'replicate'}, \code{'copy'}, \code{'shared'} or \code{NULL} if no species-specific spatial effects should be estimated.
#' @param speciesIntercept Argument to control the species intercept term. Defaults to \code{TRUE} which creates a random intercept term, \code{FALSE} creates a fixed intercept term, and \code{NULL} removes the intercept term.
#' @param speciesEnvironment Argument to control the species environmental term. Defaults to \code{"stack"} which creates species specific effects. Setting the argument to \code{"community"} will create a community mean and variations around for each covariate, and setting the argument to \code{"shared"} will create shared responses to the covariates. Note that if \code{covariateFormula} is provided in the \code{Formulas} argument, then this argument cannot be set to \code{"community"}. 
#' @param pointCovariates The non-spatial covariates to be included in the integrated model (for example, in the field of ecology the distance to the nearest road or time spent sampling could be considered). These covariates must be included in the same data object as the points.
#' @param Offset Name of the offset variable (class \code{character}) in the datasets. Defaults to \code{NULL}; if the argument is non-\code{NULL}, the variable name needs to be standardized across datasets (but does not need to be included in all datasets). The offset variable will be transformed onto the log-scale in the integrated model.
#' @param pointsIntercept Logical argument: should the points be modeled with intercepts. Defaults to \code{TRUE}.  Note that if this argument is non-\code{NULL} and \code{pointsIntercepts} is missing, \code{pointsIntercepts} will be set to \code{FALSE}.
#' @param pointsSpatial Argument to determine whether the spatial field is shared between the datasets, or if each dataset has its own unique spatial field. The datasets may share a spatial field with \pkg{INLA}'s "copy" feature if the argument is set to \code{copy}. May take on the values: \code{"shared"}, \code{"individual"}, \code{"copy"}, \code{"correlate"} or \code{NULL} if no spatial field is required for the model. Defaults to \code{"copy"}.
#' @param responseCounts Name of the response variable in the counts/abundance datasets. This variable name needs to be standardized across all counts datasets used in the integrated model. Defaults to \code{'counts'}.
#' @param responsePA Name of the response variable (class \code{character}) in the presence absence/detection non-detection datasets. This variable name needs to be standardized across all present absence datasets. Defaults to \code{'present'}.
#' @param trialsPA Name of the trials response variable (class \code{character}) for the presence absence datasets. Defaults to \code{NULL}.
#' @param speciesName Name of the species variable name (class \code{character}). Specifying this argument turns the model into a stacked species distribution model, and calculates covariate values for the individual species, as well as a species group model in the shared spatial field. Defaults to \code{NULL}. Note that if this argument is non-\code{NULL} and \code{pointsIntercepts} is missing, \code{pointsIntercepts} will be set to \code{FALSE}.
#' @param temporalName Name of the temporal variable (class \code{character}) in the model. This variable is required to be in all the datasets. Defaults to \code{NULL}.
#' @param Formulas A named list with two objects. The first one, \code{covariateFormula}, is a formula for the covariates and their transformations for the distribution part of the model. Defaults to \code{NULL} which includes all covariates specified in \code{spatialCovariates} into the model. The second, \code{biasFormula}, specifies which covariates are used for the PO datasets. Defaults to \code{NULL} which includes no covariates for the PO datasets.
#' 
#' @return A \code{\link{specifySpecies}} object (class \code{R6}). Use \code{?specifySpecies} to get a comprehensive description of the slot functions associated with this object.
#' 
#' @note The idea with this function is to describe the full model: that is, all the covariates and spatial effects will appear in all the formulas for the datasets and species.
#' If some of these terms should not be included in certain observation models in the integrated model, they can be thinned out using the \code{.$updateFormula} function.
#' Note: the point covariate will only be included in the formulas for where they are present in a given dataset, and so these terms do not need to be thinned out if they are not required by certain observation models.
#' 
#' @import methods
#' 
#' @examples 
#'  
#'  \dontrun{
#'  if (requireNamespace('INLA')) {
#'  
#'  ##REDO WITH OTHER DATA
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  
#'  
#'  proj <- "+proj=longlat +ellps=WGS84"
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set base model up
#'  baseModel <- startSpecies(data, Mesh = mesh,
#'                            Projection = proj, responsePA = 'Present',
#'                            speciesName = 'speciesName')
#'  
#'  #Print summary
#'  baseModel
#'  
#'  #Set up model with dataset specific spatial fields
#'  
#'  indSpat <- startSpecies(data, Mesh = mesh,
#'                          Projection = proj, pointsSpatial = 'individual', 
#'                          responsePA = 'Present', speciesName = 'speciesName')
#'                      
#'  #Model with offset variable
#'  offSet <- startSpecies(data, Mesh = mesh,
#'                      Projection = proj, 
#'                      Offset = 'area', 
#'                      responsePA = 'Present',
#'                      speciesName = 'speciesName')
#'                      
#'  #Non-random effects for the species
#'  speciesInt <- startSpecies(data, Mesh = mesh,
#'                         Projection = proj, 
#'                         speciesIntercept = FALSE,
#'                         responsePA = 'Present',
#'                         speciesName = 'speciesName')
#'                         
#' #Turn off species level field
#'
#'  speciesInt <- startSpecies(data, Mesh = mesh,
#'                         Projection = proj, 
#'                         speciesSpatial = NULL,
#'                         responsePA = 'Present',
#'                         speciesName = 'speciesName')
#'                      
#'  }
#'  }
#' 
#' @export

startSpecies <- function(..., spatialCovariates = NULL, 
                      Projection, Mesh, 
                      speciesSpatial = 'replicate', 
                      speciesIntercept = TRUE, 
                      speciesEnvironment = 'stack', 
                      speciesName,
                      IPS = NULL, Boundary = NULL, pointCovariates = NULL,
                      Offset = NULL, pointsIntercept = TRUE, 
                      pointsSpatial = 'copy', 
                      responseCounts = 'counts', 
                      responsePA = 'present', 
                      trialsPA = NULL, temporalName = NULL,
                      Formulas = list(covariateFormula = NULL,
                                      biasFormula = NULL)) {
  
  if (!is.null(Boundary) && !inherits(Boundary, c('sf', 'sfc'))) stop('Boundary needs to be an sf object.')
  
  if (!inherits(Projection, 'character')) stop('Projection needs to be a character object.')
  
  if (!inherits(Mesh, 'fm_mesh_2d')) stop('Mesh needs to be a fm_mesh_2d object.')
  
  if (missing(speciesName)) stop('speciesName needs to be provided.')
  
  if (!is.null(pointsSpatial)) {
    
    if (length(pointsSpatial) != 1 || !pointsSpatial %in% c('shared', 'individual', 'copy', 'correlate')) stop('PointsSpatial needs to be one of: "shared", "copy", "individual", "correlate" or NULL.')
    
  } 
  
  if (!is.null(speciesSpatial)) {
    
    if (length(speciesSpatial) != 1 || !speciesSpatial %in% c('shared', 'copy', 'replicate')) stop('speciesSpatial needs to be one of: "shared", "copy", "replicate" or NULL.')
    
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
  
  if (any(!names(Formulas) %in% c('covariateFormula', 'biasFormula'))) stop('Formulas must be a named list containing at least one of covariateFormula or biasFormula.')
  
  if (is.logical(speciesEnvironment)) {
    
    lifecycle::deprecate_warn(when = '2.1.4', 
   what = "PointedSDMs::startSpecies(speciesEnvironment = 'Must now be one of community, stack or shared. Setting speciesEnvironment to TRUE will create a community model and setting it to FALSE will create a stacked model')",always = TRUE)
    
    if (speciesEnvironment) speciesEnvironment <- 'community'
    else 'stacked'
  } else 
    if (!speciesEnvironment %in% c('community', 'stack', 'shared') | length(speciesEnvironment) != 1) stop('speciesEnvironment must be one of "community", "stacked" or "shared".')
  
  if (!is.null(Formulas$covariateFormula) & speciesEnvironment == 'community') {
    
    warning ('speciesEnvironment cannot be "community" if covariateFormula is provided. Setting speciesEnvironment to "stack".')
    speciesEnvironment <- 'stack'
    
  }
  
  bruData <- specifySpecies$new(data = dataPoints, projection = Projection,
                             Inlamesh = Mesh, initialnames = initialnames,
                             responsecounts = responseCounts,
                             responsepa = responsePA,
                             pointcovariates = pointCovariates,
                             trialspa = trialsPA,
                             spatial = pointsSpatial,
                             intercepts = pointsIntercept,
                             spatialcovariates = spatialCovariates,
                             boundary = Boundary,
                             speciesindependent = FALSE, 
                             speciesintercept = speciesIntercept,
                             speciesname = speciesName, 
                             speciesenvironment = speciesEnvironment, 
                             speciesspatial = speciesSpatial,
                             ips = IPS,
                             temporal = temporalName,
                             temporalmodel = temporalModel,
                             offset = Offset,
                             copymodel = copyModel,
                             formulas = Formulas)
  
  bruData
  
  
  
}
