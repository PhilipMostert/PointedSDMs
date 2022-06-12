#' @title \emph{blockedCV}: run spatial blocked cross-validation on the integrated model.
#' 
#' @description This function is used to perform spatial blocked cross-validation with regards to model selection for the integrated model. It does so by leaving out a block of data in the full model, running a model with the remaining data, and then calculating the deviance information criteria (DIC) as a score of model fit.
#' @param data An object produced by \code{\link{intModel}}. Requires the slot function, \code{.$spatialBlock} to be run first in order to specify how the data in the model is blocked.
#' @param options A list of \pkg{INLA} or \pkg{inlabru} options to be used in the model. Defaults to \code{list()}.
#' 
#' @import inlabru
#' @import stats

#' 
#' @examples 
#' 
#'\dontrun{
#'  if(requireNamespace('INLA')) {
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  proj <- CRS("+proj=longlat +ellps=WGS84")
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj)
#'  
#'  #Set up spatial block
#'  organizedData$spatialBlock(k = 2, rows = 2, cols = 1)
#'  
#'  #Run spatial block cross-validation
#'  blocked <- blockedCV(organizedData)
#'  
#'  #Print summary
#'  blocked
#'    
#'  }
#'}
#' 
#' @return An object of class \code{blockedCV}, which is essentially a list of DIC values obtained from each iteration of the model.
#' 
#' @export
#' 
#' 
#' 
blockedCV <- function(data, options = list()) {
  
  #How do we do this?
   #Should we make data a list of data files;
   # or should we make another argument for thinned formulas to test based on the full model?
  
  if (!inherits(data, 'dataSDM')) stop('data needs to be a dataSDM object.')
  
  if (is.null(data$.__enclos_env__$private$INLAmesh)) stop('An inla.mesh object is required before any model is run.')
  
  if (!data$.__enclos_env__$private$blockedCV) stop('Please use ".$spatialBlock" before using this function.')
  
  data2ENV(data = data, env = environment())
  
  deviance <- list()
  
  block_index <- lapply(unlist(data$.__enclos_env__$private$modelData), function(x) x@data[,'.__block_index__'])
  
  for (fold in unique(unlist(block_index))) {
    
    ##Maybe make block_index in dataSDM as a list such that we can see which datasets are not in block i to easily remove them.
     #Get formula terms only after likelihood construction, and then thin components from there.
     #And also for the whole, control.family thing
    
    
    trainData <- lapply(data$.__enclos_env__$private$modelData, function(data) {
      
      lapply(data, function(x) {
        
        x[x$.__block_index__ != fold,]
        
      })
      
      
    })
    
    trainLiks <- do.call(inlabru::like_list,
                 makeLhoods(data = trainData,
                 formula = data$.__enclos_env__$private$Formulas,
                 family = data$.__enclos_env__$private$Family,
                 mesh = data$.__enclos_env__$private$INLAmesh,
                 ips = data$.__enclos_env__$private$IPS,
                 paresp = data$.__enclos_env__$private$responsePA,
                 ntrialsvar = data$.__enclos_env__$private$trialsPA,
                 markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                 speciesname = data$.__enclos_env__$private$speciesName,
                 speciesindex = data$.__enclos_env__$private$speciesIndex))
      
    formula_terms <- unique(unlist(lapply(trainLiks, function(x) {
      
      if (!is.null(x$include_components)) x$include_components
      else labels(terms(x$formula))
      
    })))
    
    comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
    
    comp_keep <- comp_terms %in% formula_terms
    
    thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))

    foldOptions <- data$.__enclos_env__$private$optionsINLA
    
    fold_ind <- unique(unlist(block_index))[unique(unlist(block_index)) != fold]
    
    foldOptions$control.family <- foldOptions$control.family[sapply(unlist(data$.__enclos_env__$private$modelData), 
                                                                    function(x) any(fold_ind %in% x@data[, '.__block_index__']))]

    optionsTrain <- append(options, foldOptions)
    
    ##Calculate DIC for just this model?
    trainedModel <- inlabru::bru(components = thinnedComponents,
                                 trainLiks,
                                 options = optionsTrain)
    
    ## -log(intensity)
    ## add an offset argument...
    deviance[[paste0('DIC_fold_', fold)]] <- trainedModel$dic
    
    }
  
  
  
  comps <- formula(paste0(' ~ ', paste0(gsub('\\(.*$', '', data$.__enclos_env__$private$Components), collapse = ' + ')))
  deviance <- append(deviance, list(Formula = comps))
  class(deviance) <- c('blockedCV', 'list')
  deviance

}


#' Export class blockedCV
#' 
#' @export

setClass('blockedCV')

#' Print for blockedCV
#' 
#' @export print.blockedCV

#' Export print.blockedCV
#' @title Print function for \code{blockedCV}.
#' @param x A blockedCV object.
#' @param ... Unused argument.
#' 
#' @exportS3Method 

print.blockedCV <- function(x, ...) {
  
  cat('Spatial block cross-validation score:')
  cat('\n\n')
  cat('Formula: ')

  cat(deparse1(x$Formula))
  cat('\n\n')
  x$Formula <- NULL
  mean.deviance <- sapply(x, function(y) y$mean.deviance)
  p.eff <- sapply(x, function(y) y$p.eff)
  dic <- sapply(x, function(y) y$dic)
  
  dataobj <- data.frame(mean.deviance = mean.deviance,
                        p.eff = p.eff,
                        dic = dic)
  
  row.names(dataobj) <- paste0('fold ', 1:nrow(dataobj))
  
  print.data.frame(dataobj)
  
  cat('\nmean DIC score: ')
  cat(mean(dataobj$dic))


}
