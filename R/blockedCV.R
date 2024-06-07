#' @title \emph{blockedCV}: run spatial blocked cross-validation on the integrated model.
#' 
#' @description This function is used to perform spatial blocked cross-validation with regards to model selection for the integrated model. It does so by leaving out a block of data in the full model, running a model with the remaining data, and then calculating the deviance information criteria (DIC) as a score of model fit.
#' @param data An object produced by either \code{\link{startISDM}} of \code{\link{startSpecies}}. Requires the slot function, \code{.$spatialBlock} to be run first in order to specify how the data in the model is blocked.
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
#'  proj <- "+proj=longlat +ellps=WGS84"
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- startISDM(data, Mesh = mesh,
#'                             responsePA = 'Present',
#'                             Projection = proj)
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
blockedCV <- function(data, options = list()) {
  
  #How do we do this?
   #Should we make data a list of data files;
   # or should we make another argument for thinned formulas to test based on the full model?
  
  if (!inherits(data, 'dataSDM') && !inherits(data, 'specifySpecies') && !inherits(data, 'specifyISDM')) stop('data needs to be a dataSDM object.')
  
  if (is.null(data$.__enclos_env__$private$INLAmesh)) stop('An inla.mesh object is required before any model is run.')
  
  if (!data$.__enclos_env__$private$blockedCV) stop('Please use ".$spatialBlock" before using this function.')
  
  data2ENV(data = data, env = environment())
  
  deviance <- list()
  
  block_index <- lapply(unlist(data$.__enclos_env__$private$modelData, recursive = FALSE), function(x) data.frame(x)[, '.__block_index__'])
  
  if (!is.null(data$.__enclos_env__$private$temporalName)) {
    
    numTime <- length(unique(unlist(data$.__enclos_env__$private$temporalVars)))
    
    newIPS <- rep(list(data$.__enclos_env__$private$IPS), numTime)
    
    newIPS <- do.call(rbind, newIPS)
    
    newIPS[, data$.__enclos_env__$private$temporalName] <- rep(1:numTime, each = nrow(data$.__enclos_env__$private$IPS))
    
    newIPS <- st_transform(newIPS, data$.__enclos_env__$private$Projection)
    
    data$.__enclos_env__$private$IPS <- newIPS
    
  }
  
  for (fold in unique(unlist(block_index))) {
    
    ##Maybe make block_index in dataSDM as a list such that we can see which datasets are not in block i to easily remove them.
     #Get formula terms only after likelihood construction, and then thin components from there.
     #And also for the whole, control.family thing
    
    
    trainData <- lapply(data$.__enclos_env__$private$modelData, function(data) {
      
      lapply(data, function(x) {
        
        data <- x[x$.__block_index__ != fold,]
        if (nrow(data) > 0) data
        else NULL
        
      })
      
      
    })
    
    ##Check if all copy + bias Main in here
    
    testData <- lapply(data$.__enclos_env__$private$modelData, function(data) {
      
      lapply(data, function(x) {
        
        data <- x[x$.__block_index__ == fold,]
        if (nrow(data) > 0) data
        else NULL
        
      })
      
      
    })
    
    trainLiks <- do.call(inlabru::like_list,
                 makeLhoods(data = trainData,
                 formula = data$.__enclos_env__$private$Formulas,
                 family = data$.__enclos_env__$private$Family,
                 mesh = data$.__enclos_env__$private$INLAmesh,
                 ips = data$.__enclos_env__$private$IPS,
                 samplers = data$.__enclos_env__$private$Samplers,
                 paresp = data$.__enclos_env__$private$responsePA,
                 ntrialsvar = data$.__enclos_env__$private$trialsPA,
                 markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                 speciesname = data$.__enclos_env__$private$speciesName,
                 speciesindex = data$.__enclos_env__$private$speciesIndex))
      
    formula_terms <- unique(unlist(lapply(trainLiks, function(x) {
      
      if (!identical(unlist(x$used), character(0))) unlist(x$used)
      else labels(terms(x$formula))
      
    })))
    
    comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
    
    comp_keep <- comp_terms %in% formula_terms
    
    if (!all(comp_keep)) {
      
      if (data$.__enclos_env__$private$Spatial != 'shared' | data$.__enclos_env__$private$biasCopy) {
        
        Main <- grepl('_spatial', data$.__enclos_env__$private$Components) & grepl('_field', data$.__enclos_env__$private$Components)
        if (!any(Main)) Main <- 'NOTMAIN'
        else Main <- sub("\\(.*", "", comp_terms[Main])
        
        MainBias <- grepl('_biasField', data$.__enclos_env__$private$Components) & grepl('_bias_field', data$.__enclos_env__$private$Components)
        if (!any(MainBias)) MainBias <- 'NOTALLMAINBIAS'
        else MainBias <- sub("\\(.*", "", comp_terms[MainBias])
        
        whichMissing <- names(trainData)[sapply(unlist(trainData, recursive = F), is.null)]
        warning('More than 2 datasets missing from the block with either pointsSpatial = "copy" or copyModel = TRUE for the bias field.\n Will choose the first available dataset to copy on.')
        if(paste0(whichMissing[1],'_biasField') == MainBias | paste0(whichMissing[1],'_spatial') == Main) {
          
          thinnedComponents <- reduceComps(componentsOld =  formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components, collapse = ' + '))),
                                            pointsCopy = ifelse(data$.__enclos_env__$private$Spatial == 'copy', 
                                                                TRUE, FALSE),
                                            biasCopy = data$.__enclos_env__$private$biasCopy,
                                            datasetName = whichMissing[1],
                                            reducedTerms = comp_terms[comp_keep])
          
        } else  thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
        
        
      } else thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
      
      
    }
    else thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))

    foldOptions <- data$.__enclos_env__$private$optionsINLA
    
    fold_ind <- unique(unlist(block_index))[unique(unlist(block_index)) != fold]
    
    foldOptions$control.family <- foldOptions$control.family[sapply(unlist(data$.__enclos_env__$private$modelData, recursive = FALSE), 
                                                                    function(x) any(fold_ind %in% data.frame(x)[, '.__block_index__']))]

    optionsTrain <- append(options, foldOptions)
    
    ##Calculate DIC for just this model?
    trainedModel <- try(inlabru::bru(components = thinnedComponents,
                                 trainLiks,
                                 options = optionsTrain))
    
    ## -log(intensity)
    ## add an offset argument...
    if (inherits(trainedModel, 'try-error')) {
      
      warning('Model failed for a block. Please change your block layout to ensure all datasets are included in all blocks')
      deviance[[paste0('DIC_fold_', fold)]] <- NA
    }
    else deviance[[paste0('DIC_fold_', fold)]] <- trainedModel$dic
    
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
  cat(mean(dataobj$dic, na.rm = TRUE))


}
