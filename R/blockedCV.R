#' @title Run blocked cross-validation.
#' 
#' @param data A bruSDM data file to be used in the integrated model.
#' @param options A list of INLA options used in the model. Defaults to \code{list()}.
#' 
#' @export
#' 
blockedCV <- function(data, options = list()) {
  
  if (!inherits(data, 'dataSDM')) stop('data needs to be a dataSDM object.')
  
  if (!data$.__enclos_env__$private$blockedCV) stop('Please use ".$spatialBlock" before using this function.')
  
  data2ENV(data = data, env = environment())
  
  #Need to subset data based on value of block_id
  #Remove all species not in index
  #Remember to delete all empty likelihoods
  #Remember to delete all un-used components (ie similar to datasetOut so make function)
  #Run model on these data
  #Predict using the left out data
  #calculate score
  
  #Need another index for number of blocks
  
  block_index <- lapply(unlist(data$.__enclos_env__$private$modelData), function(x) x@data[,'.__block_index__'])
  
  for (fold in unique(unlist(block_index))) {
    
    ##Maybe make block_index in dataSDM as a list such that we can see which datasets are not in block i to easily remove them.
     #Get formula terms only after likelihood construction, and then thin components from there.
     #And also for the whole, control.family thing
    
    
    trainData <- lapply(dataObj$.__enclos_env__$private$modelData, function(data) {
      
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

    trainedModel <- inlabru::bru(components = thinnedComponents,
                                 trainLiks,
                                 options = options)
    
    test <- do.call(rbind.SpatialPoints,
            lapply(unlist(data$.__enclos_env__$private$modelData, recursive = TRUE), function (x, idx) {
              
              x[x$.__block_index__ == idx, ]
                
                }, idx = fold))
    
    
    test_formula <- formula(paste('~ ', paste(gsub('\\(.*$', '', data$.__enclos_env__$private$Components[comp_keep]), collapse = ' + ')))

    predictTest <- predict(object = trainedModel, data = test, formula = test_formula)
    
    stop(return(predictTest))
    
    }
  
  
  
}