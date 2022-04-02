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
  
  block_index <- lapply(unlist(data$.__enclos_env__$private$modelData), function(x) x@data[,'block_index'])
  
  for (fold in unique(unlist(block_index))) {
    
    ##Maybe make block_index in dataSDM as a list such that we can see which datasets are not in block i to easily remove them.
     #Get formula terms only after likelihood construction, and then thin components from there.
     #And also for the whole, control.family thing
    
    trainLiks <- do.call(inlabru::like_list,
                 makeLhoods(data = data$.__enclos_env__$private$modelData,
                 formula = data$.__enclos_env__$private$Formulas,
                 family = data$.__enclos_env__$private$Family,
                 mesh = data$.__enclos_env__$private$INLAmesh,
                 ips = data$.__enclos_env__$private$IPS,
                 paresp = data$.__enclos_env__$private$responsePA,
                 ntrialsvar = data$.__enclos_env__$private$trialsPA,
                 markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                 speciesname = data$.__enclos_env__$private$speciesName,
                 speciesindex = data$.__enclos_env__$private$speciesIndex,
                 filter = block_index))
      ##Get formulas here
       #Then thin components
    trainedModel <- inlabru::bru(components = thinnedComponents,
                                 trainLiks,
                                 options = thinedOptions)
    
    test <- do.call(rbind.SpatialPoints,
            lapply(unlist(data$.__encos_enc__$private$modelData, recursive = TRUE), function (x) {
        
                x[x$block_index == block_index, ]}))
    
    
    predictTest <- predict(object = trainedModel, data = test, formula = ~(predictor))
    #...  
    
    }
  
  
  formula_terms <- unique(unlist(lapply(data$.__enclos_env__$private$modelData, function(x) {
    
    if (is.null(x$include_components))  attributes(terms(x$formula))[['term.labels']]
    else x$include_components
    
  })))
  
  comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
  
  comp_keep <- comp_terms %in% formula_terms
  
  componentsJoint <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
  
  ##in case there are duplicates, will it cause an error??
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  
}