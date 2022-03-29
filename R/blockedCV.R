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
  
  block_index <- lapply(data$.__enclos_env__$private$modelData, function(x) x$data@data[,'block_index'])
  
  for (fold in unique(unlist(block_index))) {
    
    for (dataset in names(data$.__enclos_env__$private$modelData)) {
      
      # if family == 'cp'
       # index family by removing all the data points not in the fold
       # go into data@data[,'BRU_aggregate'] and count the number of TRUEs
       # Then go into response_data and change BRU_response_cp to that sum
      ##Need to also remove ipoints not in fold::
      
      
      #train: all data not fold
      train <- zz
      #test: all data fold
      test <- xx  
      
    }
    
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