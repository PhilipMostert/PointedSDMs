#' Function which calculate the difference in covariate values between a full bru_sdm model and a model with one dataset left out
#' 
#' @param model Model of class bru_sdm run with all datasets
#' @param dataset Dataset to leave out
#' @param predictions Will new models be used for predictions. If \code{TRUE} returns marginals and bru_info in model. Defaults to \code{FALSE}. 

leave_one_out <- function(model, dataset,
                          predictions = FALSE) {
  
  if (!inherits(model, 'bru_sdm')) stop('Model needs to be of class "bru_sdm".')
  
  if (!all(dataset%in%model[['sources_of_information']])) stop('Dataset provided not run in initial model.')

  if (length(unique(model[['sources_of_information']])) == 1) stop('Model was only run with one dataset. Leaving a dataset out is not possible.')
    
  model_results <- list()
  
  for (dataname in dataset) {
  
  index <- !model[['sources_of_information']]%in%dataname
  
  reduced_options <- model[['bru_sdm_options']]
  
  reduced_options <- model[['bru_sdm_options']]$control.family[index]
  
  if (!predictions) {
  
  reduced_options$control.compute <- list(return.marginals = FALSE)
  
  }
  model_reduced <- bru(components = model$components,
                       model$bru_info$lhoods[index],
                       options = reduced_options)
  
  
  if (!predictions) {
  
  model_reduced[['bru_info']] <- NULL
  model_reduced[['call']] <- NULL
  
  }
  
  model_reduced[['data_type']] <- model[['data_type']][index]
  model_reduced[['dataset_names']] <- model[['dataset_names']][index]
  model_reduced[['multinom_vars']] <- model[['multinom_vars']]
  
  class(model_reduced) <- c('bru_sdm',class(model_reduced))
  
  model_results[[paste0('Leaving_out_',dataname)]] <- model_reduced
  
  }
  
  
  return(model_results)

  }