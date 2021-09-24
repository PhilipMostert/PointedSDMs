#' Function which calculate the difference in covariate values between a full bru_sdm model and a model with one dataset left out
#' 
#' @param model Model of class bru_sdm run with all datasets.
#' @param dataset Datasets to leave out.
#' @param predictions Will new models be used for predictions. If \code{TRUE} returns marginals and bru_info in model. Defaults to \code{FALSE}. 
#' 
#' @export

leave_one_out <- function(model, dataset,
                          predictions = FALSE) {
  
  if (!inherits(model, 'bru_sdm')) stop('Model needs to be of class "bru_sdm".')
  
  if (!all(dataset%in%model[['sources_of_information']])) stop('Dataset provided not run in initial model.')
  
  if (length(unique(model[['sources_of_information']])) == 1) stop('Model was only run with one dataset. Leaving a dataset out is not possible.')
  
  model_results <- list()
  
  model_differences <- list()
  
    
  if (!is.null(model$spatial_covariates_used)) {
      
  for (names in model$spatial_covariates_used) {
        
  assign(names, model$bru_info$model$effects[[names]]$env[[names]])  
        
  }
      
  if (!is.null(model$bru_info$model$effects[[1]]$env[['spdemodel']])) {
        
  assign('spdemodel', model$bru_info$model$effects[[1]]$env[['spdemodel']])
        
  } else spdemodel <- NULL
      
  if (is.list(spdemodel[[1]])) {
        
  for (name in model[['spatial_datasets']]) {
          
  assign(paste0(name,'_spde'),spdemodel[[name]])  
          
  }
        
  }    
    
  validation_results <- list()  
    
  }
  else validation_results <- list()
  
  for (dataname in dataset) {
    
  index <- !model[['sources_of_information']]%in%dataname
  reduced_options <- list()
  
  #for (option in names(model[['bru_sdm_options']])) {
  #
  #reduced_options[[option]] <- model[['bru_sdm_options']][[option]][index] 
  #  
  #}
  reduced_options <- model[['bru_sdm_options']]
  reduced_options$control.family <- model[['bru_sdm_options']]$control.family[index]
    
  if (!predictions) {
      
  reduced_options$control.compute <- list(return.marginals = FALSE)
      
  }
  ##is it update formula that crashes the model?  
  reduced_components <- update(model$components, paste0(' ~ . -',dataname,'_intercept(1)'))
  reduced_components <- update(reduced_components, paste0(' ~ . - ',
                                                          dataname,'_spde(main = coordinates, model = spdemodel)'))
 
  if (!is.null(model[['spatial_datasets']])) {
    
  if (!any(names(model$bru_info$lhoods)[!names(model$bru_info$lhoods)%in%dataname])%in%model[['spatial_datasets']]) {
    
  reduced_components <- update(reduced_components, paste0(' ~ . -','shared_spatial(main = coordinates, model = spdemodel)'))    
  
  }  
    
  }
  
  if (!is.null(model$marks_used)) {
  
  if (any(names(model$marks_used) == dataname)) {   
       
  for (i in 1:length(model$marks_used)) {
          
  reduced_components <- update(reduced_components, paste0(' ~ . -', dataname, '_',
                                                          model$marks_used[i],'_intercept(1)'))
  reduced_components <- update(reduced_components, paste0(' ~ . -', dataname, '_',
                                                          model$marks_used[i],'_spde(main = coordinates, model = spdemodel)'))
  reduced_components <- update(reduced_components, paste0(' ~ . -', dataname, '_',
                                                          model$marks_used[i],'_spde(main = coordinates, copy = \"', dataname,
                                                          '_spde\", hyper = list(beta = list(fixed = FALSE)))'))

  }
    
  dataset_marks <- model$marks_used[dataname]  
    
  for (mark in dataset_marks) {
      
  if (sum(model$marks_used == mark) == 1) {
        
  reduced_components <- gsub(paste0('[+] ', mark,'*\\(.*?\\) *'), '', reduced_components, perl = T)
  reduced_components <- gsub(paste0('[+] ', mark,'_phi*\\(.*?\\) *'), '', reduced_components, perl = T)   
        
  }
      
  }
      
  }
      
  }

  model_reduced <- inlabru::bru(components = formula(reduced_components),
                                model$bru_info$lhoods[index],
                                options = reduced_options)
   
  model_reduced[['components']] <- reduced_components
    
  if (!predictions) {
      
  model_reduced[['bru_info']] <- NULL
  model_reduced[['call']] <- NULL
      
  }
    
  model_reduced[['data_type']] <- model[['data_type']][index]
  model_reduced[['dataset_names']] <- model[['dataset_names']][index]
  model_reduced[['multinom_vars']] <- model[['multinom_vars']]
    
  class(model_reduced) <- c('bru_sdm',class(model_reduced))
    
  model_results[[paste0('Leaving_out_',dataname)]] <- model_reduced
    
  var_names <- row.names(model_reduced$summary.fixed)[row.names(model_reduced$summary.fixed)%in%row.names(model$summary.fixed)]
  model_differences[[paste0('Leaving_out_',dataname)]] <- model_reduced$summary.fixed[var_names,] - model$summary.fixed[var_names,]
    
  if (predictions) {
      
  #if (model$bru_info$lhoods[[dataname]]$family == 'cp') {
  #
  #reduced_link <- 'default'
  #
  #}
  #else {
  #
  #reduced_link <- 'cloglog'
  #
  #}
  reduced_lik <- model$bru_info$lhoods
  
  for (data in names(model$bru_info$lhoods)[!index]) { 
  
  train <- predict(model_reduced, data = model$bru_info$lhoods[[data]]$data, formula = eval(parse(text = paste0('~ (',paste(model$spatial_covariates_used, collapse = ' + '),')'))))  
  
  reduced_lik[[data]]$data@data['offset'] <- train$mean
  reduced_lik[[data]]$include_components <- c(reduced_lik[data]$include_components, 'offset')
  

  }
  #reduced_lik <- model$bru_info$lhoods[index]
  #reduced_lik[['offset']] <- train_lik
  #Add other options later...
  #reduced_options[['control.family']][[length(reduced_options$control.family) + 1]] <- list(link = reduced_link)
  
  offset_components <- update(model$components, ~ . + offset)
  
  reduced_mlik <- inlabru::bru(offset_components, reduced_lik,
                               options = model$bru_sdm_options)

  validation_results[[dataname]] <- model$mlik[1] - reduced_mlik$mlik[1]
      ##Need to add a bunch of loss functions here
      ## i.e Add RMSE SEL ...
  }
    
    ##Add some sort of cross validation here
    ## To predict the left out dataset and then
    ## Produce some score.
  }
  
  attributes(model_results)[['differences']] <- model_differences
  #names(attributes(model_results)[['differences']]) <- dataset
  
  if (predictions) {
    
  attributes(model_results)[['validation_results']] <- unlist(validation_results)
    
  }
  else attributes(model_results)['validation_results'] <- NULL
  
  #names(model_results) <- paste0('Leaving_out_',dataset)
  class(model_results) <- c('bru_sdm_leave_one_out', 'list')
  
  return(model_results)
  
}
