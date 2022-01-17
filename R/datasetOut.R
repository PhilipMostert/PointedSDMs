#' Function which calculate the difference in covariate values between a full bru_sdm model and a model with one dataset left out
#' 
#' @param model Model of class bru_sdm run with all datasets.
#' @param dataset Datasets to leave out.
#' @param predictions Will new models be used for predictions. If \code{TRUE} returns marginals and bru_info in model. Defaults to \code{FALSE}. 
#' 
#' @export

datasetOut <- function(model, dataset,
                       predictions = FALSE) {
  
  if (!inherits(model, 'bruSDM')) stop('Model needs to be of class "bru_sdm".')
  
  if (!all(dataset%in%model[['source']])) stop('Dataset provided not run in initial model.')
  
  if (length(unique(model[['source']])) == 1) stop('Model was only run with one dataset. Leaving a dataset out is not possible.')

  model_results <- list()
  
  model_differences <- list()
  
  if (!is.null(model[['spatCovs']][['name']])) {
    
    for (names in model[['spatCovs']][['name']]) {
      
      if (!is.null(model[['species']][['speciesIn']])) {
     
      for (species in unique(unlist(model[['species']][['speciesIn']]))) {
       
        assign(paste0(species,'_',names), model$bru_info$model$effects[[paste0(species,'_',names)]]$env[[paste0(species,'_',names)]])
        
      }
        
      }
      else assign(names, model$bru_info$model$effects[[names]]$env[[names]])  
      
    }
    
    if (!is.null(model$bru_info$model$effects[[1]]$env[['spdeModel']])) {
      ##Also assign the speciesModel and the marksModel
      assign('spdeModel', model$bru_info$model$effects[[1]]$env[['spdeModel']])
      
    } else spdeModel <- NULL

    if (is.list(spdeModel[[1]])) {
      ##TODO if we are doing
      for (name in model[['spatial_datasets']]) {
        
        assign(paste0(name,'_spde'),spdemodel[[name]])  
        
      }
      
    }    
    ##Why do I have two validation results
    validation_results <- list()  
    
  }
  else validation_results <- list()
 
  for (dataname in dataset) {
    
    index <- !model[['source']]%in%dataname
    reduced_options <- list()

    reduced_options <- model[['optionsJoint']]
    reduced_options$control.family <- model[['optionsJoint']]$control.family[index]
    
    if (!predictions) {
      
      reduced_options$control.compute <- list(return.marginals = FALSE)
      
    }
    
    reduced_components <- update(model$components, paste0(' ~ . -',dataname,'_intercept(1)'))
    
    ##Need to add this later...
    reduced_components <- update(reduced_components, paste0(' ~ . - ',
                                                            dataname,'_spde(main = coordinates, model = spdemodel)'))
    ##Also need to add biasField #pref something like spatia_datasets...
    
    #if (!is.null(model[['spatial_datasets']])) {
      
    #  if (!dataname%in%model[['spatial_datasets']]) {
        
    #    reduced_components <- update(reduced_components, paste0(' ~ . -','shared_spatial(main = coordinates, model = spdemodel)'))    
        
    #  }  
      
  #  }
    ##Add a new marks used
    if (!is.null(model$marks_used)) {
      
      if (any(names(model$marks_used) == dataname)) {   
        
        for (i in 1:length(model$marks_used)) {
          ##Still need to add markintercepts to the model...
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
            ##Not sure why sum == 1???
            reduced_components <- gsub(paste0('[+] ', mark,'*\\(.*?\\) *'), '', reduced_components, perl = T)
            reduced_components <- gsub(paste0('[+] ', mark,'_phi*\\(.*?\\) *'), '', reduced_components, perl = T)   
            
          }
          
        }
        
      }
      
    }

    if (!is.null(model[['species']][['speciesIn']])) {
      ##Species in all but left out data
      reduced_species <- unlist(unique(model[['species']][['speciesIn']][!names(model[['species']][['speciesIn']])%in%dataname]))
      #species in left out data
      species_dataset <- unlist(unique(model[['species']][['speciesIn']][dataname]))
      #Need diff ## REDO
      species_rm <- intersect(reduced_species, species_dataset)
      species_rm <- species_dataset[!species_dataset %in% reduced_species]

      if (!identical(species_rm, 'character(0)')) {
        
        if (any(!species_dataset%in%reduced_species)) {

          for (species in unique(species_dataset[!species_dataset%in%reduced_species])) {
          
            species_covs <- as.vector(outer(paste0(species,'_'), model[['spatCovs']][['name']], FUN = 'paste0'))

            for (i in 1:length(species_covs)) {
   
              reduced_components <- update(reduced_components, paste('~ . -',  paste0(species_covs[i],'(main = ', species_covs[i], ', model = "linear")')))  
              reduced_components <- update(reduced_components, paste('~ . -',  paste0(species_covs[i],'(main = ', species_covs[i], ', model = "factor_full")')))  
              reduced_components <- update(reduced_components, paste('~ . -',  paste0(species_covs[i],'(main = ', species_covs, ', model = "factor_contrast")')))  
            
              }
     
            reduced_components <- update(reduced_components, paste('~ . -', paste0(species,'_intercept(1)')))
            
          }
          n_species <- as.character(max(as.numeric(unlist(model[['species']][['speciesIn']]))))
          n_species_reduced <- as.character(max(as.numeric(reduced_species)))
          
          #reduced_components <- update(reduced_components, paste('~ . -', paste0(attributes(model)$Species,'_spde(main = coordinates, model = spdeModel, group = ',attributes(model)$Species, ', ngroup = ', n_species,', control.group = ', list(attributes(model)$Speciesmodel), ')')))
          #reduced_components <- update(reduced_components, paste('~ . +', paste0(attributes(model)$Species,'_spde(main = coordinates, model = spdeModel, group = ',attributes(model)$Species, ', ngroup = ', n_species_reduced,', control.group = ', list(attributes(model)$Speciesmodel), ')')))
          
        }
        
      }
      
    }
    stop(return(reduced_components))
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
      
      reduced_lik <- model$bru_info$lhoods
      
      for (data in names(model$bru_info$lhoods)[!index]) { 
        
        train <- predict(model_reduced, data = model$bru_info$lhoods[[data]]$data, formula = eval(parse(text = paste0('~ (',paste(model$spatial_covariates_used, collapse = ' + '),')'))))  
        
        reduced_lik[[data]]$data@data['offset'] <- train$mean
        reduced_lik[[data]]$include_components <- c(reduced_lik[data]$include_components, 'offset')
        
        
      }
      
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
  class(model_results) <- c('datasetOut', 'list')
  
  return(model_results)
  
}

#' Export class bru_sdm_leave_one_out
#' 
#' @export

setClass('datasetOut')

#' Print for bru_sdm_leave_one_out
#' 
#' @export print.datasetOut

#' Export print.bru_sdm_leave_one_out
#' 
#' @exportS3Method 

print.datasetOut <- function(x, ...) {
  
  for (name in 1:length(x)) {
    
    cat('Changes in fixed values by leaving out', gsub("Leaving_out_.*?","",paste0(names(attributes(x)$differences)[name]),':'))
    cat('\n\n')
    
    print(attributes(x)$differences[[name]])
    
    if (!is.null(attributes(x)$validation_results)) {
      cat('\n')
      
      cat('Leave-one-out cross-validation score:', attributes(x)$validation_results[name])
      
      cat('\n\n')
      
    }
  } 
}