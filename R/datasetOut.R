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
        
      if (model$spatial$species) assign(paste0(species,'_field'), model$bru_info$model$effects[[paste0(species,'_field')]]$env[[paste0(species,'_field')]])  
        
      }
        
      }
      else assign(names, model$bru_info$model$effects[[names]]$env[[names]])  
      
    }
    
  }
    
  
  if (model$spatial$points) {
      ##Also assign the speciesModel and the marksModel
      assign('shared_field', model$bru_info$model$effects[[1]]$env[['spdeModel']])
      
    } else spdeModel <- NULL
    
    if (!is.null(model$biasData)) {
      
      for (data in model$biasData) { 
        
      assign(paste0(data, '_bias_field'), model$bru_info$model$effects[[paste0(data,'_bias_field')]]$env[[paste0(data,'_bias_field')]])
        
      }
      
    }
  
  if (model$spatial$marks) {
    
    for (mark in unlist(unique(model$marks$marksIn))) {
      
      assign(paste0(mark,'_field'), model$bru_info$model$effects[[paste0(mark,'_field')]]$env[[paste0(mark,'_field')]])  
      
    }
    
  }

    validation_results <- list()  
    
  
 
  if (!is.null(model$species$speciesVar)) {
    # if species spatial ...
    for (species in unique(unlist(model$species$speciesIn))) {
      
      assign(paste0(species,'_field'), model$bru_info$model$effects[[paste0(species,'_field')]]$env[[paste0(species,'_field')]])
      
    }
    
    #assign('speciesModel', model$bru_info$model$effects[[paste0(model$species$speciesVar,'_spatial')]]$env$speciesModel)
  }
    
  if (model$spatial$temporalSpatial) assign('temporal_field', model$bru_info$model$effects$temporal_field$env$temporal_field)
  
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
    #reduced_components <- update(reduced_components, paste0(' ~ . -', dataname,'_bias_field(main = coordinates, model = biasField)'))
    reduced_components <- update(reduced_components, paste0(' ~ . -', dataname,'_biasField(main = coordinates, model = ', paste0(dataname,'_bias_field)')))
    
    ##Also need to add biasField #pref something like spatia_datasets...
    
    #if (!is.null(model[['spatial_datasets']])) {
    
    #  if (!dataname%in%model[['spatial_datasets']]) {
    
    #    reduced_components <- update(reduced_components, paste0(' ~ . -','shared_spatial(main = coordinates, model = spdemodel)'))    
    
    #  }  
    
    #  }
    ##Add a new marks used
    
    ##Redo
    
    if (!all(is.na(unlist(model$marks$marksIn)))) {
      
      marksOut <- model$marks$marksIn[[dataname]]
      marksIn <- unique(unlist(model$marks$marksIn[!names(model$marks$marksIn)%in%dataname]))
      
      marksRM <- marksOut[!marksOut %in% marksIn]
      
      if (!identical(marksRM, 'charachter(0)')) {
        ##Still need to add markintercepts to the model...
        reduced_components <- update(reduced_components, paste0(' ~ . -', marksRM,
                                                                '_intercept(1)'))
        reduced_components <- update(reduced_components, paste0(' ~ . -', marksRM,
                                                                '_spatial(main = coordinates, model =', paste0(marksRM,'_field)')))
        reduced_components <- update(reduced_components, paste0(' ~ . -', marksRM,
                                                                paste0('_phi(main =', marksRM, '_phi, model = "iid", initial = -10, fixed = TRUE)')))
        reduced_components <- update(reduced_components, paste0(' ~ . -', marksRM,
                                                                paste0('(main =', marksRM, ', model = "iid", constr = FALSE, fixed = TRUE)')))
        
        
      }
      
      
    }
    
    if (!is.null(model[['species']][['speciesIn']])) {
      ##Species in all but left out data
      reduced_species <- unlist(unique(model[['species']][['speciesIn']][!names(model[['species']][['speciesIn']])%in%dataname]))
      #species in left out data
      species_dataset <- unlist(unique(model[['species']][['speciesIn']][dataname]))
      #Need diff ## REDO
      #species_rm <- intersect(reduced_species, species_dataset)
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
            
            reduced_components <- update(reduced_components, paste('~ . -', paste0(species,'_spatial(main = coordinates, model = ', paste0(species,'_field)'))))
            reduced_components <- update(reduced_components, paste('~ . -', paste0(species,'_intercept(1)')))
            
          }

        }
        
      }
      
    }
    
    model_reduced <- inlabru::bru(components = formula(reduced_components),
                                  model$bru_info$lhoods[index],
                                  options = reduced_options)
    
    model_reduced[['componentsJoint']] <- reduced_components
    model_reduced[['optionsJoint']] <- reduced_options
    model_reduced[['spatCovs']] <- model[['spatCovs']]
    model_reduced[['species']] <- list(speciesIn = model[['species']][['speciesIn']][!names(model[['species']][['speciesIn']]) %in% dataset],
                                       speciesVar = model[['species']][['speciesVar']])
    
    model_reduced[['dataType']] <- model[['dataType']][index]
    model_reduced[['source']] <- model[['source']][index]
    model_reduced[['marks']] <- list(marksIn = model[['marks']][['marksIn']][!names(model[['marks']][['marksIn']]) %in% dataset],
                                 multinomVars = model[['marks']][['multinomVars']])

    if (!predictions) {
      
      model_reduced[['bru_info']] <- NULL
      model_reduced[['call']] <- NULL
      
    }
    

    
    class(model_reduced) <- c('bru_sdm',class(model_reduced))
    
    model_results[[paste0('Leaving_out_',dataname)]] <- model_reduced
    
    var_names <- row.names(model_reduced$summary.fixed)[row.names(model_reduced$summary.fixed)%in%row.names(model$summary.fixed)]
    model_differences[[paste0('Leaving_out_',dataname)]] <- model_reduced$summary.fixed[var_names,] - model$summary.fixed[var_names,]
    
    
    if (predictions) {
      
      reduced_lik <- model$bru_info$lhoods
      
      for (data in names(model$bru_info$lhoods)[!index]) { 
        
        if (!is.null(model[['species']][['speciesIn']])) covs <- as.vector(outer(paste0(reduced_species,'_'), model$spatCovs$name, FUN = 'paste0'))
        else covs <- model$spatCovs$name
        
        train <- predict(model_reduced, data = model$bru_info$lhoods[[data]]$data, formula = eval(parse(text = paste0('~ (',paste(covs, collapse = ' + '),')'))))  
        
        reduced_lik[[data]]$data@data['offset'] <- train$mean
        #reduced_lik[[data]]$formula <- update(reduced_lik[[data]]$formula, ~ . + offset)
        reduced_lik[[data]]$include_components <- c(reduced_lik[[data]]$include_components, 'offset')
        
      }
     
      #for (data in names(model$bru_info$lhoods[index])) {
    #  
    #  reduced_lik[[data]]$data@data['offset'] <- 0
    #  reduced_lik[[data]]$formula <- update(reduced_lik[[data]]$formula, ~ . + offset)

#      }
      offset_components <- update(model$componentsJoint, ~ . + offset)

      reduced_mlik <- inlabru::bru(offset_components, reduced_lik,
                                   options = model$optionsJoint)
      
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
  
  model_results
  
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
    
    cat('Changes in fixed values by leaving out', gsub("Leaving_out_.*?","",paste0(names(attributes(x)$differences)[name],':')))
    cat('\n\n')
    
    print(attributes(x)$differences[[name]])
    
    if (!is.null(attributes(x)$validation_results)) {
      cat('\n')
      
      cat('Leave-one-out cross-validation score:', attributes(x)$validation_results[name])
      
      cat('\n\n')
      
    }
  } 
}