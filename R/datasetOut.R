#' @title \emph{datasetOut}: function that removes a dataset out of the main model, and calculates some cross-validation score.
#' 
#' @description This function calculates the difference in covariate values between a full integrated model and a model with one dataset left out, as well as some cross-validation score, which is used to obtain a score of the relative importance of the dataset in the full model. The score is calculated as follows:
#' \enumerate{
#' 
#'   \item Running a new model with one less dataset (from the main model) -- resulting in a reduced model,
#'   \item predicting the intensity function at the locations of the left-out dataset with the reduced model,
#'   \item using the predicted values as an offset in a new model,
#'   \item finding the difference between the marginal-likelihood of the main model (ie the model with all the datasets considered) and the marginal-likelihood of the offset model.
#' 
#' }
#' 
#' 
#' @param model Model of class modISDM run with multiple datasets.
#' @param dataset Names of the datasets to leave out. If missing, will run for all datasets used in the full model.
#' @param predictions Will new models be used for predictions. If \code{TRUE} returns marginals and bru_info in model. Defaults to \code{TRUE}. 
#' 
#' @return A list of inlabru models with the specified dataset left out. If predictions is \code{FALSE}, these objects will be missing their \code{bru_info} and \code{call} lists.
#' 
#' @import inlabru
#' 
#' @export
#' 
#' @examples 
#' 
#'\dontrun{
#'  if (requireNamespace('INLA')) {
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
#'                             Projection = proj, 
#'                             responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- fitISDM(organizedData,
#'               options = list(control.inla = list(int.strategy = 'eb')))
#'    
#'   #Choose dataset to leave out
#'   eBirdOut <- datasetOut(modelRun, dataset = 'eBird')
#'   
#'   #Print datasetOut summary
#'   eBirdOut
#'   
#'   
#'    
#'  }
#'}

datasetOut <- function(model, dataset,
                       predictions = TRUE) {
  
  if (!inherits(model, 'bruSDM') && !inherits(model, 'modISDM')) stop('model needs to be either a bruSDM or modISDM  object.')
  
  if (missing(dataset)) dataset <- unique(model[['source']])
  
  if (!all(dataset%in%model[['source']])) stop('Dataset provided not run in initial model.')
  
  if (length(unique(model[['source']])) == 1) stop('Model was only run with one dataset. Leaving a dataset out is not possible.')

  model_results <- list()
  
  model_differences <- list()
  ##re-write all of this with data2env
  
  if (!is.null(model$spatial$points)) {
      
      if (model$spatial$points %in% c('shared', 'correlate')) assign('shared_field', model$bru_info$model$effects[[1]]$env[['shared_field']])
      else {
        
        for (data in unique(model$source)) {
          
          assign(paste0(data, '_field'),  model$bru_info$model$effects[[paste0(data,'_spatial')]]$env[[paste0(data,'_field')]])
          
        }
        
      }
      
    } else spdeModel <- NULL
    
    if (!is.null(model$biasData)) { #if shared 
      
      for (data in model$biasData$Fields) { 
        
      assign(paste0(data, '_bias_field'), model$bru_info$model$effects[[paste0(data,'_biasField')]]$env[[paste0(data,'_bias_field')]])
        
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
    
  #if (model$spatial$temporalSpatial) assign('temporal_field', model$bru_info$model$effects$temporal_field$env$temporal_field)
  
  for (dataname in dataset) {
    
    index <- !model[['source']]%in%dataname
    reduced_options <- list()

    reduced_options <- model[['optionsJoint']]
    reduced_options$control.family <- model[['optionsJoint']]$control.family[index]
    
    if (!predictions) {
      
      reduced_options$control.compute <- list(return.marginals = FALSE)
      
    }

    
    reduced_terms <- unique(unlist(lapply(model$bru_info$lhoods[index], function(x) {
      
      if (!identical(unlist(x$used), character(0))) unlist(x$used)
      else labels(terms(x$formula))
      
    })))
    
    reduced_components <- reduceComps(componentsOld = model$componentsJoint,
                          pointsCopy = ifelse(model$spatial$points == 'copy', 
                                      TRUE, FALSE),
                          biasCopy = model$biasData$Copy,
                          datasetName = dataname,
                          reducedTerms = reduced_terms)
    
    model_reduced <- inlabru::bru(components = reduced_components,
                                  model$bru_info$lhoods[index],
                                  options = reduced_options)
  
    model_reduced[['componentsJoint']] <- reduced_components
    model_reduced[['optionsJoint']] <- reduced_options
    model_reduced[['spatCovs']] <- model[['spatCovs']]
    model_reduced[['spatial']] <- model[['spatial']]
    model_reduced[['species']] <- list(speciesIn = model[['species']][['speciesIn']][!names(model[['species']][['speciesIn']]) %in% dataset],
                                       speciesVar = model[['species']][['speciesVar']],
                                       speciesEffects = list(Intercepts = model[['species']][['speciesEffects']][['Intercepts']],
                                                             Environmental = model[['species']][['speciesEffects']][['Environmental']]))
    model_reduced[['species']][['speciesTable']] <- model[['species']][['speciesTable']][model$species$speciesTable[,2] %in% unique(unlist(model_reduced$species$speciesIn)),]
    
    model_reduced[['dataType']] <- model[['dataType']][index]
    model_reduced[['source']] <- model[['source']][index]
    model_reduced[['marks']] <- list(marksIn = model[['marks']][['marksIn']][!names(model[['marks']][['marksIn']]) %in% dataset],
                                 multinomVars = model[['marks']][['multinomVars']])

    if (!predictions) {
      
      model_reduced[['bru_info']] <- NULL
      model_reduced[['call']] <- NULL
      
    }
    

    
    class(model_reduced) <- c('bruSDM',class(model_reduced))
    
    model_results[[paste0('Leaving_out_',dataname)]] <- model_reduced
    
    var_names <- row.names(model_reduced$summary.fixed)[row.names(model_reduced$summary.fixed)%in%row.names(model$summary.fixed)]
    model_differences[[paste0('Leaving_out_',dataname)]] <- model_reduced$summary.fixed[var_names,] - model$summary.fixed[var_names,]
    
    
    if (predictions) {
      
      reduced_lik <- model$bru_info$lhoods
      
      for (data in names(model$bru_info$lhoods)[!index]) { 
        
        #if (!is.null(model[['species']][['speciesIn']])) covs <- as.vector(outer(paste0(reduced_species,'_'), model$spatCovs$name, FUN = 'paste0'))
        #else covs <- model$spatCovs$name
        train <- predict(model_reduced, data = model$bru_info$lhoods[[data]]$data, formula = eval(parse(text = paste0('~(',paste(gsub('\\(.*$', '', labels(terms(reduced_components))), collapse = ' + '),')'))))  
       
        reduced_lik[[data]]$data[,'offset'] <- train$predictions$mean
        #reduced_lik[[data]]$formula <- update(reduced_lik[[data]]$formula, ~ . + offset)
        reduced_lik[[data]]$used$effect <- c(reduced_lik[[data]]$used$effect, 'offset')
        
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
#' @title Generic print function for \code{datasetOut}.
#' @param x datasetOut object.
#' @param ... Unused argument.
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