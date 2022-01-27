#' Export class predict_bru_sdm
#' 
#' @export 

setClass('bruSDM_predict')

#' Predict for bru_sdm
#' @param object A \code{bru_sdm object}.
#' @param data Data containing points of the map with which to predict on. May be \code{NULL} if one of \code{mesh} or \code{mask} is \code{NULL}.
#' @param formula Formula to predict. May be \code{NULL} if other arguments: \code{covariates}, \code{spatial}, \code{intercept} are not \code{NULL}.
#' @param mesh An inla.mesh object
#' @param mask A mask of the study background. Defaults to \code{NULL}.
#' @param species Create plot of species. Note: \code{speciesname} in the \code{organize_data} function needs to be specified first.
#' @param covariates Name of covariates to predict.
#' @param spatial Include spatial effects in prediction.
#' @param intercept Include intercept in prediction.
#' @param interceptnames Names of the datasets to include intercept term. Defaults to \code{NULL}.
#' @param biasfield Include bias field in prediction.
#' @param biasnames Names of the datasets to include bias term. Defaults to \code{NULL}.
#' @param predictor Should all terms run in the linear predictor be included. Defaults to \code{FALSE}.
#' @param fun Function used to predict. Set to \code{'linear'} if effects on the linear scale are desired.
#' @param ... Additional arguments used by the inlabru \code{predict} function.
#' 
#' @method predict bruSDM
#' @export
#' 

predict.bruSDM <- function(model, data = NULL, formula = NULL, mesh = NULL, 
                           mask = NULL, species = FALSE,
                           covariates = NULL, spatial = TRUE,
                           intercept = FALSE, interceptnames = NULL,
                           bisafield = FALSE, biasnames = NULL, predictor = FALSE,
                           fun = 'exp', ...) {

  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an inla.mesh object is required.")
  
  speciesin <- factor((unlist(model[['species']][['speciesIn']])))
  
  if (!all(covariates%in%row.names(model$summary.fixed)) && !all(as.vector(outer(paste0(speciesin,'_'), covariates, FUN = 'paste0'))%in%row.names(model$summary.fixed))) stop("Covariates provided not in model.")
  
  if (is.null(formula) & !spatial & !intercept & is.null(covariates)) stop("Please provide at least one of spatial, intercept or covariates.")
  
  if (species & is.null(model[['species']]['speciesIn'])) stop('Predictions for species was selected but no species were included in model. Please provide a species name in bruSDM with the argument "species".')
  
  if (!species) {
    
    if (is.null(formula) && !intercept && !spatial && is.null(covariates)) stop("Please provide either a formula or components of a formula to be predicted.")
    
  }
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- pixels(mesh, mask = mask)
      
    }   
    else data <- pixels(mesh)
  }
 

  if (is.null(formula)) {
    
    int <- list()
    
    class(model) <- c('bru','inla','iinla')
    
    if (species) {
      
      numeric_species <- as.numeric(speciesin)
      
      species_variable <- model[['species']][['speciesVar']]
      
      species_data <- data.frame(seq_len(max(numeric_species)))
      names(species_data) <- species_variable
     
      speciesData <- inlabru::cprod(data, data.frame(species_data))
      names(speciesData@data) <- c(species_variable, 'weight')
    
      if (is.null(fun) | fun == 'linear') {fun <- ''}
      
      ##How do we do this?? species covs are now: speciesname_covariatename??
      if (!is.null(covariates)) {
      
        speciescovs <- as.vector(outer(paste0(unique(as.character(speciesin)),'_'), model[['spatCovs']][['name']], FUN = 'paste0'))

        species_formula <- paste0(fun,'(',paste(paste0(species_variable,'_spatial'), paste(speciescovs, collapse = ' + '), sep = ' + '),')')
        
      }
      else  species_formula <- paste0(fun,'(',paste(paste0(species_variable,'_spatial')),')')

      int <- predict(model, speciesData, ~ data.frame(species_variable = eval(parse(text = species_variable)), formula = eval(parse(text = species_formula))))
      
      int[[species_variable]] <- as.character(rep(unique(speciesin[order(speciesin)])), length(data)) ##order?
      
      int <- list(int)
      names(int) <- 'Species predictions'
      class(int) <- c('bruSDM_predict', class(int))
      
      return(int)
      
    }
    else {
      
      if (predictor) formula_components <- c(row.names(model$summary.fixed), names(model$sumary.random))
        
      else{
        if (spatial) {
          
          if (!'shared_spatial' %in% names(model$summary.random)) stop('Model run without spatial effects. Please specify Spatial = TRUE in bruSDM.')
          else spatial_obj <- 'shared_spatial'
            
          } 
        else spatial_obj <- NULL
        
        if (intercept) {
          
          if (is.null(interceptnames)) interceptnames <- unique(model$source)
          
          if (!all(paste0(interceptnames,'_intercept')%in%row.names(model$summary.fixed))) stop('Either dataset name is incorrect or bru_sdm model run without intercepts.')
          else intercept_obj <- paste0(interceptnames,'_intercept')
            
        } 
        else intercept_obj <- NULL
        
        if (!is.null(covariates)) {
          
          if (!identical(as.character(speciesin), character(0))) covariates <- as.vector(outer(paste0(unique(as.character(speciesin)),'_'), model[['spatCovs']][['name']], FUN = 'paste0'))
          
        }
        
        if (biasfield) {
          
          if (is.null(biasnames)) bias_obj <- names(model$summary.random)[grepl('_bias_field$', names(model$summary.ramdom))]
          else bias_obj <- paste0(biasnames, '_bias_field')
          
          if (!all(bis_obj %in% names(model$summary.random))) stop('Either no bias field has been used or an incorrect dataset name was given.')
          
        }
        else bias_obj <- NULL
        
        }
      
        formula_components <- c(covariates, spatial_obj, intercept_obj, bias_obj)
      
        if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial, intercept or biasfield.')
        if (is.null(fun) | fun == 'linear') {fun <- ''}
  
        formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))

        #int[[i]] <- predict(model, data = data, formula = formula, ...)
        int <- predict(model, data = data, formula = formula, ...)
        int <- list(int)
        names(int) <- 'predictions'
      }
      
      #names(int) <- datasetNames
      class(int) <- c('bruSDM_predict', class(int))
      return(int)
    
  }
  
  else {
    
    class(model) <- c('bru','inla','iinla')
    int <- predict(model, data = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    class(int) <- c('bruSDM_predict', class(int))
    
    return(int)
    
  }
  
}

print.bruSDM_predict <- function(x, ...) {

    
    cat('Summary of predicted data:')
    cat('\n\n')
    print(summary(x[['predictions']]@data))
    cat('\n\n')
    
  
}

#' Plot for preduct_bru_sdm
#' 
#' @exportS3Method 

plot.bruSDM_predict <- function(x, plotall = TRUE,
                                datasettoplot = NULL,
                                whattoplot = c('mean'),
                                colours = NULL,
                                cols = NULL,
                                layout = NULL,
                                plot = TRUE,
                                ...) {
  ## Add species plot::
  #if (!plotall & is.null(datasettoplot)) stop('Please provide a list of datasets to plot or set plotall to TRUE.')
  
  #if (any(!datasettoplot%in%names(x))) stop('Dataset name to plot not provided in prediction object.')
  
  #if (plotall) datasettoplot <- names(x)
  
  if (any(!whattoplot%in%c("mean", "sd", "q0.025", "median","q0.975",
                           "smin", "smax", "cv", "var" ))) stop('Whattoplot is not a valid variable to plot')
  
  
  if (length(x) == 1 & names(x) == 'Species predictions') {
    
    if (length(whattoplot) > 1) stop('Please only plot one variable at a time for species plots.')
    
    species_name_var <- names(x[[1]]@data)[sapply(x[[1]]@data, class) == 'character']
    
    species_title <- ggtitle('Plot of the species predictions')
    
    if (!is.null(colours)) {
      
      plot_colours <- scale_fill_gradientn(colours = rev(brewer.pal(9,colours)),
                                           limits = range(x[['Species predictions']]@data[,stat]))
      
    }
    else plot_colours <- NULL
    
    plot_grid <- ggplot() + gg(x[[1]], aes_string(fill = whattoplot)) + facet_grid(~ eval(parse(text = species_name_var))) + plot_colours + species_title
    return(plot_grid)
  }
  
  if (!plot) {
    
    all_plots <- list()
    prediction_list <- list()
  }
  datasettoplot <- 'predictions'
  for (plotname in datasettoplot) {
    
    plot_list <- list()
    
    for (stat in whattoplot) {
      
      #title <- ggtitle(paste('Plot of',stat,'for',plotname))
      title <- ggtitle('Plot of predictions')
      
      prediction <- gg(x[[plotname]], aes_string(fill = stat))
      
      if (!plot) prediction_list[[stat]] <- prediction
      
      if (!is.null(colours)) {
        
        plot_colours <- scale_fill_gradientn(colours = rev(brewer.pal(9,colours)),
                                             limits = range(x[[plotname]]@data[,stat]))
        
      }
      else plot_colours <- NULL
      
      
      plot_list[[stat]] <- ggplot() + prediction + title + plot_colours
      
    }
    
    if (plot) {
      if (is.null(cols)) cols <- length(whattoplot)
      plot_grid <- inlabru::multiplot(plotlist = plot_list, cols = cols, layout = layout)
      
    }
    else{
      
      all_plots[[plotname]] <- prediction_list
      
    }
    
  }
  
  if (!plot) return(all_plots)
  
}