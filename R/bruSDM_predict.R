#' Predict for bru_sdm
#' @param object A \code{bru_sdm object}.
#' @param data Data containing points of the map with which to predict on. May be \code{NULL} if one of \code{mesh} or \code{mask} is \code{NULL}.
#' @param formula Formula to predict. May be \code{NULL} if other arguments: \code{covariates}, \code{spatial}, \code{intercepts} are not \code{NULL}.
#' @param mesh An inla.mesh object
#' @param mask A mask of the study background. Defaults to \code{NULL}.
#' @param covariates Name of covariates to predict.
#' @param temporal Make predictions for the temporal component of the model.
#' @param spatial Include spatial effects in prediction.
#' @param intercepts Include intercept in prediction.
#' @param datasets Names of the datasets to include intercept term. Default of \code{NULL} results in all datasets being predicted.
#' @param species Names of the species to predict. Default of \code{NULL} results in all species being predicted.
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
                           mask = NULL, temporal = FALSE, covariates = NULL, spatial = FALSE,
                           intercepts = FALSE, datasets = NULL, species = NULL,
                           biasfield = FALSE, biasnames = NULL, predictor = FALSE,
                           fun = 'exp', ...) {
  
  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an inla.mesh object is required.")
  
  ## if datasets !is.null but at least one not in model stop
  # else datasets <- all datasets in the model.
  ##Need another defense: only one of predictor; biasfield; temporal
  if (sum(predictor, temporal, biasfield) > 1) stop('Please only choose one of: predictor, temporal and biasfield.')
  ## if non-null biasfields ## if no bias fields in stop: if biasnames not in biasfields stop
  if (biasfield && spatial) stop('Please choose one of biasfield and spatial.')
  
  if (!is.null(unlist(model[['species']][['speciesIn']]))) {
    
    speciespreds <- TRUE
    
    if (is.null(species)) speciesin <- unique(unlist(model[['species']][['speciesIn']]))
    if (!all(species %in% unique(unlist(model[['species']][['speciesIn']])))) stop('Species provided not in model.')
    
  }
  else {
    
    speciespreds <- FALSE
    if (intercepts) intercept_terms <- paste0(datasets, '_intercept')
    
  }
  
  if (!intercepts) intercept_terms <- NULL
  
  if (!all(covariates%in%row.names(model$summary.fixed)) && !all(as.vector(outer(paste0(unlist(model[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))%in%row.names(model$summary.fixed))) stop("Covariates provided not in model.")
  
  if (is.null(formula) && !intercepts && !spatial && is.null(covariates) && !temporal && !biasfield) stop("Please provide either a formula or components of a formula to be predicted.")
  
  if (temporal && is.null(model$temporal$temporalVar)) stop('Temporal is set to TRUE but no temporal component found in the model.')
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- pixels(mesh, mask = mask)
      
    }   
    else data <- pixels(mesh)
  }
  
  
  if (is.null(formula)) {
    
    int <- list()
    
    class(model) <- c('bru','inla','iinla')
    
    if (is.null(fun) | fun == 'linear') fun <- ''
    
    if (temporal) {
      
      numeric_time <- order(as.numeric(unique(unlist(model$temporal$temporalIn))))
      time_variable <- model$temporal$temporalVar
      
      time_data <- data.frame(seq_len(max(numeric_time)))
      names(time_data) <- time_variable
      
      timeData <- inlabru::cprod(data, data.frame(time_data))
      names(timeData@data) <- c(time_variable, 'weight')  
      
      if (!is.null(covariates)) covariates <- as.vector(outer(paste0(unlist(model[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))
      if (intercepts) intercept_terms <- paste0(unlist(model[['species']][['speciesIn']]), '_intercept')
      
      time_formula <- paste0(fun,'(',paste(covariates, intercept_terms, 'shared_spatial', collapse = ' + '),')')

      int[['temporalPredictions']] <- predict(model, timeData, ~ data.frame(time_variable = eval(parse(text = time_variable)), formula = eval(parse(text = time_formula))))
      int[['temporalPredictions']] <- int[['temporalPredictions']][,!names(int[['temporalPredictions']]@data) %in% time_variable]
      class(int) <- c('bruSDM_predict', class(int))
      
      return(int)  
      
    }
    
    if (biasfield) {
      
      if (is.null(biasnames)) biasnames <- model$biasData
      
      if (!all(paste0(biasnames,'_biasField') %in% names(model$summary.random))) stop('Either no bias field has been used or an incorrect dataset name was given.')
      
      int[['biasFields']] <- vector(mode = 'list', length = length(biasnames))
      
      for (bias in biasnames) {
        
        formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(paste0(bias,'_biasField'),')')))
        int[[1]][[bias]] <- predict(model, data = data, formula = formula, ...)
        
      }
      
      names(int[[1]]) <- paste0(biasnames, '_biasField')
      
      class(int) <- c('bruSDM_predict', class(int))
      return(int) 
      
    }
    
    if (speciespreds && ! predictor) {
      
      int[['speciesPredictions']] <- vector(mode = 'list', length(unlist(unique(model$species$speciesIn))))
      names(int[['speciesPredictions']]) <- speciesin
      
      for (spec in speciesin) {
        
        if (!is.null(covariates)) species_covs <- paste0(spec, '_', covariates)
        else species_covs <- NULL
        
        if (intercepts) species_int <- paste0(spec,'_intercept')
        else species_int <- NULL
        
        if (spatial) species_spat <- paste0(spec,'_spatial')
        else species_spat <- NULL
        
        species_formula <- formula(paste0('~', fun, '(', paste0(c(species_covs, species_int, species_spat), collapse = ' + '),')'))
        
        int[[1]][[spec]] <- predict(model, data, formula = species_formula, ...)
      
        }
      
      class(int) <- c('bruSDM_predict', class(int))
      return(int) 
      
      
    }
    
    if (spatial) {
      
      if (!'shared_spatial' %in% names(model$summary.random)) stop('Model run without spatial effects. Please specify Spatial = TRUE in bruSDM.')
      else spatial_obj <- 'shared_spatial'
      
    } 
    else spatial_obj <- NULL
    
    if (predictor) formula_components <- c(row.names(model$summary.fixed), names(model$summary.random))
    else formula_components <- c(covariates, intercept_terms, spatial_obj)
    
    if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial, intercepts or biasfield.')
    
    formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
    
    #int[[i]] <- predict(model, data = data, formula = formula, ...)
    int <- predict(model, data = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    
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
    if (names(x)[[1]] == 'speciesPredictions') {
      
      for (species in names(x[[1]])) {
        
        cat('Predictions for', paste0(species,':'))
        print(summary(x[[1]][[species]]@data))
        cat('\n')
         
        }
      
    }
    else
      if (names(x)[[1]] == 'temporalPredictions') {
        
        cat('Predictions for the temporal variable:')
        cat('\n\n')
        print(summary(x[[1]]@data))
        
      }
    else print(summary(x[['predictions']]@data))
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
  
  
  if (length(x) == 1 & names(x) == 'speciesPredictions') {
    
    if (length(whattoplot) > 1) stop('Please only plot one variable at a time for species plots.')
    
    species_title <- ggtitle('Plot of the species predictions')
    
    if (!is.null(colours)) {
      
      plot_colours <- scale_fill_gradientn(colours = rev(brewer.pal(9,colours)),
                                           limits = range(x[['Species predictions']]@data[,stat]))
      
    }
    else plot_colours <- NULL
    
    all_plots <- list()
    
    for (species in names(x$speciesPredictions)) {
      
      all_plots[[species]] <- ggplot() + gg(x$speciesPredictions[[species]], aes_string(fill = whattoplot)) + ggtitle(paste('Plot of predictions for', species)) + plot_colours
      
    }

    if (plot) {
      
      plots <- inlabru::multiplot(plotlist = all_plots, cols = length(all_plots), layout = layout)
      return(plots)
      
      }
    else return(all_plots)
    
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