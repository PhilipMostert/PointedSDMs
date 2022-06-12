
#' Export class predict_bru_sdm
#' 
#' @export 

setClass('bruSDM_predict')

#' Predict for bru_sdm
#' @title Generic predict function for \code{bru_SDM} objects.
#' @description Predict function for the object produced by \code{\link{runModel}}. Should act identically to \pkg{inlabru}'s generic predict function if wanted, but has additional arguments to help predict certain components created by the model. This is needed since \code{\link{intModel}} creates variable names which might not be directly apparent to the user.
#' @param object A \code{bru_sdm} objects.
#' @param data Data containing points of the map with which to predict on. May be \code{NULL} if one of \code{mesh} or \code{mask} is \code{NULL}.
#' @param formula Formula to predict. May be \code{NULL} if other arguments: \code{covariates}, \code{spatial}, \code{intercepts} are not \code{NULL}.
#' @param mesh An \code{inla.mesh} object.
#' @param mask A mask of the study background. Defaults to \code{NULL}.
#' @param covariates Name of covariates to predict.
#' @param temporal Make predictions for the temporal component of the model.
#' @param spatial Logical: include spatial effects in prediction. Defaults to \code{FALSE}.
#' @param intercepts Logical: include intercept terms in prediction. Defaults to \code{FALSE}.
#' @param datasets Names of the datasets to include intercept and spatial term.
#' @param species Names of the species to predict. Default of \code{NULL} results in all species being predicted.
#' @param biasfield Logical include bias field in prediction. Defaults to \code{FALSE}.
#' @param biasnames Names of the datasets to include bias term. Defaults to \code{NULL}. Note: the chosen dataset needs to be run with a bias field first; this can be done using \code{.$addBias} with the object produced by \code{\link{intModel}}.
#' @param predictor Should all terms run in the linear predictor be included in the predictions. Defaults to \code{FALSE}.
#' @param fun Function used to predict. Set to \code{'linear'} if effects on the linear scale are desired.
#' @param ... Additional arguments used by the inlabru \code{predict} function.
#' 
#' @method predict bruSDM
#' @rdname predict
#' @return A list of inlabru predict objects.
#' @export
#' 
#' @examples 
#'
#' \dontrun{
#'  
#'  if (requireNamespace('INLA')) {
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  proj <- CRS("+proj=longlat +ellps=WGS84")
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- runModel(organizedData, options = list(control.inla = list(int.strategy = 'eb')))
#'    
#'   #Predict spatial field on linear scale
#'   predictions <- predict(modelRun, mesh = mesh, spatial = TRUE, fun = 'linear')
#'    
#'  }
#'}
#' 

predict.bruSDM <- function(object, data = NULL, formula = NULL, mesh = NULL, 
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
  
  if (is.null(datasets)) datasets <- unique(object$source)
  
  if (!is.null(unlist(object[['species']][['speciesIn']]))) {
    
    speciespreds <- TRUE
    
    if (is.null(species)) speciesin <- unique(unlist(object[['species']][['speciesIn']]))
    if (!all(species %in% unique(unlist(object[['species']][['speciesIn']])))) stop('Species provided not in model.')
    
  }
  else {
    
    speciespreds <- FALSE
    if (intercepts) intercept_terms <- paste0(datasets, '_intercept')
    
  }
  
  if (!intercepts) intercept_terms <- NULL
  
  if (!all(covariates%in%row.names(object$summary.fixed)) && !all(as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))%in%row.names(object$summary.fixed))) stop("Covariates provided not in model.")
  
  if (is.null(formula) && !intercepts && !spatial && is.null(covariates) && !temporal && !biasfield && !predictor) stop("Please provide either a formula or components of a formula to be predicted.")
  
  if (temporal && is.null(object$temporal$temporalVar)) stop('temporal is set to TRUE but no temporal component found in the model.')
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- inlabru::pixels(mesh, mask = mask)
      
    }   
    else data <- inlabru::pixels(mesh)
  }
  
  
  if (is.null(formula)) {
    
    int <- list()
    
    class(object) <- c('bru','inla','iinla')
    
    if (is.null(fun) | fun == 'linear') fun <- ''
    
    if (temporal) {
      
      numeric_time <- order(as.numeric(unique(unlist(object$temporal$temporalIn))))
      time_variable <- object$temporal$temporalVar
      
      time_data <- data.frame(seq_len(max(numeric_time)))
      names(time_data) <- time_variable
      
      timeData <- inlabru::cprod(data, data.frame(time_data))
      names(timeData@data) <- c(time_variable, 'weight')  
      
      if (!is.null(covariates) && speciespreds) covariates <- as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))
      if (intercepts) intercept_terms <- paste0(unlist(object[['species']][['speciesIn']]), '_intercept')
      
      time_formula <- paste0(fun,'(',paste0(c(covariates, intercept_terms, 'shared_spatial'), collapse = ' + '),')')

      formula = formula(paste('~',paste0('data.frame(', time_variable,' = ', time_variable, ',formula =', time_formula,')')))

      #int[['temporalPredictions']] <- predict(object, timeData, ~ data.frame(..temporal_variable_index.. = eval(parse(text = time_variable)), formula = eval(parse(text = time_formula))))
      int[['temporalPredictions']] <- predict(object, timeData, formula)
      
      #int[['temporalPredictions']] <- int[['temporalPredictions']][,!names(int[['temporalPredictions']]@data) %in% time_variable]
    
      class(int) <- c('bruSDM_predict', class(int))
      
      return(int)  
      
    }
    
    if (biasfield) {
      
      if (is.null(biasnames)) biasnames <- object$biasData
      
      if (!all(paste0(biasnames,'_biasField') %in% names(object$summary.random))) stop('Either no bias field has been used or an incorrect dataset name was given.')
      
      for (bias in biasnames) {
        
        formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(paste0(bias,'_biasField'),')')))
        int[['biasFields']][[bias]] <- predict(object, data = data, formula = formula, ...)
        
      }
      
      class(int) <- c('bruSDM_predict', class(int))
      return(int) 
      
    }
    
    if (speciespreds && ! predictor) {
      
      int[['speciesPredictions']] <- vector(mode = 'list', length(unlist(unique(object$species$speciesIn))))
      names(int[['speciesPredictions']]) <- speciesin
      
      for (spec in speciesin) {
        
        if (!is.null(covariates)) species_covs <- paste0(spec, '_', covariates)
        else species_covs <- NULL
        
        if (intercepts) species_int <- paste0(spec,'_intercept')
        else species_int <- NULL
        
        if (spatial) species_spat <- paste0(spec,'_spatial')
        else species_spat <- NULL
        
        species_formula <- formula(paste0('~', fun, '(', paste0(c(species_covs, species_int, species_spat), collapse = ' + '),')'))
        
        int[[1]][[spec]] <- predict(object, data, formula = species_formula, ...)
      
        }
      
      class(int) <- c('bruSDM_predict', class(int))
      return(int) 
      
      
    }
    
    if (spatial) {
      
      if ('shared_spatial' %in% names(object$summary.random))  spatial_obj <- 'shared_spatial'
      else 
        if (!all(paste0(datasets,'_spatial') %in% names(object$summary.random))) stop('Spatial effects not provided in intModel.')
      else spatial_obj <- paste0(datasets, '_spatial')
      
    } 
    else spatial_obj <- NULL
    
    if (predictor) formula_components <- c(row.names(object$summary.fixed), names(object$summary.random))
    else formula_components <- c(covariates, intercept_terms, spatial_obj)
    
    if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial, intercepts or biasfield.')
    
    formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
    
    #int[[i]] <- predict(object, data = data, formula = formula, ...)
    int <- predict(object, data = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    
    class(int) <- c('bruSDM_predict', class(int))
    return(int)
    
  }
  else {
    
    class(object) <- c('bru','inla','iinla')
    int <- predict(object, data = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    class(int) <- c('bruSDM_predict', class(int))
    
    return(int)
    
  }
  
  
  
}

#' @title Generic print function for \code{bru_sdm_predict}.
#' @param x bruSDM_predict object
#' @param ... Not used.

print.bruSDM_predict <- function(x, ...) {

    
    cat('Summary of predicted data:')
    cat('\n\n')
    if (names(x)[[1]] == 'speciesPredictions') {
      
      for (species in names(x[[1]])) {
        
        cat('Predictions for', paste0(species,':'))
        cat('\n')
        print(summary(x[[1]][[species]]@data))
        cat('\n')
         
        }
      
    }
    else
      if (names(x)[[1]] == 'temporalPredictions') {
        
        cat('Predictions for the temporal variable:')
        cat('\n')
        print(summary(x[[1]]@data))
        
      }
    else
      if (names(x)[[1]] == 'biasFields') {
        
        for(bias in names(x[['biasFields']])) {
          
          cat('Predictions of the bias field for', paste0(bias,':'))
          cat('\n')
          print(summary(x[[1]][[bias]]@data))
          cat('\n')
          
        }
        
      }
    else print(summary(x[['predictions']]@data))
    cat('\n\n')
    
  
}

#' Plot for preduct_bru_sdm
#' @title Generic plot function for \code{predict_bru_sdm}.
#' @param x A bruSDM_predict object.
#' @param whattoplot One of the following statistics to plot: "mean", "sd", "q0.025", "median","q0.975", "smin", "smax", "cv", "var" 
#' @param cols Number of columns required for the plotting. Used by inlabru's multiplot function.
#' @param layout Layout of the plots. Used by inlabru's multiplot function.
#' @param colourLow Colour for the low values in the predictions (see ?scale_colour_gradient from \code{ggplot2}). Defaults to \code{NULL}. If non-\code{NULL}, \code{colourHigh} is required.
#' @param colourHigh Colour for the high values in the predictions (see ?scale_colour_gradient from \code{ggplot2}). Defaults to \code{NULL}. If non-\code{NULL}, \code{colourLow} is required.
#' @param plot Should the plots be printed, defaults to \code{TRUE}. If \code{FALSE} will  produce a list of ggplot objects.
#' @param ... Argument not used
#' @return A ggplot2 object.
#' 
#' @method plot bruSDM_predict
#' @rdname plot
#' 
#' @exportS3Method
#' 
#' @examples 
#' \dontrun{
#'  
#'  if (requireNamespace('INLA')) {
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  proj <- CRS("+proj=longlat +ellps=WGS84")
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- runModel(organizedData, options = list(control.inla = list(int.strategy = 'eb')))
#'    
#'   #Predict spatial field on linear scale
#'   predictions <- predict(modelRun, mesh = mesh, spatial = TRUE, fun = 'linear')
#'    
#'   #Make generic plot of predictions
#'   plot(predictions, colourHigh = 'red', colourLow = 'orange')
#'  
#'  }
#'}
#' 

plot.bruSDM_predict <- function(x,
                                whattoplot = c('mean'),
                                cols = NULL,
                                layout = NULL,
                                colourLow = NULL,
                                colourHigh = NULL,
                                plot = TRUE,
                                ...) {
  ## Add species plot::
  #if (!plotall & is.null(datasettoplot)) stop('Please provide a list of datasets to plot or set plotall to TRUE.')
  
  #if (any(!datasettoplot%in%names(x))) stop('Dataset name to plot not provided in prediction object.')
  
  #if (plotall) datasettoplot <- names(x)
  
  if (any(!whattoplot%in%c("mean", "sd", "q0.025", "median","q0.975",
                           "smin", "smax", "cv", "var" ))) stop('Whattoplot is not a valid variable to plot')
  
  if (!is.null(colourHigh) && is.null(colourLow) || !is.null(colourLow) & is.null(colourHigh)) stop ('Both colourLow and colourHigh are requried to be non-NULL.')
  
  if (!is.null(colourLow) && !is.null(colourHigh)) colours <- ggplot2::scale_fill_gradient(low = colourLow, high = colourHigh)
  else colours <- NULL
  
  if (length(x) == 1 && names(x) %in% 'temporalPredictions') {
    
    if (length(whattoplot) > 1) stop('Please only plot one variable at a time for species plots.')
    
    
    temporalName <- names(x[[1]]@data)[!names(x[[1]]@data) %in% c('weight',  'mean', 'sd', 'q0.025', 'median', 'q0.975', 'smin', 'smax', 'cv','var')]
    class(x[[1]]@data[,temporalName]) <- 'character'
    names(x[[1]]@data)[names(x[[1]]@data) == temporalName] <- '..temporal_variable_index..'

    ##Would be nice to get full temporal variable names in here ...
    plot_grid <- ggplot() + inlabru::gg(x[[1]], aes_string(fill = whattoplot)) + facet_grid(~ ..temporal_variable_index..) + ggtitle('Plot of the temporal predictions')
    return(plot_grid)
    
  }
  else 
    if (length(x) == 1 && names(x) %in% c('speciesPredictions', 'biasFields')) {
    
    nameObj <- names(x)
    
    if (length(whattoplot) > 1) stop('Please only plot one variable at a time for species plots.')
    
    
    all_plots <- list()
    
    for (object in names(x[[nameObj]])) {
      
      if (nameObj ==  'speciesPredictions') title <- ggtitle(paste('Plot of predictions for', object))
      else title <- ggtitle(paste('Plot of bias field for', object))

      all_plots[[object]] <- ggplot() + inlabru::gg(x[[nameObj]][[object]], aes_string(fill = whattoplot)) + title + colours
      
    }

    if (plot) {
      
      return(inlabru::multiplot(plotlist = all_plots, cols = length(all_plots), layout = layout))
      
      
      }
    else {
      
      return(all_plots)
      
    }

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
      
      prediction <- inlabru::gg(x[[plotname]], aes_string(fill = stat))
      
      if (!plot) prediction_list[[stat]] <- ggplot() + prediction

      
      plot_list[[stat]] <- ggplot() + prediction + title + colours
      
    }
    
    if (plot) {
      if (is.null(cols)) cols <- length(whattoplot)
      plot_grid <- inlabru::multiplot(plotlist = plot_list, cols = cols, layout = layout)
      
    }
    else{
      
      all_plots[[plotname]] <- prediction_list
      return(all_plots)
      
    }
    
  }
  

}
