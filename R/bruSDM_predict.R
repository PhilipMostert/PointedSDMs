
#' Export class predict_bru_sdm
#' 
#' @export 

setClass('bruSDM_predict')

#' Predict for bru_sdm
#' @title Generic predict function for \code{bru_SDM} objects.
#' @description Predict function for the object produced by \code{\link{fitISDM}}. Should act identically to \pkg{inlabru}'s generic predict function if wanted, but has additional arguments to help predict certain components created by the model. This is needed since \code{\link{intModel}} creates variable names which might not be directly apparent to the user.
#' @param object A \code{bru_sdm} objects.
#' @param data Data containing points of the map with which to predict on. May be \code{NULL} if one of \code{mesh} or \code{mask} is \code{NULL}.
#' @param formula Formula to predict. May be \code{NULL} if other arguments: \code{covariates}, \code{spatial}, \code{intercepts} are not \code{NULL}.
#' @param mesh An \code{fm_mesh_2d} object.
#' @param mask A mask of the study background. Defaults to \code{NULL}.
#' @param covariates Name of covariates to predict.
#' @param temporal Make predictions for the temporal component of the model.
#' @param spatial Logical: include spatial effects in prediction. Defaults to \code{FALSE}.
#' @param intercepts Logical: include intercept terms in prediction. Defaults to \code{FALSE}.
#' @param datasets Names of the datasets to include intercept and spatial term.
#' @param marks Names of the marks to include intercept and spatial term.
#' @param species Names of the species to predict. Default of \code{NULL} results in all species being predicted.
#' @param biasfield Logical include bias field in prediction. Defaults to \code{FALSE}.
#' @param biasnames Names of the datasets to include bias term. Defaults to \code{NULL}. Note: the chosen dataset needs to be run with a bias field first; this can be done using \code{.$addBias} with the object produced by \code{\link{intModel}}.
#' @param predictor Should all terms (except the bias terms) included in the linear predictor be used in the predictions. Defaults to \code{FALSE}.
#' @param fun Function used to predict. Set to \code{'linear'} if effects on the linear scale are desired.
#' @param format Class of the data for which to predict on. Must be one of \code{'sp'}, \code{'sf'} or \code{'terra'}. Defaults to \code{'sf'}.
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
#'  proj <- "+proj=longlat +ellps=WGS84"
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- fitISDM(organizedData, options = list(control.inla = list(int.strategy = 'eb')))
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
                           marks = NULL, biasfield = FALSE, biasnames = NULL, predictor = FALSE,
                           fun = 'linear', format = 'sf', ...) {
  
  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an fm_mesh_2d object is required.")
  
  ## if datasets !is.null but at least one not in model stop
  # else datasets <- all datasets in the model.
  ##Need another defense: only one of predictor; biasfield; temporal
  if (sum(predictor, temporal) > 1 || sum(predictor, biasfield) > 1) stop('You cannot combine predictor with either "temporal" or "biasfield".')
  ## if non-null biasfields ## if no bias fields in stop: if biasnames not in biasfields stop
  if (biasfield && spatial) stop('Please choose one of biasfield and spatial.')
  
  if (!is.null(marks)) {

    if (!all(marks %in% unlist(object[['marks']][['marksIn']]))) stop('Marks provided not in model.')
    
  }
  
  if (is.null(datasets)) datasets <- unique(object$source)

  
  if (!is.null(unlist(object[['species']][['speciesIn']]))) speciespreds <- TRUE
  else speciespreds <- FALSE
    
    if (predictor) {
      
      intercepts <- TRUE
      covariates <- object$spatCovs$name
      spatial <- TRUE
      
    }
    
    if (intercepts) {
      
      if (!is.null(marks)) {
        
        marks_intercepts <- paste0(marks,'_intercept')
        
        if (identical(rownames(object$summary.fixed)[rownames(object$summary.fixed) %in% mark_intercepts], character(0))) mark_intercepts <- NULL
        
        
      }
      
      if (!is.null(datasets)) {
        
        intercept_terms <- paste0(datasets, '_intercept')
        
        if (identical(rownames(object$summary.fixed)[rownames(object$summary.fixed) %in% intercept_terms], character(0))) intercept_terms <- NULL
        
      }
      
      
    } 
    else {
      
      intercept_terms <- NULL
      marks_intercepts <- NULL
      
    }
    
    if (is.null(species)) speciesin <- unique(unlist(object[['species']][['speciesIn']]))
    else speciesin <- species
    if (!all(species %in% unique(unlist(object[['species']][['speciesIn']])))) stop('Species provided not in model.')
    
    
  
  if (is.null(object$spatCovs$covariateFormula)) {
    
  if (!all(covariates%in%row.names(object$summary.fixed)) && !all(as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))%in%row.names(object$summary.fixed))) stop("Covariates provided not in model.")

  }
  
  if (is.null(formula) && !intercepts && !spatial && is.null(covariates) && !temporal && !biasfield && !predictor) stop("Please provide either a formula or components of a formula to be predicted.")
  
  if (temporal && is.null(object$temporal$temporalVar)) stop('temporal is set to TRUE but no temporal component found in the model.')
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- fmesher::fm_pixels(mesh, mask = mask, format = format)
      
    }   
    else data <- fmesher::fm_int(mesh, format = format)
  }
  
  if (speciespreds) {
  
  if (object[['species']][['speciesEffects']][['Intercepts']]) {
    
    data <- fmesher::fm_cprod(data, data.frame(speciesIndexREMOVE = 1:length(unique(unlist(object$species$speciesIn)))))
    names(data)[names(data) == 'speciesIndexREMOVE'] <- object[['species']][['speciesVar']]
    
  }
    
  if (object$spatial$species == 'replicate') {
      
  if (!object[['species']][['speciesVar']] %in% names(data)) data <- fmesher::fm_cprod(data, data.frame(speciesSpatialGroup = 1:length(unique(unlist(object$species$speciesIn)))))
  else data$speciesSpatialGroup <- data[[object[['species']][['speciesVar']]]]
    
  }
    
  }
  
  if (!any(names(data) %in% object$spatCovs$name)) {
    
    for (spatCov in object$spatCovs$name) {
      
      if (!is.null(object$species$speciesIn) && object$species$speciesEffects$Environmental) covIndex <- paste0(unique(unlist(object$species$speciesIn)), '_', spatCov)
      else covIndex <- spatCov
      
      data[, covIndex] <- inlabru::eval_spatial(where =  data, 
                                                data = get('spatialcovariates', 
                                                envir = object$spatCov$env)[spatCov],
                                                layer = spatCov)
      
      if (any(is.na( data[, covIndex]))) {
        
        for (indFix in covIndex) {
        
        data[[indFix]] <- inlabru::bru_fill_missing(where =  data, 
                                                      data = get('spatialcovariates', 
                                                             envir = object$spatCov$env)[spatCov],
                                                      layer = spatCov,
                                                      values = data[[indFix]])
        
        }
        
      }
        
      
    }

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
      
      timeData <- fmesher::fm_cprod(data, data.frame(time_data))
      names(timeData)[!names(timeData) %in% c('geometry', '.block')] <- c(time_variable, 'weight')  
      
      #bias
      
      if (biasfield) {
        
        if (is.null(biasnames)) biasnames <- object$biasData
        
        biasnames <- paste0(biasnames,'_biasField')
        
        if (!all(biasnames %in% names(object$summary.random))) stop('Either no bias field has been used or an incorrect dataset name was given.')
        
        
      } else biasnames <- NULL
      
      if (!is.null(covariates) && speciespreds) covariates <- as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))
      if (intercepts) {
        
      if (!is.null(object$species$speciesIn)) intercept_terms <- paste0(unlist(object[['species']][['speciesIn']]), '_intercept')
      else intercept_terms <- paste0(object$source, '_intercept')
        
      }
      
      if ('shared_spatial' %in% names(object$summary.random))  spatial_obj <- 'shared_spatial'
      else
      if (!all(paste0(datasets,'_spatial') %in% names(object$summary.random))) stop('Spatial effects not provided in intModel.')
      else spatial_obj <- paste0(datasets, '_spatial')
      
      time_formula <- paste0(fun,'(',paste0(c(covariates, biasnames, intercept_terms, spatial_obj, intercepts_species), collapse = ' + '),')')

      formula <- formula(paste('~',paste0('data.frame(', time_variable,' = ', time_variable, ',formula =', time_formula,')')))

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
        int[['biasFields']][[bias]] <- predict(object, newdata = data, formula = formula, ...)
        
      }
      
      class(int) <- c('bruSDM_predict', class(int))
      return(int) 
      
    }
    
    if (speciespreds) {
      
      int[['speciesPredictions']] <- vector(mode = 'list', length(speciesin))
      names(int[['speciesPredictions']]) <- speciesin
      
      for (spec in speciesin) {
        
        if (!is.null(covariates)) {
          
          if (object[['species']][['speciesEffects']][['Environmental']]) species_covs <- paste0(spec, '_', covariates)
          else species_covs <- covariates
          
        }
        else species_covs <- NULL
        
        if (intercepts) {
          
          if (!object[['species']][['speciesEffects']][['Intercepts']]) {
            
            intercepts_species <- paste0(spec,'_intercept')
            
          } else if (object[['species']][['speciesEffects']][['Intercepts']]) {
              
            intercepts_species <- paste0(object[['species']][['speciesVar']], '_intercepts')
            
            } else intercepts_species <- NULL
        
        }
        else intercepts_species <- NULL
        
        if (spatial) {
          ##Fix this
          allSpat <- c(paste0(spec, '_spatial'),
                       paste0(spec, '_', names(object$dataType), '_spatial'),
                       'speciesShared')
          
          species_spat <- allSpat[allSpat %in% names(object$summary.random)]
          
          
          }
        else species_spat <- NULL
        
        species_formula <- formula(paste0('~', fun, '(', paste0(c(species_spat, species_covs, intercepts_species, intercept_terms), collapse = ' + '),')'))
        
        if (any(c('speciesSpatialGroup', object[['species']][['speciesVar']]) %in% names(data))) {
          
          spInd <- object$species$speciesTable[object$species$speciesTable$species == spec,]$index
          which <- names(data)[names(data) %in% c('speciesSpatialGroup', object[['species']][['speciesVar']]) ][1]
          dataSP <- data[data[[which]] == spInd,]
          
          
        }
        else dataSP <- data
        
        int[[1]][[spec]] <- predict(object, dataSP, formula = species_formula, ...)
      
        }
      
      class(int) <- c('bruSDM_predict', class(int))
      return(int) 
      
      
    }
    
    if (spatial) {
      
      if (!is.null(marks)) {
        
        marks_spatial <- paste0(marks,'_spatial')
        spatial_obj <- NULL
      }
      else marks_spatial <- NULL
      
      if (is.null(marks_spatial)) {
      
      if ('shared_spatial' %in% names(object$summary.random))  spatial_obj <- 'shared_spatial'
      else
        if (object$spatial$points == 'copy') spatial_obj <- paste0(object$source[1], '_spatial')
      else
        if (!all(paste0(datasets,'_spatial') %in% names(object$summary.random))) stop('Spatial effects not provided in intModel.')
      else spatial_obj <- paste0(datasets, '_spatial')
      
      } 
      }
    else {
      
      spatial_obj <- NULL
      marks_spatial <- NULL
      
    }
    
    if (predictor) formula_components <- c(row.names(object$summary.fixed), names(object$summary.random)[!names(object$summary.random) %in% paste0(object[['source']], '_biasField')])
    else formula_components <- c(covariates, intercept_terms, spatial_obj, marks_spatial, marks_intercepts)
    
    if (!is.null(object$spatCovs$biasFormula)) formula_components <- formula_components[!formula_components %in% c('Bias__Effects__Comps', paste0(unique(object$species$speciesIn),'_Bias__Effects__Comp'))]
    
    if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial, intercepts or biasfield.')
    
    formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
    
    #int[[i]] <- predict(object, data = data, formula = formula, ...)
    int <- predict(object, newdata = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    
    class(int) <- c('bruSDM_predict', class(int))
    return(int)
    
  }
  else {
    
    class(object) <- c('bru','inla','iinla')
    int <- predict(object, newdata = data, formula = formula, ...)
    if (any(c('speciesSpatialGroup', object[['species']][['speciesVar']]) %in% names(int))) {
      int[['..speciesPlotVar..']] <- NA
      which <- names(int)[names(int) %in% c('speciesSpatialGroup', object[['species']][['speciesVar']])][1]
      
      for (spec in speciesin) {
      
      spInd <- object$species$speciesTable[object$species$speciesTable$species == spec,]$index
      try(int[int[[which]] == spInd,][[which]] <- spec, silent = TRUE)
      try(int[int[[which]] == spec,][['..speciesPlotVar..']] <- spec, silent = TRUE)
      
      }
      
    }
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
        print(summary(data.frame(x[[1]][[species]])[!names(data.frame(x[[1]][[species]])) %in% c('geometry', 'coords.x1', 'coords.x2')]))
        cat('\n')
         
        }
      
    }
    else
      if (names(x)[[1]] == 'temporalPredictions') {
        
        cat('Predictions for the temporal variable:')
        cat('\n')
        print(summary(data.frame(x[[1]])[!names(data.frame(x[[1]])) %in% c('geometry', 'coords.x1', 'coords.x2')]))
        
      }
    else
      if (names(x)[[1]] == 'biasFields') {
        
        for(bias in names(x[['biasFields']])) {
          
          cat('Predictions of the bias field for', paste0(bias,':'))
          cat('\n')
          print(summary(data.frame(x[[1]])[[bias]][!names(data.frame(x[[1]][[bias]])) %in% c('geometry', 'coords.x1', 'coords.x2')]))
          cat('\n')
          
        }
        
      }
    else print(summary(x[['predictions']]@data))
    cat('\n\n')
    
  
}

#' Plot for predict_bru_sdm
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
#'  proj <- "+proj=longlat +ellps=WGS84"
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- fitISDM(organizedData, options = list(control.inla = list(int.strategy = 'eb')))
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
    
    
    temporalName <- names(x[[1]])[!names(x[[1]]) %in% c(".block", 'geometry', 'weight', 'mean', 'sd', 'q0.025', 'median', 'q0.975', 'q0.5', 'smin', 'smax', 'cv','mean.mc_std_err', 'sd.mc_std_err')]
    #class(x[[1]][,temporalName]) <- 'character'
    #names(x[[1]]@data)[names(x[[1]]@data) == temporalName] <- '..temporal_variable_index..'
    x[[1]]$..temporal_variable_index.. <- as.character(data.frame(x[[1]])[, temporalName])

    if (inherits(x[[1]], 'sf')) plot_obj <- geom_sf(data = x[[1]], aes(col = .data[[whattoplot]]))
    else plot_obj <- inlabru::gg(x[[1]], aes(col = .data[[whattoplot]]))
    
    ##Would be nice to get full temporal variable names in here ...
    plot_grid <- ggplot() + plot_obj + facet_wrap(~ ..temporal_variable_index..) + ggtitle('Plot of the temporal predictions')
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
      
      if (inherits(x[[nameObj]][[object]], 'sf')) plot_obj <- geom_sf(data = x[[nameObj]][[object]], aes(col = .data[[whattoplot]]))
      else plot_obj <- inlabru::gg(x[[nameObj]][[object]], aes(col = .data[[whattoplot]]))

      all_plots[[object]] <- ggplot() + plot_obj + title + colours
      
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
      
      if ('..speciesPlotVar..' %in% names(x[[plotname]])) speciesTerm <- facet_wrap(~..speciesPlotVar..)
      else speciesTerm <- NULL
      
      if (inherits(x[[plotname]], 'sf')) prediction <- geom_sf(data = x[[plotname]], aes(col = .data[[stat]]))
      else prediction <- inlabru::gg(x[[plotname]], aes(col = .data[[stat]]))
      
      if (!plot) prediction_list[[stat]] <- ggplot() + prediction + speciesTerm

      
      plot_list[[stat]] <- ggplot() + prediction + title + colours + speciesTerm
      
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
