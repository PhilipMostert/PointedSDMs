
#' Export class predict_modMarks
#' 
#' @export 

setClass('modMarks_predict')

#' Predict for modMarks
#' @title Generic predict function for \code{modMarks} objects.
#' @description Predict function for the object produced by \code{\link{fitISDM}}. Should act identically to \pkg{inlabru}'s generic predict function if wanted, but has additional arguments to help predict certain components created by the model. This is needed since \code{\link{startMarks}} creates variable names which might not be directly apparent to the user.
#' @param object A \code{modMarks} object.
#' @param data Data containing points of the map with which to predict on. May be \code{NULL} if one of \code{mesh} or \code{mask} is \code{NULL}.
#' @param formula Formula to predict. May be \code{NULL} if other arguments: \code{covariates}, \code{spatial}, \code{intercepts} are not \code{NULL}.
#' @param mesh An \code{fm_mesh_2d} object.
#' @param mask A mask of the study background. Defaults to \code{NULL}.
#' @param covariates Name of covariates to predict.
#' @param spatial Logical: include spatial effects in prediction. Defaults to \code{FALSE}.
#' @param intercepts Logical: include intercept terms in prediction. Defaults to \code{FALSE}.
#' @param datasets Names of the datasets to include intercept and spatial term.
#' @param marks Names of the marks to include intercept and spatial term.
#' @param bias Logical include bias field in prediction. Defaults to \code{FALSE}.
#' @param biasnames Names of the datasets to include bias term. Defaults to \code{NULL}. Note: the chosen dataset needs to be run with a bias field first; this can be done using \code{.$addBias} with the object produced by \code{\link{startMarks}}.
#' @param predictor Should all terms (except the bias terms) included in the linear predictor be used in the predictions. Defaults to \code{FALSE}.
#' @param fun Function used to predict. Set to \code{'linear'} if effects on the linear scale are desired.
#' @param ... Additional arguments used by the inlabru \code{predict} function.
#' 
#' @method predict modMarks
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
#'  data <- lapply(data, function(x) {x$mark = runif(nrow(x));x})
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- startMarks(data, Mesh = mesh, markNames = 'mark',
#'                              markFamily = 'gaussian',
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- fitISDM(organizedData, options = list(control.inla = list(int.strategy = 'eb',
#'                                                     diagonal = 1)))
#'    
#'   #Predict spatial field on linear scale
#'   predictions <- predict(modelRun, mesh = mesh, marks = 'mark', fun = 'linear')
#'    
#'  }
#'}
#' 

predict.modMarks <- function(object, data = NULL, formula = NULL, mesh = NULL, 
                             mask = NULL, covariates = NULL, spatial = FALSE,
                             intercepts = FALSE, datasets = NULL, 
                             marks = NULL, bias = FALSE, biasnames = NULL, predictor = FALSE,
                             fun = 'linear', ...) {
  
  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an fm_mesh_2d object is required.")
  
  ## if non-null biasfields ## if no bias fields in stop: if biasnames not in biasfields stop
  if (bias && spatial) stop('Please choose one of bias and spatial.')
  
  if (!is.null(marks)) {
    
    if (!all(marks %in% unlist(object[['marks']][['marksIn']]))) stop('Marks provided not in model.')
    
  }
  
  if (is.null(datasets)) datasets <- unique(object$source)
  
  if (predictor) {
    
    intercepts <- TRUE
    covariates <- object$spatCovs$name
    
    if (!is.null(object$covariateFormula)) covariates <- c(covariates, 'Fixed__Effects__Comps')
    
    if (!is.null(object$spatCovs$biasFormula)) covariates <- covariates[!covariates %in% labels(terms(object$spatCovs$biasFormula))]
    
    if (is.character(object$spatial$points)) spatial <- TRUE
    else {
      
      if(!object$spatial$points) spatial <- FALSE
      else spatial <- TRUE
      
    }
    
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
  
  if (!is.null(object$temporal$temporalVar)) temporal <- TRUE
  else temporal <- FALSE
  
  if (is.null(object$spatCovs$covariateFormula)) {
    
    if (!all(covariates%in%row.names(object$summary.fixed)) && !all(as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))%in%row.names(object$summary.fixed))) stop("Covariates provided not in model.")
    
  }
  
  if (is.null(formula) && is.null(marks) && !intercepts && !spatial && is.null(covariates) && !bias && !predictor) stop("Please provide either a formula or components of a formula to be predicted.")
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- fmesher::fm_pixels(mesh, mask = mask)
      
    }   
    else data <- fmesher::fm_int(mesh)
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
      ##Is this needed?
      numeric_time <- order(as.numeric(unique(unlist(object$temporal$temporalIn))))
      time_variable <- object$temporal$temporalVar
      
      time_data <- data.frame(seq_len(max(numeric_time)))
      names(time_data) <- time_variable
      
      data <- fmesher::fm_cprod(data, data.frame(time_data))
      data$.__plot__index__ <- data[[time_variable]]
      
    }
    
    if (bias) {
      
      if (!is.null(object$biasData$Fields)) {
        
        if (is.null(biasnames)) biasnames <- object$biasData$Fields
        biasnames <- paste0(biasnames,'_biasField')
        
      } else biasnames <- NULL
      
      #paste0 specieshere
      if ('Bias__Effects__Comps' %in% names(object$summary.random)) biasnames <- c('Bias__Effects__Comps', biasnames)
      else biasnames <- c(biasnames, NULL)
      
      if (!all(biasnames %in% names(object$summary.random))) stop('Either no bias field has been used or an incorrect dataset name was given.')
      
      if (temporal) { 
        
        formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(biasnames,')')))
        int[['temporalBiasFields']] <- predict(object, newdata = data, formula = formula, ...)
        
      }
      else {
        
        for (bias in biasnames) {
          
          formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(bias,')')))
          int[['biasFields']][[bias]] <- predict(object, newdata = data, formula = formula, ...)
          
        }
      }
      
      class(int) <- c('modMarks_predict', class(int))
      return(int) 
      
    }
    
      if (!is.null(marks)) {
        
        marks_spatial <- paste0(marks,'_spatial')
        marks_copy_spat <- names(object$summary.random)[names(object$summary.random) %in% apply(expand.grid(datasets, paste0(marks, '_spatial')), 1, paste, collapse="_")]
      }
      else {
        
        marks_spatial <- NULL
        marks_copy_spat <- NULL
      }
    
    if (spatial) {
        
        if ('shared_spatial' %in% names(object$summary.random))  spatial_obj <- 'shared_spatial'
        else
          if (object$spatial$points == 'copy') spatial_obj <- paste0(object$source[1], '_spatial')
          else
            if (!all(paste0(datasets,'_spatial') %in% names(object$summary.random))) stop('Spatial effects not provided in startMarks')
          else spatial_obj <- paste0(datasets, '_spatial')
          
          if (object$spatial$points == 'correlate') {
            
            if (any(object$dataType == "Present absence")) data$._dataset_index_var_. <- which(object$dataType == "Present absence")[1]
            else data$._dataset_index_var_. <- 1
            
          }
          
    }
    else spatial_obj <- NULL

    
    if (predictor) formula_components <- c(row.names(object$summary.fixed), names(object$summary.random)[!names(object$summary.random) %in% paste0(object[['source']], '_biasField')])
    else formula_components <- c(covariates, intercept_terms, spatial_obj, marks_spatial, marks_intercepts)
    
    if (!is.null(object$spatCovs$biasFormula)) formula_components <- formula_components[!formula_components %in% c('Bias__Effects__Comps', paste0(unique(object$species$speciesIn),'_Bias__Effects__Comp'))]
    
    if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial, intercepts or bias.')
    
    formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
    
    #int[[i]] <- predict(object, data = data, formula = formula, ...)
    int <- predict(object, newdata = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    
    class(int) <- c('modMarks_predict', class(int))
    return(int)
    
  }
  else {
    
    class(object) <- c('bru','inla','iinla')
    int <- predict(object, newdata = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    class(int) <- c('modMarks_predict', class(int))
    
    return(int)
    
  }
  
  
  
}

#' Plot for modMarks_predict
#' @title Generic plot function for \code{modMarks_predict}.
#' @param x A bruSDM_predict object.
#' @param variable One of the following statistics to plot: "mean", "sd", "q0.025", "median","q0.975", "smin", "smax", "cv", "var" 
#' @param plot Should the plots be printed, defaults to \code{TRUE}. If \code{FALSE} will  produce a list of ggplot objects.
#' @param ... Argument not used
#' @return A ggplot2 object.
#' 
#' @method plot modMarks_predict
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
#'  organizedData <- startMarks(data, Mesh = mesh,
#'                              Projection = proj, responsePA = 'Present',
#'                              markNames = 'speciesName', 
#'                              markFamily = 'multinomial')
#'  
#'   ##Run the model
#'  modelRun <- fitISDM(organizedData, options = list(control.inla = list(int.strategy = 'eb',
#'                                                                      diagonal = 1)))
#'    
#'   #Predict spatial field on linear scale
#'   predictions <- predict(modelRun, mesh = mesh, spatial = TRUE, fun = 'linear')
#'    
#'   #Make generic plot of predictions
#'   plot(predictions)
#'  
#'  }
#'}
#' 

plot.modMarks_predict <- function(x,
                                  variable = 'mean',
                                  plot = TRUE,
                                  ...) {
  
  if (any(!variable%in%c("mean", "sd", "q0.025", "median","q0.975",
                         "smin", "smax", "cv", "var" ))) stop('variable is not a valid variable to plot')
  
  if (length(x) == 1 && '.__plot__index__' %in% names(x[[1]])) {
    
    if (length(variable) > 1) stop('Please only plot one variable at a time for species plots.')
    
    ##Need to create a new var called ..temporal_variable_index.. which is the tempVar
    #temporalName <- names(x[[1]])[!names(x[[1]]) %in% c(".block", 'geometry', 'weight', 'mean', 'sd', 'q0.025', 'median', 'q0.975', 'q0.5', 'smin', 'smax', 'cv','mean.mc_std_err', 'sd.mc_std_err')]
    temporalName <- '.__plot__index__'
    
    x[[1]]$..temporal_variable_index.. <- as.character(data.frame(x[[1]])[, temporalName])
    
    if (inherits(x[[1]], 'sf')) plot_obj <- geom_sf(data = x[[1]], aes(col = .data[[variable]]))
    else plot_obj <- inlabru::gg(x[[1]], aes(col = .data[[variable]]))
    
    plot_grid <- ggplot() + plot_obj + facet_wrap(~ ..temporal_variable_index..) + ggtitle('Plot of the temporal predictions')
    return(plot_grid)
    
  }
  
  
  if (!plot) {
    
    all_plots <- list()
    prediction_list <- list()
    
  }
  
  if (names(x) == 'biasFields') {
    
    biasPlot <- TRUE
    namesBias <- names(x[[1]])
    x <- unlist(x, recursive = FALSE, use.names = FALSE)
    names(x) <- namesBias
  } else biasPlot <- FALSE
  
  ##If plots for the individual fields??
  
  datasettoplot <- names(x)
  
  plot_list <- list()
  
  for (plotname in datasettoplot) {
    
    for (stat in variable) {
      
      if (biasPlot) title <- ggtitle(paste('Predictions of', stat, 'for', plotname, 'bias field'))
      else title <- ggtitle(paste('Predictions of', stat))
      
      prediction <- inlabru::gg(x[[plotname]], aes(col = .data[[stat]]))
      
      if (!plot) prediction_list[[plotname]][[stat]] <- ggplot() + prediction + title
      
      plot_list[[plotname]][[stat]] <- ggplot() + prediction + title
      
    }
    
  }
  
  if (plot) {
    
    rows <- length(variable)
    cols <- length(datasettoplot)
    
    ind <- max(rows, cols)      
    plot_list <- unlist(plot_list, recursive = FALSE)
    
    plot_grid <- inlabru::multiplot(plotlist = plot_list, cols = ind)
    
  }
  else{
    
    all_plots[[plotname]] <- prediction_list
    return(all_plots)
    
  }
  
  
}
