
#' Export class predict_modSpecies
#' 
#' @export 

setClass('modSpecies_predict')

#' Predict for modSpecies
#' @title Generic predict function for \code{bru_SDM} objects.
#' @description Predict function for the object produced by \code{\link{fitISDM}}. Should act identically to \pkg{inlabru}'s generic predict function if wanted, but has additional arguments to help predict certain components created by the model. This is needed since \code{\link{startSpecies}} creates variable names which might not be directly apparent to the user.
#' @param object A \code{modSpecies} object.
#' @param data Data containing points of the map with which to predict on. May be \code{NULL} if one of \code{mesh} or \code{mask} is \code{NULL}.
#' @param formula Formula to predict. May be \code{NULL} if other arguments: \code{covariates}, \code{spatial}, \code{intercepts} are not \code{NULL}.
#' @param mesh An \code{inla.mesh} object.
#' @param mask A mask of the study background. Defaults to \code{NULL}.
#' @param covariates Name of covariates to predict.
#' @param spatial Logical: include spatial effects in prediction. Defaults to \code{FALSE}.
#' @param intercepts Logical: include intercept terms in prediction. Defaults to \code{FALSE}.
#' @param datasets Names of the datasets to include intercept and spatial term.
#' @param species Names of the species to predict. Default of \code{NULL} results in all species being predicted.
#' @param bias Logical include bias field in prediction. Defaults to \code{FALSE}.
#' @param biasnames Names of the datasets to include bias term. Defaults to \code{NULL}. Note: the chosen dataset needs to be run with a bias field first; this can be done using \code{.$addBias} with the object produced by \code{\link{intModel}}.
#' @param predictor Should all terms (except the bias terms) included in the linear predictor be used in the predictions. Defaults to \code{FALSE}.
#' @param fun Function used to predict. Set to \code{'linear'} if effects on the linear scale are desired.
#' @param ... Additional arguments used by the inlabru \code{predict} function.
#' 
#' @method predict modSpecies
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
#'  organizedData <- startSpecies(data, Mesh = mesh, speciesName = 'speciesName',
#'                              Projection = proj, responsePA = 'Present')
#'  
#'   ##Run the model
#'   modelRun <- fitISDM(organizedData, options = list(control.inla = list(int.strategy = 'eb',
#'                                                     diagonal = 1)))
#'    
#'   #Predict spatial field on linear scale
#'   predictions <- predict(modelRun, mesh = mesh, spatial = TRUE, fun = 'linear')
#'    
#'  }
#'}
#' 

predict.modSpecies <- function(object, data = NULL, formula = NULL, mesh = NULL, 
                           mask = NULL, covariates = NULL, spatial = FALSE,
                           intercepts = FALSE, datasets = NULL, species,
                           bias = FALSE, biasnames = NULL, predictor = FALSE,
                           fun = 'linear', ...) {
  
  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an inla.mesh object is required.")
  
  if (sum(predictor, bias) > 1) stop('You cannot combine predictor with bias.')
  
  #Why can't you do both here?
  if (bias && spatial) stop('Please choose one of bias and spatial.')
  
  if (is.null(datasets)) datasets <- unique(object$source)
  
  if (!missing(species))   {
    
    if (!all(species %in% unique(unlist(object[['species']][['speciesIn']])))) stop('Species provided not in model.')
   
    species <- sort(species)
     
  }
  else species <-  unique(unlist(object$species$speciesIn))
  
  if(is.null(species)) speciespreds <- FALSE
  else speciespreds <- TRUE
  
  
  if (predictor) {
    
    intercepts <- TRUE
    covariates <- object$spatCovs$name
    if (!is.null(object$covariateFormula)) covariates <- c(covariates, paste0(species,'_Fixed__Effects__Comps'))
    
    if (!is.null(object$spatCovs$biasFormula)) covariates <- covariates[!covariates %in% labels(terms(object$spatCovs$biasFormula))]
    
    spatial <- TRUE
    
  }
  
  if (intercepts) {
    
    if (!is.null(datasets)) {
      
      intercept_terms <- paste0(datasets, '_intercept')
      
      if (identical(rownames(object$summary.fixed)[rownames(object$summary.fixed) %in% intercept_terms], character(0))) intercept_terms <- NULL
      
    } else intercept_terms <- NULL
    
    
  } else intercept_terms <- NULL
  
  if (is.null(object$spatCovs$covariateFormula)) {
    
    covsEst <- c(row.names(object$summary.fixed), names(object$summary.random)[names(object$summary.random) %in% c(object$spatCovs$name, 
               as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), object$spatCovs$name, FUN = 'paste0')))])
    
    if (!all(covariates%in%covsEst) && !all(as.vector(outer(paste0(unlist(object[['species']][['speciesIn']]),'_'), covariates, FUN = 'paste0'))%in%covsEst)) stop("Covariates provided not in model.")
    
  }
  
  if (is.null(formula) && !intercepts && !spatial && is.null(covariates) && !bias && !predictor) stop("Please provide either a formula or components of a formula to be predicted.")
  
  if (!is.null(object$temporal$temporalVar)) temporal <- TRUE
  else temporal <- FALSE
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- inlabru::fm_pixels(mesh, mask = mask)
      
    }   
    else data <- inlabru::fm_int(mesh)
  }
  
  if (speciespreds) {
    
    if (object[['species']][['speciesEffects']][['Intercepts']]) {
      
      data <- fm_cprod(data, data.frame(speciesIndexREMOVE = 1:length(unique(unlist(object$species$speciesIn)))))
      names(data)[names(data) == 'speciesIndexREMOVE'] <- object[['species']][['speciesVar']]
      
    }
    
    if (object$spatial$species == 'replicate') {
      
      if (!object[['species']][['speciesVar']] %in% names(data)) data <- fm_cprod(data, data.frame(speciesSpatialGroup = 1:length(unique(unlist(object$species$speciesIn)))))
      else data$speciesSpatialGroup <- data[[object[['species']][['speciesVar']]]]
      
    }
    
    if (!object$species$speciesEffects$Intercepts) {
      
      data <- fm_cprod(data, data.frame(temp_species_index_var = 1:length(unique(unlist(object$species$speciesIn)))))
      names(data)[names(data) == 'temp_species_index_var'] <- object[['species']][['speciesVar']]
      
    }
    
  }
  
  if (!any(names(data) %in% object$spatCovs$name)) {
    
    for (spatCov in object$spatCovs$name) {
      
      if (!is.null(object$spatCovs$biasFormula)) {
        
        if (spatCov %in% labels(terms(object$spatCovs$biasFormula))) covIndex <- spatCov
        else 
          if (object$species$speciesEffects$Environmental) covIndex <- paste0(unique(unlist(object$species$speciesIn)), '_', spatCov)
          else covIndex <- spatCov
      }
      else
        if (object$species$speciesEffects$Environmental) covIndex <- paste0(unique(unlist(object$species$speciesIn)), '_', spatCov)
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
      
      data <- inlabru::fm_cprod(data, data.frame(time_data))
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
      
      class(int) <- c('modSpecies_predict', class(int))
      return(int) 
      
    }
    
    if (speciespreds) {
      
      int[['speciesPredictions']] <- vector(mode = 'list', length(species))
      names(int[['speciesPredictions']]) <- species
      
      speciesEff <- vector(mode = 'list', length = length(species))
      names(speciesEff) <- species
      
      for (spec in species) {
        
        if (!is.null(covariates)) {
          
          if (paste0(spec, '_Fixed__Effects__Comps') %in% names(object$summary.random)) species_covs <- paste(spec, '_Fixed__Effects__Comps')
          
          
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
        speciesEff[spec] <-  paste(spec, '=', fun, '(0 +',paste0(c(species_covs, species_spat, intercepts_species, intercept_terms), collapse = '+'), ')')
        
      }
      
      .__speciesFormulas.__ <- paste(do.call(paste0, list(speciesEff, sep = ';')), collapse = '')
      
      .__speciesEval.__ <- paste('Predictions = list(', paste(species,'=',species, collapse = ' , '),')')
      
      .__thin.__ <- paste0(paste(paste0(species, '[!1:length(',species,') %in% seq(', 1:length(species),',length(',species,'),', length(species), ')] <- FALSE'), collapse=';'),';')
      
      
      predictionFormula <- paste('{',
                                 .__speciesFormulas.__,
                                 .__thin.__,
                                 .__speciesEval.__ ,'}')
      
        int <- predict(object, data, formula = parse(text = predictionFormula), ...)
        
        int <- list(mapply(function(x, seq) {
          
          pred <- x[x[[object$species$speciesVar]] == seq,]
          list(pred)
          
          
        }, int, seq = 1:length(int)))
        
      names(int) <- 'speciesPredictions'  
      
      class(int) <- c('modSpecies_predict', class(int))
      return(int) 
      
      
    }
    
    if (spatial) {
      
        
        if ('shared_spatial' %in% names(object$summary.random))  spatial_obj <- 'shared_spatial'
        else
          if (object$spatial$points == 'copy') spatial_obj <- paste0(object$source[1], '_spatial')
          else
            if (!all(paste0(datasets,'_spatial') %in% names(object$summary.random))) stop('Spatial effects not provided in intModel.')
          else spatial_obj <- paste0(datasets, '_spatial')
          }
         else spatial_obj <- NULL

    
    if (predictor) formula_components <- c(row.names(object$summary.fixed), names(object$summary.random)[!names(object$summary.random) %in% paste0(object[['source']], '_biasField')])
    else formula_components <- c(covariates, intercept_terms, spatial_obj)
    
    if (!is.null(object$spatCovs$biasFormula)) formula_components <- formula_components[!formula_components %in% c('Bias__Effects__Comps', paste0(unique(object$species$speciesIn),'_Bias__Effects__Comp'))]
    
    if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial, intercepts or bias.')
    
    formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
    
    #int[[i]] <- predict(object, data = data, formula = formula, ...)
    int <- predict(object, newdata = data, formula = formula, ...)
    int <- list(int)
    names(int) <- 'predictions'
    
    class(int) <- c('modSpecies_predict', class(int))
    return(int)
    
  }
  else {
    
    ##Fix this?
    class(object) <- c('bru','inla','iinla')
    int <- predict(object, newdata = data, formula = formula, ...)
    if (any(c('speciesSpatialGroup', object[['species']][['speciesVar']]) %in% names(int))) {
      int[['..speciesPlotVar..']] <- NA
      which <- names(int)[names(int) %in% c('speciesSpatialGroup', object[['species']][['speciesVar']])][1]
      
      for (spec in species) {
        
        spInd <- object$species$speciesTable[object$species$speciesTable$species == spec,]$index
        try(int[int[[which]] == spInd,][[which]] <- spec, silent = TRUE)
        try(int[int[[which]] == spec,][['..speciesPlotVar..']] <- spec, silent = TRUE)
        
      }
      
    }
    int <- list(int)
    names(int) <- 'predictionsFormula'
    class(int) <- c('modSpecies_predict', class(int))
    
    return(int)
    
  }
  
  
  
}

#' Plot for modSpecies_predict
#' @title Generic plot function for \code{modSpecies_predict}.
#' @param x A modSpecies_predict object.
#' @param variable One of the following statistics to plot: "mean", "sd", "q0.025", "median","q0.975", "smin", "smax", "cv", "var" 
#' @param plot Should the plots be printed, defaults to \code{TRUE}. If \code{FALSE} will  produce a list of ggplot objects.
#' @param ... Argument not used
#' @return A ggplot2 object.
#' 
#' @method plot modSpecies_predict
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

plot.modSpecies_predict <- function(x,
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
  
  if (names(x)[1] == 'biasFields') {
    
    biasPlot <- TRUE
    namesBias <- names(x[[1]])
    x <- unlist(x, recursive = FALSE, use.names = FALSE)
    names(x) <- namesBias
  } else biasPlot <- FALSE
  
  if (names(x)[] == 'speciesPredictions') {
    
    speciesPlot <- TRUE
    namesSpecies <- names(x[[1]])
    x <- unlist(x, recursive = FALSE, use.names = FALSE)
    names(x) <- namesSpecies
  } else speciesPlot <- FALSE
  ##If plots for the individual fields??
  
  datasettoplot <- names(x)
  
  plot_list <- list()
  
  for (plotname in datasettoplot) {
    
    for (stat in variable) {
      
      if (biasPlot) title <- ggtitle(paste('Predictions of', stat, 'for', plotname, 'bias field'))
      else 
        if (speciesPlot) title <-  ggtitle(paste('Predictions of', stat, 'for', plotname))
      else title <- ggtitle(paste('Predictions of', stat, 'for', plotname))
      
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
