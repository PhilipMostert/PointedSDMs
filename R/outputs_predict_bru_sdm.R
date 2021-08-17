#' Outputs predict_bru_sdm
#' 

predict.bru_sdm <- function(object, data = NULL, formula = NULL, mesh = NULL, 
                            mask = NULL, datasetstopredict = NULL, 
                            covariates = NULL, spatial = TRUE,
                            intercept = FALSE, fun = 'exp',
                            n.samples = 100, ...) {
  
  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an inla.mesh object is required.")
  
  if (!all(covariates%in%row.names(object$summary.fixed))) stop("Covariates provided not in model.")
  
  if (is.null(formula) & !spatial & !intercept & is.null(covariates)) stop("Please provide at least one of spatial, intercept or covariates.")
  
  if (is.null(formula) & is.null(datasetstopredict)) stop("Please provide either a formula or a dataset included in the bru_sdm model to be predicted.")
  
  #if (predictmark & is.null(markstopredict)) stop("Marks prediction is chosen but no marks to predict given.")
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- pixels(mesh, mask = mask)
      
    } 
    
    else data <- pixels(mesh)
    
  }
  
  if (is.null(formula)) {
    
    int <- list()
    
    class(object) <- c('bru','inla','iinla')
    
    for(i in 1:length(datasetstopredict)) {
      
      if (spatial) {
        
        if (!paste0(datasetstopredict[[i]],'_spde')%in%names(object$summary.random)) stop('Either dataset name is incorrect or bru_sdm model run without spatial effects.')
        else spatial_obj <- paste0(datasetstopredict[[i]],'_spde')
        
      } 
      else spatial_obj <- NULL
      
      if (intercept) {
        
        if (!paste0(datasetstopredict[[i]],'_intercept')%in%row.names(object$summary.fixed)) stop('Either dataset name is incorrect or bru_sdm model run without intercepts.')
        else intercept_obj <- paste0(datasetstopredict[[i]],'_intercept')
        
      } 
      else intercept_obj <- NULL
      
      formula_components <- c(covariates, spatial_obj, intercept_obj)
      
      if (is.null(fun) | fun == 'linear') {fun <- ''}
      
      formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
      
      
      int[[i]] <- predict(object, data = data, formula = formula, n.samples = n.samples, ...)
      
    }
    
    names(int) <- datasetstopredict
    class(int) <- c('predict_bru_sdm', class(int))
    return(int)
    
  }
  
  else {
    
    class(object) <- c('bru','inla','iinla')
    int <- predict(object, data = data, formula = formula, n.samples = n.samples, ...)
    int <- list(int)
    class(int) <- c('bru_sdm_predict',class(int))
    return(int)
    
  }
  
}

plot.predict_bru_sdm <- function(x, plotall = TRUE,
                                 datasettoplot = NULL,
                                 whattoplot = c('mean','var'),
                                 colours = NULL,
                                 cols = length(whattoplot),
                                 layout = NULL,
                                 plot = TRUE,
                                 ...) {
  
  if (!plotall & is.null(datasettoplot)) stop('Please provide a list of datasets to plot or set plotall to TRUE.')
  
  if (any(!datasettoplot%in%names(x))) stop('Dataset name to plot not provided in prediction object.')
  
  if (plotall) datasettoplot <- names(x)
  
  if (any(!whattoplot%in%c("mean", "sd", "q0.025", "median","q0.975",
                           "smin", "smax", "cv", "var" ))) stop('Whattoplot is not a valid ...')
  
  
  if (!plot) {
    
    all_plots <- list()
    prediction_list <- list()
  }
  
  for (plotname in datasettoplot) {
    
    plot_list <- list()
    
    for (stat in whattoplot) {
      
      title <- ggtitle(paste('Plot of',stat,'for',plotname))
      
      prediction <- gg(x[[plotname]], aes_string(fill = stat))
      
      if (!plot) prediction_list[[stat]] <- prediction
      
      if (!is.null(colours)) {
        
        plot_colours <- scale_fill_gradientn(colours = rev(brewer.pal(9,colours)),
                                             limits = range(x[[plotname]]@data[,stat]))
        
      }
      else plot_colours <- NULL
      
      plot_list[[stat]] <- ggplot() + prediction + title  + plot_colours
      
    }
    
    if (plot) {
      
      plot_grid <- inlabru::multiplot(plotlist = plot_list, cols = cols, layout = layout)
      
    }
    else{
      
      all_plots[[plotname]] <- prediction_list
      
    }
    
  }
  
  if (!plot) return(all_plots)
  
}

