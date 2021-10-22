setClass('predict_bru_sdm')

#' Predict for bru_sdm
#' 
#' @export predict.bru_sdm

predict.bru_sdm <- function(object, data = NULL, formula = NULL, mesh = NULL, 
                            mask = NULL, datasetstopredict = NULL, species = FALSE,
                            covariates = NULL, spatial = TRUE,
                            intercept = FALSE, fun = 'exp',
                            n.samples = 100, ...) {
  
  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an inla.mesh object is required.")
  
  if (!all(covariates%in%row.names(object$summary.fixed))) stop("Covariates provided not in model.")
  
  if (is.null(formula) & !spatial & !intercept & is.null(covariates)) stop("Please provide at least one of spatial, intercept or covariates.")
  
  if (species & is.null(object[['species_in']])) stop('Predictions for species was selected but no species were included in model. Please run bru_sdm with specieseffects = TRUE.')
  
  if (!species) {
  
  if (is.null(formula) & is.null(datasetstopredict)) stop("Please provide either a formula, and species or a dataset included in the bru_sdm model to be predicted.")
  
  }
    
  if (is.null(data)) {
    
  if (!is.null(mask)) {
      
  data <- pixels(mesh, mask = mask)
      
  } 
    
  else data <- pixels(mesh)
    
  }
  
  if (is.null(formula)) {
    
  int <- list()
    
  class(object) <- c('bru','inla','iinla')
  
  if (species) {
    
  numeric_species <- as.numeric(unlist(object[['species_in']]))
  
  species_variable <- attributes(object)$Species
  
  species_data <- data.frame(seq_len(max(numeric_species)))
  names(species_data) <- species_variable
  
  pixels <- inlabru::cprod(data, species_data)
  
  if (is.null(fun) | fun == 'linear') {fun <- ''}
  
  if (!is.null(covariates)) {
  species_formula <- paste0(fun,'(',paste(paste0(species_variable,'_spde'), covariates, sep = ' + '),')')
  
  }
  else  species_formula <- paste0(fun,'(',paste(paste0(species_variable,'_spde')),')')
  
  int <- predict(object, pixels, ~ data.frame(species_variable = eval(parse(text = species_variable))  ,formula = eval(parse(text = species_formula))))
  
  int[[species_variable]] <- as.character(rep(unique(unlist(object[['species_in']])), length(data)))
  
  int <- list(int)
  names(int) <- 'Species predictions'
  class(int) <- c('predict_bru_sdm', class(int))
  
  return(int)
    
  }
  else {  
  
  for(i in 1:length(datasetstopredict)) {
      
  if (spatial) {
      
  if(any(grepl('shared_spatial',object[['components']], fixed = TRUE))) {
        
  if (datasetstopredict[[i]]%in%object[['spatial_datasets']]) {
        
  spatial_obj <- 'shared_spatial'  
        
  }
  else spatial_obj <- NULL
        
  }
  else       
      
  if (!datasetstopredict[[i]]%in%object[['spatial_datasets']]) {
        
  stop('Either dataset name is incorrect or bru_sdm model run without spatial effects.')
    
  }
  else spatial_obj <- paste0(datasetstopredict[[i]],'_spde')
        
  } 
  else spatial_obj <- NULL
      
  if (intercept) {
        
  if (!paste0(datasetstopredict[[i]],'_intercept')%in%row.names(object$summary.fixed)) stop('Either dataset name is incorrect or bru_sdm model run without intercepts.')
  else intercept_obj <- paste0(datasetstopredict[[i]],'_intercept')
        
  } 
  else intercept_obj <- NULL
      
  formula_components <- c(covariates, spatial_obj, intercept_obj)
  if (all(is.null(formula_components))) stop('Please specify at least one of: covariates, spatial or intercept.')
  if (is.null(fun) | fun == 'linear') {fun <- ''}
      
  formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
      
  int[[i]] <- predict(object, data = data, formula = formula, n.samples = n.samples, ...)
      
  }
    
  names(int) <- datasetstopredict
  class(int) <- c('predict_bru_sdm', class(int))
  return(int)
    
  }
  
  }
  
  else {
    
  class(object) <- c('bru','inla','iinla')
  int <- predict(object, data = data, formula = formula, n.samples = n.samples, ...)

  return(int)
    
  }
  
  }

print.predict_bru_sdm <- function(x, ...) {
  
  
  for(name in names(x)) {
    
    cat('Summary of predicted data for', paste0(name,':'))
    cat('\n\n')
    print(summary(x[[name]]@data))
    cat('\n\n')
    
    
  }
  
  
  
}

#' Plot for preduct_bru_sdm
#' 
#' @export plot.predict_bru_sdm

plot.predict_bru_sdm <- function(x, plotall = TRUE,
                                 datasettoplot = NULL,
                                 whattoplot = c('mean'),
                                 colours = NULL,
                                 cols = length(whattoplot),
                                 layout = NULL,
                                 plot = TRUE,
                                 ...) {
  ## Add species plot::
  if (!plotall & is.null(datasettoplot)) stop('Please provide a list of datasets to plot or set plotall to TRUE.')
  
  if (any(!datasettoplot%in%names(x))) stop('Dataset name to plot not provided in prediction object.')
  
  if (plotall) datasettoplot <- names(x)
  
  if (any(!whattoplot%in%c("mean", "sd", "q0.025", "median","q0.975",
                           "smin", "smax", "cv", "var" ))) stop('Whattoplot is not a valid variable to plot')
  
  if (!is.null(colours)) {
    
    plot_colours <- scale_fill_gradientn(colours = rev(brewer.pal(9,colours)),
                                         limits = range(x[[plotname]]@data[,stat]))
    
  }
  else plot_colours <- NULL
  
  
  if (length(x) == 1 & names(x) == 'Species predictions') {
    
  if (length(whattoplot) > 1) stop('Please only plot one variable at a time for species plots.')
    
    species_name_var <- names(x[[1]]@data)[sapply(x[[1]]@data, class) == 'character']
    
    species_title <- ggtitle('Plot of the species predictions')
    
    plot_grid <- ggplot() + gg(x[[1]], aes_string(fill = whattoplot)) + facet_grid(~ eval(parse(text = species_name_var))) + plot_colours + species_title
    return(plot_grid)
  }
  
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

