#' Bru_sdm output functions
#' 


print.bru_sdm = function(x,...) {
  print(summary(x))
}

summary.bru_sdm = function(x,...) {
  #cat('----bru_sdm summary STILL IN DEVELOPMENT----\n\n')
  
  cat("Summary of 'bru_sdm' object:\n\n")
  cat(paste0("inlabru version: ", x$bru_info$inlabru_version, "\n"))
  cat(paste0("INLA version: ", x$bru_info$INLA_version, "\n\n"))
  
  cat('Types of data modelled:\n')
  names_data = data.frame(x[['data_type']][!duplicated(names(x[['data_type']]))])
  names(names_data) = c('                              ')
  print(names_data)
  cat('\n')
  
  ##Make this a for loop to print out
   # all 'mulitnom_vars' summaries
  if(!is.null(x[['species_in_model']])){
    #cat('Species in model:\n\n')
    #cat(paste0(x[['species_in_model']],collapse = ", "))
    #cat('\n\n')
    cat('Summary of species:\n\n')
    species <- x$summary.random$species
    names(species)[1] <- 'species'
    print(species[,1:7])#, digits = 3, row.names = FALSE)
    cat('\n\n')
  }
  
  if(!is.null(x[['multinom_vars']])){
    cat('Summary of multinomial variables:\n\n')
    for (i in x[['multinom_vars']]) {
      
      variable <- x$summary.random[[i]]
      names(variable)[1] <- i
      print(variable[,1:7], row.names = FALSE, digits = 3)
      cat('\n\n')
      
    }
    }
  
  ##Fix naming of datasets in bru_sdm file
  if(length(x[['model_residuals']]) > 0){
  cat('Summary of residuals:\n\n')
  summary_residuals = sapply(x[['model_residuals']], summary, digits = 3)
  print(summary_residuals)
  cat('\n\n')
  }
  
  #Add Likelihoods with: Family, predictor and which dataset it came from
  class(x) = 'inla'
  x$call = NULL
  summary(x)
  
}

predict.bru_sdm <- function(object, data = NULL, formula = NULL, mesh = NULL, 
                            mask = NULL, datasettopredict = NULL, 
                            covariates = NULL, marktopredict = NULL, 
                            spatial = TRUE, intercept = FALSE,
                            fun = 'exp', n.samples = 100, ...) {

  if (is.null(data) & is.null(mesh)) stop("Either data covering the entire study region or an inla.mesh object is required.")
  
  if (!all(covariates%in%row.names(object$summary.fixed))) stop("Covariates provided not in model.")
  
  if (is.null(formula) & !spatial & !intercept & is.null(covariates)) stop("Please provide at least one of spatial, intercept or covariates.")
  
  if (is.null(formula) & is.null(datasettopredict)) stop("Please provide either a formula or a dataset included in the bru_sdm model to be predicted.")
  
  #if (predictmark & is.null(markstopredict)) stop("Marks prediction is chosen but no marks to predict given.")
  
  if (is.null(data)) {
    
    if (!is.null(mask)) {
      
      data <- pixels(mesh, mask = mask)
      
    } 
    else data <- pixels(mesh)
    
  }
  
  datasettopredict <- substitute(datasettopredict)
  
  if (!is.null(marktopredict)) {dataset <- paste0(datasettopredict,'_',marktopredict)}
  
  if (is.null(formula)) {
    
    if (spatial) {
      
      if (!paste0(datasettopredict,'_spde')%in%names(object$summary.random)) stop('Either dataset name is incorrect or bru_sdm model run without spatial effects.')
      else spatial <- paste0(datasettopredict,'_spde')
      
    } 
    else spatial <- NULL
    
    if (intercept) {
      
      if (!is.null(marktopredict)) {
        
        if (!paste0(marktopredict,'_intercept')%in%row.names(object$summary.fixed)) stop('Either dataset name is incorrect or bru_sdm model run without intercepts.')
        else intercept <- paste0(marktopredict,'_intercept')
        
      }
      
      else {
        
        if (!paste0(datasettopredict,'_intercept')%in%row.names(object$summary.fixed)) stop('Either dataset name is incorrect or bru_sdm model run without intercepts.')
        else intercept <- paste0(datasettopredict,'_intercept')
        
      }
      
    } 
    else intercept <- NULL
    
    formula_components <- c(covariates, spatial, intercept)
    
    if (is.null(fun)) {fun <- ''}
    
    formula <- as.formula(paste0('~ ',as.character(fun),'(',paste(formula_components, collapse = ' + '),')'))
    
  }
  
  class(object) <- c('bru','inla','iinla')
  
  int <- predict(object, data = data, formula = formula, n.samples = n.samples, ...)
  return(int)
  
}

