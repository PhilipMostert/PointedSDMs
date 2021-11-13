#' Function to make formulas for the points
#' 
#' @param response The response variable names of the points.
#' @param dataset The names of the datasets.
#' @param covariates The names of the covariates used in the model.
#' @param pointsspatial Are spatial effects used for the points.
#' @param pointsintercept Are intercepts used for the points.
#' @param sharedspatial Are the spatial effects shared across the datasets.
#' @param spatiadatasets Which datasets should have spatial effects.
#' @param covariatesbydataset Are there specific covariates for specific datasets.
#' @param species The name of the species effects variable.
#' @param speciesdataset If species, what species are in the datasets.
#' @param pointcovariates Additional covariates to attach to the points.
#' 
#' @export

formula_maker <- function(response,
                          dataset,
                          covariates,
                          pointsspatial,
                          pointsintercept,
                          sharedspatial,
                          spatialdatasets,
                          covariatesbydataset,
                          species,
                          speciesdataset,
                          pointcovariates) {
  ##Create component_maker after this::
  formula <- list()
  
  for (i in 1:length(response)) {
  
  if (!is.null(covariatesbydataset)) {
    
  if (dataset[i]%in%names(covariatesbydataset)) {
      
  covs <- covariatesbydataset[[dataset[i]]]
      
  } else {
      
  covs <- covariates
      
  }
    
  } else covs <- covariates
  
  if (is.null(covariates)) covs <- '.'
  
  formula[[i]] <- formula(paste0(c(response[i],'~', paste(covs, collapse = ' + ')),collapse = " ")) 
  
  if (!is.null(species)) {
    
  species_in <- unique(speciesdataset[[i]])
    
  }
  
  if (covs == '.') covs <- NULL
  
  if (pointsintercept) {
    
  if (!is.null(species)) {
      
  if (!is.null(covs)) {
        
  if (length(unique(unlist(speciesdataset))) > 1) {  
          
  formula[[i]] <- update(formula[[i]], paste0(' ~ . +', paste(paste0(species_in,'_intercept'), collapse = ' + ')))
          
  }
  else formula[[i]] <- update(formula[[i]], paste0(' ~ . + intercept'))
        
  }
  else {
        
  resp <- as.character(formula[[i]])[2]
        
  if (length(unique(unlist(speciesdataset))) > 1) { 
          
  formula[[i]] <- formula(paste(resp, ' ~ ', paste0(species_in,'_intercept',collapse = ' + ')))
          
  }
  else formula[[i]] <- formula(paste(resp, ' ~ ', paste0('intercept',collapse = ' + ')))
        
  }
      
  }
  else
  if (!is.null(covs))  { 
  
  formula[[i]] <- update(formula[[i]], paste0(' ~ . +', paste0(dataset[i],'_intercept'), collapse = ' + '))
        
  }
  else {
      
  resp <- as.character(formula[[i]])[2]
      
  formula[[i]] <- update(formula[[i]], paste(resp, ' ~ ', paste0(dataset[i],'_intercept', collapse = ' + ')))
      
  }
    
  }
  
  if (!is.null(species)) {
    
  if (length(unique(unlist(speciesdataset))) > 1) {
    
  species_covs <- apply(expand.grid(paste0(species_in,'_'),covs), MARGIN = 1, FUN = paste0,collapse='')
  
  if (!identical(species_covs, character(0))) {
        
  #for(i in 1:length(species_covs)) {
  formula[[i]] <- update(formula[[i]], paste('~ . +', paste(species_covs, collapse = ' + ')))
          
  #}
        
  }

  }

  } 

  if (pointsspatial) {
    
  if (sharedspatial) {
      
  if (!is.null(spatialdatasets)) {
        
  if (dataset[[i]]%in%spatialdatasets) {
          
  if (is.null(covs) & ! pointsintercept) {
            
  resp <- as.character(formula[[i]])[2]
            
  formula[[i]] <- formula(paste(resp, ' ~ +', paste0('~ . +','shared_spatial')))
            
  }   
  else formula[[i]] <- update(formula[[i]], paste0('~ . +','shared_spatial'))
          
  }  
        
  }
  else
  if (is.null(covs) & !pointsintercept) {
          
  resp <- as.character(formula[[i]])[2]
          
  formula[[i]] <- formula(paste(resp, ' ~ +', paste0('~ . +','shared_spatial')))
          
  }    
  else formula[[i]] <- update(formula[[i]], paste0('~ . +','shared_spatial'))
      
  }  
  else  
  if (!is.null(spatialdatasets)) {
        
  if (dataset[[i]]%in%spatialdatasets) {
          
  if (is.null(covs) & !pointsintercept) {
            
  resp <- as.character(formula[[i]])[2]
            
  formula[[i]] <- formula(paste0(resp, ' ~ ',dataset[[i]],'_spde'))
            
  }   
  else formula[[i]] <- update(formula[[i]], paste0('~ . +',dataset[[i]],'_spde'))
          
  }
    
  }  
  else 
  if (is.null(covs) & !pointsintercept) {
        
  resp <- as.character(formula[[i]])[2]
        
  formula[[i]] <- formula(paste0(resp, ' ~ ', dataset[[i]],'_spde'))
        
  }  
  else formula[[i]] <- update(formula[[i]], paste0('~ . +',dataset[[i]],'_spde'))
    
  }

  if (!is.null(species)) {
    
  formula[[i]] <- update(formula[[i]], paste0(' ~ . + ', species, '_spde'))
    
  }
  
  if (!is.null(pointcovariates)) {
    
  formula[[i]] <- update(formula[[i]], paste(' ~ . +', paste(pointcovariates_incl, collapse = ' + ')))  
    
  }
  
  }
  
  return(formula)
  
}
  
