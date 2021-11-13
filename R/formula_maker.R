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


#' Function to make formulas for the marks
#' @param responsemarks Names of the response variables for the marks.
#' @param namesdatasetmarks Names of the datasets of the marks. In the form: 'dataset','_','mark'.
#' @param marksname The name of the marks.
#' @param covariates The names of the covariates used in the model.
#' @param covariatesbydataset Are there specific covariates for specific datasets.
#' @param marksspatial Are spatial effects considered for the marks.
#' @param sharedspatial Are the spatial effects shared across the marks.
#' @param spatialdatasets Which datasets should have spatial effects.
#' @param marksintercept Are intercepts run for the marks.
#' @param marktype What is the marks type.
#' @param markphi Phi variable for mark if multinomial.
#' 
#' @export

mark_formula_maker <- function(responsemarks,
                               namesdatasetmarks,
                               marksname,
                               covariates,
                               covariatesbydataset,
                               marksspatial,
                               sharedspatial,
                               spatialdatasets,
                               marksintercept,
                               marktype,
                               markphi) {
  
  formula_marks <- list()
  
  for (i in 1:length(responsemarks)) {
    
  if (!is.null(covariatesbydataset)) {
      
  if (gsub('_.*$',"",namesdatasetmarks)%in%names(covariatesbydataset)) {
        
  markscovs <- covariatesbydataset[[gsub('_.*$',"",namesdatasetmarks[[i]])]]  
        
  } 
  else markscovs <- covariates
      
  } 
  else markscovs <- covariates
    
  if (is.null(covariates)) markscovs <- '.'
    
  formula_marks[[i]] <- formula(paste0(c(responsemarks[i],'~',paste(markscovs, collapse = ' + ')),collapse = " "))
    
  if (markscovs == '.') markscovs <- NULL
    
  if (marksspatial) {
      
  if (sharedspatial) {
        
  if (!is.null(spatialdatasets)) {
          
  if (gsub('_.*$',"",namesdatasetmarks[[i]])%in%spatialdatasets) {
            
  if (is.null(markscovs)) {
              
  formula_marks[[i]] <- formula(paste(responsemarks[i], ' ~ + shared_spatial'))  
              
  }
  else formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +",'shared_spatial'))
            
  }   
          
  }
  else
  if (is.null(markscovs)) {
            
  formula_marks[[i]] <- formula(paste(responsemarks[i], ' ~ + shared_spatial'))  
            
  }
  else formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +",'shared_spatial'))
        
  }
  else
  if (!is.null(spatialdatasets)) {
          
  if (gsub('_.*$',"",namesdatasetmarks[[i]])%in%spatialdatasets) {
            
  if (is.null(markscovs)) {
              
  formula(paste(responsemarks[i], ' ~ + ',paste0(namesdatasetmarks[[i]],'_spde')))   
              
  }
  else formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +",namesdatasetmarks[[i]],'_spde'))
            
  }
          
  }  
  else
  if (is.null(markscovs)) {
          
  formula(paste(responsemarks[i], ' ~ + ',paste0(namesdatasetmarks[[i]],'_spde')))   
          
  }  
  else formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +",namesdatasetmarks[[i]],'_spde'))
      
  }
  
  if (marksintercept) {
      
  if (marktype != 'Multinomial mark') {
        
  if (is.null(markscovs) & !marksspatial) {
          
  formula_marks[[i]] <- formula(paste(responsemarks[i], ' ~ + ',paste0(namesdatasetmarks[[i]],'_intercept')))  
          
  }
  else formula_marks[[i]] <- update(formula_marks[[i]],paste0('. ~ . +', paste0(namesdatasetmarks[[i]],'_intercept'), collapse = ' + '))
        
  }
      
  }
    
  if (marktype == 'Multinomial mark') {
      
  if (is.null(markscovs) & !marksspatial) {
        
  formula_marks[[i]] <- formula(paste0(paste0(namesdatasetmarks[[i]],'_response'), ' ~ +',  paste(marksname, markphi, sep = ' + ')))
        
  }
      
  formula_marks[[i]] <- update(formula_marks[[i]], paste0(paste0(namesdatasetmarks[[i]],'_response'), ' ~ . + ', paste(marksname, markphi, sep = ' + ')))
      
  }
    
  }
  
  return(formula_marks)
  
}
