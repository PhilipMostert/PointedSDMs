#' Export bru_sdm class
#' 
#' @export

setClass('bru_sdm')

#' Print method for bru_sdm
#' 
#' @exportS3Method 


print.bru_sdm <- function(x, ...) {
  
  print(summary(x))
  
}

#' Summary for bru_sdm
#' 
#' @exportS3Method 

summary.bru_sdm <- function(x,...) {
  #cat('----bru_sdm summary STILL IN DEVELOPMENT----\n\n')
  
  cat("Summary of 'bru_sdm' object:\n\n")
  cat(paste0("inlabru version: ", x$bru_info$inlabru_version, "\n"))
  cat(paste0("INLA version: ", x$bru_info$INLA_version, "\n\n"))
  
  cat('Types of data modelled:\n')
  names_data = data.frame(x[['data_type']][!duplicated(names(x[['data_type']]))])
  names(names_data) = c('                              ')
  print(names_data)
  cat('\n')
  
  if(!is.null(x[['multinom_vars']])){
    cat('Summary of multinomial variables:\n\n')
    for (i in x[['multinom_vars']]) {
      if (i%in%names(x$summary.random)) { 
        variable <- x$summary.random[[i]]
        names(variable)[1] <- i
        print(variable[,1:7], row.names = FALSE, digits = 3)
        cat('\n\n')
      } else {
        
        cat('Variable',i,'not run in model')
        cat('\n\n')
        
      }
      
    }
  }
  
  
  if(any(x[['spatial_covariates_used']]%in%names(x[['summary.random']]))) {
    
    cat('Summary of categorical variables:\n\n')
    
    for (cov in x[['spatial_covariates_used']]) {
    if (cov %in%names(x[['summary.random']])) {
    cat(cov)
    cat('\n')
    print(x[['summary.random']][[cov]], row.names = FALSE, digits = 3)

    cat('\n')
    
    }
    }
  }
  
  if (!is.null(x[['species_in']])) {
  
  if (length(unique(unlist(x[['species_in']]))) > 1) {
  
  cat('Summary of the fixed effects for the species:')
  cat('\n\n')
  
  for (species in as.character(unique(unlist(x$species_in)))) {
    
  cat('Summary for:', species)
  cat('\n')
  print.data.frame(x[['summary.fixed']][grepl(paste0('\\<',species,'_'), row.names(x[['summary.fixed']])),])    
    
  cat('\n')
  
      
  }
  
  class(x) = 'inla'
  x$call = NULL
  x$summary.fixed = NULL
  summary(x)

  }
  else {
    
  class(x) = 'inla'
  x$call = NULL
  summary(x)  
  
  }
  }
  
  else {
    
    class(x) = 'inla'
    x$call = NULL
    summary(x)
    
  }
  
}
