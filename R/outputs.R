#' Bru_sdm output functions
#' 

print.bru_sdm <- function(x, ...) {
  
  print(summary(x))
  
}

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
