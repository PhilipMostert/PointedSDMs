
#' Export bru_sdm class
#' 
#' @export

setClass('bruSDM')

#' Print method for bru_sdm
#' @param x bruSDM object.
#' @param ... Un used argument.
#' 
#' @exportS3Method 

print.bruSDM <- function(x, ...) {
  
  print(summary(x))
  
}

#' Summary for bru_sdm
#' @param x bruSDM object.
#' @param ... Un used argument.
#' @rdname summary
#' @export

summary.bruSDM <- function(x, ...) {
  #cat('----bru_sdm summary STILL IN DEVELOPMENT----\n\n')
  
  cat("Summary of 'bruSDM' object:\n\n")
  cat(paste0("inlabru version: ", x$bru_info$inlabru_version, "\n"))
  cat(paste0("INLA version: ", x$bru_info$INLA_version, "\n\n"))
  
  cat('Types of data modelled:\n')
  names_data = data.frame(x[['dataType']][!duplicated(names(x[['dataType']]))])
  names(names_data) = c('                              ')
  print(names_data)
  cat('\n')
  
  if(!is.null(x[['multinomVars']])){
    cat('Summary of multinomial variables:\n\n')
    for (i in x[['multinomVars']][order(x[['multinomVars']])]) {
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
  
  if (!is.null(x[['species']][['speciesIn']])) {
      
      cat('Summary of the fixed effects for the species:')
      cat('\n\n')
      
      for (species in as.character(unique(unlist(x[['species']]['speciesIn'])))) {
        
        cat('Summary for', paste0(species,':'))
        cat('\n')
        print.data.frame(x[['summary.fixed']][grepl(paste0('\\<',species,'_'), row.names(x[['summary.fixed']])),])    
        
        cat('\n')
        
        
      }
      
      class(x) = 'inla'
      x$call = NULL
      x$summary.fixed = NULL
      summary(x, ...)
      
  }
  
  else {
    
    class(x) = 'inla'
    x$call = NULL
    summary(x, ...)
    
  }
  
}