
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
#' @param object bruSDM object.
#' @param ... Un used argument
#' @rdname summary
#' @exportS3Method 

summary.bruSDM <- function(object, ...) {
  #cat('----bru_sdm summary STILL IN DEVELOPMENT----\n\n')
  
  cat("Summary of 'bruSDM' object:\n\n")
  cat(paste0("inlabru version: ", object$bru_info$inlabru_version, "\n"))
  cat(paste0("INLA version: ", object$bru_info$INLA_version, "\n\n"))
  
  cat('Types of data modelled:\n')
  names_data = data.frame(object[['dataType']][!duplicated(names(x[['dataType']]))])
  names(names_data) = c('                              ')
  print(names_data)
  cat('\n')
  
  if(!is.null(object[['multinomVars']])){
    cat('Summary of multinomial variables:\n\n')
    for (i in x[['multinomVars']][order(x[['multinomVars']])]) {
      if (i%in%names(object$summary.random)) { 
        variable <- object$summary.random[[i]]
        names(variable)[1] <- i
        print(variable[,1:7], row.names = FALSE, digits = 3)
        cat('\n\n')
      } else {
        
        cat('Variable',i,'not run in model')
        cat('\n\n')
        
      }
      
    }
  }
  
  if (!is.null(object[['species']][['speciesIn']])) {
      
      cat('Summary of the fixed effects for the species:')
      cat('\n\n')
      
      for (species in as.character(unique(unlist(object[['species']]['speciesIn'])))) {
        
        cat('Summary for', paste0(species,':'))
        cat('\n')
        print.data.frame(x[['summary.fixed']][grepl(paste0('\\<',species,'_'), row.names(object[['summary.fixed']])),])    
        
        cat('\n')
        
        
      }
      
      class(object) = 'inla'
      object$call = NULL
      object$summary.fixed = NULL
      summary(object)
      
  }
  
  else {
    
    class(object) = 'inla'
    object$call = NULL
    summary(object)
    
  }
  
}