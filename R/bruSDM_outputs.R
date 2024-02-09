
#' Export bru_sdm class
#' 
#' @export

setClass('bruSDM')

#' Print method for bru_sdm
#' @title Generic print function for \code{bruSDM}.
#' @param x bruSDM object.
#' @param ... Un used argument.
#' 
#' @exportS3Method 

print.bruSDM <- function(x, ...) {
  
  print(summary(x))
  
}

#' Summary for bru_sdm
#' @title Generic summary function for \code{bruSDM}.
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
  names_data = data.frame(object[['dataType']][!duplicated(names(object[['dataType']]))])
  names(names_data) = c('                              ')
  print(names_data)
  cat('\n')
  
  if(!is.null(object[['multinomVars']])){
    cat('Summary of multinomial variables:\n\n')
    for (i in object[['multinomVars']][order(object[['multinomVars']])]) {
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
  
  if (!is.null(object[['species']][['speciesIn']]) 
      && any(unlist(object[['species']][['speciesEffects']]))) {
      
      cat('Summary of the fixed effects for the species:')
      cat('\n\n')
      
      for (species in as.character(unique(unlist(object[['species']]['speciesIn'])))) {
        
        cat('Summary for', paste0(species,':'))
        cat('\n')
        
        if (!is.null(object$spatCovs$covariateFormula)) {
          
          speciesCovs <- object$summary.random[[paste0(species, '_Fixed__Effects__Comps')]]
          row.names(speciesCovs) <- speciesCovs$ID
          speciesCovs <- speciesCovs[, -c(1)]
          
        } else speciesCovs <- data.frame()
        
        if (!is.null(object$spatCovs$biasFormula)) {
          
          biasCovs <- object$summary.random[[paste0(species, '_Bias__Effects__Comps')]]
          row.names(biasCovs) <- biasCovs$ID
          biasCovs <- biasCovs[, -c(1)]
          
        } else biasCovs <- data.frame()
        
        if (any(paste0(species, '_', object$spatCovs$name) %in% names(object$summary.random))) {
          
          factorCovs <- do.call(rbind, object$summary.random[paste0(species, '_', object$spatCovs$name)])
          row.names(factorCovs) <- paste0(species, '_', factorCovs$ID)
          factorCovs$ID <- NULL
        } else factorCovs <- data.frame()
        
        if(object$species$speciesEffects$Intercepts) {
          
          interceptTerms <- object$summary.random[[paste0(object$species$speciesVar, '_intercepts')]]
          interceptTerms <- interceptTerms[row.names(interceptTerms) == species,]
          row.names(interceptTerms) <- paste0(row.names(interceptTerms), '_random_intercept')
          interceptTerms$ID <- NULL
          
        } else interceptTerms <- data.frame()
        print.data.frame(rbind(object[['summary.fixed']][grepl(paste0('\\<',species,'_'), row.names(object[['summary.fixed']])),], speciesCovs, factorCovs, interceptTerms, biasCovs))   
        
        cat('\n')
        
        
      }
      
      class(object) = 'inla'
      object$call = NULL
      object$summary.fixed = NULL
      summary(object)
      
  }
  
  else {
    
    if (!is.null(object$spatCovs$covariateFormula)) {
      
        fixedEffects <- object$summary.random[['Fixed__Effects__Comps']]
        row.names(fixedEffects) <- fixedEffects$ID
        fixedEffects <- fixedEffects[,-c(1)]
        object$summary.fixed <- rbind(object$summary.fixed, fixedEffects)
        
    }
    
    if (!is.null(object$spatCovs$biasFormula)) {
      
      biasComps <- object$summary.random[['Bias__Effects__Comps']]
      row.names(biasComps) <- biasComps$ID
      biasComps <- biasComps[, -c(1)]
      biasComps <- biasComps 
      object$summary.fixed <- rbind(object$summary.fixed, biasComps)
      
    }

    class(object) = 'inla'
    object$call = NULL
    summary(object)
    
  }
  
}