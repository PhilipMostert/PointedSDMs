#' @title \emph{makeFormulaComps}: function to make components for the covariate and bias Formulas.
#' @description An internal function used to make the formula components required for inlabru.
#' @param form The formula which needs to be changed into a component.
#' @param species Logical indicating if species occur in the model.
#' @param speciesnames The names of the species occurring in the model.
#' @param type What type of component is being created: one for the covariate formula or the bias formula.
#' @return A vector of components required for inlabru.

makeFormulaComps <- function(form, species, speciesnames, type) {
  
  ##CHANGE PREDICT FUNCTION WITH ID FOR DIFFERENT TERMS
  
  terms <- labels(terms(form))
  
  if (type == 'Bias') frontPart <- 'Bias__Effects__Comps(main = '
  else frontPart <- 'Fixed__Effects__Comps(main = '
  
  if (!species) {
    
  terms <- paste(terms, collapse = ' + ')
  
  newFormula <- formula(paste('~ ', terms, '-1 '))
  
  newComps <- paste0(frontPart, paste(newFormula, collapse = ''),', model = "fixed")')
  
  newComps
  
  }
  
  else {
    
    newTerms <- list()
    newList <- list()
    
    for (species in speciesnames) {
      
      speciespart <- paste0('\\1',species,'_\\2')
      
      for (cov in terms) {
        
        covpart <- paste0("(.*)(",cov,')')
        
        newParts <- sub(covpart, speciespart, terms)
        newParts <- newParts[!newParts %in% terms]
        
        if (all(newParts %in% terms)) newTerms[[species]][[cov]] <- NULL ##Assume there won't be some half parts?
        else newTerms[[species]][[cov]] <- newParts
        
      }
      
      speciesFormula <- formula(paste('~', paste(unique(unlist(newTerms[[species]])), collapse = ' + '), '-1'))
      
      
      newList[[species]] <- paste0(species, '_', frontPart, paste(speciesFormula, collapse = ''),', model = "fixed")')
      
    }
    
    unlist(newList)
    
    
  }
  
}