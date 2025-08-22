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
  
  if (!species %in% c('stack', 'community')) {
    
  terms <- paste(terms, collapse = ' + ')
  
  newFormula <- formula(paste('~ ', terms, '-1 '))
  
  newComps <- paste0(frontPart, paste(newFormula, collapse = ''),', model = "fixed")')
  
  newComps
  
  }
  
  else {
    
    #newTerms <- list()
    newList <- list()
    
    for (species in speciesnames) {
      
      terms <- paste(terms, collapse = ' + ')
      
      newFormula <- formula(paste('~ ', terms, '-1 '))
      
      newList[[species]] <- paste0(species, '_', frontPart, paste(newFormula, collapse = ''),', model = "fixed")')

    }
    
    unlist(newList)
    
    
  }
  
}