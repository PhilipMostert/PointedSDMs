#' @title removeFormula: Function to remove term from a formula.
#' @description formulaChanger: Internal function used to remove formula terms.
#' @param formulaRemove The formula with all the components to remove.
#' @param oldFormula The formula which needs to change.

removeFormula <- function(formulaRemove,
                           oldFormula) {
  
  updated_formula <- update(formula(paste('~ ', paste0(oldFormula, collapse = ' + '))), formulaRemove)
  
  termsIn <- labels(terms(updated_formula))
  
  termsIn <- termsIn[!termsIn %in% oldFormula]
  
  updated_formula <- all.vars(updated_formula)[all.vars(updated_formula) != '.']
  
  updated_formula
  

}