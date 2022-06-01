#' @title \emph{nameChanger}: function to change a variable name.
#' @description An internal function used to change the name of a variable.
#' @param data A list of datasets.
#' @param oldName The old variable name.
#' @param newName The new variable name.
#' @return A list of data.frame or spatial objects with the name of the variable changes.

nameChanger <- function(data, oldName, newName) {
  
  lapply(data, function(dat) {
    
    if (inherits(dat,'data.frame')) {
      
      if (oldName %in% names(dat)) {
        
        names(dat)[names(dat) == oldName] <- newName
        dat
      }
      else dat

      
    }
    else {
      
      if (oldName %in% names(dat@data)) {
        
        names(dat@data)[names(dat@data) == oldName] <- newName
        dat
        
      }
      else dat
      
    }
    
    
  })
  
  
}