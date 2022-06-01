#' @title \emph{changeCoords}: function used to change coordinate names.
#' @description An internal function used to change the coordinate names for datasets.
#' @param data A list of datasets.
#' @param oldcoords The old coordinate names.
#' @param newcoords The new coordinate names.
#' @return A list of data.frame or spatial objects with the coordinate names changed.

changeCoords <- function(data, oldcoords, newcoords) {
  
  lapply(data, function(dat) {
    
    if (inherits(dat,'data.frame')) {
      
      if (!all(oldcoords %in% names(dat))) stop('Coordinate names specified not in dataset.')
      else {
        
        names(dat)[names(dat) == oldcoords] <- newcoords
        dat
        
      }
      
    }
    else {
    
      if (!all(oldcoords %in% colnames(dat@coords))) stop('Coordinate names specified not in dataset.')
      else {
        
        colnames(dat@coords) <- newcoords
        dat
      
      }
      
    }
    
    
  })
  
}