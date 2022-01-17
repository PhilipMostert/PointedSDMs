#' @description An internal function used to change the coordinate names for datasets.
#' @param data A list of datasets.
#' @param oldcoords The old coordinate names.
#' @param newcoords The new coordinate names.

changeCoords <- function(data, oldcoords, newcoords) {
  
  lapply(data, function(dat) {
    
    if (class(dat) == 'data.frame') {
      
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