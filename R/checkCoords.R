#' @title \emph{checkCoords}: function used to check coordinate names.
#' @description An internal function used to check if all the coordinates are the same.
#' @param data A list of datasets.
#' @param coords A vector of length 2 of the coordinate names.
#' @return A logical variable.

checkCoords <- function(data, coords) {
  
  coords_in <- sapply(data, function(dat) {
    if (inherits(dat,'data.frame')) 
      coords %in% names(dat)
    else if (inherits(dat, "Spatial")) {
      x_coord <- colnames(dat@coords)[1]
      y_coord <- colnames(dat@coords)[2]
      coords %in% c(x_coord, y_coord)
    }
  })
  if (!all(coords_in)) OK <- FALSE
  else OK <- TRUE
  
  OK
  
}