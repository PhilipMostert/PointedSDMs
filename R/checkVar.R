#' @title \emph{checkVar}: Function used to check variable names.
#' @description Internal function used to check if variable is in dataset.
#' 
#' @param data A list of datasets.
#' @param var A variable name.
#' 
#' @return A logical variable

checkVar <- function(data, var) {
  
  var_in <- sapply(data, function(dat) {
    if (inherits(dat,'data.frame')) 
      var %in% names(dat)
    else if (inherits(dat, "Spatial")) {
   
      var %in% names(dat@data)
      
    }
  })
  if (!all(var_in)) OK <- FALSE
  else OK <- TRUE
  
  OK
  
}