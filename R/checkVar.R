#' @description Internal function used to check if variable is in dataset.
#' 
#' @param data A list of datasets.
#' @param var A variable name.

checkVar <- function(data, var) {
  
  var_in <- sapply(data, function(dat) {
    if (class(dat) == "data.frame") 
      var %in% names(dat)
    else if (inherits(dat, "Spatial")) {
   
      var %in% names(dat@data)
      
    }
  })
  if (!all(var_in)) OK <- FALSE
  else OK <- TRUE
  
  OK
  
}