setClass('bru_sdm_leave_one_out')

#' Print for bru_sdm_leave_one_out
#' 
#' @export print.bru_sdm_leave_one_out

print.bru_sdm_leave_one_out <- function(x, ...) {
  
  for (name in 1:length(x)) {
    
    cat('Changes in fixed values by leaving out', paste0(names(attributes(x)$differences)[name],':'))
    cat('\n\n')
    
    print(attributes(x)$differences[[name]])
    
    if (!is.null(attributes(x)$validation_results)) {
      cat('\n')
      
      cat('Leave one out cross-validation score:', attributes(x)$validation_results[name])
      
   }
    
    
    cat('\n\n')
    
  }
} 
