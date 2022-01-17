makeMarks <- function(pointdata, markname, coords, speciesname,
                      proj, datasetname, marktrialname, markfamily) {
  #Don't run point covariates on marks.
    ##Select mark and get class
  
  #pointdata == data.frame object
  
  classMark <- class(pointdata[,markname])
  
  if (classMark == 'factor' || classMark == 'character') {
    
    ##First question:: how do marks work with PA data or even countdata????
     #For now assume that all multinomial marks have a response of 1?
    
    pointdata[, paste0(markname, '_response')] <- 1
    pointdata[, paste0(markname, '_phi')] <- 1
    
    markSP <- sp::SpatialPointsDataFrame(coords = pointdata[, names(pointdata) %in% coords],
                                         data = pointdata[, c(speciesname, markname,
                                                              paste0(markname, '_response'), paste0(markname, '_phi'))],
                                         proj4string = proj)
    
    markNtrials <- NULL
    markFamily <- 'poisson'
    markDatatype <- 'Multinomial mark'
    
    
  }
  else {
  
  if (!is.null(marktrialname)) {
   
    if (!marktrialname %in% names(pointdata)) marktrialname <- NULL
  
   }
  markSP <- sp::SpatialPointsDataFrame(coords = pointdata[, names(pointdata) %in% coords],
                                       data = data.frame(pointdata[, c(markname, marktrialname, speciesname)]),
                                       proj4string = proj)
  
  if (ncol(markSP@data) == 1) names(markSP@data) <- markname
  if (is.null(names(markname))) {
    
  warning(paste('No family given to mark ', markname,'. Will assume family is Gaussian.'))  
  markFamily = 'gaussian'  
  }
  else markFamily <- markfamily[markname]
  if (markFamily == 'binomial') {
    
    if (marktrialname %in% names(markSP@data)) markNtrials <- markSP@data[, marktrialname]
    else markNtrials <- NULL
  }
  else markNtrials <- NULL
  
  markDatatype <- paste(gsub("^(\\w)(\\w+)", 
                             "\\U\\1\\L\\2",
                             markFamily,
                             perl = TRUE), 'mark')
    
  }
  
  return(list(markData = markSP, markNtrials = markNtrials,
              markFamily = markFamily,
              markDatatype = markDatatype))
  
}