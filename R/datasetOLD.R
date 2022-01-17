##Thing to do for tomorrow::
#Create a vector of all possible vars in the data object.
# i.e what point covs are in the data (and for multinom all the phi etc.)
# and then make an includeComps function
# with arguments Intercept, spatial, pointcovs, spcovs.


##Change this::
#Dont have the makeMarks function since all we are doing is data duplication:
#Rather create a list of all the marks within each dataset and recycle data in likelihood construction
#And then name appropriatly in lik construction paste0(mark_dataset)
#Obvs still need to create all phi vars given multinomial mark

makeData <- function(datapoints, datanames, coords, proj, pointcovnames,
                     paresp, countsresp, trialname, speciesname,
                     marks, marktrialname, markfamily) {
  
  if (paresp != 'poresp' || countsresp != 'poresp') poresp <- 'poresp'
  else
    if (paresp != 'respPO' || countsresp != 'respPO') poresp <- 'respPO'
    else poresp <- 'POresponse'
    
    if (length(datapoints) != length(datanames)) stop('Number of dataset names needs to equal length of datasets.')
    
    dataOrganized <- list()
    Ntrials <- list()
    Family <- list()
    dataType <- list()
    
    marksData <- list()
    marksNtrials <- list()
    marksFamily <- list()
    marksType <- list()
    
    for (dat in 1:length(datapoints)) {
      
      datasetname <- datanames[dat]
      dtSub <- datapoints[[dat]]
      data <- as.data.frame(dtSub)
      data_vars <- names(data)
      
      if (!is.null(marks)) marksInData <- data_vars[data_vars %in% marks]
      else marksInData <- NULL
      if (length(marksInData) > 0) {
        
        for (mark in 1:length(marksInData)) {
          if (length(marksData) == 0) index <- 1  
          else index <- length(marksData) + 1  
          markObtained <- makeMarks(pointdata = data, markname = marksInData[mark],
                                    coords = coords, speciesname = speciesname,
                                    datasetname = datanames, marktrialname = marktrialname,
                                    markfamily = markfamily, proj = proj)
          marksData[[index]] <- markObtained[['markData']]
          marksNtrials[index] <- markObtained['markNtrials']
          marksFamily[index] <- markObtained['markFamily']
          marksType[index] <- markObtained['markDatatype']
          
          names(marksData)[[index]] <- paste0(datasetname,'_',marksInData[mark])
          
        }
        
      }
      
      if (paresp %in% data_vars) {
        
        if (!trialname %in% data_vars) subtrialname <- NULL
        else subtrialname <- trialname
        
        datSP <- sp::SpatialPointsDataFrame(coords = data[,coords], 
                                            data = data.frame(data[, c(paresp, subtrialname, pointcovnames)]),
                                            proj4string = proj)
        
        if (ncol(datSP@data) == 1) names(datSP@data) <- paresp
        
        if (!is.null(trialname)) {
          
          if (trialname %in% data_vars) Ntrials[dat] <- data.frame(datSP@data[, trialname])[, 1]
          else Ntrials[dat] <- NULL
          
        }
        else Ntrials[dat] <- NULL
        
        Family[dat] <- "binomial"
        dataType[dat] <- "Present absence"
        dataOrganized[[dat]] <- datSP
        
      }
      else 
        if (countsresp %in% data_vars) {
          
          datSP <- sp::SpatialPointsDataFrame(coords = data[, coords],
                                              data = data.frame(data[, c(countsresp, pointcovnames)]),
                                              proj4string = proj)
          
          if (ncol(datSP@data) == 1) names(datSP@data) = countsresp
          
          Ntrials[dat] <- NULL
          Family[dat] <- "poisson"
          dataType[dat] <- "Count data"
          dataOrganized[[dat]] <- datSP
          
        }
      else {
        
        data[, poresp] <- rep(1, nrow(data))
        
        data_vars <- c(data_vars, poresp)
        
        datSP <- sp::SpatialPointsDataFrame(coords = data[, coords], 
                                            data = data.frame(data[, c(poresp,pointcovnames)]),
                                            proj4string = proj)
        
        if (ncol(datSP@data) == 1) names(datSP@data) <- poresp
        Ntrials[dat] <- NULL
        Family[dat] <- "cp"
        dataType[dat] <- "Present only"
        dataOrganized[[dat]] <- datSP
        
      }
    }
    
    names(dataOrganized) <- datanames
    
    return(list(data = dataOrganized,
                Ntrials = Ntrials,
                Family = Family,
                dataType = dataType,
                marksData = marksData,
                marksNtrials = marksNtrials,
                marksFamily = marksFamily,
                marksType = marksType))
    
}