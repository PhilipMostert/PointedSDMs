#' Internal function used to standardize datasets, as well as assign metadata.
#' @description Internal function used to assist in structuring the data.
#' @param datapoints A list of datasets as either sf, data.frame or SpatialPoints objects.
#' @param datanames A vector of the names of the datasets.
#' @param coords Names of the coordinates used in the model.
#' @param proj The projection reference system used in the model.
#' @param marktrialname Name of the trial variable used by the binomial marks.
#' @param paresp Name of the response variable used by the presence absence datasets.
#' @param countsresp Name of the response variable used by the counts data.
#' @param trialname Name of the trial variable used by the presence absence datasets.
#' @param speciesname Name of the species name variable.
#' @param temporalvar Name of the temporal variable.
#' @param marks Name of the marks considered in the model.
#' @param pointcovnames Name of the point covariates used in the model.
#' @param markfamily A vector describing the distribution of the marks.
#' @param offsetname Name of the offset column in the datasets.
#' 
#' @return A list of relevant metadata
#' 
#' @export

dataSet <- function(datapoints, datanames, coords, proj, pointcovnames,
                    paresp, countsresp, trialname, speciesname,
                    marks, marktrialname, markfamily, temporalvar, offsetname) {
  
  if (paresp != 'poresp' || countsresp != 'poresp') poresp <- 'poresp'
  else
    if (paresp != 'respPO' || countsresp != 'respPO') poresp <- 'respPO'
    else poresp <- 'POresponse'
    
    if (length(datapoints) != length(datanames)) stop('Number of dataset names needs to equal length of datasets.')
    
    if (!is.null(marks)) {
      
      if (length(marks) != length(markfamily)) stop('Number of marks needs to equal the number of mark families.')
      
    }
    ##Things to consider for changes here...
    #When making likelihoods it goes in the order:
     #Dataset -> species -> process
    #Note that there can only be 2 Ntrials in a given dataset (one for points; one for marks)
    #So is it worth creating an Ntrials list here?
    
    #What I think I need is:
     #A vector of what processes are in each dataset
     #ie points response + marks
     #Name the list by dataset
    
    #Keep the get family for multinomial marks here

    dataOrganized <- vector(mode = 'list', length = length(datapoints))
    names(dataOrganized) <- datanames
    
    #Won't need Ntrials
    #Ntrials <- list()
    #Keep family -- make it similar to the prcoesses in lisr
    Family <- list()
    #Keep dataType
    dataType <- c()
    numObs <- c()
    Marks <- list()
    marksType <- list()
    #Keep processIn
    #Don't actually need? Just use the names from the formula?
    #processIn <- list()
    #Keep multinomVars
    multinomVars <- c()
    varsIn <- list()
    
    for (dat in 1:length(datapoints)) {
      
      datasetname <- datanames[dat]
      
      if (inherits(datapoints[[dat]], 'sf')) {
        
        st_geometry(datapoints[[dat]]) <- 'geometry'
        coordsSF <- sf::st_coordinates(datapoints[[dat]])
        colnames(coordsSF) <- coords
        datapoints[[dat]][, coords] <- coordsSF
        
      }
      
      dtSub <- datapoints[[dat]]
      data <- as.data.frame(dtSub)
      data_vars <- names(data)
      
      varsin <- data_vars[data_vars %in% c(pointcovnames, offsetname)]
      
      if (identical(varsin, character(0))) varsin <- NULL
      
      marksin <- marks[marks %in% data_vars]
   
      if (identical(marksin, character(0))) marksin <- NULL
        
      if (!is.null(marktrialname)) {
      
        if (!marktrialname %in% data_vars) MTrialssub <- NULL
      
        else MTrialssub <- marktrialname
      }
      else MTrialssub <- NULL
      
      if (!is.null(marksin)) {

        markstype <- paste0(gsub("^(\\w)(\\w+)", 
                                 "\\U\\1\\L\\2",
                                 markfamily[marks %in% marksin],
                                 perl = TRUE),' mark')
   
        if (length(marksin) == 1) classMarks <- class(data[, marksin])
        else classMarks <- sapply(data[, marksin], class)

      if (any(classMarks == 'character' | classMarks == 'factor')) {

        characterMarks <- marksin[classMarks %in% c('character', 'factor')]
      
        for (mark in characterMarks) {
          ##Need a list to say which marks are multinomial
          #To add to the formulas.
        data[, paste0(mark, '_phi')] <- rep(1,nrow(data))
        
        if (paresp %in% data_vars) charResp <- data[, paresp]
        else 
          if (countsresp %in% data_vars) charResp <- data[,countsresp]
          else charResp <- rep(1,nrow(data))
        
        data[, paste0(mark,'_response')] <- charResp#rep(1,nrow(data))
        
        phiVars <- paste0(characterMarks, '_phi')
        responseVars <- paste0(characterMarks,'_response')
    
        #marksin <- c(marksin, paste0(mark, '_phi'), paste0(mark,'_response')) ##Do I want these phi's here??
        
        markfamily[match(mark, marksin)] <- 'poisson'
        markstype[match(mark, marksin)] <- 'Multinomial mark'
        
        multinomVars <- c(multinomVars, characterMarks)
   
      }

      }
        
      else {
        
      characterMarks <- NULL
      phiVars <- NULL
      responseVars <- NULL
        
      }
      
      } 
      else {
        
      classMarks <- NULL
      #markfamily <- NULL
      markstype <- NULL
      phiVars <- NULL
      responseVars <- NULL
      
      }
      
      if (!is.null(speciesname)) {
        
        speciesIndexVAR <- paste0(speciesname, 'INDEX_VAR')
        data[, speciesIndexVAR] <- data[, speciesname]
        
      } else speciesIndexVAR <- NULL
      
      if (paresp %in% data_vars) {
      if (!is.null(trialname)) {  
        if (!trialname %in% data_vars) subtrialname <- NULL
        else subtrialname <- trialname
      } 
      else subtrialname <- NULL
      
        #datSP <- sp::SpatialPointsDataFrame(coords = data[,coords], 
        #                                    data = data.frame(data[, c(paresp, subtrialname, temporalvar,
        #                                                               marksin, MTrialssub, speciesname,
        #                                                               phiVars, responseVars,varsin)]),
        #                                    proj4string = proj)
        datSP <- sf::st_as_sf(x = data.frame(data[, c(paresp, subtrialname, temporalvar,
                              marksin, MTrialssub, speciesname, coords, speciesIndexVAR,
                              phiVars, responseVars,varsin)]),
                              coords = coords,
                              crs = proj)
        
        if (ncol(datSP[names(datSP) != 'geometry']) == 1) names(datSP[names(datSP) != 'geometry']) <- paresp
        
        #if (!is.null(trialname)) {
        #  
        #  if (trialname %in% data_vars) Ntrials[[dat]] <- data.frame(datSP@data[, trialname])[, 1]
        #  else Ntrials[[dat]] <- rep(1, nrow(datSP@coords))

          
        #}
        #else Ntrials[[dat]] <- rep(1, nrow(datSP@coords))

        family <- 'binomial'
        #Family[dat] <- "binomial"
        datatype <- "Present absence"
        dataOrganized[[dat]][[1]] <- datSP
        
      }
      else 
        if (countsresp %in% data_vars) {
          
          #datSP <- sp::SpatialPointsDataFrame(coords = data[, coords],
          #                                    data = data.frame(data[, c(countsresp, marksin, temporalvar,
          #                                                               speciesname, MTrialssub,
          #                                                               phiVars, responseVars, varsin)]),
          #                                    proj4string = proj)
          
          datSP <- sf::st_as_sf(x = data.frame(data[, c(countsresp, marksin, temporalvar,
                                speciesname, MTrialssub, coords, speciesIndexVAR,
                                phiVars, responseVars, varsin)]),
                                coords = coords,
                                crs = proj)
          
          if (ncol(datSP[names(datSP) != 'geometry']) == 1) names(datSP[names(datSP) != 'geometry']) <- countsresp
          
          #Ntrials[[dat]] <- rep(1, nrow(datSP@coords))
          family <- 'poisson'
          #Family[dat] <- "poisson"
          datatype <- "Count data"
          dataOrganized[[dat]][[1]] <- datSP
          
        }
        else {
        
        family <- 'cp'

        data[, poresp] <- rep(1, nrow(data))
        
        data_vars <- c(data_vars, poresp)
    
        #datSP <- sp::SpatialPointsDataFrame(coords = data[, coords], 
        #                                    data = data.frame(data[, c(poresp, marksin, temporalvar,
        #                                                               speciesname, MTrialssub,
        #                                                               phiVars, responseVars, varsin)]),
        #                                    proj4string = proj)
        
        datSP <- sf::st_as_sf(x = data.frame(data[, c(poresp, marksin, temporalvar,
                              speciesname, MTrialssub, coords, speciesIndexVAR,
                              phiVars, responseVars, varsin)]),
                              coords = coords,
                              crs = proj)
        
        if (ncol(datSP[names(datSP) != 'geometry']) == 1) names(datSP[names(datSP) != 'geometry']) <-poresp
        #Ntrials[[dat]] <- rep(1, nrow(datSP@coords))
      
        #Family[dat] <- "cp"
        datatype <- "Present only"
        dataOrganized[[dat]][[1]] <- datSP
       
        }
      
      
      #multinomVars[[dat]] <- characterMarks
      #processIn <- c(pointsresp, marksin)
      
      #if (is.null(marksin)) marksin <- NA
      #if (is.null(varsin)) varsin <- NA
      
      #marksIn[[dat]] <- list(marksin)
      #marksClass[[dat]] <- list(classMarks)
      #marksFamily[[dat]] <- markfamily[marks %in% marksIn[[dat]][[1]]]
      #marksType[[dat]] <- list(markstype)
      #if (!is.null(MTrialssub)) marksNtrials[dat] <- datSP@data[,MTrialssub]
      #else marksNtrials[[dat]] <- rep(1, nrow(datSP@coords))
      numobs <- nrow(datSP)
      numObs[dat] <- numobs
      varsIn[[dat]] <- list(varsin)
      Family[[dat]] <- c(family, markfamily[marks %in% marksin])
      
      if (!is.null(marks)) {
        
        Marks[[dat]] <-  marks[marks %in% marksin]
        if (identical(Marks[[dat]], character(0))) Marks[[dat]] <- NA
        
        if (is.null(markstype)) markstype <- NA
        marksType[[dat]] <- markstype
        names(marksType[[dat]]) <- Marks[[dat]]
      
      }
      else {
      
        Marks[[dat]] <- NA
        marksType[[dat]] <- NA
      
      }


      dataType[dat] <- datatype
    
    }

    names(dataOrganized) <- datanames
    names(Family) <- datanames
    names(dataType) <- datanames
    names(Marks) <- datanames
    names(marksType) <- datanames
    #names(dataClass) <- datanames
    #names(processIn) <- datanames
    #names(marksIn) <- datanames
    names(varsIn) <- datanames
    
    #if any identical to character 0, set to NULL
    #Maybe split dataType into pointsType and MarksType; makes it less dodge
    
    #return(list(Data = dataOrganized,
    #            Ntrials = Ntrials,
    #            Family = Family,
    #            dataType = dataType,
    #            marksClass = marksClass,
    #            marksNtrials = marksNtrials,
    #            marksFamily = marksFamily,
    #            marksIn = marksIn,
    #            marksType = marksType,
    #            varsIn = varsIn))
    
    ##New return

    return(list(Data = dataOrganized,
                Family = Family,
                dataType = dataType,
                varsIn = varsIn,
                Marks = Marks,
                marksType = marksType,
                multinomVars = multinomVars,
                numObs = numObs))
    
}
