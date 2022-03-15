#' @title R6 class for creating a dataSDM object to be used in runModel
#' @description A data object containing the data and the relevant information about the integrated model.
#' @export
#' @importFrom R6 R6Class

dataSDM <- R6::R6Class(classname = 'dataSDM', lock_objects = FALSE, cloneable = FALSE)

dataSDM$set('private', 'Projection', NULL)
dataSDM$set('private', 'Coordinates', NULL)
dataSDM$set('private', 'responseCounts', NULL)
dataSDM$set('private', 'responsePA', NULL)
dataSDM$set('private', 'trialsPA', NULL)
dataSDM$set('private', 'trialsMarks', NULL)
dataSDM$set('private', 'markFamily', NULL)
dataSDM$set('private', 'speciesName', NULL)
dataSDM$set('private', 'speciesIn', NULL)

dataSDM$set('private', 'Boundary', NULL)
dataSDM$set('private', 'INLAmesh', NULL)
dataSDM$set('private', 'Components', NULL)
#dataSDM$set('private', 'spatialModel', NULL)
dataSDM$set('private', 'markNames', NULL)
dataSDM$set('private', 'pointCovariates', NULL)
dataSDM$set('private', 'ptcovsClass', NULL)
dataSDM$set('private', 'initialnames', NULL)

dataSDM$set('private', 'temporalName', NULL)
dataSDM$set('private', 'temporalVars', NULL)
dataSDM$set('private', 'temporalModel', NULL)
dataSDM$set('private', 'speciesSpatial', TRUE)

dataSDM$set('private', 'modelData', list())
dataSDM$set('private', 'pointsField', list())
##Make speciesField a named list for each species
 # if no field given for a specific species: then inla.spde2.matern()
#dataSDM$set('private', 'speciesField', list())
#dataSDM$set('private', 'biasField', list())
#dataSDM$set('private', 'marksField', list())

dataSDM$set('private', 'spatcovsObj', NULL)
dataSDM$set('private', 'spatcovsNames', NULL)
dataSDM$set('private', 'spatcovsEnv', NULL)
dataSDM$set('private', 'spatcovsClass', NULL)
dataSDM$set('private', 'dataSource', NULL)

dataSDM$set('private', 'Spatial', TRUE)
dataSDM$set('private', 'marksSpatial', TRUE)
dataSDM$set('private', 'Intercepts', TRUE)
dataSDM$set('private', 'marksIntercepts', TRUE)
dataSDM$set('private', 'IPS', NULL)
dataSDM$set('private', 'multinomVars', NULL)
dataSDM$set('private', 'printSummary', NULL)
dataSDM$set('private', 'multinomIndex', list())
dataSDM$set('private', 'optionsINLA', list())

#' @description Initialize function for dataSDM: used to store some compulsory arguments.
#' @param coordinates A vector of length 2 containing the names of the coordinates.
#' @param projection The projection of the data.
#' @param Inlamesh An inla.mesh object.
#' @param initialnames The names of the datasets if data is passed through bruSDM.
#' @param responsecounts The name of the response variable for the count data.
#' @param responsepa The name of the response variable for the presence absence data.
#' @param marksnames The names of the marks contained in the data.
#' @param marksfamily The statistical family of the marks.
#' @param pointcovariates Names of the additional, non-spatial covariates describing the points.
#' @param trialspa The name of the trials variable for the presence absence datasets.
#' @param trialsmarks The name of the trials variable for the binomial marks datasets.
#' @param speciesname Name of the species variable used in the data.
#' @param marksspatial Should spatial fields be included for the marks
#' @param spatial Logical argument describing if spatial effects should be included.
#' @param intercepts Logical argument describing if intercepts should be included in the model.
#' @param spatialcovariates Spatial covariates object used in the model.
#' @param marksintercept Logical argument describing if the marks should have interceptes.
#' @param boundary A polygon map of the study area.
#' @param ips Integration points and their respective weights to be used in the model.
#' @param temporal Name of the temporal variable used in the model.

dataSDM$set('public', 'initialize', function(coordinates, projection, Inlamesh, initialnames,
                                             responsecounts, responsepa, 
                                             marksnames, marksfamily, pointcovariates,
                                             trialspa, trialsmarks, speciesname, marksspatial,
                                             spatial, intercepts, spatialcovariates, marksintercepts,
                                             boundary, ips, temporal, temporalmodel, speciesspatial) {
  
  if (missing(coordinates)) stop('Coordinates need to be given.')
  if (missing(projection)) stop('projection needs to be given.')
  if (missing(Inlamesh)) stop('Mesh needs to be given.')

  if (class(Inlamesh) != 'inla.mesh') stop('Mesh needs to be an inla.mesh object.')
  
  if (class(projection)[1] != 'CRS') stop('Projection needs to be a CRS object.')
  
  if (length(coordinates) != 2) stop('Coordinates needs to be a vector of length 2 containing the coordinate names.')
  
  if (coordinates[1] == coordinates[2]) stop('Coordinates need to be unique values.')
  
  if (any(missing(responsecounts), missing(responsepa)) ||
      any(is.null(responsecounts), is.null(responsepa))) stop('At least one of responseCounts and responsePA are NULL. Please provide both.')
  
  private$responseCounts <- responsecounts
  private$responsePA <- responsepa
  
  if (!missing(trialspa)) private$trialsPA <- trialspa
  if (!missing(trialsmarks)) private$trialsMarks <- trialsmarks
  
  if (!missing(speciesname)) private$speciesName <- speciesname
  
  if (!missing(temporal)) private$temporalName <- temporal
  private$temporalModel <- temporalmodel
  
  if (!missing(initialnames)) private$initialnames <- initialnames
  if (!missing(boundary)) private$Boundary <- boundary
  
  private$markNames <- marksnames
  private$markFamily <- marksfamily
  
  private$pointCovariates <- pointcovariates
  
  if (!is.null(spatialcovariates)) private$spatialCovariates(spatialcovariates)
  
  if (!is.null(ips)) private$IPS <- ips
  else {
    
    if (is.null(boundary)) private$IPS <- inlabru::ipoints(samplers = boundary, domain = Inlamesh)
    else private$IPS <- inlabru::ipoints(domain = Inlamesh)
    
  }
  
  private$Spatial <- spatial
  private$marksSpatial <- marksspatial
  private$Intercepts <- intercepts
  private$marksIntercepts <- marksintercepts
  
  private$speciesSpatial <- speciesspatial
  
  #if (!private$Spatial && private$markSpatial) warning('Spatial has been set to FALSE but marksSpatial is TRUE. Spatial effects for the marks will still be run.')
  
  private$Coordinates <- coordinates
  private$Projection <- projection
  private$INLAmesh <- Inlamesh
  invisible(self)
  
})

#' @description Prints the datasets, their data type and the number of observations, as well as the marks and their respective families.

dataSDM$set('public', 'print', function(...) {
  
  if (length(private$modelData) == 0) cat('No data found. Please add data using the `$.addData` function.')
  else {
    cat('Summary of dataSDM data file:\n\n')
    ## Find present absence dataset
    if (length(names(private$printSummary$Type)[private$printSummary$Type == 'Present absence']) > 0) {
      cat('Summary of present absence datasets:\n\n')
      dataIn <- data.frame(c('-----', names(private$printSummary$Type)[private$printSummary$Type == 'Present absence']),
                           c('',rep('|  ---  |', length(private$printSummary$Type[private$printSummary$Type == 'Present absence']))),
                           c('------------------', private$printSummary$numObs[private$printSummary$Type == 'Present absence']))
      names(dataIn) <- c('Name:','','# of observations:')
      print.data.frame(dataIn[,1:3], row.names = FALSE, right = FALSE)
      cat('\n')
      
    }
    
    if (length(names(private$printSummary$Type)[private$printSummary$Type == 'Present only']) > 0) {
      cat('Summary of present only datasets:\n\n')
      dataIn <- data.frame(c('-----', names(private$printSummary$Type)[private$printSummary$Type == 'Present only']),
                           c('',rep('|  ---  |', length(private$printSummary$Type[private$printSummary$Type == 'Present only']))),
                           c('------------------', private$printSummary$numObs[private$printSummary$Type == 'Present only']))
      names(dataIn) <- c('Name:','','# of observations:')
      print.data.frame(dataIn[,1:3], row.names = FALSE, right = FALSE)
      cat('\n')
      
    }
    
    if (length(names(private$printSummary$Type)[private$printSummary$Type == 'Count data']) > 0)  {
      cat('Summary of count datasets:\n\n')
      dataIn <- data.frame(c('-----', names(private$printSummary$Type)[private$printSummary$Type == 'Count data']),
                           c('',rep('|  ---  |', length(private$printSummary$Type[private$printSummary$Type == 'Count data']))),
                           c('------------------', private$printSummary$numObs[private$printSummary$Type == 'Count data']))
      names(dataIn) <- c('Name:','','# of observations:')
      print.data.frame(dataIn[,1:3], row.names = FALSE, right = FALSE)
      cat('\n')
      
    }
    
    if (!is.null(private$markNames)) {
      
      cat('Marks included:\n\n')
      mark_data <- na.omit(data.frame(c('-----',unlist(private$printSummary$Marks)),
                                      c('',rep('|  ---  |', length(unlist(private$printSummary$Marks)))),
                                      c('-----',unlist(private$printSummary$marksType))))
      #mark_data <- mark_data[!duplicated(lapply(x@Mark_data, function(dat) attributes(dat)$mark_name)),]
      mark_data <- unique(mark_data)
      names(mark_data) <- c('Name:','', 'Type:')
      print.data.frame(mark_data[,1:3], row.names = FALSE, right = FALSE)
      cat('\n')
      
    }
  }
  
})

#' @description Makes a plot of the region as well as the points.
#' @param Datasets Name of the datasets to plot.
#' @param Species Should species be plotted as well? Defaults to \code{FALSE}.

dataSDM$set('public', 'plot', function(Datasets, Species = FALSE, ...) {
  
  if (length(private$modelData) == 0) stop('Please provide data before running the plot function.')
  
  if (missing(Datasets)) Datasets <- unique(private$dataSource)
  
  if (!all(Datasets %in% private$dataSource)) stop('Dataset provided not provided to the object.') 
  
  if (Species && is.null(private$speciesName)) stop('speciesName in bruSDM required before plotting species.')
  
  ##Get data
  points <- list()
  
  for (data in Datasets) {
    
    index <- which(private$dataSource == data)
    
    if (!is.null(private$markNames)) index <- index[!endsWith(names(private$modelData[index]), paste0('_', private$markNames))]
    ##if species then dataset_species_response
    # else paste dataset_response ## but also only need point response -- NOT MARKS
    
    #Probably want to create one df object with another variable called dataset placeholder or something
    #Then colour by species or dataset or whatever
    
    #Also need to create a boundary of sorts... either if Boundary is non null; else can make from mesh...
    #Maybe even allow maps if SpatialPolygon provided...

    for (i in 1:length(index)) {
  
    idx <- index[i]  
    points[[data]][[i]] <- private$modelData[[idx]]$data[, names(private$modelData[[idx]]$data) %in% c(private$speciesName, private$responseCounts,
                                                                                                       private$responsePA,'BRU_aggregate')]

   if (!Species) points[[data]][[i]]@data[,'..Dataset_placeholder_var..'] <- rep(data, nrow( points[[data]][[i]]@data))
    
    }
    
  #points[[data]] <- do.call(rbind.SpatialPointsDataFrame, points[[data]])
      
  }

  plotData <- do.call(rbind.SpatialPointsDataFrame, unlist(points))
  ## Make boundary
  stop(return(plotData))
  
  
})

#' @description Function used to add additional, and possibly heterogeneous data to the dataSDM object.
#' @param ... The datasets added into the model: should be a data.frame or a SpatialPoints* object.
#' @param responseCounts The name of the response variable for the counts data.
#' @param responsePA The name of the response variable for the presence absence data.
#' @param trialsPA The name of the trials variable for the presence absence data.
#' @param markNames The names of the marks found in the data.
#' @param markFamily The associated distributions of the marks.
#' @param pointCovariates The additional, non-spatial covariates describing the data.
#' @param trialsMarks The name of the trials variable for the binomial marks.
#' @param speciesName The name of the species variable included in the data.
#' @param Coordinates A vector of length 2 describing the names of the coordinates of the data.
#' @param pointsField An inla.spde model describing the random field for the points.
#' @param speciesField An inla.spde model describing the random field for the species.
#' @param marksField An inla.spde model describing the random field for the marks.

dataSDM$set('public', 'addData', function(..., responseCounts, responsePA, trialsPA,
                                          markNames, markFamily, pointCovariates,
                                          trialsMarks, speciesName, temporalName,
                                          Coordinates, pointsField,
                                          speciesField,
                                          marksField) {
  pointData <- dataOrganize$new()

  dataPoints <- list(...)
  
  if (missing(markNames)) markNames <- private$markNames
  else {
    
    newmarksLens <- length(markNames)
    private$markNames <- markNames <- unique(c(private$markNames, markNames))
    
  }
  
  if (length(dataPoints) == 0) stop('Please provide data in the ... argument.')
  
  datasetClass <- unlist(lapply(dataPoints, class))
  
  if (length(datasetClass) == 1 && datasetClass == "list") {
    
    dataNames <- NULL
    dataPoints <- unlist(dataPoints)
    datasetClass <- lapply(dataPoints, class)
    dataList <- TRUE
    
  }
  else dataList <- FALSE
  
  if (any(!unlist(datasetClass) %in% c("SpatialPointsDataFrame", "SpatialPoints", "data.frame"))) stop("Datasets need to be either a SpatialPoints* object or a data frame.")
  
  if (!is.null(private$initialnames)) dataNames <- private$initialnames
  else
    if (dataList) {
      
      if (is.null(dataNames)) {
        
        dataNames <- setdiff(gsub("list[(]|[)]", "", as.character(match.call(expand.dots = TRUE))), 
                             gsub("list[(]|[)]", "", as.character(match.call(expand.dots = FALSE))))
        dataNames <- unlist(strsplit(x = dataNames, split = ", "))
        
        if (length(dataNames) != length(dataPoints)) {
          
          warning("Issues with naming the datasets from a list. Will create generic dataset names.", 
                  immediate. = FALSE)
          
          cat("\n")
          dataNames <- paste0("dataset_", seq_len(length(dataList)))
          
        }
        
      }
      
    }
  else {
    
    dataNames <- setdiff(as.character(match.call(expand.dots = TRUE)), 
                         as.character(match.call(expand.dots = FALSE)))
  }
  
  if (!missing(pointsField)) {
    
    if (!is.null(pointsField)) self$spatialFields$sharedField <- pointsField #private$pointsField <- pointsField
    else self$spatialFields$sharedField <- INLA::inla.spde2.matern(mesh = private$INLAmesh) #private$pointsField <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
    
  }
  
  if (!is.null(markNames)) {
    
    if (private$marksSpatial) {
      
      self$spatialFields$markFields <- vector(mode = 'list', length = length(markNames))
      names(self$spatialFields$markFields) <- markNames
  
      if (!is.null(marksField)) self$spatialFields$markFields[1:length(self$spatialFields$markFields)] <- list(marksField)
      else self$spatialFields$markFields[1:length(self$spatialFields$markFields)] <- list(INLA::inla.spde2.matern(mesh = private$INLAmesh))
      
    
        
      
    }
  }
  
  if (missing(responseCounts)) responseCounts <- private$responseCounts
  else {
    
    if (responseCounts != private$responseCounts) {
      
      dataPoints <- nameChanger(data = dataPoints, oldName = responseCounts,
                                newName = private$responseCounts)
      responseCounts <-  private$responseCounts
      
    }
    
  }
  
  if (missing(responsePA)) responsePA <- private$responsePA
  else {
    
    if (responsePA != private$responsePA) {
      
      dataPoints <- nameChanger(data = dataPoints, oldName = responsePA,
                                newName = private$responsePA)
      responsePA <-  private$responsePA
    }
    
  }
  
  if (missing(trialsPA)) trialsPA <- private$trialsPA
  else {
    
    if (!is.null(trialsPA) && trialsPA != private$trialsPA) {
      
      dataPoints <- nameChanger(data = dataPoints, oldName = trialsPA,
                                newName = private$trialsPA)
      trialsPA <-  private$trialsPA
    }
    
  }
  
  if (!is.null(private$speciesName)) {
    
    if (missing(speciesName)) speciesName <- private$speciesName
    else {
      
      if (speciesName != private$speciesName) {
        
        dataPoints <- nameChanger(data = dataPoints, oldName = speciesName,
                                  newName = private$speciesName)
        speciesName <-  private$speciesName
        
      }
      
    }
    
  } else speciesName <- NULL
  
  if (!is.null(private$temporalName)) {
    
    if (missing(temporalName)) temporalName <- private$temporalName
    else {
      
      if (temporalName != private$temporalName) {
        
        dataPoints <- nameChanger(data = dataPoints, oldName = temporalName,
                                  newName = private$temporalName)
        temporalName <-  private$temporalName
        
      }
      
    }
    
  } else temporalName <- NULL
  
  
  if (missing(trialsMarks)) trialsMarks <- private$trialsMarks
  else {
    
    if (!is.null(private$trialsMarks) && trialsMarks != private$trialsMarks) {
      
      dataPoints <- nameChanger(data = dataPoints, oldName = trialsMarks,
                                newName = private$trialsMarks)
      trialsMarks <-  private$trialsMarks
    }
    
  }
  
  if (!is.null(markNames)) {
    
    if (!missing(markFamily)) {
      
      markFamily <- private$markFamily <- c(private$markFamily, markFamily)
      
    }
    else
      if (missing(markFamily) && !is.null(private$markFamily)) {
        if (length(markNames) != length(private$markFamily)) {
          
          warning('Mark families not given. Will assume marks as gaussian.')
          markFamily <- private$markFamily <- c(private$markFamily, rep('gaussian', length(newmarksLens)))
          
        } else markFamily <- private$markFamily
      }  
    else 
      if (missing(markFamily) && is.null(private$markFamily)) {
        
        warning('Mark families not given. Will assume marks as gaussian.')
        private$markFamily <- markFamily <- rep('gaussian', length(markNames))
        
      }
    else 
      if (length(markNames) != length(markFamily)) stop('markFamily needs to be the same length as markNames.')
    else private$markFamily <- markFamily <- unique(c(private$markFamily, markFamily))
    
  } else markFamily <- NULL
  
  if (missing(pointCovariates)) pointCovariates <- private$pointCovariates
  else private$pointCovariates <- pointCovariates <- unique(c(private$pointCovariates, pointCovariates))
  
  if (!is.null(private$pointCovariates) && class(pointCovariates) != 'character') stop('pointCovariates is required to be a vector containing the names of the covariates found in the datasets.')
  
  if (!missing(Coordinates) && Coordinates != private$Coordinates) {
    
    dataPoints <- changeCoords(data = dataPoints, oldcoords = Coordinates,
                               newcoords = private$Coordinates)
    
  } 
  else { 
    
    coordsOK <- checkCoords(data = dataPoints,
                            coords = private$Coordinates)
    
    if (!coordsOK) stop('At least one dataset does not have the correct coordinates. Please try again.')  
    
  }
  
  if (!is.null(private$speciesName)) {
    
    speciesOK <- checkVar(data = dataPoints,
                          var = private$speciesName)
    
    if (!speciesOK) stop('The species variable name is required to be present in all the datasets.')
    
  }
  
  if (!is.null(private$temporalName)) {
    
    timeOK <- checkVar(data = dataPoints,
                          var = private$temporalName)
    
    if (!timeOK) stop('The temporal variable name is required to be present in all the datasets.')
    
    #self$spatialFields$temporalField <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
    
  }
  
  pointData$makeData(datapoints = dataPoints, datanames = dataNames,
                     coords = private$Coordinates, proj = private$Projection,
                     countsresp = responseCounts, paresp = responsePA,
                     trialname = trialsPA, speciesname = speciesName,
                     marktrialname = trialsMarks, temporalvar = private$temporalName,
                     marks = markNames, markfamily = markFamily,
                     pointcovnames = pointCovariates)
  
  if (is.null(private$printSummary))  private$printSummary <- list(Type = unlist(pointData$dataType),
                                                                   numObs = unlist(pointData$numObs),
                                                                   Marks = pointData$Marks,
                                                                   marksType = pointData$marksType)
  
  else {
    
    private$printSummary <- list(Type = c(private$printSummary$Type, unlist(pointData$dataType)),
                                 numObs = c(private$printSummary$numObs, unlist(pointData$numObs)),
                                 Marks = c(private$printSummary$Marks, pointData$Marks),
                                 marksType = c(private$printSummary$marksType, pointData$marksType))
    
  }
  
  if (!is.null(speciesName)) {
    
    pointData$makeSpecies(speciesname = speciesName) 
    
    if (is.null(private$speciesIn)) private$speciesIn <- pointData$SpeciesInData
    else private$speciesIn <- c(private$speciesIn, pointData$SpeciesInData)
    
    #ADD argument common field for species
    self$spatialFields$speciesFields <- vector(mode = 'list', length = length(unique(unlist(private$speciesIn))))
    names(self$spatialFields$speciesFields) <- unique(unlist(private$speciesIn))
    
    if (private$speciesSpatial) {
    
    if (!is.null(speciesField)) self$spatialFields$speciesFields[1:length(self$spatialFields$speciesField)] <- list(speciesField)
    else self$spatialFields$speciesFields[1:length(self$spatialFields$speciesField)] <- list(INLA::inla.spde2.matern(mesh = private$INLAmesh))
      
    
    }
  }
  
  if (is.null(private$dataSource)) private$dataSource <- unlist(as.vector(pointData$dataSource))
  else private$dataSource <- c(private$dataSource, unlist(as.vector(pointData$dataSource)))
  
  
  ##Add here that if markSpatial then add mark_spatial
  #Also add markModel in the initial call.
  pointData$makeFormulas(spatcovs = private$spatcovsNames, speciesname = speciesName, temporalname = private$temporalName,
                         paresp = responsePA, countresp = responseCounts, marksspatial = private$marksSpatial,
                         marks = markNames, spatial = private$Spatial, 
                         intercept = private$Intercepts, markintercept = private$marksIntercepts, speciesspatial = private$speciesSpatial)
  
  if (!is.null(private$temporalName)) {
    
    pointData$makeMultinom(multinomVars = private$temporalName,
                           return = 'time', oldVars = NULL)
    
    private$temporalVars <- pointData$timeIndex #??
    
  }
  
  if (is.null(private$multinomVars)) {
    
    if (length(pointData$multinomVars) != 0) private$multinomVars <- multinomNames <- pointData$multinomVars
    else multinomNames <- NULL
    
  } else private$multinomVars <- multinomNames <- unique(c(private$multinomVars, pointData$multinomVars))
  
  if (!is.null(multinomNames)) {
    
    if (length(private$multinomIndex) != 0) oldVars <- private$multinomIndex
    else oldVars <- NULL
    
    pointData$makeMultinom(multinomVars = multinomNames, 
                           return = 'marks',
                           oldVars = oldVars)
    
    ##Will need to create a new function which can update the indexing on multinoVars if new data is added...
    for (i in names(pointData$multinomIndex)) {
      
      if (i %in% names(private$multinomIndex)) private$multinomIndex[[i]] <- c(private$multinomIndex[[i]],unique(unlist(pointData$multinomIndex[[i]])))[!is.na(c(private$multinomIndex[[i]],unique(unlist(pointData$multinomIndex[[i]]))))]
      else private$multinomIndex[[i]] <- unique(unlist(pointData$multinomIndex[[i]]))[!is.na(unique(unlist(pointData$multinomIndex[[i]])))]
      
    }
    
  }
  
  if (is.null(private$Components)) {
    
    #Add if markSpatial then ...
    #Then assign the relevant spatial terms in runModel
    private$Components <- pointData$makeComponents(spatial = private$Spatial, intercepts = private$Intercepts,
                                                   marks = markNames, datanames = dataNames, speciesname = speciesName,
                                                   multinomnames = multinomNames, pointcovariates = pointCovariates,
                                                   covariatenames = private$spatcovsNames, 
                                                   covariateclass = private$spatcovsClass,
                                                   marksspatial = private$marksSpatial,
                                                   marksintercept = private$marksIntercepts,
                                                   temporalname = private$temporalName,
                                                   #speciesspatial = private$speciesField,
                                                   numtime = length(unique(unlist(private$temporalVars))),
                                                   temporalmodel = private$temporalModel,
                                                   speciesspatial = private$speciesSpatial)
    
  }
  else {
    ##Add if species; then remove first species_spatial
    
    if (!is.null(private$speciesName) && private$Spatial) {
      ##Select species_spatial with lowest ngroup...
      ##or maybe just do ngroups again ...
      private$Components <- private$Components[!grepl(paste0('^',private$speciesName,'_spatial'), private$Components)]
      
    }
    
    newComponents <- pointData$makeComponents(spatial = private$Spatial, intercepts = private$Intercepts,
                                              marks = markNames, datanames = dataNames, speciesname = speciesName,
                                              multinomnames = multinomNames, pointcovariates = pointCovariates,
                                              covariatenames = private$spatcovsNames, 
                                              covariateclass = private$spatcovsClass,
                                              marksspatial = private$marksSpatial,
                                              marksintercept = private$marksIntercepts,
                                              #speciesspatial = private$speciesField,
                                              numtime = length(unique(unlist(private$temporalVars))),
                                              temporalmodel = private$temporalModel,
                                              temporalname = private$temporalName,
                                              speciesspatial = private$speciesSpatial)
    
    
    private$Components <- union(private$Components, newComponents)
    
    
  }
  
  
  if (!is.null(private$pointCovariates)) {
    
    datMatrix <- as.data.frame(matrix(0, nrow = nrow(private$IPS@coords), ncol = length(private$pointCovariates)))
    names(datMatrix) <- private$pointCovariates
    private$IPS@data <- cbind(private$IPS@data, datMatrix)
    
  }
  
  pointData$makeLhoods(mesh = private$INLAmesh,
                       ips = private$IPS, paresp = responsePA,
                       ntrialsvar = trialsPA,
                       markstrialsvar = trialsMarks,
                       speciesname = speciesName)
  
  if (!is.null(private$optionsINLA[['control.family']])) {
    
    index1 <- length(private$optionsINLA[['control.family']]) + 1
    index2 <- length(private$optionsINLA[['control.family']]) + length(unlist(pointData$Family))
    
  }
  else {
    
    index1 <- 1
    index2 <- length(pointData$Data)
    
  }
  
  familyIndex <- c(rep(NA, length(private$modelData)), sapply(pointData$Data, function(x) x$family))
  
  for (i in index1:index2) {
    
    if (familyIndex[i] == 'binomial') {
      
      private$optionsINLA[['control.family']][[i]] <- list(link = 'cloglog')
      
    }
    else private$optionsINLA[['control.family']][[i]] <- list(link = 'log')
    
  }
  
  
  if (length(private$modelData) == 0) {
    
    private$modelData <- pointData$Data
    
  }
  else {
    
    private$modelData <- append(private$modelData, pointData$Data)
    
  }
  
})

#' @description Function used to add or change spatial covariates in the model.
#' @param spatialCovariates A SpatialPixelsDataFrame or Raster* object describing covariates at each spatial point.

dataSDM$set('private', 'spatialCovariates', function(spatialCovariates) {
  
  if (missing(spatialCovariates)) stop('Please add spatialCovariates as a Raster* or SpatialPixelsDataFrame object.')
  
  objName <- as.character(match.call())[2]
  
  if (!objName %in% names(parent.frame())) {
    
    if (!objName %in% names(globalenv())) stop('Spatial covariates object not available.')
    else spatcovsEnv <- globalenv()
    
  } 
  else spatcovsEnv <- parent.frame()
  
  if (!class(spatialCovariates) %in% c('RasterLayer', 'RasterBrick',
                                       'RasterStack',
                                       'SpatialPixelsDataFrame')) stop('The spatial Covariates need to be a Raster* object or a SpatialPixelsDataFrame.')
  
  spatcovsIncl <- names(spatialCovariates)
  
  if (!inherits(spatialCovariates, 'Spatial')) {
    #This won't work... will have to convert in runModel
    objSpat <- as(spatialCovariates, 'SpatialPixelsDataFrame')
    covsClass <- sapply(objSpat@data, class)
    
  } else covsClass <- sapply(spatialCovariates@data, class)
  
  if (is.null(private$ptcovsClass))   private$ptcovsClass <- covsClass
  else private$ptcovsClass <- c(private$ptcovsClass, covsClass) #correct? ## maybe even do this by names...
  
  covsClass <- ifelse(covsClass == 'factor', ifelse(private$Intercepts, 'factor_contrast', 'factor_full'), 'linear')
  
  if (length(private$modelData) > 0) {
    
    if (is.null(private$spatcovsObj)) {
      
      #  newForms <- sapply(private$modelData, function(x) {
      
      #    update(x$formula, paste('~ . +', paste(spatcovsIncl, collapse = ' + ')))
      
      #  })
      
      for (form in 1:length(private$modelData)) {
        ##If is.null(private$speciesName) ...
        private$modelData[[form]][['include_components']] <- c(private$modelData[[form]][['include_components']], spatcovsIncl)
        
      }
      
      private$Components <- c(private$Components, newComps)
      
      ## Add all the new spat covs names to formulas ie update func can I sapply it all?
      ## Add all new spat covs to components ... 
    }
    else {
      
      covsKeep <- spatcovsIncl[spatcovsIncl %in% private$spatcovsNames]
      
      if (identical(covsKeep, 'character(0)')) covsKeep <- NULL
      
      covsOut <- private$spatcovsNames[!private$spatcovsNames %in% spatcovsIncl]
      
      if (identical(covsOut, 'character(0)')) covsOut <- NULL
      
      ## Remove all covs if covsOut not NULL
      ## Add all new covs if covsKeep not NULL
      ## Do same for components ...
      
    }
  }
  
  private$spatcovsObj <- objName
  private$spatcovsNames <- spatcovsIncl
  private$spatcovsEnv <- spatcovsEnv
  private$spatcovsClass <- covsClass
  
})

#' @description Function used to add additional bias fields to different data processes.
#' @param datasetNames A vector of datasets to add the bias fields to.
#' @param allPO Should the bias fields be added to all the presence only datasets. Defaults to \code{NULL}.
#' @param biasField An inla.spde model descrbing the random bias field.


## if we are doing spatial fields in a public list, then don't need biasField?
dataSDM$set('public', 'addBias', function(datasetNames = NULL,
                                          allPO = FALSE,
                                          biasField = NULL) {
  
  if (allPO) datasetNames <- names(private$printSummary)[private$printSummary == 'Present Only']
  else
    if (is.null(datasetNames)) stop('Dataset names need to be given.')
  
  if (!all(datasetNames %in% private$dataSource)) stop('Dataset provided not available.')
  
  for (dat in datasetNames) {
    
    index <- which(private$dataSource == dat)
    
    for (lik in index) {
      
      #private$modelData[[lik]]$formula <- update(private$modelData[[lik]]$formula, paste0(' ~ . + ', dat,'_bias_field'))
      private$modelData[[lik]]$include_components <- c(private$modelData[[lik]]$include_components, paste0(dat, '_biasField'))
    }
    
    if (is.null(biasField)) self$spatialFields$biasFields[[dat]] <- inla.spde2.matern(mesh = private$INLAmesh)
    else self$spatialFields$biasFields[[dat]] <- biasField
    
  }
  
  
  
  ##Add a way for separate Fields per dataset.
  #Shouldn't be too difficult -- have already done it before.
  
  
  #Should I copy the bias fields for the marks?
  private$Components <- c(private$Components, paste0(datasetNames ,'_biasField(main = coordinates, model = ', datasetNames, '_bias_field)'))
  ##Things to do here:
  #Go into the liks of PO datasets and add the biasfield
  #Go inth the components and add the bias field component
  #Then store the bias field somewhere else. If field NULL, then add generic bias field
  
  
})

#' @description Function to update formulas for a given process.
#' @param datasetName Name of the dataset to change the formula for.
#' @param speciesName Name of the species to change the formula for.
#' @param markName Name of the mark to change the formula for.
#' @param Formula An updated formula to give to the process.
#' @param allDataset Logical argument: if \code{TRUE} changes the formulas for all processes in a dataset.
#' @param keepSpatial Logical argument: should the spatial effects remain in the formula. Defaults to \code{TRUE}.
#' @param keepIntercepts Logical argument: should the intercepts remain in the formula. Defaults to \code{TRUE}.
#' @param newFormula Completely change the formula for a process. Note: all terms need to be correctly specified here. ## TO ADD

dataSDM$set('public', 'updateFormula', function(datasetName = NULL, speciesName = NULL,
                                                markName = NULL, Formula, allDataset = FALSE,
                                                keepSpatial = TRUE, keepIntercepts = TRUE,
                                                newFormula, ...) {
  ## Will need to update this such that if keepSpatial & species then paste0(species,_spatial)
  
  if (all(is.null(datasetName), is.null(speciesName), is.null(markName))) stop ('At least one of: datasetName, speciesName or markName needs to be specified.')
  
  if (!is.null(speciesName) && is.null(private$speciesName)) stop ('Species are given but none are present in the model. Please specify species in the model with "speciesName" in bruSDM.')
  
  if (is.null(speciesName) && !is.null(private$speciesName)) speciesName <- unlist(private$speciesIn[datasetName])
  
  if (!datasetName %in% private$dataSource) stop ('Dataset name provided not in model.')
  
  if (!missing(Formula) && !missing(newFormula)) stop ('Please provide only one of Formula and newFormula. \n Use Formula to update the current formula; use newFormula to make a completely new formula.')
  
  if (!missing(Formula) && class(Formula) != 'formula') stop ('Formula must be of class "formula".')
  
  if (!missing(newFormula) && class(newFormula) != 'formula') stop('newFormula must be of class "formula".')
  
  if (!missing(Formula) && length(as.character(Formula)) == 3) stop ("Please remove the response variable of the formula.")
  
  if (!missing(newFormula) && length(as.character(newFormula)) == 3) stop ("Please remove the response variable of the formula.")
  
  if (allDataset && is.null(datasetName)) stop ('Please provide a dataset name in conjunction with allDataset.')
  
  if (!is.null(markName)) {
    
    if (!markName %in% private$markNames) stop ('Mark provided not in model.')
    
  }
  
  if (keepSpatial) {
    
    if (!private$Spatial && is.null(private$markNames)) keepSpatial <- FALSE
    else
      if (!private$Spatial && !private$marksSpatial) keepSpatial <- FALSE
      
  }
  
  if (!is.null(speciesName) && !is.null(markName)) {
    
    speciesName <- list()
    
    for (dataset in datasetName) {
      
      if (!is.na(private$printSummary$Marks[dataset])) {
        
        if (any(markName %in% private$printSummary$Marks[dataset])) speciesName[dataset] <- rep(private$speciesIn[datasetName], each = length(sum(markName %in% private$printSummary$Marks[dataset])))
        else speciesName[dataset] <- private$speciesIn[dataset]
        
      } else speciesName[dataset] <- private$speciesIn[dataset]
      
    } 
    
    speciesInd <- unlist(speciesName)  
    
  } else speciesInd <- speciesName
  
  if (allDataset) name_index <- names(private$modelData)[startsWith(private$modelData, paste0(datasetName, '_'))]
  
  else{
    
    if (!is.null(markName)) process_index <- paste0('_', markName)
    else process_index <- paste0('_', c('coordinates', private$responsePA, private$responseCounts))
    
    if (!is.null(private$speciesName)) {
      
      if (!is.null(speciesName)) species_index <- paste0('_', speciesName)
      else species_index <- NULL
      
      name_index <- apply(expand.grid(datasetName, species_index), MARGIN = 1, FUN = paste0,collapse='')
      name_index <- apply(expand.grid(name_index, process_index), MARGIN = 1, FUN = paste0,collapse='')
      
      name_index <- name_index[name_index %in% names(private$modelData)]
      
      if (identical(name_index, character(0))) stop('Species name provided not in dataset given.')
      
    }
    else {
      
      name_index <-  apply(expand.grid(datasetName, process_index), MARGIN = 1, FUN = paste0,collapse='')  
      name_index <- name_index[name_index %in% names(private$modelData)]
      
    }
    
  }
  
  if (missing(Formula) && missing(newFormula)) {
    
    get_formulas <- lapply(private$modelData[name_index], function(x) list(formula = x$formula,
                                                                           components = x$include_components))
    
    get_formulas <- lapply(get_formulas, function(x) {
      
      if (as.character(x$formula)[3] == '.') update(x$formula, paste('~', paste0(x$components, collapse = '+'))) 
      else x$formula
      
    })
    ##Does this return the names of the formulas??
    return(get_formulas)
    
  }
  else
    if (!missing(Formula)) {
    
    #if (length(as.character(Formula)) == 2 && as.character(Formula)[2] == '.') formula_terms <- c()
    #else formula_terms <- attributes(terms(formula))[['term.labels']]
    
    index_species <- 0
    
    
    for (dataset in name_index) {
      
      index_species <- index_species + 1
      
      if (!is.null(markName)) {
        
        if (!is.null(private$speciesName)) {
          
          if (dataset == paste0(datasetName, '_', markName, '_', speciesInd[index_species])) mark_p <- TRUE
          else mark_p <- FALSE
          
        }
        else {
          
          if (dataset == paste0(datasetName, '_', markName)) mark_p <- TRUE
          else mark_p <- FALSE
        } 
      }
      else mark_p <- FALSE
      
      formula_update <- Formula ## This needs to be the formula
      
      if (!keepSpatial) {
        ##Here update formula_update to - shared_spatial or mark_spatial or whatever
        if (mark_p) {
          
          if (private$marksSpatial) formula_update <- update(formula_update, formula(paste('~ . -', paste0(markName, '_spatial'))))
          
        }
        else {
          
          if (private$Spatial) formula_update <- update(formula_update, formula(~ . - shared_spatial))
          
        }
        
      }
      
      if (!keepIntercepts) {
        ##Here update formula update to - intercepts of whatever
        if (mark_p) {
          
          if (private$marksIntercepts) formula_update <- update(formula_update, formula(paste('~ . -', paste0(markName, '_intercept'))))
          
        }
        
        else {
          
          if (!is.null(private$speciesName)) formula_update <- update(formula_update, formula(paste('~ . -', paste0(speciesInd[index_species], '_intercept'))))
          else formula_update <-update(formula_update, formula(paste('~ . -', paste0(datasetName, '_intercept'))))
          
          
        }
        
      }
      
      if (!is.null(speciesName)) {
        
        covs_in <- all.vars(Formula)[all.vars(Formula) %in% private$spatcovsNames]
        
        if (!identical(covs_in, character(0))) {
          #Change this to an update formula side
          # but this is more difficult since what if there is a plus?
          # need to do a strsplit; get all the covariate terms; and change them to species_cov
          char_formula <- unlist(strsplit(as.character(Formula), split = ' '))
          
          char_formula[char_formula == covs_in] <- paste0(speciesInd[index_species], '_', covs_in)
          
          formula_update <- formula(paste(char_formula, collapse = ' '))
          
        }
        
      }
      
      updated_formula <- update(formula(paste('~ ', paste0(private$modelData[[dataset]]$include_components, collapse = ' + '))), formula_update)
      updated_formula <- all.vars(updated_formula)[all.vars(updated_formula) != '.']
      ##This will be removed once we move the like construction to runModel.
      private$modelData[[dataset]]$include_components <- updated_formula
      
    }
    
    }
  else {
    
    for (dataset in name_index) { 
      
      
      ## Get old formula:
      cat('Old formula for', paste0(dataset, ': '))
      if (length(all.vars(private$modelData[[dataset]]$formula)) == 2) {
        
        oldForm <- update(private$modelData[[dataset]]$formula, paste(' ~ ', paste(private$modelData[[dataset]]$include_components, collapse = ' + ')))
        
        print(oldForm)
        cat('New formula: ')
        ## Maybe I should do a check: if cov in newFormula then paste0(species, _ , covariate) ## how else does it work within a for loop
        newForm <- update(private$modelData[[dataset]]$formula, newFormula)
        print(newFormula)
        cat('\n')
        
      }
      else {
        
        print(private$modelData[[dataset]]$formula)
        cat('\n')
        cat('New formula: ')
        newForm <- update(paste(as.character(private$modelData[[dataset]]$formula), '~ '), newFormula)
      }
      
      private$modelData[[dataset]]$formula <- newForm
      private$modelData[[dataset]]$include_components <- c()
      
    }
      

  }
  
})

#' @description Function to add custom components to the integrated modeling.
#' @param component Component to add to the integrated model.
#' @param datasetName Names of the dataset to add the component to.
#' @param speciesName Name of the species to add the component to.
#' @param markName Name of the mark to add the component to.
#' @param allDatasets Add the component to all the datasets.
#' @param ... Any additional objects associated with the component, such as an inla.mesh.
#' 

dataSDM$set('public', 'changeComponents', function(addComponent, removeComponent) {
  
  terms <- gsub('\\(.*$', '', private$Components)
  
  if (!missing(addComponent)) {
   
   if (gsub('\\(.*$', '', addComponent) %in% terms) private$Components <- private$Components[! terms %in% gsub('\\(.*$', '', addComponent)]
   
   private$Components <- c(private$Components, addComponent)
    
  } 
  
  if (!missing(removeComponent)) private$Components <- private$Components[!terms%in%removeComponent]
    
  
  
  cat('Components:')
  cat('\n')
  componentsJoint <- formula(paste('~ - 1 +', paste(private$Components, collapse = ' + ')))
  componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
  
  print(componentsJoint)
  
  
})

#' @description Function used to account for preferential sampling in the modeling framework.
#' 
#' 
#' 

dataSDM$set('public', 'samplingBias', function(...) {
  
  stop('Do later...')
  
  # Things to consider::
   # Do we treat every point as a sampling location???
   # If all we have are points of species 
    # ie reflecting only the location of the species
    # then we are essentially just duplicating the data?
     # Can we just use "expert maps"
      # ie use PA data to infer where the sampling locations are?  
})

#' @description Function to spatially block the datasets
#' @param k Number of cross-validation folds.
#' @param rows Integer value by which the area is divided into latitudinal bins.
#' @param cols Integer value by which the area is divided into longitudinal bins.
#' @param plot Plot the cross-validation folds. Defaults to \code{FALSE}.
#' 
#' 
dataSDM$set('public', 'spatialBlock', function(k, rows, cols, plot = FALSE, ...) {
  
  ##Replace datasets
  blocks <- blockCV::spatialBlock(speciesData = do.call(rbind.SpatialPointsDataFrame, lapply(private$modelData, function(x) x$data)),
                                  k = k, rows = rows, cols = cols, selection = 'random',
                                  verbose = FALSE)

  folds <- blocks$blocks$folds
  
  blocksPoly <- list(sapply(1:(rows * cols), function(s) SpatialPolygons(blocks$blocks@polygons[s], proj4string = private$Projection)))
  
  blocked_data <- list()
  in_where <- list()
  
  ##Integration points?
  
  for (data in names(private$modelData)) {
    
    in_where[[data]] <- lapply(1:(rows * cols), function(i) !is.na(over(private$modelData[[data]]$data, blocksPoly[[1]][[i]])))
    
    for (i in 1:(rows * cols)) {
      
      blocked_data[[data]][[i]] <- private$modelData[[data]]$data[in_where[[data]][[i]], ]
      blocked_data[[data]][[i]]$block_index <- as.character(folds[i])
      
    }
    
    private$modelData[[data]]$data <- do.call(rbind.SpatialPointsDataFrame, blocked_data[[data]])
    
  }
  
  if (plot) {
    
    # if boundary is.null
    
    loc <- private$INLAmesh$loc
    segm <- private$INLAmesh$segm$bnd
    
    coords <- na.omit(data.frame(loc[t(cbind(segm$idx[,, drop=FALSE], NA)), 1],
                                 loc[t(cbind(segm$idx[,, drop=FALSE], NA)), 2]))

    Polys <- Polygon(coords = coords)
    Polys <- Polygons(srl = list(Polys), ID = 'id')
    SpatPolys <- SpatialPolygons(list(Polys), proj4string = private$Projection)
    
    all_data <- do.call(rbind.SpatialPointsDataFrame, lapply(private$modelData, function(x) x$data))
    
    blocks$plots + gg(all_data, aes(col = block_index)) + gg(SpatPolys)
    ## Need block block$plots + gg(species_data) + gg(boundary), which we need to get from the mesh
    
    
  }
  
})

## Need to change all the spatialFields to self$spatialFields and then the relevent sublist?
dataSDM$set('public', 'spatialFields', list(sharedField = list(),
                                            speciesFields = list(),
                                            markFields = list(),
                                            biasFields = list()))

