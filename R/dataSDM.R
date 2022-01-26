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
dataSDM$set('private', 'spatialModel', NULL)
dataSDM$set('private', 'markNames', NULL)
dataSDM$set('private', 'pointCovariates', NULL)
dataSDM$set('private', 'ptcovsClass', NULL)
dataSDM$set('private', 'initialnames', NULL)

dataSDM$set('private', 'modelData', list())
dataSDM$set('private', 'pointsField', list())
dataSDM$set('private', 'speciesField', list())
dataSDM$set('private', 'biasField', NULL)
dataSDM$set('private', 'marksField', list())

dataSDM$set('private', 'spatcovsObj', NULL)
dataSDM$set('private', 'spatcovsNames', NULL)
dataSDM$set('private', 'spatcovsEnv', NULL)
dataSDM$set('private', 'spatcovsClass', NULL)
dataSDM$set('private', 'dataSource', NULL)

dataSDM$set('private', 'Spatial', TRUE)
dataSDM$set('private', 'Intercepts', TRUE)
dataSDM$set('private', 'IPS', NULL)
dataSDM$set('private', 'multinomVars', NULL)
dataSDM$set('private', 'printSummary', NULL)
dataSDM$set('private', 'multinomIndex', list())
dataSDM$set('private', 'optionsINLA', list())

#' @description Initialize function for dataSDM: used to store some compulsoury arguments.
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
#' @param spatial Logical argument describing if spatial effects should be included.
#' @param intercepts Logical argument describing if intercepts should be included in the model.
#' @param spatialcovariates Spatial covariates object used in the model.
#' @param boundary A polygon map of the study area.
#' @param ips Integration points and their respective weights to be used in the model.

dataSDM$set('public', 'initialize', function(coordinates, projection, Inlamesh, initialnames,
                                             responsecounts, responsepa, 
                                             marksnames, marksfamily, pointcovariates,
                                             trialspa, trialsmarks, speciesname,
                                             spatial, intercepts, spatialcovariates,
                                             boundary, ips) {
  
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
  
  if (!missing(initialnames)) private$initialnames <- initialnames
  if (!missing(boundary)) private$Boundary <- boundary

  private$markNames <- marksnames
  private$markFamily <- marksfamily

  private$pointCovariates <- pointcovariates

  if (!is.null(spatialcovariates)) self$spatialCovariates(spatialcovariates)

  if (!is.null(ips)) private$IPS <- ips
  else {
    
    if (is.null(boundary)) private$IPS <- inlabru::ipoints(samplers = boundary, domain = Inlamesh)
    else private$IPS <- inlabru::ipoints(domain = Inlamesh)
    
  }
  
  private$Spatial <- spatial
  private$Intercepts <- intercepts

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

#' @description Makes a plot of the region as well as the points. (DO FOR BOTH DATASETS AND SPECIES.)

dataSDM$set('public', 'plot', function(...) {
  
  
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
                                          trialsMarks, speciesName,
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
  
  if (!is.null(private$initialnames) && length(private$modelData) == 0) dataNames <- private$initialnames
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
          dataNames <- paste0("dataset_", seq_len(length(datasets)))
          
        }
        
      }
      
    }
  else {
    
    dataNames <- setdiff(as.character(match.call(expand.dots = TRUE)), 
                         as.character(match.call(expand.dots = FALSE)))
  }
  
  if (!missing(pointsField)) {
    
    if (!is.null(pointsField)) private$pointsField <- pointsField
    else private$pointsField <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
    
  }
  
  if (!is.null(private$speciesName)) {
    
    if (!missing(speciesField)) {
      
      if (!is.null(speciesField)) private$speciesField <- speciesField
      else private$speciesField <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
      
    }
  }
  
  if (!is.null(markNames)) {
    
    if (!missing(marksField)) {
      
      if (!is.null(marksField)) private$marksField <- marksField
      else private$marksField <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
      
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
  
  pointData$makeData(datapoints = dataPoints, datanames = dataNames,
                     coords = private$Coordinates, proj = private$Projection,
                     countsresp = responseCounts, paresp = responsePA,
                     trialname = trialsPA, speciesname = speciesName,
                     marktrialname = trialsMarks,
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
    
    pointData$makeSpecies(speciesname = speciesName) ## update here is wll
    
    if (is.null(private$speciesIn)) private$speciesIn <- pointData$SpeciesInData
    else private$speciesIn <- c(private$speciesIn, pointData$SpeciesInData) 
    
  }
  
  if (is.null(private$dataSource)) private$dataSource <- unlist(pointData$dataSource)
  else private$dataSource <- c(private$dataSource, unlist(pointData$dataSource))
  
  
  ##Add here that if markSpatial then add mark_spatial
  #Also add markModel in the initial call.
  pointData$makeFormulas(spatcovs = private$spatcovsNames, speciesname = speciesName,
                         paresp = responsePA, countresp = responseCounts,
                         marks = markNames, spatial = private$Spatial, intercept = private$Intercepts)
  
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
                                                   #speciesspatial = private$speciesField,
                                                   numspecies = length(unique(unlist(private$speciesIn))))
    
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
                                              #speciesspatial = private$speciesField,
                                              numspecies = length(unique(unlist(private$speciesIn))))
    
    
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

dataSDM$set('public', 'spatialCovariates', function(spatialCovariates) {
  
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
        
        private$modelData[[form]][['include_components']] <- c(private$modelData[[form]][['include']], spatcovsIncl)
        
      }
      
      if (!is.null(private$speciesName)) {
        
        speciesCovs <- apply(expand.grid(paste0(unique(unlist(private$speciesIn)),'_'), spatcovsIncl), MARGIN = 1, FUN = paste0,collapse='')
        newComps <- paste0(speciesCovs, '(main = ', spatcovsIncl, ', model = \"',covsClass,'\")', collapse = ' + ')
        
      }
      else newComps <- paste0(spatcovsIncl, '(main = ', spatcovsIncl, ', model = \"',covsClass,'\")', collapse = ' + ')
      
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

dataSDM$set('public', 'addBias', function(datasetNames = NULL,
                                          allPO = FALSE,
                                          biasField = NULL) {
  
  if (allPO) datasetNames <- names(private$printSummary)[private$printSummary == 'Present Only']
  else
    if (is.null(datasetNames)) stop('Dataset names need to be given.')
  
  for (dat in datasetNames) {
    
    index <- which(private$dataSource == dat)
    
    for (lik in index) {
      
      #private$modelData[[lik]]$formula <- update(private$modelData[[lik]]$formula, paste0(' ~ . + ', dat,'_bias_field'))
      private$modelData[[lik]]$include_components <- c(private$modelData[[lik]]$include_components, paste0(dat, '_bias_field'))
    }
    
  }
  
  if (is.null(biasField)) private$biasField <- inla.spde2.matern(mesh = private$INLAmesh)
  else private$biasField <- biasField
  
  ##Add a way for separate Fields per dataset.
  #Shouldn't be too difficult -- have already done it before.
  
  
  #Should I copy the bias fields for the marks?
  private$Components <- c(private$Components, paste0(datasetNames ,'_bias_field(main = coordinates, model = biasField)'))
  ##Things to do here:
  #Go into the liks of PO datasets and add the biasfield
  #Go inth the components and add the bias field component
  #Then store the bias field somewhere else. If field NULL, then add generic bias field
  
  
})

#' @description Function to add INLA options to the model.
#' @param ... A list of INLA options.

#Obsolete function? Maybe just add options in runModel
#dataSDM$set('public', 'addOptions', function(...) {
#  
#  private$optionsINLA <- append(private$optionsINLA, ...)
#  
#})


dataSDM$set('public', 'speciesFormula', function(...) {
  
  
})

dataSDM$set('public', 'addComponents', function(...) {
  
  
})


