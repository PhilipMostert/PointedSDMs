#' @title R6 class for creating a dataSDM object to be used in runModel
#' @description A data object containing the data and the relevant information about the integrated model.
#' @export
#' @importFrom R6 R6Class
#' 

dataSDM <- R6::R6Class(classname = 'dataSDM', lock_objects = FALSE, cloneable = FALSE, public = list(
  
  
  #' @description Initialize function for dataSDM: used to store some compulsory arguments. Please refer to the wrapper function, \code{bruSDM} for creating new dataSDM objects.
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
  
  initialize = function(coordinates, projection, Inlamesh, initialnames,
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
  }
  ,
  #' @description Prints the datasets, their data type and the number of observations, as well as the marks and their respective families.
  
  print = function(...) {
    
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
    
  }
  ,
  #' @description Makes a plot of the region as well as the points.
  #' @param Datasets Name of the datasets to plot.
  #' @param Species Should species be plotted as well? Defaults to \code{FALSE}.
  #' @return A ggplot object.
  
  plot = function(Datasets, Species = FALSE, ...) {
    
    if (length(private$modelData) == 0) stop('Please provide data before running the plot function.')
    
    if (missing(Datasets)) Datasets <- unique(private$dataSource)
    
    if (!all(Datasets %in% private$dataSource)) stop('Dataset provided not provided to the object.') 
    
    if (Species && is.null(private$speciesName)) stop('speciesName in bruSDM required before plotting species.')
    
    ##Get data
    points <- list()
    
    for (data in Datasets) {
      
      index <- which(private$dataSource == data)
      
      if (!is.null(private$markNames)) {
        
       if (!is.null(private$speciesName)) index <- unique(private$speciesIn[[data]])
       else index <- 1
      #index <- index[!endsWith(names(private$modelData[index]), paste0('_', private$markNames))]
        
      } 

      for (i in 1:length(index)) {
        
        points[[data]][[i]] <- private$modelData[[data]][[i]][, names(private$modelData[[data]][[i]]) %in% c(private$speciesName, private$responseCounts,
                                                                                                                   private$responsePA,'BRU_aggregate')]

        if ('BRU_aggregate' %in% names(points[[data]][[i]])) points[[data]][[i]] <- points[[data]][[i]][points[[data]][[i]]$BRU_aggregate,]
   
        if (!Species) points[[data]][[i]]@data[,'..Dataset_placeholder_var..'] <- rep(data, nrow(points[[data]][[i]]@data))
       
        
      }
      
      #points[[data]] <- do.call(rbind.SpatialPointsDataFrame, points[[data]])
      
    }
    
    plotData <- do.call(rbind.SpatialPointsDataFrame, lapply(unlist(points), function(x) x[, names(x) %in% c('..Dataset_placeholder_var..', private$speciesName)]))
    
    bound <- private$polyfromMesh()
    
    if (Species) {
      
      ## Need to add the species in here
      
      plotData@data[, private$speciesName] <- unlist(lapply(1:length(plotData[,private$speciesName]), function(x,y,z) y[z[x]], y = unlist(private$speciesIn), z= plotData@data[,private$speciesName]))
      
      colOption <- gg(plotData, aes(col = eval(parse(text = private$speciesName))))
      
      ggplot() + colOption + gg(bound) + guides(col = guide_legend(title = 'Species Name')) 
      
    }
    else { 
      
      colOption <- gg(plotData, aes(col = eval(parse(text = '..Dataset_placeholder_var..'))))
      ggplot() + colOption + gg(bound) + guides(col = guide_legend(title = 'Dataset Name')) 
      
    }
    
    
    
    
  }
  ,
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

  addData = function(..., responseCounts, responsePA, trialsPA,
                     markNames, markFamily, pointCovariates,
                     trialsMarks, speciesName, temporalName,
                     Coordinates) {
    pointData <- dataOrganize$new()
    
    dataPoints <- list(...)
    
    if (missing(markNames)) markNames <- private$markNames # do same for species
    else {
      
      newmarksLens <- length(markNames)
      private$markNames <- markNames <- unique(c(private$markNames, markNames))
      
    }
    
    if (!is.null(private$speciesName)) {
      #Test
      speciesOld <- unique(unlist(private$speciesIn))
      
      speciesNew <- unique(unlist(lapply(dataPoints, function(x) {
        
        if (inherits(x, 'Spatial')) x@data[,private$speciesName]
        else x[,private$speciesName]
        
        
      })))
      
      speciesIn <- c(speciesOld, speciesNew)
      
      
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
    
    if (private$Spatial) {
      
      if (is.null(self$spatialFields$sharedField[['sharedField']])) self$spatialFields$sharedField[['sharedField']] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
      
    }
    
    if (!is.null(markNames)) {
      
      if (private$marksSpatial) {
        #re do this such that only add new marks to self$spatialFields$markFields
        if (!all(markNames %in% names(self$spatialFields$markFields))) {
        
        new_marks <- vector(mode = 'list', length = sum(!markNames %in% names(self$spatialFields$markFields)))
        names(new_marks) <- markNames[!markNames %in% names(self$spatialFields$markFields)]
        
        self$spatialFields$markFields <- append(self$spatialFields$markFields, new_marks)
        #self$spatialFields$markFields <- vector(mode = 'list', length = length(markNames))
        #names(self$spatialFields$markFields) <- markNames
        
        }
        ## re do this such that if is non null, don't touch
        if (any(unlist(lapply(self$spatialFields$markFields, is.null)))) {
          
          for (mark in names(self$spatialFields$markFields)) {
            
            if (is.null(self$spatialFields$markFields[[mark]])) self$spatialFields$markFields[[mark]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
            
          }
          
        }
        
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
        
        if (!all(speciesIn %in% names(self$spatialFields$speciesFields))) {
          ##Not species Name
          new_species <- vector(mode = 'list', length = sum(!speciesIn %in% names(self$spatialFields$speciesFields)))
          names(new_species) <- speciesIn[!speciesIn %in% names(self$spatialFields$speciesFields)]
          
          self$spatialFields$speciesFields <- append(self$spatialFields$speciesFields, new_species)
          #self$spatialFields$markFields <- vector(mode = 'list', length = length(markNames))
          #names(self$spatialFields$markFields) <- markNames
          
        }
        
        if (any(unlist(lapply(self$spatialFields$speciesFields, is.null)))) {
          
          for (species in names(self$spatialFields$speciesFields)) {
            
            if (is.null(self$spatialFields$speciesFields[[species]])) self$spatialFields$speciesFields[[species]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
            
          }
          
        }
        
        
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
    
    private$Formulas <- pointData$Formulas
    
    #pointData$makeLhoods(mesh = private$INLAmesh,
    #                     ips = private$IPS, paresp = responsePA,
    #                     ntrialsvar = trialsPA,
    #                     markstrialsvar = trialsMarks,
    #                     speciesname = speciesName)
    
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
    
  }
  ,
  #' @description Function used to add additional bias fields to different data processes.
  #' @param datasetNames A vector of datasets to add the bias fields to.
  #' @param allPO Should the bias fields be added to all the presence only datasets. Defaults to \code{NULL}.
  #' @param biasField An inla.spde model descrbing the random bias field.
  
  addBias = function(datasetNames = NULL,
                     allPO = FALSE,
                     biasField = NULL) {
    
    if (allPO) datasetNames <- names(private$printSummary)[private$printSummary == 'Present Only']
    else
      if (is.null(datasetNames)) stop('Dataset names need to be given.')
    
    if (!all(datasetNames %in% private$dataSource)) stop('Dataset provided not available.')
    
    for (dat in datasetNames) {
      
      #index <- which(private$dataSource == dat)
      
      for (lik in 1:length(private$Formulas[[dat]])) {
        
        #private$modelData[[lik]]$formula <- update(private$modelData[[lik]]$formula, paste0(' ~ . + ', dat,'_bias_field'))
        private$Formulas[[dat]][[lik]][[1]]$RHS <- c(private$Formulas[[dat]][[lik]][[1]]$RHS, paste0(dat, '_biasField'))
        #private$modelData[[lik]]$include_components <- c(private$modelData[[lik]]$include_components, paste0(dat, '_biasField'))
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
    
    
  }
  ,
  #' @description Function to update formulas for a given process.
  #' @param datasetName Name of the dataset to change the formula for.
  #' @param Points Logical: should the formula be changed for the points. Defaults to \code{TRUE}. If \code{FALSE}, then markNames needs to be non-null.
  #' @param speciesName Name of the species to change the formula for.
  #' @param markName Name of the mark to change the formula for.
  #' @param Formula An updated formula to give to the process.
  #' @param allProcesses Logical argument: if \code{TRUE} changes the formulas for all processes in a dataset.
  #' @param newFormula Completely change the formula for a process -- primarily used to add non-linear components into the formula. Note: all terms need to be correctly specified here.
  
  updateFormula = function(datasetName = NULL, Points = TRUE, speciesName = NULL,
                           markName = NULL, Formula, allProcesses = FALSE,
                           newFormula, ...) {
    ## Will need to update this such that if keepSpatial & species then paste0(species,_spatial)
    
    if (all(is.null(datasetName), is.null(speciesName), is.null(markName))) stop ('At least one of: datasetName, speciesName, markName or allProcesses needs to be specified.')
    
    if (!is.null(speciesName) && is.null(private$speciesName)) stop ('Species are given but none are present in the model. Please specify species in the model with "speciesName" in bruSDM.')
    
    if (is.null(speciesName) && !is.null(private$speciesName)) speciesName <- unlist(private$speciesIn[datasetName])
    
    if (!Points && is.null(markName)) stop ('markNames cannot be non-null if Points is FALSE.')
    
    if (length(datasetName) != 1) stop ('Please only provide one dataset name.')
    
    if (!datasetName %in% private$dataSource) stop ('Dataset name provided not in model.')
    
    if (!missing(Formula) && !missing(newFormula)) stop ('Please provide only one of Formula and newFormula. \n Use Formula to update the current formula; use newFormula to make a completely new formula.')
    
    if (!missing(Formula) && class(Formula) != 'formula') stop ('Formula must be of class "formula".')
    
    if (!missing(newFormula) && class(newFormula) != 'formula') stop('newFormula must be of class "formula".')
    
    if (!missing(Formula) && length(as.character(Formula)) == 3) stop ("Please remove the response variable of the formula.")
    
    if (allProcesses && is.null(datasetName)) stop ('Please provide a dataset name in conjunction with allProcesses')
    
    if (!is.null(markName)) {
      
      if (!markName %in% private$markNames) stop ('Mark provided not in model.')
      
    }
    
  #  if (!is.null(speciesName) && !is.null(markName)) {
      
    #  speciesName <- list()
      
    #  for (dataset in datasetName) {
        
    #    if (!is.na(private$printSummary$Marks[dataset])) {
          
    #      if (any(markName %in% private$printSummary$Marks[dataset])) speciesName[dataset] <- rep(private$speciesIn[datasetName], each = length(sum(markName %in% private$printSummary$Marks[dataset])))
    #      else speciesName[dataset] <- private$speciesIn[dataset]
          
    #    } else speciesName[dataset] <- private$speciesIn[dataset]
        
    #  } 
      
    #  speciesInd <- unlist(speciesName)  
      
    #} else speciesInd <- speciesName
    
    if (allProcesses) {
      
      name_index <- names(private$Formulas[[datasetName]]) #name index should be species or dataset name
      index2 <- seq_along(private$Formulas[[datasetName]][[1]]) ## Since they should all be the same ...
      
    }
    
    else {
      
    if (!is.null(private$speciesName)) {
      
      if (!is.null(speciesName))  {

        if (!all(speciesName %in% names(private$Formulas[[datasetName]]))) stop('Species provided not available in the dataset.') 
        else name_index <- speciesName
      
      }
      else name_index <- names(private$Formulas[[datasetName]])
      
      
      
    } else name_index <- datasetName
    
    if (Points) index2 <- 1
    else index2 <- NULL
    
    if (!is.null(markName)) {
      
      if (!markName %in% names(private$Formulas[[datasetName]][[name_index[1]]])) stop('markName not provided in datasetName.')
      else index2 <- c(index2, which(names(private$Formulas[[datasetName]][[name_index[1]]]) %in% markName)) #Since they should all be the same ...
      
      
    }


    }
    
    if (missing(Formula) && missing(newFormula)) {
    
      get_formulas <- list()
      for (i in index2) {
        
      get_formulas[[i]] <- lapply(private$Formulas[[datasetName]][name_index], function(x) list(formula = x[[i]]$LHS,
                                                                             components = x[[i]]$RHS))
      
      }

      get_formulas <- lapply(unlist(get_formulas, recursive = FALSE), function(x) {
        
        if (as.character(x$formula)[3] == '.') update(x$formula, paste('~', paste0(x$components, collapse = '+'))) 
        else x$formula
        
      })
      ##Does this return the names of the formulas??
      return(get_formulas)
      
    }
    else
      if (!missing(Formula)) {
        
        for (dataset in name_index) {
          
          for (process in index2) {
            
            formula_update <- Formula
            
            if (!is.null(private$speciesName)) {
              
              covs_in <- all.vars(Formula)[all.vars(Formula) %in% private$spatcovsNames]
       
              ##What happens if people do species_var?
              if (!identical(covs_in, character(0))) {
                
                char_formula <- unlist(strsplit(as.character(Formula), split = ' '))
                
                char_formula[char_formula == covs_in] <- paste0(dataset, '_', covs_in)
           
                formula_update <- formula(paste(char_formula, collapse = ' '))
              
              }
              
            }
            
            updated_formula <- update(formula(paste('~ ', paste0(private$Formulas[[datasetName]][[dataset]][[process]]$RHS, collapse = ' + '))), formula_update)
            updated_formula <- all.vars(updated_formula)[all.vars(updated_formula) != '.']

            ##This will be removed once we move the like construction to runModel.
            private$Formulas[[datasetName]][[dataset]][[process]]$RHS <- updated_formula
            
          }
          
        }
        
      }
    else {
      
      for (dataset in name_index) { 
        
        for (process in index2) {
        ##Change all of these...
        ## Get old formula:
        #cat('Old formula for', paste0(dataset, ': '))
        if (length(all.vars(private$Formulas[[datasetName]][[dataset]][[process]]$LHS)) == 2) {
        
          oldForm <- update(private$Formulas[[datasetName]][[dataset]][[process]]$LHS, formula(paste(' ~ ',paste0(private$Formulas[[datasetName]][[dataset]][[process]]$RHS, collapse = ' + '))))
          
          #print(oldForm)
          #cat('New formula: ')
          
          ## Maybe I should do a check: if cov in newFormula then paste0(species, _ , covariate) ## how else does it work within a for loop
          newForm <- update(oldForm, newFormula)
          #newForm <- update(private$modelData[[dataset]]$formula, newFormula)
          #print(newForm)
          #cat('\n')
          
        }
        else {
          
          #print(private$modelData[[dataset]]$formula)
          #cat('\n')
          #cat('New formula: ')
          newForm <- update(paste(update(private$Formulas[[datasetName]][[dataset]][[process]]$LHS), '~ '), newFormula)
        }
        #Change
          
        private$Formulas[[datasetName]][[dataset]][[process]]$LHS <- newForm
        private$Formulas[[datasetName]][[dataset]][[process]]$RHS <- c()
        
      
      
      
    }
    }
    }
    }
  ,
  #' @description Function to add custom components to the integrated modeling.
  #' @param addComponent Component to add to the model.
  #' @param removeComponent Component (or name of a component) present in the model which should be removed.
  
  changeComponents = function(addComponent, removeComponent, print = TRUE) {
    
    terms <- gsub('\\(.*$', '', private$Components)
    
    if (!missing(addComponent)) {
      
      if (any(gsub('\\(.*$', '', addComponent) %in% terms)) private$Components <- private$Components[! terms %in% gsub('\\(.*$', '', addComponent)]
      
      private$Components <- c(private$Components, addComponent)
      
    } 
    
    if (!missing(removeComponent)) private$Components <- private$Components[!terms%in%removeComponent]
    
    componentsJoint <- formula(paste('~ - 1 +', paste(private$Components, collapse = ' + ')))
    componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
    
    if (print) {
      
      cat('Components:')
      cat('\n')
      print(componentsJoint)
      
    }
    
    
  }
  ,
  #' @description Function to change priors for the fixed (and possibly random) effects of the model.
  #' @param effect Name of the fixed effect covariate to change the prior for.
  #' @param species Name of the species for which the prior should change. Defaults to \code{NULL} which will change the prior for all species added to the model.
  #' @param dataset Name of the dataset for which the prior of the intercept should change (if fixedEffect = 'intercet'). Defaults to \code{NULL} which will change the prior effect of the intercepts for all the datasets in the model.
  #' @param mean.linear Mean value for the prior of the fixed effect. Defaults to \code{0}.
  #' @param prec.linear Precision value for the prior of the fixed effect. Defaults to \code{0.001}.

  priorsFixed = function(effect, species = NULL, dataset = NULL,
                         mean.linear = 0, prec.linear = 0.001) {
    
    
    if (!missing(effect)) {
      
      if (effect == 'intercept') {
        
        intTRUE <- TRUE
        
        if (is.null(private$speciesName)) {
          
          if (!private$Intercepts) stop('Fixed effect is given as "intercept", but intercepts have been turned off in bruSDM.')
          
          if (is.null(dataset)) effect <- paste0(unique(private$dataSource),'_intercept')
          
        } 
        else {
          
          if (is.null(species)) effect <- paste0(unique(unlist(private$speciesIn)), '_intercept')  #this won't work, unless we run a for loop...
          else effect <- paste0(species, '_intercept')
          
          
        }
        
      }
      else {
        
        intTRUE <- FALSE
        
        if (!effect %in% c(private$spatcovsNames, private$pointCovariates)) stop('Fixed effect provided not present in the model. Please add covariates using the "spatialCovariates" or "pointCovariates" argument in bruSDM.')
        
        
      } 
      
      if (effect %in% private$spatcovsNames) cov_class <- private$spatcovsClass[effect]
      else cov_class <- 'linear'
      
      if (!is.null(private$speciesName)) {
        
        if (!is.null(species)) {
          
          if (!species %in% unlist(private$speciesIn)) stop('Species given is not available in the model.')
          
          effect <- paste0(species, '_', effect) #this won't work, unless we run a for loop...
          
        }
        else effect <- paste0(unique(unlist(private$speciesIn)), '_', effect)
        
      }
      
      if (intTRUE) newComponent <- paste0(effect, '(1, mean.linear = ', mean.linear, ', prec.linear = ', prec.linear,' )')
      else newComponent <- paste0(effect,'(main = ', effect, ', model = \"', cov_class, '\", mean.linear = ', mean.linear, ', prec.linear = ', prec.linear, ')')
      
      self$changeComponents(addComponent = newComponent, print = FALSE)
      
    }
    
  }
  ,
  #' @description Function to specify the random fields in the model using PC priors for the parameters.
  #' 
  #' @param sharedSpatial Logical: specify the shared spatial field in the model. Defaults to \code{FALSE}.
  #' @param species Name of the species's spatial field to be specified.
  #' @param mark Name of the mark's spatial field to be specified.
  #' @param bias Name of the dataset's bias field to be specified.
  #' @param pc Logical: should the Matern model be specified with pc priors. Defaults to \code{TRUE}. 
  #' @param remove Logical: should the spatial field be removed. Requires one of sharedSpatial, species, mark or bias to be non-missing.
  #' @param ... Additional arguments used by INLA's \code{inla.spde2.pcmatern} or \code{inla.spde2.matern} function.
  
  specifySpatial = function(sharedSpatial = FALSE, 
                            species, mark,
                            bias, pc = TRUE,
                            remove = FALSE, ...) {
    
    if (all(!sharedSpatial && missing(species)  && missing(mark)  &&  missing(bias))) stop('At least one of sharedSpatial, dataset, species or mark needs to be provided.')
    
    if (sum(sharedSpatial, !missing(species), !missing(mark), !missing(bias)) != 1) stop('Please only choose one of sharedSpatial, species, mark or bias.')
    
    if (remove && sum(sharedSpatial, !missing(species), !missing(mark), !missing(bias)) !=1) stop('Please choose one of sharedSpatial, species, mark or bias to remove.')
    
    if (sharedSpatial) {
      
      if (!private$Spatial) stop('Shared spatial field not included in the model. Please use pointsSpatial = TRUE in bruSDM.')
      
      field_type <- 'sharedField'
      if (!remove) index <- 'sharedField'
      else index <- 'shared_spatial'
      
    }
    
    if (!missing(species)) {
      
      if (is.null(unlist(private$speciesIn))) stop('Species name provided but no species present in the model.') 
      
      if (!species %in% unlist(private$speciesIn)) stop('Species name provided is not currently in the model.')
      
      field_type <- 'speciesFields'
      if (!remove) index <- species
      else index <- paste0(species, '_spatial')
      
    }
    
    if (!missing(mark)) {
      
      if (is.null(unlist(private$markNames))) stop('Mark name provided but no marks present in the model.') 
      
      if (!mark %in% unlist(private$markNames)) stop('Mark name provided is not currently in the model.')
      
      if (!remove) field_type <- 'markFields'
      else index <- paste0(mark, '_spatial')
      
    } 
    
    if (!missing(bias)) {
      
      if (!bias %in% names(self$spatialFields$biasFields)) stop('Dataset name provided does not have a bias field. Please use ".$biasField()" beforehand.')
      
      field_type <- 'biasFields'
      if (!remove) index <- bias
      else index <- paste0(bias, 'biasField')
      
    }
    
    #if (missing(prior.range) || missing(prior.sigma)) stop('Both prior.range and prior.sigma need to be spefied.')
    
    if (!remove) {
      
      if (pc) model <- INLA::inla.spde2.pcmatern(mesh = private$INLAmesh, ...)
      else model <- INLA::inla.spde2.matern(mesh = private$INLAmesh, ...)
      
      for (field in index) {
        
        self$spatialFields[[field_type]][[field]] <- model
        
        
      }  
      
    }
    else {
      
      #do I need to remove from self$spatialFields[[index?]]  
      
      self$changeComponents(removeComponent = index, print = FALSE)
      
      for (data in unique(private$dataSource)) {  
        
        for (term in index) {
          
          self$updateFormula(datasetName = data, allProcesses = TRUE, Formula = formula(paste(' ~ . -', term)))  
          
        } 
      }
      
    }
    
  }
  
  ))

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
dataSDM$set('private', 'blockedCV', FALSE)
dataSDM$set('private', 'Formulas', list())

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

dataSDM$set('private', 'polyfromMesh', function(...) {
  
  loc <- private$INLAmesh$loc
  segm <- private$INLAmesh$segm$int
  
  coords <- na.omit(data.frame(loc[t(cbind(segm$idx[,, drop=FALSE], NA)), 1],
                               loc[t(cbind(segm$idx[,, drop=FALSE], NA)), 2]))
  
  Polys <- Polygon(coords = coords)
  Polys <- Polygons(srl = list(Polys), ID = 'id')
  SpatPolys <- SpatialPolygons(list(Polys), proj4string = private$Projection)
  SpatPolys
  
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
#' @param ... Extra arguments used by blockCV's spatialBlock
#' 
#' 
dataSDM$set('public', 'spatialBlock', function(k, rows, cols, plot = FALSE, ...) {
   
  #stop('Need to completely re do')
  #The easiest thing to do may be to do likelihood construction in runModel and blockedCV
  #So need to move all the objects from dataOrganize to here
  #And then be careful with regards to indexing for all the slot functions here...
  
  #Functions this would effect:
   #...$plot()
   #...$addData()
   #...$addBias()
   #...$updateFormula()
   #runModel()
   #dataOrganize()
  
  blocks <- R.devices::suppressGraphics(blockCV::spatialBlock(speciesData = do.call(rbind.SpatialPoints, lapply(private$modelData, function(x) x$data)),
                                                              k = k, rows = rows, cols = cols, selection = 'random',
                                                              verbose = FALSE, progress = FALSE, ...))
  
  folds <- blocks$blocks$folds
  
  blocksPoly <- list(sapply(1:(rows * cols), function(s) SpatialPolygons(blocks$blocks@polygons[s], proj4string = private$Projection)))
  
  blocked_data <- list()
  in_where <- list()
  
  ##Integration points?
  
  for (data in names(private$modelData)) {
    
    in_where[[data]] <- lapply(1:(rows * cols), function(i) !is.na(over(private$modelData[[data]]$data, blocksPoly[[1]][[i]])))
    
    for (i in 1:(rows * cols)) {
      
      blocked_data[[data]][[i]] <- private$modelData[[data]]$data[in_where[[data]][[i]], ]
      
      if (nrow(blocked_data[[data]][[i]]) !=0) blocked_data[[data]][[i]]$block_index <- as.character(folds[i])
      
    }
    
    private$modelData[[data]]$data <- do.call(rbind.SpatialPointsDataFrame, blocked_data[[data]])
    
  }
  
  private$blockedCV <- TRUE 
  
  if (plot) {
    
    spatPolys <- private$polyfromMesh()
    
    all_data <- do.call(rbind.SpatialPointsDataFrame, lapply(private$modelData, function(x) {
      
      if ('BRU_aggregate' %in% names(x$data)) ob <- x$data[x$data$BRU_aggregate, 'block_index']
      else ob <- x$data[,'block_index']
      
      ob
      
    }))
    
    ggplot() + gg(blocks$blocks) + blocks$plot$layers[[2]] +
      gg(all_data, aes(col = block_index)) +
      gg(blocks$blocks) +
      gg(spatPolys) +
      ggtitle('Plot of the blocked data') +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Need block block$plots + gg(species_data) + gg(boundary), which we need to get from the mesh
    
  }
  
})

## Need to change all the spatialFields to self$spatialFields and then the relevent sublist?
dataSDM$set('public', 'spatialFields', list(sharedField = list(),
                                            speciesFields = list(),
                                            markFields = list(),
                                            biasFields = list()))

