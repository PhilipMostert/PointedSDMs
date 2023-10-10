#' @title R6 class for creating a \code{dataSDM} object.
#' @description A data object containing the data and the relevant information about the integrated model. The function \code{\link{intModel}} acts as a wrapper in creating one of these objects. The output of this object has additional functions within the object which allow for further specification and customization of the integrated model.
#' @export
#' @importFrom R6 R6Class
#' 

dataSDM <- R6::R6Class(classname = 'dataSDM', lock_objects = FALSE, cloneable = FALSE, public = list(
  
  #' @description Prints the datasets, their data type and the number of observations, as well as the marks and their respective families.
  #' @param ... Not used.
  #' @import stats
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
  #' @description Makes a plot of the points surrounded by the boundary of the region where they were collected. The points may either be plotted based on which dataset they come from, or which species group they are part of (if \code{speciesName} is non-\code{NULL} in \code{\link{intModel}}).
  #' @param datasetNames Name of the datasets to plot. If this argument is missing, the function will plot all the data available to the model.
  #' @param Species Logical: should the points be plotted based on the species name. Defaults to \code{FALSE}.
  #' @param Boundary Logical: should a boundary (created using the \code{Mesh} object) be used in the plot. Defaults to \code{TRUE}.
  #' @param ... Not used.
  #' @return A ggplot object.
  #' @import ggplot2
  #' @examples
  #' 
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  library(ggplot2)
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                              Projection = proj, responsePA = 'Present')
  #'  
  #'   #Create plot of data
  #'   organizedData$plot()
  #' 
  #' }
  
  plot = function(datasetNames, Species = FALSE, 
                  Boundary = TRUE, ...) {
    
    if (length(private$modelData) == 0) stop('Please provide data before running the plot function.')
    
    if (missing(datasetNames)) datasetNames <- unique(private$dataSource)
    
    if (!all(datasetNames %in% private$dataSource)) stop('datasetNames provided not provided to the object.') 
    
    if (Species && is.null(private$speciesName)) stop('speciesName in intModel required before plotting species.')
    
    ##Get data
    points <- list()
    
    for (data in datasetNames) {
      
      index <- which(private$dataSource == data)
      
      if (!is.null(private$markNames)) {
        
        if (!is.null(private$speciesName)) index <- unique(private$speciesIn[[data]])
        else index <- 1
        #index <- index[!endsWith(names(private$modelData[index]), paste0('_', private$markNames))]
        
      } 
      
      for (i in 1:length(index)) {
        
        points[[data]][[i]] <- private$modelData[[data]][[i]][, names(private$modelData[[data]][[i]]) %in% c(private$temporalName, private$speciesName, 'geometry')]
       
        if ('BRU_aggregate' %in% names(points[[data]][[i]])) points[[data]][[i]] <- points[[data]][[i]][points[[data]][[i]]$BRU_aggregate,]
        
        if (!Species) points[[data]][[i]][,'..Dataset_placeholder_var..'] <- rep(data, nrow(points[[data]][[i]]))
        
        
      }
      
      #points[[data]] <- do.call(rbind.SpatialPointsDataFrame, points[[data]])
      
    }
    
    plotData <- do.call(rbind, unlist(points, recursive = F)) #plotData <- do.call(rbind, lapply(unlist(points), function(x) x[, names(x) %in% c('..Dataset_placeholder_var..', private$speciesName, private$temporalName)]))
    
    if (Boundary) {
      
      if (!is.null(private$Boundary)) bound <- geom_sf(data = sf::st_boundary(private$Boundary))
      else {
        
        bound <- try(geom_sf(data = sf::st_boundary(private$polyfromMesh())), silent = TRUE)
        
        if (inherits(bound, 'try-error')) {
          
          warning('Could not make a polygon from the mesh, polygon will be switched off.')
          bound <- NULL
          
        }
        
      }
      
    }
    else bound <- NULL
    
    if (!is.null(private$temporalName)) {
      
      if (Species) {
       
        plotData[, private$speciesName] <- unlist(private$speciesIndex) 
        
        colOption <- geom_sf(data = plotData, aes(col = eval(parse(text = private$speciesName))))#gg(plotData, aes(col = eval(parse(text = private$speciesName))))
        
        ggplot() + colOption + bound + guides(col = guide_legend(title = 'Species Name')) + facet_wrap(formula(paste('~', private$temporalName)))
        
      }
      else {
        
        colOption <- geom_sf(data = plotData, aes(col = eval(parse(text = '..Dataset_placeholder_var..'))))
  
        ggplot() + colOption + bound + guides(col = guide_legend(title = 'Dataset Name')) + facet_wrap(formula(paste('~', private$temporalName)))
        
      }
      
    } else {
    
    if (Species) {
      
      ## Need to add the species in here
      
      plotData[, private$speciesName] <- unlist(private$speciesIndex) #unlist(lapply(1:length(plotData[,private$speciesName]), function(x,y,z) y[z[x]], y = unlist(private$speciesIn), z= plotData@data[,private$speciesName]))
      
      colOption <- geom_sf(data = plotData, aes(col = eval(parse(text = private$speciesName))))
      
      ggplot() + colOption + bound + guides(col = guide_legend(title = 'Species Name')) 
      
    }
    else { 
      
      colOption <- geom_sf(data = plotData, aes(col = eval(parse(text = '..Dataset_placeholder_var..'))))
      
      ggplot() + colOption + bound + guides(col = guide_legend(title = 'Dataset Name')) 
      
    }
    
    }
    
    
  }
  ,
  #' @description Function used to add additional datasets to the \code{dataSDM} object. This function should be used if the user would like to add any additional datasets to the integrated model, but do not have the same standardized variable names as those added initially with \code{\link{intModel}}. Use \code{?intModel} for a more comprehensive description on what each argument in this function is used for.
  #' @param ... The datasets to be added to the integrated model: should be either \code{sf}, \code{data.frame} or \code{SpatialPoints*} objects, or a list of objects from these classes.
  #' @param responseCounts The name of the response variable for the counts data.
  #' @param responsePA The name of the response variable for the presence absence data.
  #' @param trialsPA The name of the trials variable for the presence absence data.
  #' @param markNames The names of the marks found in the data.
  #' @param markFamily The associated distributions of the marks.
  #' @param pointCovariates The additional, non-spatial covariates describing the data.
  #' @param trialsMarks The name of the trials variable for the binomial marks.
  #' @param speciesName The name of the species variable included in the data. Used to make a stacked species distribution model.
  #' @param temporalName The name of the temporal variable in the datasets.
  #' @param Coordinates A vector of length 2 describing the names of the coordinates of the data.
  #' @param Offset Name of the offset column in the dataset
  #' 
  #' @note The arguments of this function may be missing (ie not provided) if they have already been specified in \code{\link{intModel}}, and do not need changing. Therefore this function is useful if there are some variable names not standardized across the datasets; this function will thus standardize the variable names to those provided initially in \code{\link{intModel}}.
  #' 
  #' @import methods
  #' @import sf
  #' 
  #' @examples
  #'  
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  
  #'  #Only select eBird data
  #'  ebird <- SolitaryTinamou$datasets$eBird
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(ebird, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                              Projection = proj)
  #'                              
  #'  #Print summary
  #'  organizedData
  #'  
  #'  #Add new dataset
  #'  Parks = SolitaryTinamou$datasets$Parks
  #'  organizedData$addData(Parks, responsePA = 'Present')
  #'  
  #'  #Print summary
  #'  organizedData
  #'    
  #'  }
  addData = function(..., responseCounts, responsePA, trialsPA,
                     markNames, markFamily, pointCovariates,
                     trialsMarks, speciesName, temporalName,
                     Coordinates, Offset) {
    
    pointData <- dataOrganize$new()
    
    dataPoints <- list(...)
    
    if (missing(markNames)) markNames <- private$markNames # do same for species
    else {
      
      newmarksLens <- length(markNames)
      private$markNames <- markNames <- unique(c(private$markNames, markNames))
      
    }
    
    if (length(dataPoints) == 0) stop('Please provide data in the ... argument.')
    
    datasetClass <- unlist(lapply(dataPoints, class))
    
    if (length(datasetClass) == 1 && datasetClass == "list") {
      
      dataNames <- NULL
      dataPoints <- unlist(dataPoints, recursive = FALSE)
      datasetClass <- lapply(dataPoints, class)
      dataList <- TRUE
      
    }
    else dataList <- FALSE
    
    if (!all(unlist(datasetClass) %in% c("SpatialPointsDataFrame", "SpatialPoints", "data.frame", 'sf',
                                         'tbl', 'tbl_df'))) stop("Datasets need to be either a SpatialPoints* object, sf or a data frame.")
    
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
            private$initialNames <- dataNames
          }
          
        }
        
      }
    else {
      
      dataNames <- setdiff(as.character(match.call(expand.dots = TRUE)), 
                           as.character(match.call(expand.dots = FALSE)))
      private$initialNames <- dataNames
    
      }
    
    if (!is.null(private$speciesName)) {
      
      speciesOK <- checkVar(data = dataPoints,
                            var = private$speciesName)
      
      if (!speciesOK) stop('The species variable name is required to be present in all the datasets.')
      #Test
      speciesOld <- unique(unlist(private$speciesIn))
      
      speciesNew <- unique(unlist(lapply(dataPoints, function(x) {
        
        if (inherits(x, 'Spatial')) x@data[,private$speciesName]
        else data.frame(x)[,private$speciesName]
        
        
      })))
      
      speciesIn <- c(speciesOld, speciesNew)
      
      
    }
    
    if (!is.null(private$Spatial)) {
      
      if (private$Spatial == 'shared' && is.null(self$spatialFields$sharedField[['sharedField']])) self$spatialFields$sharedField[['sharedField']] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
      
      else 
        if (private$Spatial == 'copy') {
          
          mainName <- dataNames[[1]]
          
          if (is.null(self$spatialFields$datasetFields[[mainName]])) self$spatialFields$datasetFields[[mainName]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
          
          } 
      else {
        
        for (dataset in dataNames) {
          
          if (is.null(self$spatialFields$datasetFields[[dataset]])) self$spatialFields$datasetFields[[dataset]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
          
        }
        
      }
      
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
    
    if (missing(Offset)) Offset <- private$Offset
    else {
      
      if (!is.null(Offset) && any(Offset != private$Offset)) {
        
        dataPoints <- nameChanger(data = dataPoints, oldName = Offset,
                                  newName = private$Offset)
        Offset <-  private$Offset
        
      }
      
    }
    
    if (!is.null(Offset) && any(Offset %in% 'offset')) stop('The offset variable cannot be called "offset". Please choose a new name.')
    
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
      
      namesIn <- unlist(sapply(dataPoints, function(x) names(x)))
      if (length(private$modelData) != 0) namesIn <- c(namesIn,  unlist(sapply(unlist(private$modelData, recursive = FALSE), function(x) names(x)[names(x) != 'geometry'])))
        
      if (!all(markNames %in% unlist(namesIn))) stop('At least one mark specified is not present in the datasets, please check again.')
        
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
    
    if (!is.null(private$pointCovariates) && !inherits(pointCovariates, 'character')) stop('pointCovariates is required to be a vector containing the names of the covariates found in the datasets.')
    
    if (!missing(Coordinates) && Coordinates != private$Coordinates) {
      
      dataPoints <- changeCoords(data = dataPoints, oldcoords = Coordinates,
                                 newcoords = private$Coordinates)
      
    } 
    else { 
      
      coordsOK <- checkCoords(data = dataPoints,
                              coords = private$Coordinates)
      
      if (!coordsOK) stop('At least one dataset does not have the correct coordinates. Please try again.')  
      
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
                       pointcovnames = pointCovariates, offsetname = private$Offset)
    
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
      self$spatialFields$speciesFields <- vector(mode = 'list', length = length(unique(unlist(private$speciesIn))) * length(names(private$speciesIn)))
      
      names(self$spatialFields$speciesFields) <- do.call(paste0, expand.grid(paste0(unique(unlist(private$speciesIn)),'_'), names(private$speciesIn)))
     
      if (!is.null(private$speciesSpatial)) {
        
        if (!all(do.call(paste0, expand.grid(paste0(speciesIn,'_'), names(private$speciesIn))) %in% names(self$spatialFields$speciesFields))) {
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
      
      private$temporalVars <- pointData$timeIndex
      
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
                                                     speciesspatial = private$speciesSpatial,
                                                     offsetname = private$Offset,
                                                     copymodel = private$copyModel)
      
    }
    else {
      ##Add if species; then remove first species_spatial
      
      #if (!is.null(private$speciesName) && private$speciesSpatial) {
      #  ##Select species_spatial with lowest ngroup...
      #  ##or maybe just do ngroups again ...
      #  private$Components <- private$Components[!grepl(paste0('^',private$speciesName,'_spatial'), private$Components)]
      #  
      #}
      
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
                                                speciesspatial = private$speciesSpatial,
                                                offsetname = private$Offset,
                                                copymodel = private$copyModel)
      
      
      private$Components <- union(private$Components, newComponents)
      
      
    }
    
    if (!is.null(c(private$Offset, private$pointCovariates))) {
      
      datMatrix <- as.data.frame(matrix(NA, nrow = nrow(private$IPS), ncol = length(c(private$Offset, private$pointCovariates))))
      names(datMatrix) <- c(private$pointCovariates, private$Offset)
      private$IPS <- cbind(private$IPS, datMatrix)
      
    }
    if (length(private$Formulas) == 0)  private$Formulas <- pointData$Formulas
    else private$Formulas <- append(private$Formulas, pointData$Formulas)
    
    if (length(private$Family) == 0) private$Family <- pointData$Family
    else private$Family <- append(private$Family, pointData$Family)
    
    if (length(private$speciesIndex) == 0) private$speciesIndex <- pointData$speciesIndex
    else private$speciesIndex <- append(private$speciesIndex, pointData$speciesIndex)
    
    #pointData$makeLhoods(mesh = private$INLAmesh,
    #                     ips = private$IPS, paresp = responsePA,
    #                     ntrialsvar = trialsPA,
    #                     markstrialsvar = trialsMarks,
    #                     speciesname = speciesName)
    
    if (!is.null(private$speciesName)) newFamily <- mapply(function(family, number) rep(family , times = number), family = private$Family,
                                                           number = lapply(private$speciesIn, length))
    
    else newFamily <- pointData$Family
    
    if (!is.null(private$optionsINLA[['control.family']])) {
      
      index1 <- length(private$optionsINLA[['control.family']]) + 1
      index2 <- length(private$optionsINLA[['control.family']]) + length(unlist(newFamily))
      
    }
    else {
      
      index1 <- 1
      index2 <- length(unlist(newFamily))
      
    }
    ##Rather just re do this with the changeLink function ... 
    familyIndex <- c(rep(NA, length(private$optionsINLA[['control.family']])), unlist(newFamily)) #Make this the length rep(NA, length(private$inlaOptions$link whatervers))
    
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
    
    if (!is.null(private$initialnames)) private$initialnames <- NULL
    
  }
  ,
  #' @description Function used to add additional spatial fields (called \emph{bias fields}) to a selected dataset present in the integrated model. \emph{Bias fields} are typically used to account for sampling biases in opportunistic citizen science data in the absence of any covariate to do such.
  #' @param datasetNames A vector of dataset names (class \code{character}) for which a bias field needs to be added to. If \code{NULL} (default), then \code{allPO} has to be \code{TRUE}.
  #' @param allPO Logical: should a bias field be added to all datasets classified as presence only in the integrated model. Defaults to \code{FALSE}.
  #' @param biasField An \code{inla.spde} object used to describe the bias field. Defaults to \code{NULL} which uses \code{\link[INLA]{inla.spde2.matern}} to create a Matern model for the field.
  #' @param temporalModel List of model specifications given to the control.group argument in the time effect component. Defaults to \code{list(model = 'ar1')}; see \code{\link[INLA]{control.group}} from the \pkg{INLA} package for more details. \code{temporalName} needs to be specified in \code{intModel} prior.
  #' @examples
  #'  
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                              Projection = proj, responsePA = 'Present')
  #'  
  #' #Add bias field to eBird records
  #' organizedData$addBias(datasetNames = 'eBird')
  #' 
  #' }
  addBias = function(datasetNames = NULL,
                     allPO = FALSE,
                     biasField = NULL,
                     temporalModel = list(model = 'ar1')) {
    
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
      
      if (is.null(biasField)) self$spatialFields$biasFields[[dat]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
      else self$spatialFields$biasFields[[dat]] <- biasField
      
    }
    
    
    
    ##Add a way for separate Fields per dataset.
    #Shouldn't be too difficult -- have already done it before.
    
    
    #Should I copy the bias fields for the marks?
    
    if (!is.null(private$temporalName)) {

      temporalModel <- deparse1(temporalModel)
      
      private$Components <- c(private$Components, paste0(datasetNames ,'_biasField(main = geometry, model = ', datasetNames, '_bias_field, group = ', private$temporalName, ', ngroup = ', length(unique(unlist(private$temporalVars))),', control.group = ', temporalModel,')'))
      
      
    }
    else {
      
      private$Components <- c(private$Components, paste0(datasetNames ,'_biasField(main = geometry, model = ', datasetNames, '_bias_field)'))
      
    }
    
    
    ##Things to do here:
    #Go into the liks of PO datasets and add the biasfield
    #Go inth the components and add the bias field component
    #Then store the bias field somewhere else. If field NULL, then add generic bias field
    
    
  }
  ,
  #' @description Function used to update the formula for a selected observation model. The function is designed to work similarly to the generic \code{update} formula, and should be used to thin terms out of a process from the full model specified in \code{\link{intModel}}. The function also allows the user to add their own formula to the model, such that they can include non-linear components in the model. The function can also be used to print out the formula for a process by not specifying the \code{Formula} or \code{newFormula} arguments.
  #' @param datasetName Name of the dataset (class \code{character}) for which the formula needs to be changed.
  #' @param Points Logical: should the formula be changed for the points (or otherwise, a marked process). Defaults to \code{TRUE}. If \code{FALSE}, then \code{markNames} needs to be non-\code{NULL} in \code{\link{intModel}}.
  #' @param speciesName Name of the species (class \code{character}) to change the formula for. Defaults to \code{NULL}. If \code{NULL} and \code{speciesName} is non-\code{NULL} in \code{\link{intModel}}, will update the formula for all species within the dataset. Cannot be non-\code{NULL} if \code{speciesName} is \code{NULL} in \code{\link{intModel}}.
  #' @param markName Name of the mark (class \code{character}) to change the formula for. Defaults to \code{NULL}. If \code{NULL} and \code{markNames} is non-\code{NULL} in \code{\link{intModel}}, will update the formula for all marks within the dataset. Cannot be non-\code{NULL} if \code{markNames} is \code{NULL} in \code{\link{intModel}}.
  #' @param Formula An updated formula to give to the process. The syntax provided for the formula in this argument should be identical to the formula specification as in base \strong{R}. Should be used to thin terms out of a formula but could be used to add terms as well. If adding new terms not specified in \code{intModel}, remember to add the associated component using \code{.$addComponents} as well.
  #' @param allProcesses Logical argument: if \code{TRUE} changes the formulas for all of the processes in a dataset. Defaults to \code{FALSE}.
  #' @param newFormula Completely change the formula for a process -- primarily used to add non-linear components into the formula. Note: all terms need to be correctly specified here.
  #' 
  #' @import methods
  #' @import stats
  #' 
  #' @examples
  #' 
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  Forest <- SolitaryTinamou$covariates$Forest
  #'  
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present',
  #'                            pointsSpatial = 'individual')
  #' 
  #'  #Remove Forest from eBird
  #'  organizedData$updateFormula(datasetName = 'eBird', Formula = ~ . - Forest)
  #'  
  #'  #Add some scaling to Forest for Parks
  #'  organizedData$updateFormula('Parks', newFormula = ~ I(. +(Forest+1e-6)*scaling))
  #'  
  #'  #Now dd scaling to components
  #'  organizedData$changeComponents(addComponent = 'scaling') 
  #'  
  #' }
  #' 
  #' @return If \code{Formula} and \code{newFormula} are missing, will print out the formula for the specified processes. 
  #' 
  updateFormula = function(datasetName = NULL, Points = TRUE, speciesName = NULL,
                           markName = NULL, Formula, allProcesses = FALSE,
                           newFormula) {
    ## Will need to update this such that if keepSpatial & species then paste0(species,_spatial)
    
    if (all(is.null(datasetName), is.null(speciesName), is.null(markName))) stop ('At least one of: datasetName, speciesName, markName or allProcesses needs to be specified.')
    
    if (!is.null(speciesName) && is.null(private$speciesName)) stop ('Species are given but none are present in the model. Please specify species in the model with "speciesName" in intModel.')
    
    if (is.null(speciesName) && !is.null(private$speciesName)) speciesName <- unlist(private$speciesIn[datasetName])
    
    if (!Points && is.null(markName)) stop ('markNames cannot be non-null if Points is FALSE.')
    
    if (is.null(datasetName) && !is.null(markName)) {
      
      if (length(private$Formulas) > 1) stop('Please supply a dataset name with the mark.')
      else datasetName <- names(private$Formulas)
      
    }
    
    if (length(datasetName) != 1) stop ('Please only provide one dataset name.')
    
    if (!datasetName %in% private$dataSource) stop ('Dataset name provided not in model.')
    
    if (!missing(Formula) && !missing(newFormula)) stop ('Please provide only one of Formula and newFormula. \n Use Formula to update the current formula; use newFormula to make a completely new formula.')
    
    if (!missing(Formula) && !inherits(Formula, 'formula')) stop ('Formula must be of class "formula".')
    
    if (!missing(newFormula) && !inherits(newFormula, 'formula')) stop('newFormula must be of class "formula".')
    
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
        else index2 <- c(which(names(private$Formulas[[datasetName]][[name_index[1]]]) %in% markName)) #Since they should all be the same ...
        
        
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
  #' @description Function to add and specify custom components to model, which are required by \pkg{inlabru}. The main purpose of the function is to re-specify or completely change components already in the model, however the user can also add completely new components to the model as well. In this case, the components need to be added to the correct formulas in the model using the \code{.$updateFormula} function. If \code{addComponent} and \code{removeComponent} are both missing, the function will print out the components to be supplied to \pkg{inlabru}'s \code{\link{bru}} function.
  #' @param addComponent Component to add to the integrated model. Note that if the user is re-specifying a component already present in the model, they do not need to remove the old component using \code{removeComponent}.
  #' @param removeComponent Component (or just the name of a component) present in the model which should be removed.
  #' @param print Logical: should the updated components be printed. Defaults to \code{TRUE}.
  #' @examples 
  #' \dontrun{
  #' 
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  Forest <- SolitaryTinamou$covariates$Forest
  #'  
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present')
  #' 
  #'  #Remove Forest from components
  #'  organizedData$changeComponents(removeComponent = 'Forest')
  #'
  #' }
  #' 
  #' }
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
  #' @param Effect Name of the fixed effect covariate to change the prior for. Can take on \code{'intercept'}, which will change the specification for an intercept (specified by one of \code{species} or \code{datasetName}).
  #' @param Species Name of the species (class \code{character}) for which the prior should change. Defaults to \code{NULL} which will change the prior for all species added to the model.
  #' @param datasetName Name of the dataset for which the prior of the intercept should change (if fixedEffect = 'intercept'). Defaults to \code{NULL} which will change the prior effect of the intercepts for all the datasets in the model.
  #' @param mean.linear Mean value for the prior of the fixed effect. Defaults to \code{0}.
  #' @param prec.linear Precision value for the prior of the fixed effect. Defaults to \code{0.001}.
  #' @examples
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  Forest <- terra::rast(
  #'  system.file(
  #'  'extdata/SolitaryTinamouCovariates.tif', 
  #'  package = "PointedSDMs"))$Forest
  #'  
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present',
  #'                            pointsSpatial = 'individual')
  #' 
  #'  #Add prior to Forest
  #'  organizedData$priorsFixed(Effect = 'Forest', mean.linear = 2, prec.linear = 0.1)
  #'
  #' }
  
  priorsFixed = function(Effect, Species = NULL, datasetName = NULL,
                         mean.linear = 0, prec.linear = 0.001) {
    
    
    if (!missing(Effect)) {
      
      if (Effect == 'intercept') {
        
        intTRUE <- TRUE
        
        if (is.null(private$speciesName)) {
          
          if (!private$Intercepts) stop('Fixed effect is given as "intercept", but intercepts have been turned off in intModel.')
          
          if (is.null(datasetName)) Effect <- paste0(unique(private$dataSource),'_intercept')
          
        } 
        else {
          
          if (is.null(Species)) Effect <- paste0(unique(unlist(private$speciesIn)), '_intercept')  #this won't work, unless we run a for loop...
          else Effect <- paste0(Species, '_intercept')
          
          
        }
        
      }
      else {
        
        intTRUE <- FALSE
        
        if (!Effect %in% c(private$spatcovsNames, private$pointCovariates)) stop('Fixed effect provided not present in the model. Please add covariates using the "spatialCovariates" or "pointCovariates" argument in intModel.')
        
        
      } 
      
      if (Effect %in% private$spatcovsNames) cov_class <- private$spatcovsClass[Effect]
      else cov_class <- 'linear'
      
      if (!is.null(private$speciesName)) {
        
        if (!is.null(Species)) {
          
          if (!Species %in% unlist(private$speciesIn)) stop('Species given is not available in the model.')
          
          Effect <- paste0(Species, '_', Effect) #this won't work, unless we run a for loop...
          
        }
        else Effect <- paste0(unique(unlist(private$speciesIn)), '_', Effect)
        
      }
      
      if (intTRUE) newComponent <- paste0(Effect, '(1, mean.linear = ', mean.linear, ', prec.linear = ', prec.linear,' )')
      else newComponent <- paste0(Effect,'(main = ', Effect, ', model = \"', cov_class, '\", mean.linear = ', mean.linear, ', prec.linear = ', prec.linear, ')')
      
      self$changeComponents(addComponent = newComponent, print = FALSE)
      
    }
    
  }
  ,
  #' @description Function to specify random fields in the model using penalizing complexity (PC) priors for the parameters.
  #' 
  #' @param sharedSpatial Logical: specify the shared spatial field in the model. Requires \code{pointsSpatial == 'shared'} in \code{\link{intModel}}. Defaults to \code{FALSE}.
  #' @param datasetName Name of which of the datasets' spatial fields to be specified. Requires \code{pointsSpatial = 'individual'} in \code{\link{intModel}}.
  #' @param Species Name of which of the species' spatial field to be specified. Requires \code{speciesName} to be non-\code{NULL} in \code{\link{intModel}}.
  #' @param Mark Name of which of the mark' spatial field to be specified. Requires \code{markNames} to be non-\code{NULL} in \code{\link{intModel}}.
  #' @param Bias Name of the dataset for which the bias field to be specified.
  #' @param Copy Name of the dataset to create a copy using \code{\link[INLA]'s copy feature for a mark or bias field. Defaults to \code{NULL} which does not copy any field. Requires one of \code{Mark} or \code{Bias} to be specified.
  #' @param PC Logical: should the Matern model be specified with pc priors. Defaults to \code{TRUE}, which uses \code{\link[INLA]{inla.spde2.pcmatern}} to specify the model; otherwise uses \code{\link[INLA]{inla.spde2.matern}}.
  #' @param Remove Logical: should the chosen spatial field be removed. Requires one of \code{sharedSpatial}, \code{species}, \code{mark} or \code{bias} to be non-missing, which chooses which field to remove.
  #' @param ... Additional arguments used by \pkg{INLA}'s \code{\link[INLA]{inla.spde2.pcmatern}} or \code{\link[INLA]{inla.spde2.matern}} function, dependent on the value of \code{PC}.
  #' 
  #' 
  #' @examples 
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  Forest <- SolitaryTinamou$covariates$Forest
  #'  
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present')
  #' 
  #'  #Specify the shared spatial field
  #'  organizedData$specifySpatial(sharedSpatial = TRUE, PC = TRUE, 
  #'                        prior.range = c(1,0.001),
  #'                        prior.sigma = c(1,0.001))
  #'
  #' } 
  
  specifySpatial = function(sharedSpatial = FALSE,
                            datasetName,
                            Species, Mark,
                            Bias, PC = TRUE,
                            Copy = NULL,
                            Remove = FALSE, ...) {
    
    if (all(!sharedSpatial && missing(datasetName) && missing(Species)  && missing(Mark)  &&  missing(Bias))) stop('At least one of sharedSpatial, datasetName, dataset, Species or Mark needs to be provided.')
    
    if (sum(sharedSpatial, !missing(datasetName), !missing(Species), !missing(Mark), !missing(Bias)) != 1) stop('Please only choose one of sharedSpatial, datasetName, Species, Mark or Bias.')
    
    if (Remove && sum(sharedSpatial, !missing(datasetName), !missing(Species), !missing(Mark), !missing(Bias)) != 1) stop('Please choose one of sharedSpatial, datasetName, Species, Mark or Bias to remove.')
    
    if (!is.null(Copy) && !Copy %in% unlist(private$dataSource)) stop('Dataset name provided is not currently in the model.')
    
    if (sharedSpatial) {
      
      if (is.null(private$Spatial)) stop('Shared spatial field not included in the model. Please use pointsSpatial = "shared" or pointsSpatial = "copy" in intModel.')
      else 
        if (private$Spatial == 'individual') stop('pointsSpatial specified as "individual" in intModel. Please specify a dataset spatial effect to specify with datasetName.')
      
      
      if (private$Spatial == 'shared') {
        
      field_type <- 'sharedField'
      if (!Remove) index <- 'sharedField'
      else index <- 'shared_spatial'
      
      } else {
        
        field_type <- 'datasetFields'
        if (!Remove) index <- private$initialNames[[1]]
        else index <- private$initialNames[[1]]
      }
      
    }
    
    if (!missing(datasetName)) {
      
      if (!datasetName %in% unlist(private$dataSource)) stop('Dataset name provided is not currently in the model.')
      
      field_type <- 'datasetFields'
      if (!Remove) index <- datasetName
      else index <- datasetName
      
    }
    
    if (!missing(Species)) {
      
      if (is.null(unlist(private$speciesIn))) stop('Species name provided but no species present in the model.') 
      
      if (!Species %in% unlist(private$speciesIn)) stop('Species name provided is not currently in the model.')
      
      field_type <- 'speciesFields'
      if (!Remove) index <- Species
      else index <- do.call(paste0, expand.grid(paste0(Species, '_'), private$dataSource))
      
    }
    
    if (!missing(Mark)) {
      
      if (is.null(unlist(private$markNames))) stop('Mark name provided but no marks present in the model.') 
      
      if (!Mark %in% unlist(private$markNames)) stop('Mark name provided is not currently in the model.')
      
      if (!is.null(Copy)) {
        
        if (Mark %in% names(private$modelData[[copy]][[1]])) stop('Mark not in dataset for provided.')
        else fieldCopied <- paste0(Mark,'_spatial(main = geometry, copy = \"', Copy, '\")')
      }
      
      field_type <- 'markFields'
      if (!Remove) index <- Mark
      else index <- paste0(Mark, '_spatial')
      
    } 
    
    if (!missing(Bias)) {
      
      if (!Bias %in% names(self$spatialFields$biasFields)) stop('Dataset name provided does not have a bias field. Please use ".$biasField()" beforehand.')
      
      field_type <- 'biasFields'
      if (!Remove) index <- Bias
      else index <- paste0(Bias, 'biasField')
      
      if (!is.null(Copy)) fieldCopied <- paste0(Bias,'_biasField(main = geometry, copy = \"', Bias, '\")')
      
    }
    
    #if (missing(prior.range) || missing(prior.sigma)) stop('Both prior.range and prior.sigma need to be spefied.')
    
    if (!is.null(Copy)) {
      
      self$changeComponents(addComponent = fieldCopied, print = FALSE)
      
    }
    
    if (!Remove) {
      
      if (PC) model <- INLA::inla.spde2.pcmatern(mesh = private$INLAmesh, ...)
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
  ,
  #' @description Function used to change the link function for a given process.
  #' @param datasetName Name of the dataset for which the link function needs to be changed.
  #' @param Species Name of the species for which the link function needs to be changed.
  #' @param Mark Name of the mark for which the link function needs to be changed.
  #' @param Link Name of the link function to add to the process. If missing, will print the link function of the specified dataset.
  #' @param ... Not used
  #' 
  #' @examples
  #' \dontrun{
  #' 
  #' #Create data object
  #' dataObj <- intModel(...)
  #' 
  #' #Print link function for a process
  #' 
  #' dataObj$changeLink(Dataset = Dataset, Species = Species)
  #' 
  #' #Change link function
  #' dataObj$changeLink(Dataset = Dataset, Species = Species,
  #'                    Link = "log")
  #' 
  #' }
  changeLink = function(datasetName, Species, Mark,
                        Link, ...) {
    
    stop('Still need to do.')
    if (missing(Dataset)) stop('Please provide a dataset name.')
    
    if (length(datasetName) > 1) stop('Please only provide one dataset at a time.')
    
    if (!datasetName %in% names(private$dataSource)) stop('Dataset name provided not in model.')
    
    if (!is.null(private$speciesName) && missing(Species)) stop('Please provide a species name.')
    
  }
  ,
  #' @description Function to spatially block the datasets, which will then be used for model cross-validation with \code{\link{blockedCV}}. See the \code{\link[blockCV]{spatialBlock}} function from \pkg{blockCV} for how the spatial blocking works and for further details on the function's arguments.
  #' @param k Integer value reflecting the number of folds to use.
  #' @param rows_cols Integer value by which the area is divided into longitudinal and latitudinal bins.
  #' @param plot Plot the cross-validation folds as well as the points across the boundary. Defaults to \code{FALSE}.
  #' @param seed Seed used by \pkg{blockCV}'s \code{\link[blockCV]{spatialBlock}} to make the spatial blocking reproducible across different models. Defaults to \code{1234}. 
  #' @param ... Additional arguments used by \pkg{blockCV}'s \code{\link[blockCV]{spatialBlock}}.
  #' 
  #' @import ggplot2
  #' @importFrom R.devices suppressGraphics
  #' @importFrom blockCV spatialBlock
  #' 
  #' @examples
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  Forest <- SolitaryTinamou$covariates$Forest
  #'  
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present',
  #'                            pointsSpatial = 'individual')
  #' 
  #'  #Specify the spatial block
  #'  organizedData$spatialBlock(k = 2, rows = 2, cols = 1, plot = FALSE)
  #'
  #' } 
  #'
  spatialBlock = function(k, rows_cols, plot = FALSE, seed = 1234, ...) {
    
    
    private$spatialBlockCall <- paste0(gsub('.*\\(', 'self$spatialBlock(', deparse(match.call())))

    #blocks <- R.devices::suppressGraphics(blockCV::spatialBlock(speciesData = do.call(sp::rbind.SpatialPoints, append(unlist(private$modelData),private$IPS)),
    #                                                            k = k, rows = rows, cols = cols, selection = 'random',
    #                                                            verbose = FALSE, progress = FALSE, seed = seed, ...))
    
    allPoints <- append(unlist(private$modelData, recursive = FALSE), list(private$IPS))
    
    allPoints <- do.call(c, lapply(allPoints, st_geometry))
    
    blocks <- R.devices::suppressGraphics(blockCV::cv_spatial(x = allPoints,
                                                              k = k, rows_cols = rows_cols, progress = FALSE, seed = seed, report = FALSE, plot = plot, ...))
    
    
    foldID <- blocks$folds_ids
    
    dataLength <- unlist(lapply(append(unlist(private$modelData, recursive = FALSE), list(private$IPS)), nrow))
                  
    for (i in 1:length(dataLength)) {
   
      if (i != length(dataLength)) {
        
        if (i == 1) start <- 0
        else start <- cumsum(dataLength)[i-1]
        
        private$modelData[[i]][[1]]$.__block_index__ <- as.character(foldID[(1 + start):cumsum(dataLength)[i]])
        
        if (any(is.na(private$modelData[[i]][[1]]$.__block_index__))) {
          
          warning('NA values found in the blocking: will remove these observations')
          private$modelData[[i]][[1]] <-  private$modelData[[i]][[1]][!is.na(private$modelData[[i]][[1]]$.__block_index__),]
          
        }
        
      }
      
      else {
        
        start <- cumsum(dataLength)[i-1]
        
        private$IPS$.__block_index__ <- as.character(foldID[(1 + start): cumsum(dataLength)[i]])
        private$IPS <- private$IPS[!is.na(private$IPS$.__block_index__), ]
        
        
      }
      
      
    }
    ##Temporary fix
    
    #folds <- blocks$blocks$folds
    
    #blocks$blocks <- as(blocks$blocks, 'Spatial')
    
    #blocksPoly <- list(sapply(1:(rows * cols), function(s) SpatialPolygons(blocks$blocks@polygons[s], proj4string = private$Projection)))
    
    #blocksPoly <- list(sapply(1:(rows * cols), function(s) blocks$blocks$geometry[s])) #, proj4string = private$Projection
    #https://github.com/r-spatial/sf/wiki/migrating read this to see how to get points over
    
    #blocked_data <- list()
    #in_where <- list()
    
    #for (data in names(private$modelData)) {
      
    #  for (process in names(private$modelData[[data]])) {
        
    #    in_where[[data]][[process]] <- lapply(1:(rows * cols), function(i) !is.na(over(private$modelData[[data]][[process]], blocksPoly[[1]][[i]])))
        
    #    for (i in 1:(rows * cols)) {
          
    #      blocked_data[[data]][[process]][[i]] <- private$modelData[[data]][[process]][in_where[[data]][[process]][[i]], ]
          
    #      if (nrow(blocked_data[[data]][[process]][[i]]) !=0) blocked_data[[data]][[process]][[i]]$.__block_index__ <- as.character(folds[i])
          
    #    }
    
    #    blocked_data[[data]][[process]] <- lapply(blocked_data[[data]][[process]], function(x) {
          
    #      row.names(x@data) <- NULL
    #      row.names(x@coords) <- NULL
    #      x
          
    #    })
        
    #    private$modelData[[data]][[process]] <- do.call(sp::rbind.SpatialPointsDataFrame, blocked_data[[data]][[process]])
        
    #  }
      
    #}
    
    #blocked_ips <- list()
    #where_ips <- lapply(1:(rows * cols), function(i) !is.na(over(private$IPS, blocksPoly[[1]][[i]])))
    
    #for (i in 1:(rows * cols)) {
      
    #  blocked_ips[[i]] <- private$IPS[where_ips[[i]], ]
      
    #  if (nrow(blocked_ips[[i]]) !=0) blocked_ips[[i]]$.__block_index__ <- as.character(folds[i])
      
    #}
    #private$IPS <- do.call(sp::rbind.SpatialPointsDataFrame, blocked_ips)
    
    if (length(private$biasData) > 0) {
      
      blocked_samplers <- list()
      in_where_samplers <- list()
      
     for (sampler in names(private$biasData)) {
       
       if (class(private$biasData[[sampler]] %in% c('SpatialPoints', 'SpatialPointsDataFrame'))) {
         
         in_where_samplers[[sampler]] <- lapply(1:(rows * cols), function(i) !is.na(over(private$biasData[[sampler]], blocksPoly[[1]][[i]])))
         
         for (i in 1:(rows * cols)) {
           
           blocked_samplers[[sampler]][[i]] <- private$biasData[[sampler]][in_where[[sampler]][[i]], ]
           
           if (nrow(blocked_samplers[[sampler]][[i]]) !=0) blocked_samplers[[sampler]][[i]]$.__block_index__ <- as.character(folds[i])
           
         }
         
         blocked_samplers[[samplers]] <- lapply(blocked_samplers[[samplers]], function(x) {
           
           row.names(x@data) <- NULL
           row.names(x@coords) <- NULL
           x
           
         })
         
         private$biasData[[samplers]] <- do.call(sp::rbind.SpatialPointsDataFrame, blocked_samplers[[samplers]])
         
       } else 
         if (class(private$biasData[[sampler]] %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame'))) {
           
           #What happens if the SpatialPolygon overlaps on two different blocks...
           #Would we have to completely have to rethink this type of data ie have a list which includes what points are in each polygon/line (overlap included.)
           
         }
       else {
         
         #Assuming this is a SpatialLines/SpatialLinesDataFrame object ...
         
       }
       
       
     }
      
      # if class(private$biasData) == 'SpatialPoints' easy...
      # if class(private$biasData) == 'SpatialPolygons' don't know...
      ## if class(private$biasData) == 'SpatialLines': then we need to LineIn <- rgeos::gIntersects(SpatialLines(msamplers@lines[i], proj4string = private$Projection),SpatialPolygons(polygons, proj4string = private$Projection))
      # With this we would probably need to somehow separate the lines and polygons; or add metadata which states where they are in the block
            #maybe something like samplers@lines[[i]]@Lines$block = block?
    }
    
    private$blockedCV <- TRUE 
    
    if (plot) {
      
      spatPolys <- geom_sf(data = sf::st_boundary(private$polyfromMesh()))
      
      all_data <- do.call(rbind, lapply(unlist(private$modelData, recursive = FALSE), function(x) {
        
        x[, '.__block_index__']
        
        
      }))

     # ggplot() + 
        #geom_sf(data = blocks$blocks) +
        #blocks$plot$layers[[2]] +
        cv_plot(blocks) +
        geom_sf(data = all_data, aes(col = .__block_index__)) + 
        spatPolys +
        labs(col = 'Block index') +
        ggtitle('Plot of the blocked data') +
        theme(plot.title = element_text(hjust = 0.5))
  
    }
    
  }
  ,
  #' @description Function to add an integration domain for the PO datasets.
  #' @param datasetName Name of the dataset for the samplers.
  #' @param Samplers A \code{Spatial*} object representing the integration domain.
  #' 
  #' @example
  #' 
  #'  if (requireNamespace('INLA')) {
  #'    
  #'  #Get Data
  #'  data("SolitaryTinamou")
  #'  proj <- "+proj=longlat +ellps=WGS84"
  #'  data <- SolitaryTinamou$datasets
  #'  mesh <- SolitaryTinamou$mesh
  #'  mesh$crs <- proj
  #'  
  #'  #Set model up
  #'  organizedData <- intModel(data, Mesh = mesh, Coordinates = c('X', 'Y'),
  #'                              Projection = proj, responsePA = 'Present')
  #'  
  #' #Add integration domain for the eBird records
  #' organizedData$addSamplers(datasetName = 'eBird', Samplers = SolitaryTinamou$region)
  #' 
  #' }
  #' 
  addSamplers = function(datasetName, Samplers) {
    
    if (!datasetName %in% private$dataSource) stop ('Dataset name provided not in model.')
    
    if (!inherits(Samplers, 'Spatial')) stop ('Samplers needs to be a Spatial* object.')
    
    Samplers@proj4string <- private$Projection
    
    private$Samplers[[datasetName]] <- Samplers
    
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
dataSDM$set('private', 'Offset', NULL)
dataSDM$set('private', 'biasData', list())

dataSDM$set('private', 'modelData', list())
dataSDM$set('private', 'blockedCV', FALSE)
dataSDM$set('private', 'Formulas', list())
dataSDM$set('private', 'Family', list())
dataSDM$set('private', 'speciesIndex', list())

dataSDM$set('private', 'spatcovsObj', NULL)
dataSDM$set('private', 'spatcovsNames', NULL)
dataSDM$set('private', 'spatcovsEnv', NULL)
dataSDM$set('private', 'spatcovsClass', NULL)
dataSDM$set('private', 'dataSource', NULL)

dataSDM$set('private', 'Spatial', 'shared')
dataSDM$set('private', 'marksSpatial', TRUE)
dataSDM$set('private', 'Intercepts', TRUE)
dataSDM$set('private', 'marksIntercepts', TRUE)
dataSDM$set('private', 'IPS', NULL)
dataSDM$set('private', 'multinomVars', NULL)
dataSDM$set('private', 'printSummary', NULL)
dataSDM$set('private', 'multinomIndex', list())
dataSDM$set('private', 'optionsINLA', list())

dataSDM$set('private', 'spatialBlockCall', NULL)
dataSDM$set('private', 'Samplers', list())
dataSDM$set('private', 'copyModel', NULL)
dataSDM$set('private', 'datasetNames', NULL)

#' @description Initialize function for dataSDM: used to store some compulsory arguments. Please refer to the wrapper function, \code{intModel} for creating new dataSDM objects.
#' @param coordinates A vector of length 2 containing the names of the coordinates.
#' @param projection The projection of the data.
#' @param Inlamesh An inla.mesh object.
#' @param initialnames The names of the datasets if data is passed through intModel.
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
#' @param marksintercept Logical argument describing if the marks should have intercepts.
#' @param boundary A polygon map of the study area.
#' @param ips Integration points and their respective weights to be used in the model.
#' @param temporal Name of the temporal variable used in the model.
#' @param offset Name of the offset column in the datasets.
#' @param copymodel List of the specifications for the hyper parameters for the \code{"copy"} model.

dataSDM$set('public', 'initialize',  function(coordinates, projection, Inlamesh, initialnames,
                                              responsecounts, responsepa, 
                                              marksnames, marksfamily, pointcovariates,
                                              trialspa, trialsmarks, speciesname, marksspatial,
                                              spatial, intercepts, spatialcovariates, marksintercepts,
                                              boundary, ips, temporal, temporalmodel, speciesspatial, offset,
                                              copymodel) {
  
  if (missing(coordinates)) stop('Coordinates need to be given.')
  if (missing(projection)) stop('projection needs to be given.')
  if (missing(Inlamesh)) stop('Mesh needs to be given.')
  
  if (!inherits(Inlamesh, 'inla.mesh')) stop('Mesh needs to be an inla.mesh object.')
  
  if (inherits(projection, 'CRS')) projection <- as(projection, 'character')
  else if (!inherits(projection, 'character')) stop('Projection needs to be a character object.')
  
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
  
  if (!is.null(copymodel)) private$copyModel <- copymodel
  
  if (!missing(initialnames)) private$initialnames <- initialnames
    
  if (!missing(boundary)) private$Boundary <- boundary
  
  if (!missing(offset)) private$Offset <- offset
  
  private$markNames <- marksnames
  private$markFamily <- marksfamily
  
  private$pointCovariates <- pointcovariates
  
  if (!is.null(spatialcovariates)) private$spatialCovariates(spatialcovariates)
  
  if (!is.null(ips)) private$IPS <- ips
  else {
    
    if (is.null(boundary)) private$IPS <- st_transform(inlabru::fm_int(samplers = boundary, domain = Inlamesh), projection)
    else private$IPS <- st_transform(inlabru::fm_int(domain = Inlamesh), projection)
    
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

dataSDM$set('private', 'polyfromMesh', function(...) {
  
  loc <- private$INLAmesh$loc
  segm <- private$INLAmesh$segm$int
  
  coords <- na.omit(data.frame(loc[t(cbind(segm$idx[,, drop=FALSE], NA)), 1],
                               loc[t(cbind(segm$idx[,, drop=FALSE], NA)), 2]))
  
  #Polys <- sp::Polygon(coords = coords)
  #Polys <- sp::Polygons(srl = list(Polys), ID = 'id')
  #SpatPolys <- sp::SpatialPolygons(list(Polys), proj4string = private$Projection)
  SpatPolys <- st_sfc(st_polygon(list(cbind(coords[,1], coords[,2]))), crs = private$Projection)
  SpatPolys
  
})   

#' @description Function used to add or change spatial covariates in the model.
#' @param spatialCovariates A SpatialPixelsDataFrame or Raster* object describing covariates at each spatial point.
#' 
#' @import methods 
#' @importFrom raster raster

dataSDM$set('private', 'spatialCovariates', function(spatialCovariates) {
  
  if (missing(spatialCovariates)) stop('Please add spatialCovariates as a Raster* or SpatialPixelsDataFrame object.')
  
  objName <- as.character(match.call())[2]
  
  if (!objName %in% names(parent.frame())) {
    
    if (!objName %in% names(globalenv())) stop('Spatial covariates object not available.')
    else spatcovsEnv <- globalenv()
    
  } 
  else spatcovsEnv <- parent.frame()
  
  if (!class(spatialCovariates) %in% c('RasterLayer', 'RasterBrick',
                                       'RasterStack', 'SpatRaster',
                                       'SpatialPixelsDataFrame')) stop('The spatial Covariates need to be a Raster* object or a SpatialPixelsDataFrame.')
  
  spatcovsIncl <- names(spatialCovariates)
  
  #if (class(spatialCovariates) %in% c('RasterLayer', 'RasterBrick', 'RasterStack')) objSpat <- terra::rast(spatialCovariates)
  
  if (inherits(spatialCovariates, 'Spatial')) covsClass <- sapply(spatialCovariates@data, class)
  else if (inherits(spatialCovariates, 'SpatRaster')) covsClass <- sapply(as.data.frame(spatialCovariates), class)
  else covsClass <- sapply(as.data.frame(terra::rast(spatialCovariates)), class)
  
  
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
#' @param datasetName Use an existing dataset already in the model to account for preferential sampling. If \code{missing}, then \code{Data} needs to be given.
#' @param Samplers Sampling locations of the species for the structured data. May come as a: \code{SpatialPolygons}, \code{SpatialLines} or \code{SpatialPoints} object. If \code{missing}, will assume the sampling locations as the locations given in the points specified with \code{datasetName}.
#' 

dataSDM$set('public', 'samplingBias', function(datasetName, Samplers) {

  stop('Not run for now.')
  
  if (!any(datasetName %in% private$dataSource)) stop('datasetName provided in the model. If this is new data, please add it using the `.addData()` function.')
  
  if (!missing(Samplers)) {
    
    if (!inherits(Samplers, 'Spatial')) stop('Data needs to be either a SpatialPoints*, SpatialLines*, or SpatialPolygons* object. ')
    if (Samplers@proj4string != private$Projection) {
      
      message('Changing the coordinate reference system to the one specified in `intModel()`.')
      Samplers@proj4string <- private$Projection
      ##Should we check that these points are contained over the mesh area?
      
    }
    
    private$biasData[[datasetName]] <- samplers

  }
  else private$biasData[[datasetName]] <- do.call(sp::rbind.SpatialPointsDataFrame, private$modelData[[datasetName]])
    
  self$changeComponents(addComponent = paste0(datasetName, '_samplers_field(main = coordinates, copy = "shared_spatial", fixed = FALSE)'), print = FALSE)
  self$changeComponents(addComponent = paste0(datasetName,'_samplers(1)'), print = FALSE)
  
  if (private$blockedCV) {
    stop("For now don't combine blockedCV and samplers...")
    message('Re-running `.$spatialBlock()` to incorporate the new data added to the model.')
    eval(parse(text = private$spatialBlockCall))
    
  }
  
  #Add something like this::
   #private$samplersSource <- c(private$samplersSource, datasetname) # and then attach the others all together...
  
  ##Things to do here:
   #When doing spatial blocking: also spatially block the sampling locations:
   #attach these likelihoods to the ones with the points
   #attach all the other relevant metadata so that loo and all that can work (ie all the sources...)
   #Add these fields to the predict + plot functions...
  
})

## Need to change all the spatialFields to self$spatialFields and then the relevent sublist?
dataSDM$set('public', 'spatialFields', list(sharedField = list(),
                                            datasetFields = list(),
                                            speciesFields = list(),
                                            markFields = list(),
                                            biasFields = list()))
