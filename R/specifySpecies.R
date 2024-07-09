#' @title R6 class for creating a \code{startSpecies} object.
#' @description A data object containing the data and the relevant information about the integrated model. The function \code{\link{startSpecies}} acts as a wrapper in creating one of these objects. The output of this object has additional functions within the object which allow for further specification and customization of the integrated model.
#' @export
#' @importFrom R6 R6Class
#' 

specifySpecies <- R6::R6Class(classname = 'specifySpecies', lock_objects = FALSE, cloneable = FALSE, public = list(
  
  #' @description Function to provide documentation for a \code{specifySpecies} object.
  #' @param ... Not used
  #' @return Documentation.

  help = function(...) {
    
    message('Documentation for specifyISDM:')
    ?PointedSDMs:::specifySpecies 
    
  }
  ,
  #' @description Prints the datasets, their data type and the number of observations, as well as the marks and their respective families.
  #' @param ... Not used.
  #' @import stats
  print = function(...) {
    
    if (length(private$modelData) == 0) cat('No data found. Please add data using the `$.addData` function.')
    else {
      cat('Summary of startSpecies data file:\n\n')
      ## Find present absence dataset
      if (length(names(private$printSummary$Type)[private$printSummary$Type == 'Present absence']) > 0) {
        cat('Summary of presence absence datasets:\n\n')
        dataIn <- data.frame(c('-----', names(private$printSummary$Type)[private$printSummary$Type == 'Present absence']),
                             c('',rep('|  ---  |', length(private$printSummary$Type[private$printSummary$Type == 'Present absence']))),
                             c('------------------', private$printSummary$numObs[private$printSummary$Type == 'Present absence']))
        names(dataIn) <- c('Name:','','# of observations:')
        print.data.frame(dataIn[,1:3], row.names = FALSE, right = FALSE)
        cat('\n')
        
      }
      
      if (length(names(private$printSummary$Type)[private$printSummary$Type == 'Present only']) > 0) {
        cat('Summary of presence only datasets:\n\n')
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
      
      cat('Species present in the model:\n')
      cat(paste0(unique(unlist(private$speciesIn)), collapse = '\n'))
      cat('\n\n')
      cat('Use .$help() to find documentation for the slot functions.')
      
    }
    
  }
  ,
  #' @description Makes a plot of the points surrounded by the boundary of the region where they were collected.
  #' @param datasetNames Name of the datasets to plot. If this argument is missing, the function will plot all the data available to the model.
  #' @param Species Logical: should the points be plotted based on the species name. Defaults to \code{TRUE}.
  #' @param Boundary Logical: should a boundary (created using the \code{Mesh} object) be used in the plot. Defaults to \code{TRUE}.
  #' @param ... Not used.
  #' @return A ggplot object.
  #' @import ggplot2
  #' @examples
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                                speciesName = 'speciesName',
  #'                                Projection = proj, 
  #'                                responsePA = 'Present')
  #'  
  #'   #Create plot of data
  #'   organizedData$plot()
  #' 
  #' }
  #' }
  
  plot = function(datasetNames, Species = TRUE, 
                  Boundary = TRUE, ...) {
    ##REDO THIS
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
  #' @description Function used to add additional spatial fields (called \emph{bias fields}) to a selected dataset present in the integrated model. \emph{Bias fields} are typically used to account for sampling biases in opportunistic citizen science data in the absence of any covariate to do such.
  #' @param datasetNames A vector of dataset names (class \code{character}) for which a bias field needs to be added to. If \code{NULL} (default), then \code{allPO} has to be \code{TRUE}.
  #' @param allPO Logical: should a bias field be added to all datasets classified as presence only in the integrated model. Defaults to \code{FALSE}.
  #' @param biasField An \code{inla.spde} object used to describe the bias field. Defaults to \code{NULL} which uses \code{\link[INLA]{inla.spde2.matern}} to create a Matern model for the field.
  #' @param shareModel Share a bias field across the datasets specified with \code{datasetNames}. Defaults to \code{FALSE}.
  #' @param copyModel Create copy models for all the of the datasets specified with either \code{datasetNames} or \code{allPO}. The first dataset in the vector will have its own spatial effect, and the other datasets will "copy" the effect with shared hyperparameters. Defaults to \code{TRUE}.
  #' @param temporalModel List of model specifications given to the control.group argument in the time effect component. Defaults to \code{list(model = 'ar1')}; see \code{\link[INLA]{control.group}} from the \pkg{INLA} package for more details. \code{temporalName} needs to be specified in \code{intModel} prior.
  #' @return A bias field to the model.
  #' @examples
  #'  \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                                speciesName = 'speciesName',
  #'                                Projection = proj, 
  #'                                responsePA = 'Present')
  #'  
  #' #Add bias field to eBird records
  #' organizedData$addBias(datasetNames = 'eBird')
  #' 
  #' }
  #' }
  addBias = function(datasetNames = NULL,
                     allPO = FALSE,
                     biasField = NULL,
                     copyModel = TRUE,
                     shareModel = FALSE,
                     temporalModel = list(model = 'ar1')) {
    
    if (allPO) datasetNames <- names(private$printSummary$Type)[private$printSummary$Type == 'Present only']
    else
      if (is.null(datasetNames)) stop('Dataset names need to be given.')
    
    if (!all(datasetNames %in% private$dataSource)) stop('Dataset provided not available.')
    
    if (copyModel && length(datasetNames) < 2) {
      
      message('Turning copyModel off since the number of datasets specified is less than 2.')
      copyModel <- FALSE
      
    }
    
    if (copyModel && shareModel) stop('Only one of copyModel and shareModel may be TRUE.')
    
    
    for (dat in datasetNames) {
      
      #index <- which(private$dataSource == dat)
      
      for (lik in 1:length(private$Formulas[[dat]])) {
        
        if (!shareModel) private$Formulas[[dat]][[lik]][[1]]$RHS <- c(private$Formulas[[dat]][[lik]][[1]]$RHS, paste0(dat, '_biasField'))
        else private$Formulas[[dat]][[lik]][[1]]$RHS <- c(private$Formulas[[dat]][[lik]][[1]]$RHS, paste0('sharedBias', '_biasField'))
        
      }
      
      
      if (!shareModel) {
        if (is.null(biasField)) self$spatialFields$biasFields[[dat]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
        else self$spatialFields$biasFields[[dat]] <- biasField
      }
      else {
        
        if (is.null(biasField)) self$spatialFields$biasFields[['sharedBias']] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
        else self$spatialFields$biasFields[['sharedBias']] <- biasField        
        
      }
      
    }
    
    
    
    ##Add a way for separate Fields per dataset.
    #Shouldn't be too difficult -- have already done it before.
    
    
    #Should I copy the bias fields for the marks?
    
    if (!is.null(private$temporalName)) {
      
      if (shareModel) private$Components <- c(private$Components, paste0('sharedBias_biasField(main = geometry, model = sharedBias_bias_field, group =', private$temporalName,', ngroup = ', length(unique(unlist(private$temporalVars))),', control.group = ', temporalModel,')'))
      
      temporalModel <- deparse1(temporalModel)
      
      private$Components <- c(private$Components, paste0(datasetNames ,'_biasField(main = geometry, model = ', datasetNames, '_bias_field, group = ', private$temporalName, ', ngroup = ', length(unique(unlist(private$temporalVars))),', control.group = ', temporalModel,')'))
      
      
    }
    else {
      
      if (!copyModel) {
        
        if (shareModel) private$Components <- c(private$Components, 'sharedBias_biasField(main = geometry, model = sharedBias_bias_field)')
        else private$Components <- c(private$Components, paste0(datasetNames ,'_biasField(main = geometry, model = ', datasetNames, '_bias_field)'))
        
      }
      
      else {
        
        biasComponents <- c() 
        
        for (bias in datasetNames) {
          
          if (match(bias, datasetNames) == 1)  biasComponents[bias] <- paste0(bias, '_biasField(main = geometry, model = ', bias
                                                                              , '_bias_field)')
          else biasComponents[bias] <- paste0(bias, '_biasField(main = geometry, copy = \"', datasetNames[1], '_biasField\", hyper = list(beta = list(fixed = FALSE)))')
          
          private$biasCopy <- TRUE
          
          
        }
        
        private$Components <- c(private$Components, biasComponents)
        
      }
      
    }
    
    
    ##Things to do here:
    #Go into the liks of PO datasets and add the biasfield
    #Go inth the components and add the bias field component
    #Then store the bias field somewhere else. If field NULL, then add generic bias field
    
    
  }
  ,
  #' @description Function used to update the formula for a selected observation model. The function is designed to work similarly to the generic \code{update} formula, and should be used to thin terms out of a process from the full model specified in \code{\link{intModel}}. The function also allows the user to add their own formula to the model, such that they can include non-linear components in the model. The function can also be used to print out the formula for a process by not specifying the \code{Formula} or \code{newFormula} arguments.
  #' @param datasetName Name of the dataset (class \code{character}) for which the formula needs to be changed.
  #' @param speciesName Name of the species for which to change a formula for. Defaults to \code{NULL} which chnages the formula for all species present in \code{datasetName}.
  #' @param processLevel Logical argument: if \code{TRUE} changes the formulas for all of the processes in a dataset. Defaults to \code{FALSE}.
  #' @param Formula An updated formula to give to the process. The syntax provided for the formula in this argument should be identical to the formula specification as in base \strong{R}. Should be used to thin terms out of a formula but could be used to add terms as well. If adding new terms not specified in \code{intModel}, remember to add the associated component using \code{.$changeComponents} as well.
  #' @param newFormula Completely change the formula for a process -- primarily used to add non-linear components into the formula. Note: all terms need to be correctly specified here.
  #' @return An updated formula.
  #' @import methods
  #' @import stats
  #' 
  #' @examples
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, speciesName = 'speciesName',
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present',
  #'                            pointsSpatial = 'individual')
  #' 
  #'  #Remove Forest from eBird
  #'  organizedData$updateFormula(datasetName = 'eBird', Formula = ~ . - Forest)
  #'  
  #'  #Add some scaling to Forest for Parks
  #'  organizedData$updateFormula(datasetName ='Parks', newFormula = ~ I(. +(Forest+1e-6)*scaling))
  #'  
  #'  #Now dd scaling to components
  #'  organizedData$changeComponents(addComponent = 'scaling') 
  #'  
  #' }
  #' }
  #' @return If \code{Formula} and \code{newFormula} are missing, will print out the formula for the specified processes. 
  #' 
  #' 
  updateFormula = function(datasetName = NULL, 
                           speciesName = NULL,
                           Formula, 
                           processLevel = FALSE,
                           newFormula) {
    
    if (all(is.null(datasetName), is.null(speciesName)) & !processLevel) stop ('At least one of: datasetName, speciesName needs to be specified.')
    
    #if (length(speciesName) > 1) stop('Please only supply one species name at a time.')
    
    if (processLevel) {
      
      datasetName <- unique(private$dataSource)
      speciesName <-  unique(unlist(private$speciesIn))
      
    }
    else {
    
    if (is.null(datasetName)) datasetName <- names(private$speciesIn)[sapply(private$speciesIn, FUN = function(x) speciesName %in% x)] #Datasets where all species occur
    
    if (is.null(speciesName)) speciesName <- unlist(private$speciesIn[datasetName])
    
    }
    #if (length(datasetName) != 1) stop ('Please only provide one dataset name.')
    
    if (!all(datasetName %in% private$dataSource)) stop ('Dataset name provided not in model.')
    
    if (!missing(Formula) && !missing(newFormula)) stop ('Please provide only one of Formula and newFormula. \n Use Formula to update the current formula; use newFormula to make a completely new formula.')
    
    if (!missing(Formula) && !inherits(Formula, 'formula')) stop ('Formula must be of class "formula".')
    
    if (!missing(newFormula) && !inherits(newFormula, 'formula')) stop('newFormula must be of class "formula".')
    
    if (!missing(Formula) && length(as.character(Formula)) == 3) stop ("Please remove the response variable of the formula.")
        
    if (!all(sapply(private$speciesIn[datasetName], function(x) any(speciesName %in% x)))) stop('Species provided not available in the dataset.') 
    else name_index <- speciesName
          
    if (missing(Formula) && missing(newFormula)) {
      
      get_formulas <- list()
      for (i in name_index) {
        
        get_formulas[[i]] <- lapply(private$Formulas[[datasetName]][[i]], function(x) list(formula = x$LHS,
                                                                                          components = x$RHS))
        
      }
      
      get_formulas <- lapply(unlist(get_formulas, recursive = FALSE), function(x) {
        
        if (as.character(x$formula)[3] == '.') update(x$formula, paste('~', paste0(x$components, collapse = '+'))) 
        else x$formula
        
      })
      
      names(get_formulas) <- name_index
      return(get_formulas)
      
    }
    else
      if (!is.null(private$covariateFormula)) {
        
        if (!missing(Formula)) {
          if (!private$speciesEnvironment) newForm <- makeFormulaComps(form = update(private$covariateFormula, Formula), species = FALSE, speciesnames = speciesName, type = 'cov')
          else newForm <- makeFormulaComps(form = update(private$covariateFormula, Formula), species = TRUE, speciesnames = speciesName, type = 'cov')
            
          private$covariateFormula <- update(private$covariateFormula, Formula)
          
        }
        else {
          
          if (!private$speciesEnvironment) newForm <- makeFormulaComps(form = newFormula, species = FALSE, speciesnames = speciesName, type = 'cov')
          else newForm <- makeFormulaComps(form = newFormula, species = TRUE, speciesnames = speciesName, type = 'cov')          
          
          private$covariateFormula <- newFormula
          
        }
        
        for (spec in 1:length(newForm)) {
        
        self$changeComponents(addComponent = newForm[spec], print = FALSE)
        
          }
      }
    else
      if (!missing(Formula)) {
        
        for (dataset in datasetName) {
          
          for (species in name_index)  {
            
            if (private$speciesEnvironment) {
                
                varsIn <- all.vars(Formula)
                varsIn <- varsIn[varsIn %in% c(private$spatcovsNames, paste0(species, '_', dataset, '_spatial'))]
                
                if (sum(paste0(species, '_', c(private$spatcovsNames, paste0(dataset, '_spatial'))) %in% varsIn) == 0) varsIn <- paste0(species, '_', varsIn)
   
                Formula2 <- as.formula(paste0('~ . -', paste0(varsIn, collapse = '+')))
            
              
            } else Formula2 <- Formula
            
            if (!is.null(private$Formulas[[dataset]][[1]][[1]]$RHS)) {
              
              if (!is.null(private$Formulas[[dataset]][[species]][[1]]$RHS)) {
          private$Formulas[[dataset]][[species]][[1]]$RHS <- removeFormula(formulaRemove =  Formula2,
                                                                     oldFormula = private$Formulas[[dataset]][[species]][[1]]$RHS)
              }
          
            }
          
          }
        }
        
      }
    else {
      
      processIn <- names(private$Formulas[[datasetName]][[1]])
      
      for (dataset in name_index) { 
        
        for (process in processIn) {
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
  #' @description Function to add and specify custom components to model, which are required by \pkg{inlabru}. The main purpose of the function is to re-specify or completely change components already in the model, however the user can also add completely new components to the model as well. In this case, the components need to be added to the correct formulas in the model using the \code{.$updateFormula} function. If \code{addComponent} and \code{removeComponent} are both missing, the function will print out the components to be supplied to \pkg{inlabru}'s \code{\link[inlabru]{bru}} function.
  #' @param addComponent Component to add to the integrated model. Note that if the user is re-specifying a component already present in the model, they do not need to remove the old component using \code{removeComponent}.
  #' @param removeComponent Component (or just the name of a component) present in the model which should be removed.
  #' @param print Logical: should the updated components be printed. Defaults to \code{TRUE}.
  #' @return An updated components list.
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
  #'  Forest <- terra::rast(
  #'  system.file(
  #'  'extdata/SolitaryTinamouCovariates.tif', 
  #'  package = "PointedSDMs"))$Forest
  #'  
  #'  
  #'  #Set model up
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                                speciesName = 'speciesName',
  #'                                spatialCovariates = Forest,
  #'                                Projection = proj, 
  #'                                responsePA = 'Present')
  #' 
  #'  #Remove Forest from components
  #'  organizedData$changeComponents(removeComponent = 'speciesSpatial')
  #'
  #' }
  #' 
  #' }
  changeComponents = function(addComponent, removeComponent, print = TRUE) {
    
    terms <- gsub('\\(.*$', '', private$Components)
    
    if (!missing(addComponent)) {
      
      if (inherits(addComponent, 'character')) addComponent <- formula(paste('~', addComponent))
      
      compTerms <- attr(terms(addComponent), 'term.labels')
      
      for (toAdd in 1:length(compTerms)) {
        
        if (any(gsub('\\(.*$', '', compTerms[toAdd]) %in% terms)) private$Components <- private$Components[! terms %in% gsub('\\(.*$', '', compTerms[toAdd])]
        
        private$Components <- c(private$Components, compTerms[toAdd])
        
      }
      
    } 
    
    if (!missing(removeComponent)) private$Components <- private$Components[!terms%in%removeComponent]
    
    componentsJoint <- formula(paste('~ - 1 +', paste(private$Components, collapse = ' + ')))
    componentsJoint <- formula(paste(paste('~ - 1 +', paste(labels(terms(componentsJoint)), collapse = ' + '))))
    
    if (print) {
      
      cat('Model components:')
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
  #' @return New priors for the fixed effects.
  #' @examples
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                            speciesName = 'speciesName',
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present',
  #'                            pointsSpatial = 'individual')
  #' 
  #'  #Add prior to Forest
  #'  organizedData$priorsFixed(Effect = 'Forest', mean.linear = 2, prec.linear = 0.1)
  #'
  #' }
  #' }
  
  priorsFixed = function(Effect, Species = NULL, datasetName = NULL,
                         mean.linear = 0, prec.linear = 0.001) {
    
  ##Add IID prior here?  
    if (!missing(Effect)) {
      
      if (Effect %in% c('intercept', 'Intercept')) {
        
        intTRUE <- TRUE
        
        if (is.null(Species)) {
          
          if (!private$Intercepts) stop('Fixed effect is given as "intercept", but intercepts have been turned off in intModel.')
          
          if (is.null(datasetName)) datasetName <- unique(private$dataSource)
          else if (!datasetName %in% unique(private$dataSource)) stop('datasetName is not the name of a dataset added to the model.')
          
          
          Effect <- paste0(datasetName,'_intercept')
          
        } 
        else {
          
            
            if (!is.null(private$speciesIntercepts)) {
              if (!private$speciesIntercepts) Effect <- paste0(Species, '_intercept')
              else stop('Species intercepts are random effects in the model.')
              
            } else stop('Species intercepts are not in the model.')

        }
        
      }
      else {
        
        intTRUE <- FALSE
        
        if (!Effect %in% c(private$spatcovsNames, private$pointCovariates)) stop('Fixed effect provided not present in the model. Please add covariates using the "spatialCovariates" or "pointCovariates" argument in intModel.')
        
        if (any(Effect %in% private$spatcovsNames)) cov_class <- private$spatcovsClass[Effect]
        else cov_class <- 'linear'
        
        if (private$speciesEnvironment) {
        
        if (!is.null(Species)) {
            
            if (!Species %in% unlist(private$speciesIn)) stop('Species given is not available in the model.')
            
            Effect <- paste0(Species, '_', Effect) #this won't work, unless we run a for loop...
            
          }
          else Effect <- paste0(unique(unlist(private$speciesIn)), '_', Effect)
          
        }

      }
      if (intTRUE) {
        
        newComponent <- c()
        
        for (eff in Effect) {
          
          newComponent[eff] <- paste0(eff, '(1, mean.linear = ', mean.linear, ', prec.linear = ', prec.linear,' )') 
          
        }
        
      }
      else newComponent <- paste0(Effect,'(main = ', Effect, ', model = \"', cov_class, '\", mean.linear = ', mean.linear, ', prec.linear = ', prec.linear, ')')
      
      for (comp in newComponent) {
        
        self$changeComponents(addComponent = comp, print = FALSE)
        
      }
      
    }
    
  }
  ,
  #' @description Function to specify random fields in the model using penalizing complexity (PC) priors for the parameters.
  #' 
  #' @param sharedSpatial Logical: specify the shared spatial field in the model. Requires \code{pointsSpatial == 'shared'} in \code{\link{intModel}}. Defaults to \code{FALSE}.
  #' @param datasetName Name of which of the datasets' spatial fields to be specified. Requires \code{pointsSpatial = 'individual'} in \code{\link{intModel}}.
  #' @param Species Name of the species to change the spatial effect for. If \code{TRUE} then changes the spatial effect for the shared species field.
  #' @param Bias Name of the dataset for which the bias field to be specified.
  #' @param PC Logical: should the Matern model be specified with pc priors. Defaults to \code{TRUE}, which uses \code{\link[INLA]{inla.spde2.pcmatern}} to specify the model; otherwise uses \code{\link[INLA]{inla.spde2.matern}}.
  #' @param Remove Logical: should the chosen spatial field be removed. Requires one of \code{sharedSpatial}, \code{species}, \code{mark} or \code{bias} to be non-missing, which chooses which field to remove.
  #' @param ... Additional arguments used by \pkg{INLA}'s \code{\link[INLA]{inla.spde2.pcmatern}} or \code{\link[INLA]{inla.spde2.matern}} function, dependent on the value of \code{PC}.
  #' @return A new model for the spatial effects.
  #' 
  #' @examples 
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                            speciesName = 'speciesName',
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present')
  #' 
  #'  #Specify the shared spatial field
  #'  organizedData$specifySpatial(sharedSpatial = TRUE, PC = TRUE, 
  #'                        prior.range = c(1,0.001),
  #'                        prior.sigma = c(1,0.001))
  #'
  #' } 
  #' }
  
  specifySpatial = function(sharedSpatial = FALSE,
                            datasetName,
                            Species,
                            Bias, PC = TRUE,
                            Remove = FALSE, ...) {
    
    if (all(!sharedSpatial && missing(datasetName) && missing(Species)  &&  missing(Bias))) stop('At least one of sharedSpatial, datasetName, dataset or Species needs to be provided.')
    
    if (sum(sharedSpatial, !missing(datasetName), !missing(Species), !missing(Bias)) != 1) stop('Please only choose one of sharedSpatial, datasetName, Species or Bias.')
    
    if (Remove && sum(sharedSpatial, !missing(datasetName), !missing(Species), !missing(Bias)) != 1) stop('Please choose one of sharedSpatial, datasetName, Species or Bias to remove.')
    
    #if (!is.null(Copy) && !Copy %in% unlist(private$dataSource)) stop('Dataset name provided is not currently in the model.')
    
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
        if (!Remove) index <- private$originalNames[[1]]
        else index <- private$originalNames[[1]]
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
      
      if (!inherits(Species, 'character')) {
      
      if (Species) Species <- unlist(private$speciesIn)
      
      }
      else
       if (!Species %in% unlist(private$speciesIn)) stop('Species name provided is not currently in the model.')

      spWhere <- sapply(private$speciesIn, function(x) Species %in% x)
      
      field_type <- 'speciesFields'
      if (!Remove) {

        if (private$speciesSpatial == 'shared' || private$speciesSpatial == 'replicate') index <- 'speciesField'
        else
          if (!private$speciesIndependent) index <- unlist(mapply(FUN = function(x, names) paste0(x,'_', names), x = private$speciesIn[spWhere], names = names(private$speciesIn)[spWhere]))#index <- do.call(paste0, expand.grid(paste0(Species, '_'), private$dataSource))
          else index <- Species
          ##Change to paste0() species_dataset_spatia,
      }
      else {
        
        if (private$speciesSpatial == 'shared') index <- 'speciesField'
        else
          if (!private$speciesIndependent) index <- unlist(mapply(FUN = function(x, names) paste0(x,'_', names, '_spatial'), x = private$speciesIn[spWhere], names = names(private$speciesIn)[spWhere]))#index <- do.call(paste0, expand.grid(paste0(Species, '_'), private$dataSource))
          else index <-  paste0(Species, '_spatial')
          
      }
      
    }
    
    if (!missing(Bias)) {
      
      if (!inherits(Bias, 'character')) {
        
      if (Bias) Bias <- names(self$spatialFields$biasFields)
      
      }
      else 
      if (!Bias %in% names(self$spatialFields$biasFields)) stop('Dataset name provided does not have a bias field. Please use ".$biasField()" beforehand.')
      
      field_type <- 'biasFields'
      if (!Remove) index <- Bias
      else index <- paste0(Bias, 'biasField')
      
      
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
          
          self$updateFormula(datasetName = data, Formula = formula(paste(' ~ . -', term)))  
          
        } 
      }
      
    }
    
  }
  ,
  #' @description Function used to change the link function for a given process.
  #' @param datasetName Name of the dataset for which the link function needs to be changed.
  #' @param Link Name of the link function to add to the process. If missing, will print the link function of the specified dataset.
  #' @return A new link function for a process.
  #' @examples
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                            speciesName = 'speciesName',
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present')
  #' 
  #'  #Specify the shared spatial field
  #'  organizedData$changeLink('Parks', 'logit')
  #'  
  #'  
  #' } 
  #' }
  changeLink = function(datasetName,
                        Link, ...) {
    
    if (missing(datasetName)) stop('Please provide a dataset name.')
    
    if (length(datasetName) > 1) stop('Please only provide one dataset at a time.')
    
    if (!datasetName %in% private$dataSource) stop('Dataset name provided not in model.')
    
    for (i in which(private$dataSource == datasetName)) {
    
    private$optionsINLA[['control.family']][[i]] <- list(link = Link)
    
    }
    
  }
  ,
  #' @description Function to spatially block the datasets, which will then be used for model cross-validation with \code{\link{blockedCV}}. See the \code{\link[blockCV]{spatialBlock}} function from \pkg{blockCV} for how the spatial blocking works and for further details on the function's arguments.
  #' @param k Integer value reflecting the number of folds to use.
  #' @param rows_cols Integer value by which the area is divided into longitudinal and latitudinal bins.
  #' @param plot Plot the cross-validation folds as well as the points across the boundary. Defaults to \code{FALSE}.
  #' @param seed Seed used by \pkg{blockCV}'s \code{\link[blockCV]{spatialBlock}} to make the spatial blocking reproducible across different models. Defaults to \code{1234}. 
  #' @param ... Additional arguments used by \pkg{blockCV}'s \code{\link[blockCV]{spatialBlock}}.
  #' @return If \code{plot = TRUE}, a plot of the grid. 
  #' @import ggplot2
  #' @importFrom R.devices suppressGraphics
  #' @importFrom blockCV spatialBlock
  #' @importFrom blockCV cv_plot
  #' 
  #' @examples
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                            speciesName = 'speciesName',
  #'                            spatialCovariates = Forest,
  #'                            Projection = proj, responsePA = 'Present',
  #'                            pointsSpatial = 'individual')
  #' 
  #'  #Specify the spatial block
  #'  organizedData$spatialBlock(k = 2, rows = 2, cols = 1, plot = FALSE)
  #'
  #' } 
  #' }
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
    
    sumIndex <- 0
    
    for (i in 1:length(private$modelData)) {
      
      if (sumIndex == 0) sumIndex <- 0
      
      if (i != length(dataLength)) {
        
        for (j in 1:length(private$modelData[[i]])) {
          
          if (sumIndex == 0) start <- 0
          else start <- cumsum(dataLength)[sumIndex]
          
          sumIndex <- sumIndex + 1
          
          private$modelData[[i]][[j]]$.__block_index__ <- as.character(foldID[(1 + start):cumsum(dataLength)[sumIndex]])
          
        }
        
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
  #' @return New samplers for a process.
  #' @example
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh, 
  #'                               speciesName = 'speciesName',
  #'                               Projection = proj, responsePA = 'Present')
  #'  
  #' #Add integration domain for the eBird records
  #' organizedData$addSamplers(datasetName = 'eBird', Samplers = SolitaryTinamou$region)
  #' 
  #' }
  #' }
  addSamplers = function(datasetName, Samplers) {
    
    if (!datasetName %in% private$dataSource) stop ('Dataset name provided not in model.')
    
    if (!inherits(Samplers, 'sf')) stop ('Samplers needs to be a sf object.')
    
    Samplers <- st_transform(Samplers, private$Projection)
    
    private$Samplers[[datasetName]] <- Samplers
    
  },
  
  #' @description Function to specify the models and priors for the random effects included in the model.
  #' @param temporalModel List of model specifications given to the control.group argument in the time effect component. Defaults to \code{list(model = 'ar1')}; see \code{\link[INLA]{control.group}} from the \pkg{INLA} package for more details.
  #' @param copyModel List of model specifications given to the hyper parameters for the \code{"copy"} model. Defaults to \code{list(beta = list(fixed = FALSE))}.
  #' @param copyBias List of model specifications given to the hyper parameters for the \code{"copy"} bias model. Defaults to \code{list(beta = list(fixed = FALSE))}.
  #' @param speciesCopy List of model specifications given to the hyper parameters for the species  \code{"copy"} model. Defaults to \code{list(beta = list(fixed = FALSE))}.
  #' @param speciesIntercepts Prior distribution for precision parameter for the random species intercept term. Defaults to \code{INLA}'s default choice.
  #' @param speciesGroup Prior distribution for the precision parameter for the iid group model. Defaults to \code{INLA}'s default.
  #'   #' @return An updated component list. 
  #' @example
  #' \dontrun{
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
  #'  organizedData <- startSpecies(data, Mesh = mesh,
  #'                            speciesName = 'speciesName',
  #'                            Projection = proj, 
  #'                            responsePA = 'Present',
  #'                            pointsSpatial = copy)
  #'  
  #' #Add integration domain for the eBird records
  #' organizedData$specifyRandom(copyModel =  list(beta = list(fixed = TRUE)))
  #' 
  #' }
  #' }
  specifyRandom = function(temporalModel = list(model = 'ar1'),
                           copyModel = list(beta = list(fixed = FALSE)),
                           copyBias =  list(beta = list(fixed = FALSE)),
                           speciesCopy = list(beta = list(fixed = FALSE)),
                           speciesIntercepts = list(prior = 'loggamma', param = c(1, 5e-5)),
                           speciesGroup = list(model = "iid", hyper = list(prec = list(prior = 'loggamma', param = c(1, 5e-5))))) {
    
    if (!is.null(private$temporalName) & !missing(temporalModel)) { 
      
      if (private$Spatial == 'shared') whichTemp <- grepl('shared_spatial\\(main = geometry', private$Components)
      else whichTemp <- grepl(paste0('group = ', private$temporalName), private$Components)
      
      private$Components[whichTemp] <- gsub(gsub('list','list\\\\', private$temporalModel),
                                            deparse1(temporalModel), private$Components[whichTemp])
      
      private$temporalModel <- deparse1(temporalModel)
    }
    
    if (!is.null(private$Spatial)) {
    
    if (private$Spatial == 'copy' & !missing(copyModel)) {
      
      whichCopied <- grepl('beta = list\\(',private$Components) & !grepl('_biasField',private$Components)
      
      private$Components[whichCopied] <- gsub(gsub('list','list\\\\', private$copyModel),
                                              deparse1(copyModel), private$Components[whichCopied])
      
      private$copyModel <- deparse1(copyModel)
      
    }
    }
    
    if (private$biasCopy & !missing(copyBias)) {
      
      whichCopied <- grepl('beta = list\\(',private$Components) & grepl('_biasField',private$Components)
      
      private$Components[whichCopied] <- gsub(gsub('list','list\\\\', private$biasCopyModel),
                                              deparse1(copyBias), private$Components[whichCopied])
      
      private$biasCopyModel <- deparse1(copyBias)
      
      
    }
    
    if (!is.null(private$speciesSpatial)) {
    
    if (private$speciesSpatial == 'copy' & !missing(speciesCopy)) {
      
      whichSpecies <- grepl(paste0(unique(unlist(private$speciesIn)),'_',collapse = '|'), private$Components)
      whichCopied <- grepl('beta = list\\(',private$Components)
      
      private$Components[whichCopied & whichSpecies] <- gsub(gsub('list','list\\\\', private$speciesCopy),
                                                             deparse1(speciesCopy), private$Components[whichCopied & whichSpecies])
      
      private$speciesCopy <- deparse1(speciesCopy)
      
      
    }
    
    if (private$speciesSpatial == 'replicate' & !missing(speciesGroup)) {
      
      whichSpatial <- grepl('speciesShared\\(main = geometry', private$Components)
      
      private$Components[whichSpatial] <- gsub(gsub('list','list\\\\', private$speciesReplicate),
                                              deparse1(speciesGroup), private$Components[whichSpatial])
      
      private$speciesReplicate <- deparse1(speciesGroup)
      
      
    }
      
    }
    
    if (!is.null(private$speciesIntercepts)) {
    
    if (private$speciesIntercepts & !missing(speciesIntercepts)) {
      
      
      whichIntercepts <- grepl(paste0(private$speciesName, '_intercepts\\(main ='), private$Components)
      
      private$Components[whichIntercepts] <- gsub(gsub('c\\(','c\\\\(', 
                                                       gsub('list','list\\\\', private$speciesIntPrior)),
                                                  deparse1(speciesIntercepts), private$Components[whichIntercepts])
      
      private$speciesIntPrior <- deparse1(speciesIntercepts)
      
      
    }
      
    }

  }
  
))

specifySpecies$set('private', 'Projection', NULL)
specifySpecies$set('private', 'Coordinates', NULL)
specifySpecies$set('private', 'responseCounts', NULL)
specifySpecies$set('private', 'responsePA', NULL)
specifySpecies$set('private', 'trialsPA', NULL)

specifySpecies$set('private', 'Boundary', NULL)
specifySpecies$set('private', 'INLAmesh', NULL)
specifySpecies$set('private', 'Components', NULL)
specifySpecies$set('private', 'pointCovariates', NULL)
specifySpecies$set('private', 'ptcovsClass', NULL)
specifySpecies$set('private', 'initialnames', NULL)

specifySpecies$set('private', 'temporalName', NULL)
specifySpecies$set('private', 'temporalVars', NULL)
specifySpecies$set('private', 'temporalModel', NULL)
specifySpecies$set('private', 'Offset', NULL)
specifySpecies$set('private', 'biasData', list())

specifySpecies$set('private', 'modelData', list())
specifySpecies$set('private', 'blockedCV', FALSE)
specifySpecies$set('private', 'Formulas', list())
specifySpecies$set('private', 'Family', list())
specifySpecies$set('private', 'biasCopy', FALSE)
specifySpecies$set('private', 'biasCopyModel', deparse1(list(beta = list(fixed = FALSE))))

specifySpecies$set('private', 'speciesIntercepts', TRUE)
specifySpecies$set('private', 'speciesEnvironment', TRUE)
specifySpecies$set('private', 'speciesSpatial', TRUE)
specifySpecies$set('private', 'speciesIndex', list())
specifySpecies$set('private', 'speciesIndependent', TRUE)
specifySpecies$set('private', 'speciesName', NULL)
specifySpecies$set('private', 'speciesIn', NULL)
specifySpecies$set('private', 'speciesCopy', deparse1(list(beta = list(fixed = FALSE))))
specifySpecies$set('private', 'speciesReplicate', deparse1(list(model = "iid")))
specifySpecies$set('private', 'speciesIntPrior', deparse1(list(prior = 'loggamma', param = c(1, 5e-5))))


specifySpecies$set('private', 'spatcovsObj', NULL)
specifySpecies$set('private', 'spatcovsNames', NULL)
specifySpecies$set('private', 'spatcovsEnv', NULL)
specifySpecies$set('private', 'spatcovsClass', NULL)
specifySpecies$set('private', 'dataSource', NULL)

specifySpecies$set('private', 'Spatial', 'shared')
specifySpecies$set('private', 'Intercepts', TRUE)
specifySpecies$set('private', 'IPS', NULL)
specifySpecies$set('private', 'multinomVars', NULL)
specifySpecies$set('private', 'printSummary', NULL)
specifySpecies$set('private', 'multinomIndex', list())
specifySpecies$set('private', 'optionsINLA', list())
specifySpecies$set('private', 'covariateFormula', NULL)
specifySpecies$set('private', 'biasFormula', NULL)

specifySpecies$set('private', 'spatialBlockCall', NULL)
specifySpecies$set('private', 'Samplers', list())
specifySpecies$set('private', 'copyModel', NULL)
specifySpecies$set('private', 'datasetNames', NULL)
specifySpecies$set('private', 'originalNames', NULL)

#' @description Initialize function for specifySpecies: used to store some compulsory arguments. Please refer to the wrapper function, \code{intModel} for creating new specifySpecies objects.
#' @param data The points of the model.
#' @param projection The projection of the data.
#' @param Inlamesh An inla.mesh object.
#' @param initialnames The names of the datasets if data is passed through intModel.
#' @param responsecounts The name of the response variable for the count data.
#' @param responsepa The name of the response variable for the presence absence data.
#' @param marksnames The names of the marks contained in the data.
#' @param marksfamily The statistical family of the marks.
#' @param pointcovariates Names of the additional, non-spatial covariates describing the points.
#' @param trialspa The name of the trials variable for the presence absence datasets.
#' @param spatial Logical argument describing if spatial effects should be included.
#' @param intercepts Logical argument describing if intercepts should be included in the model.
#' @param boundary A polygon map of the study area.
#' @param ips Integration points and their respective weights to be used in the model.
#' @param temporal Name of the temporal variable used in the model.
#' @param offset Name of the offset column in the datasets.
#' @param copymodel List of the specifications for the hyper parameters for the \code{"copy"} model.
#' @param formulas Custom formulas to add to the model. List with two objects: covariateFormula and biasFormula.
#' @return An initiated object.

specifySpecies$set('public', 'initialize',  function(data, projection, Inlamesh, initialnames,
                                                  responsecounts, responsepa, pointcovariates, speciesintercept,
                                                  trialspa, spatial, intercepts, spatialcovariates,
                                                  boundary, ips, temporal, temporalmodel, offset,
                                                  copymodel, formulas, speciesindependent,
                                                  speciesname, speciesenvironment, speciesspatial) {
  
  if (missing(projection)) stop('projection needs to be given.')
  if (missing(Inlamesh)) stop('Mesh needs to be given.')
  
  if (!inherits(Inlamesh, 'inla.mesh')) stop('Mesh needs to be an inla.mesh object.')
  
  if (!inherits(projection, 'character')) stop('Projection needs to be a character object.')
  
  if (any(missing(responsecounts), missing(responsepa)) ||
      any(is.null(responsecounts), is.null(responsepa))) stop('At least one of responseCounts and responsePA are NULL. Please provide both.')
  
  private$responseCounts <- responsecounts
  private$responsePA <- responsepa
  
  if (!missing(trialspa)) private$trialsPA <- trialspa
  
  if (!missing(temporal)) private$temporalName <- temporal
  private$temporalModel <- temporalmodel
  
  if (!is.null(copymodel)) private$copyModel <- copymodel
  
  if (!missing(initialnames)) private$initialnames <- initialnames
  
  if (!missing(boundary)) private$Boundary <- boundary
  
  if (!missing(offset)) private$Offset <- offset
  
  private$pointCovariates <- pointcovariates
  
  if (!is.null(spatialcovariates)) private$spatialCovariates(spatialcovariates)
  
  if (is.null(ips)) {
    
    if (!is.null(boundary)) ips <- st_transform(inlabru::fm_int(samplers = boundary, domain = Inlamesh), projection)
    else ips <- st_transform(inlabru::fm_int(domain = Inlamesh), projection)
    
    
  }
  
  st_geometry(ips) <- 'geometry'
  
  if (!is.null(spatial)) {
    if (spatial == 'correlate') ips <- fm_cprod(ips, data = data.frame(._dataset_index_var_. = 1:length(initialnames)))
  }  
  private$IPS <- ips
  
  private$speciesName <- speciesname
  private$speciesSpatial <- speciesspatial
  private$speciesIndependent <- speciesindependent
  private$speciesEnvironment <- speciesenvironment
  private$speciesIntercepts <- speciesintercept
  
  
  private$Spatial <- spatial
  private$Intercepts <- intercepts
  
  private$covariateFormula <- formulas$covariateFormula
  private$biasFormula <- formulas$biasFormula
  
  private$Projection <- projection
  private$INLAmesh <- Inlamesh
  
  private$addData(dataList = data, responseCounts = responsecounts, responsePA = responsepa, trialsPA = trialspa,
                  markNames = NULL, markFamily = NULL, pointCovariates = pointcovariates, dataNames = initialnames,
                  trialsMarks = NULL, speciesName = speciesname, temporalName = temporal, Offset = offset)
  
  invisible(self)
  
  
})

#' @description Function used to add additional datasets to the \code{specifySpecies} object. This function should be used if the user would like to add any additional datasets to the integrated model, but do not have the same standardized variable names as those added initially with \code{\link{intModel}}. Use \code{?intModel} for a more  <- rehensive description on what each argument in this function is used for.
#' @param dataList The datasets to be added to the integrated model: should be either \code{sf}, \code{data.frame} or \code{SpatialPoints*} objects, or a list of objects from these classes.
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
#' @param Offset Name of the offset column in the dataset.
#' @param dataNames Names of the datasets.
#' @return Data added to the object.
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
#'  organizedData <- startSpecies(ebird, 
#'                  Mesh = mesh, 
#'                  speciesName = 'speciesName',
#'                  Projection = proj)
#'                              
#'  }
#'  
specifySpecies$set('private', 'addData', function(dataList, responseCounts, responsePA, trialsPA,
                                               markNames, markFamily, pointCovariates,
                                               trialsMarks, speciesName, temporalName,
                                               Coordinates, Offset, dataNames) {
  
  pointData <- dataOrganize$new()
  
  if (!is.null(private$Spatial)) {
    
    if (private$Spatial %in% c('shared', 'correlate')) self$spatialFields$sharedField[['sharedField']] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
    
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
  
  if (!is.null(Offset) && any(Offset %in% 'offset')) stop('The offset variable cannot be called "offset". Please choose a new name.')
  
  speciesOK <- checkVar(data = dataList,
                        var = private$speciesName)
  
  if (!is.null(private$speciesSpatial)) {
    
      if (private$speciesSpatial == 'replicate') repl = TRUE
      else repl = FALSE
      
    } 
    else repl = FALSE
    
  
  if (!speciesOK) stop('The species variable name is required to be present in all the datasets.')
  
  if (!is.null(private$temporalName)) {
    
    timeOK <- checkVar(data = dataList,
                       var = private$temporalName)
    
    if (!timeOK) stop('The temporal variable name is required to be present in all the datasets.')
    
  }
  
  ##REMOVED COORDINATES NEED TO REMOVE FROM DATA ORGANIZE TOO
  
  pointData$makeData(datapoints = dataList, datanames = dataNames,
                     proj = private$Projection,
                     countsresp = responseCounts, paresp = responsePA,
                     trialname = trialsPA, speciesname = speciesName,
                     marktrialname = NULL, temporalvar = private$temporalName,
                     marks = NULL, markfamily = NULL,
                     pointcovnames = pointCovariates, offsetname = private$Offset)
  
  pointData$makeSpecies(speciesname = speciesName, repl = repl) 
  private$speciesIn <- pointData$SpeciesInData
  
  #ADD argument common field for species
  if (!private$speciesIndependent) {
    
    speciesDataNames <- unlist(mapply(FUN = function(x, names) paste0(x,'_', names), x = private$speciesIn, names = names(private$speciesIn)))

    self$spatialFields$speciesFields <- vector(mode = 'list', length = length(speciesDataNames))#length(unique(unlist(private$speciesIn))) * length(names(private$speciesIn)))
    names(self$spatialFields$speciesFields) <- speciesDataNames#do.call(paste0, expand.grid(paste0(unique(unlist(private$speciesIn)),'_'), names(private$speciesIn)))
    speciesInd <-  private$speciesIn
    
  }
  else {
    
    self$spatialFields$speciesFields <-  vector(mode = 'list', length = length(unique(unlist(private$speciesIn))))
    names(self$spatialFields$speciesFields) <- unique(unlist(private$speciesIn))
    speciesInd <- do.call(paste0, expand.grid(paste0(private$speciesIn,'_'), names(private$speciesIn)))
    
  }
  
  if (!is.null(private$speciesSpatial)) {
    
    if (private$speciesSpatial == 'shared' || private$speciesSpatial == 'replicate') {
      
      self$spatialFields$speciesFields <- list()
      self$spatialFields$speciesFields[['speciesField']] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
      
    }
    else {
      
      for (species in names(self$spatialFields$speciesFields)) {
        
        if (is.null(self$spatialFields$speciesFields[[species]])) self$spatialFields$speciesFields[[species]] <- INLA::inla.spde2.matern(mesh = private$INLAmesh)
        
      }
    }
  }
  
  private$printSummary <- list(Type = unlist(pointData$dataType),
                               numObs = unlist(pointData$numObs),
                               Marks = pointData$Marks,
                               marksType = pointData$marksType)
  
  private$dataSource <- unlist(as.vector(pointData$dataSource))
  
  pointData$makeFormulas(spatcovs = private$spatcovsNames, speciesname = speciesName, temporalname = private$temporalName,
                         paresp = responsePA, countresp = responseCounts, marksspatial = private$marksSpatial, speciesintercept = private$speciesIntercepts, 
                         marks = NULL, spatial = private$Spatial, speciesindependent = private$speciesIndependent, speciesenvironment = private$speciesEnvironment,
                         intercept = private$Intercepts, markintercept = NULL, speciesspatial = private$speciesSpatial, biasformula = private$biasFormula,
                         covariateformula = private$covariateFormula)
  
  if (!is.null(private$temporalName)) {
    
    pointData$makeMultinom(multinomVars = private$temporalName,
                           return = 'time', oldVars = NULL)
    
    private$temporalVars <- pointData$timeIndex
    
  }
  
  ##How does this work?
  ##Need to check again to see how covariates are made.
  if (is.null(private$multinomVars)) {
    
    if (length(pointData$multinomVars) != 0) private$multinomVars <- multinomNames <- pointData$multinomVars
    else multinomNames <- NULL
    
  } else private$multinomVars <- multinomNames <- unique(c(private$multinomVars, pointData$multinomVars))
  
  #CHECK THIS AGAIN
  if (!is.null(multinomNames)) {
    
    if (length(private$multinomIndex) != 0) oldVars <- private$multinomIndex
    else oldVars <- NULL
    
    pointData$makeMultinom(multinomVars = multinomNames, 
                           return = 'marks',
                           oldVars = oldVars)
    
    for (i in names(pointData$multinomIndex)) {
      
      if (i %in% names(private$multinomIndex)) private$multinomIndex[[i]] <- c(private$multinomIndex[[i]],unique(unlist(pointData$multinomIndex[[i]])))[!is.na(c(private$multinomIndex[[i]],unique(unlist(pointData$multinomIndex[[i]]))))]
      else private$multinomIndex[[i]] <- unique(unlist(pointData$multinomIndex[[i]]))[!is.na(unique(unlist(pointData$multinomIndex[[i]])))]
      
    }
    
  }
  
  private$Components <- pointData$makeComponents(spatial = private$Spatial, intercepts = private$Intercepts,
                                                 marks = markNames, datanames = dataNames, speciesname = speciesName,
                                                 multinomnames = multinomNames, pointcovariates = pointCovariates,
                                                 covariatenames = private$spatcovsNames, 
                                                 covariateclass = private$spatcovsClass,
                                                 marksspatial = NULL,
                                                 marksintercept = NULL,
                                                 temporalname = private$temporalName,
                                                 numtime = length(unique(unlist(private$temporalVars))),
                                                 temporalmodel = private$temporalModel,
                                                 speciesspatial = private$speciesSpatial,
                                                 speciesintercept =  private$speciesIntercepts,
                                                 speciesenvironment =private$speciesEnvironment,
                                                 offsetname = private$Offset,
                                                 copymodel = private$copyModel,
                                                 speciesindependent = private$speciesIndependent,
                                                 biasformula = private$biasFormula,
                                                 covariateformula = private$covariateFormula,
                                                 marksCopy = NULL)
  
  
  ##MAKE THIS A FUNCTION TOO
  if (!is.null(private$spatcovsNames)) {
    
    for (data in names(pointData$Data)) {
      
      for (species in 1:length(pointData$Data[[data]])) {
        
        for (cov in private$spatcovsNames) {
          
          if (!is.null(private$biasFormula)) {

            if (cov %in% labels(terms(private$biasFormula))) covIndex <- cov
          else 
            if (!is.null(private$speciesName) && private$speciesEnvironment) covIndex <- paste0(pointData$SpeciesInData[[data]][species],'_', cov)
            else covIndex <- cov
          }
          else
            if (!is.null(private$speciesName) && private$speciesEnvironment) covIndex <- paste0(pointData$SpeciesInData[[data]][species],'_', cov)
            else covIndex <- cov
          
          pointData$Data[[data]][[species]][[covIndex]] <- inlabru::eval_spatial(where = pointData$Data[[data]][[species]], 
                                                                                 data = get('spatialcovariates', 
                                                                                            envir = private$spatcovsEnv)[cov],
                                                                                 layer = cov)
          
          if (any(is.na(pointData$Data[[data]][[species]][[covIndex]]))) {
            
            pointData$Data[[data]][[species]][[covIndex]] <- inlabru::bru_fill_missing(where = pointData$Data[[data]][[species]], 
                                                                                       data = get('spatialcovariates', 
                                                                                                  envir = private$spatcovsEnv)[cov],
                                                                                       layer = cov,
                                                                                       values = pointData$Data[[data]][[species]][[covIndex]])
            
          }
          
          
        }
        
    
    if (!is.null(private$IPS)) {
      
      for (covIPS in private$spatcovsNames) {
     
        if (!is.null(private$biasFormula)) {
          
          if (covIPS %in% labels(terms(private$biasFormula))) covIPSindex <- covIPS
          else 
            if (!is.null(private$speciesName) && private$speciesEnvironment) covIPSindex <- paste0(pointData$SpeciesInData[[data]][species],'_', covIPS)
            else covIPSindex <- covIPS
        }
        else
          if (!is.null(private$speciesName) && private$speciesEnvironment) covIPSindex <- paste0(unique(unlist(private$speciesIn)), '_', covIPS)
          else covIPSindex <- covIPS
        
        for (covADD in covIPSindex) {
          
          private$IPS[[covADD]] <- inlabru::eval_spatial(where =  private$IPS, 
                                                         data = get('spatialcovariates', 
                                                                    envir = private$spatcovsEnv)[covIPS],
                                                         layer = covIPS
          )
          
          if (any(is.na(private$IPS[[covADD]]))) {
            
            private$IPS[[covADD]] <- inlabru::bru_fill_missing(where = private$IPS, 
                                                               data = get('spatialcovariates', 
                                                                          envir = private$spatcovsEnv)[covIPS],
                                                               layer = covIPS,
                                                               values = private$IPS[[covADD]])
            
          }
          
        }
        
        
      }
      
      
    }
        
      }
    }
    
  }
  
  
  if (!is.null(c(private$Offset, private$pointCovariates))) {
    
    datMatrix <- as.data.frame(matrix(NA, nrow = nrow(private$IPS), ncol = length(c(private$Offset, private$pointCovariates))))
    names(datMatrix) <- c(private$pointCovariates, private$Offset)
    private$IPS <- cbind(private$IPS, datMatrix)
    
  }
  
  private$Formulas <- pointData$Formulas
  private$Family <- pointData$Family
  
  private$speciesIndex <- pointData$speciesIndex

  newFamily <- mapply(function(family, number) rep(family , times = number), family = private$Family,
                        number = lapply(private$speciesIn, length))
    
    private$speciesTable <- unique(do.call(rbind, lapply(unlist(pointData$Data, recursive = F), function(x) unique(data.frame(index = data.frame(x)[,private$speciesName], species = data.frame(x)[,paste0(private$speciesName, 'INDEX_VAR')])))))
    private$speciesTable <- private$speciesTable[order(private$speciesTable$index),]
    
  
  ##REDO THIS  
  #newFamily <- pointData$Family
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
  
  private$modelData <- pointData$Data
  private$originalNames <- private$initialnames
  private$initialnames <- NULL
  
  
  
})

specifySpecies$set('private', 'polyfromMesh', function(...) {
  
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

specifySpecies$set('private', 'spatialCovariates', function(spatialCovariates) {
  
  if (missing(spatialCovariates)) stop('Please add spatialCovariates as a Raster* or SpatialPixelsDataFrame object.')
  
  objName <- as.character(match.call())[2]
  
  if (!objName %in% names(parent.frame())) {
    
    if (!objName %in% names(globalenv())) stop('Spatial covariates object not available.')
    else spatcovsEnv <- globalenv()
    
  } 
  else spatcovsEnv <- parent.frame()
  
  if (!class(spatialCovariates) %in% c('SpatRaster',
                                       'SpatialPixelsDataFrame')) stop('The spatial Covariates need to be a spatRaster object or a SpatialPixelsDataFrame.')
  
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

specifySpecies$set('public', 'samplingBias', function(datasetName, Samplers) {
  
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
specifySpecies$set('public', 'spatialFields', list(sharedField = list(),
                                                datasetFields = list(),
                                                biasFields = list(),
                                                speciesFields = list()))
