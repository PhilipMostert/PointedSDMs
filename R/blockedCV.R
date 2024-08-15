#' @title \emph{blockedCV}: run spatial blocked cross-validation on the integrated model.
#' 
#' @description This function is used to perform spatial blocked cross-validation with regards to model selection for the integrated model. It does so by leaving out a block of data in the full model, running a model with the remaining data, and then calculating the deviance information criteria (DIC) as a score of model fit.
#' @param data An object produced by either \code{\link{startISDM}} of \code{\link{startSpecies}}. Requires the slot function, \code{.$spatialBlock} to be run first in order to specify how the data in the model is blocked.
#' @param options A list of \pkg{INLA} or \pkg{inlabru} options to be used in the model. Defaults to \code{list()}.
#' @param method Which cross-validation method to perform. Must be one of \code{'DIC'} or \code{'Predict'}. If \code{'DIC'} then the DIC values for each block are obtained. If \code{'Predict'} then predictions are made on a dataset in the left out block. For this to work, please specify the argument \code{methodOptions}.
#' @param predictName Name of the dataset to predict onto if \code{method = 'Predict'}.
#' @param datasetCombs A list of vectors containing dataset combinations to include in each model run if \code{method = 'Prediction'}. If \code{NULL} then all combinations of the dataset will be estimated.
#' @import inlabru
#' @import stats
#' @importFrom utils combn
#' 
#' @examples 
#' 
#'\dontrun{
#'  if(requireNamespace('INLA')) {
#'    
#'  #Get Data
#'  data("SolitaryTinamou")
#'  proj <- "+proj=longlat +ellps=WGS84"
#'  data <- SolitaryTinamou$datasets
#'  mesh <- SolitaryTinamou$mesh
#'  mesh$crs <- proj
#'  
#'  #Set model up
#'  organizedData <- startISDM(data, Mesh = mesh,
#'                             responsePA = 'Present',
#'                             Projection = proj)
#'  
#'  #Set up spatial block
#'  organizedData$spatialBlock(k = 2, rows = 2, cols = 1)
#'  
#'  #Run spatial block cross-validation
#'  blocked <- blockedCV(organizedData)
#'  
#'  #Print summary
#'  blocked
#'    
#'  }
#'}
#' 
#' @return An object of class \code{blockedCV}, which is essentially a list of DIC values obtained from each iteration of the model.
#' 
#' @export
#' 
blockedCV <- function(data, options = list(),
                      method = 'DIC', predictName = NULL,
                      datasetCombs = NULL) {
  
  #How do we do this?
   #Should we make data a list of data files;
   # or should we make another argument for thinned formulas to test based on the full model?
  
  if (!inherits(data, 'dataSDM') && !inherits(data, 'specifySpecies') && !inherits(data, 'specifyISDM')) stop('data needs to be a dataSDM object.')
  
  if (is.null(data$.__enclos_env__$private$INLAmesh)) stop('An inla.mesh object is required before any model is run.')
  
  if (!data$.__enclos_env__$private$blockedCV) stop('Please use ".$spatialBlock" before using this function.')
  
  data2ENV(data = data, env = environment())
  
  results <- list()
  
  if (!method %in% c('DIC', 'Predict')) stop('method needs to be one of "DIC" or "Predict".')
  
  if (method == 'Predict') {
    
    if (is.null(predictName)) stop('predictName cannot be NULL if method = "Predict".')
    
    if (!predictName %in% data$.__enclos_env__$private$dataSource) stop('predictData needs to be the name of a dataset added to the model.')
    ##CHANGE
     #Remove predict dataset here, and then add it as its own
    dataIn <- unique(data$.__enclos_env__$private$dataSource)
    
    if (!is.null(datasetCombs)) {
      
      if (!all(unique(unlist(datasetCombs)) %in% dataIn)) stop('Dataset added in datasetCombs not available.')
      
    } else datasetCombs <- do.call(c, lapply(seq_along(dataIn), combn, x = dataIn, simplify = FALSE))
    
  } else datasetCombs = list(data$.__enclos_env__$private$dataSource)
  
  block_index <- lapply(unlist(data$.__enclos_env__$private$modelData, recursive = FALSE), function(x) data.frame(x)[, '.__block_index__'])
  
  if (!is.null(data$.__enclos_env__$private$temporalName)) {
    
    numTime <- length(unique(unlist(data$.__enclos_env__$private$temporalVars)))
    
    newIPS <- rep(list(data$.__enclos_env__$private$IPS), numTime)
    
    newIPS <- do.call(rbind, newIPS)
    
    newIPS[, data$.__enclos_env__$private$temporalName] <- rep(1:numTime, each = nrow(data$.__enclos_env__$private$IPS))
    
    newIPS <- st_transform(newIPS, data$.__enclos_env__$private$Projection)
    
    data$.__enclos_env__$private$IPS <- newIPS
    
  }
  
  for (dataSub in 1:length(datasetCombs)) {
  
    dataToUse <- unique(datasetCombs[[dataSub]])
    
  for (fold in unique(unlist(block_index))) {
    
    ##Maybe make block_index in dataSDM as a list such that we can see which datasets are not in block i to easily remove them.
     #Get formula terms only after likelihood construction, and then thin components from there.
     #And also for the whole, control.family thing
    
    #Subset here based on dataSub
    trainData <- lapply(data$.__enclos_env__$private$modelData[dataToUse], function(data) {
      
      lapply(data, function(x) {
        
        data <- x[x$.__block_index__ != fold,]
        if (nrow(data) > 0) data
        else NULL
        
      })
      
      
    })
    sourceIN <- which(unique(data$.__enclos_env__$private$dataSource) %in% names(trainData))
    sourcePred <- which(data$.__enclos_env__$private$dataSource %in% predictName)
    comp_terms <- gsub('\\(.*$', '', data$.__enclos_env__$private$Components)
    
    ##Check if all copy + bias Main in here
    
    #modelData[predictName]
    
    if (method == 'Predict') {
    
      testData <- lapply(data$.__enclos_env__$private$modelData[predictName], function(data) {
      
      lapply(data, function(x) {
        
        data <- x[x$.__block_index__ == fold,]
        if (nrow(data) > 0) data
        else NULL
        
      })
      
    })
      
      #if testData dosen't cover all blocks stop
      
      testLike <- do.call(inlabru::like_list,
                          makeLhoods(data = testData, #Is this an sf dataset?
                                     formula = data$.__enclos_env__$private$Formulas[predictName],
                                     family = data$.__enclos_env__$private$Family[predictName],
                                     mesh = data$.__enclos_env__$private$INLAmesh,
                                     ips = data$.__enclos_env__$private$IPS,
                                     samplers = data$.__enclos_env__$private$Samplers,
                                     paresp = data$.__enclos_env__$private$responsePA,
                                     ntrialsvar = data$.__enclos_env__$private$trialsPA,
                                     markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                                     speciesname = data$.__enclos_env__$private$speciesName,
                                     speciesindex = data$.__enclos_env__$private$speciesIndex))
      
      #Collapse into 1 sf
       #how? do.call(rbind, testData) or do.call(c, testData)
      
      #For comps: take environmental covariates: paste
       #if shared then keep
      if (!is.null(data$.__enclos_env__$private$Spatial)) {
       if (data$.__enclos_env__$private$Spatial == 'copy') spatRM <- TRUE
       else spatRM <- FALSE #What about if species copy?
      } else spatRM <- FALSE
      
      oldIn <- lapply(testLike, function(x) x$used$effect)
      
      #predComp <- comp_terms %in% formIn
      
      #if (spatRM) {
        
      #  compsChange <- data$.__enclos_env__$private$Components
      #  whichRM <- grepl(paste0(predictName, '_spatial'), compsChange)
      #  compsChange[whichRM] <- paste0(predictName, '_spatial(main = geometry, model = predict_field)')
      #  predict_field <- data$spatialFields$datasetFields[[1]]
        
      #} else compsChange <- data$.__enclos_env__$private$Components
      
      #compPreds <- formula(paste('~ - 1 +', paste(c(compsChange[predComp], 'olikhoodvar(main = olikhoodvar, model = "offset")'), collapse = ' + ')))
      
      for (lik in 1:length(testLike)) {
        
      testLike[[lik]]$used$effect <- c(paste0('testIntercept', lik), 'olikhoodvar')#c(testLike[[1]]$used$effect, 'olikhoodvar')
      
      }

    }

    trainLiks <- do.call(inlabru::like_list,
                 makeLhoods(data = trainData,
                 formula = data$.__enclos_env__$private$Formulas[sourceIN],
                 family = data$.__enclos_env__$private$Family[sourceIN],
                 mesh = data$.__enclos_env__$private$INLAmesh,
                 ips = data$.__enclos_env__$private$IPS,
                 samplers = data$.__enclos_env__$private$Samplers,
                 paresp = data$.__enclos_env__$private$responsePA,
                 ntrialsvar = data$.__enclos_env__$private$trialsPA,
                 markstrialsvar = data$.__enclos_env__$private$trialsMarks,
                 speciesname = data$.__enclos_env__$private$speciesName,
                 speciesindex = data$.__enclos_env__$private$speciesIndex))
      
    formula_terms <- unique(unlist(lapply(trainLiks, function(x) {
      
      if (!identical(unlist(x$used), character(0))) unlist(x$used)
      else labels(terms(x$formula))
      
    })))
    
    if (method == 'Predict') {
      
      if (any(grepl('_biasField', formula_terms)) & length(dataToUse) == 1) formula_terms <- formula_terms[!grepl('_biasField', formula_terms)]
    
      fams <- unlist(lapply(trainLiks, function(x) x$family)) == 'cp'

    if (!any(fams)) {
      
      ips <- data$.__enclos_env__$private$IPS
      ips$respIPS <- 0
      
      if (inherits(data, 'specifySpecies')) {
        
        if (!is.null(data$.__enclos_env__$private$speciesSpatial)) {
          
        if (data$.__enclos_env__$private$speciesSpatial == 'replicate') ips <- fm_cprod(ips, data.frame(speciesSpatialGroup = 1:max(data$.__enclos_env__$private$speciesTable$index)))
        
        }
        if (!is.null(data$.__enclos_env__$private$Intercepts)) {
          
          if (data$.__enclos_env__$private$speciesIntercepts) ips <- fm_cprod(ips, data.frame(specIntTermRem = 1:max(data$.__enclos_env__$private$speciesTable$index)))
          names(ips)[names(ips) == 'specIntTermRem'] <- data$.__enclos_env__$private$speciesName
        }
      } 
      
      ipsLike <- inlabru::like(formula = respIPS ~ .,
                                include = formula_terms, E = ips$weight,
                    family = 'poisson', data = ips)
      
      trainLiks[['ips']] <- ipsLike
      uFam <- TRUE
    } else uFam <- FALSE
      
    } else uFam <- FALSE
    
    
    comp_keep <- comp_terms %in% formula_terms
    
    if (!all(comp_keep)) {
      
      if (!is.null(data$.__enclos_env__$private$Spatial)) {
      
      if (data$.__enclos_env__$private$Spatial != 'shared' | data$.__enclos_env__$private$biasCopy) {
        
        Main <- grepl('_spatial', data$.__enclos_env__$private$Components) & grepl('_field', data$.__enclos_env__$private$Components)
        if (!any(Main)) Main <- 'NOTMAIN'
        else Main <- sub("\\(.*", "", comp_terms[Main])
        
        MainBias <- grepl('_biasField', data$.__enclos_env__$private$Components) & grepl('_bias_field', data$.__enclos_env__$private$Components)
        if (!any(MainBias)) MainBias <- 'NOTALLMAINBIAS'
        else MainBias <- sub("\\(.*", "", comp_terms[MainBias])
        
        whichMissing <- union(data$.__enclos_env__$private$dataSource[!data$.__enclos_env__$private$dataSource %in% dataToUse], 
                              names(trainData)[sapply(unlist(trainData, recursive = F), is.null)])
        
        if(paste0(whichMissing[1],'_biasField') == MainBias | paste0(whichMissing[1],'_spatial') == Main) {
        
          warning('Main dataset or more than 2 datasets missing from the block with either pointsSpatial = "copy" or copyModel = TRUE for the bias field.\n Will choose the first available dataset to copy on.')
          
          thinnedComponents <- reduceComps(componentsOld =  formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components, collapse = ' + '))),
                                            pointsCopy = ifelse(data$.__enclos_env__$private$Spatial == 'copy', 
                                                                TRUE, FALSE),
                                            biasCopy = data$.__enclos_env__$private$biasCopy,
                                            datasetName = whichMissing[1],
                                            reducedTerms = comp_terms[comp_keep])
          
        } else  thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
        
        
      } else thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
      
      } else thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))
      
    }
    else thinnedComponents <- formula(paste('~ - 1 +', paste(data$.__enclos_env__$private$Components[comp_keep], collapse = ' + ')))

    foldOptions <- data$.__enclos_env__$private$optionsINLA
    
    fold_ind <- unique(unlist(block_index))[unique(unlist(block_index)) != fold]
    
    foldOptions$control.family <- foldOptions$control.family[sapply(unlist(data$.__enclos_env__$private$modelData, recursive = FALSE), 
                                                                    function(x) any(fold_ind %in% data.frame(x)[, '.__block_index__']))]

    #Subset fold options for each family
    sourceIN <- which(data$.__enclos_env__$private$dataSource %in% names(trainData))
    foldOptions$control.family <- foldOptions$control.family[sourceIN]
    
    optionsTrain <- append(options, foldOptions)
    
    if (uFam) optionsTrain$control.family <- append(optionsTrain$control.family, list(list(link = 'log')))
    
    ##Calculate DIC for just this model?
    trainedModel <- try(inlabru::bru(components = thinnedComponents,
                                 trainLiks,
                                 options = optionsTrain))
    
    ## -log(intensity)
    ## add an offset argument...
    
    if (method == 'DIC') {
    
    if (inherits(trainedModel, 'try-error')) {
      
      warning('Model failed for a block. Please change your block layout to ensure all datasets are included in all blocks')
      results[[paste0('DIC_fold_', fold)]] <- NA
    }
    else results[[paste0('DIC_fold_', fold)]] <- trainedModel$dic
    
    }
    else {
      
      if (inherits(trainedModel, 'try-error')) {
        
        warning('Model failed for a block. Please change your block layout to ensure all datasets are included in all blocks')
        #paste_dataset here
        results[[paste(dataToUse, collapse = ' and ')]][[paste0('fold',fold)]] <- NA
      
      }
      else {
        
        if (!is.null(data$.__enclos_env__$private$biasFormula)) biasFormlabels <- c(labels(terms(data$.__enclos_env__$private$biasFormula)), 'Bias__Effects__Comps')
        else biasFormlabels <- NULL

        for(pd in 1:length(testData[[1]])) {
          
        #covInPres <- testLike[[pd]]$used$effect
        #covInPres <- covInPres[!covInPres %in% c('olikhoodvar', biasFormlabels, grep('_biasField', covInPres))] #bias
          
        #IF species only include the species in the likelihood
          #Remove bias
          covInPres <- formula_terms
          covInPres <- covInPres[!covInPres %in% biasFormlabels]
          covInPres <- covInPres[!covInPres %in% grepl('_biasField', covInPres)]
          
        if (!is.null(data$.__enclos_env__$private$speciesName)) {
          
          likeSpec <- unique(testLike[[pd]]$data[[paste0(data$.__enclos_env__$private$speciesName, 'INDEX_VAR')]])
          notSpec <- unique(unlist(data$.__enclos_env__$private$speciesIn))
          notSpec <- notSpec[notSpec != likeSpec]
          
          if (data$.__enclos_env__$private$speciesEnvironment) {
            
            specCov <- apply(expand.grid(paste0(notSpec,'_'), data$.__enclos_env__$private$spatcovsNames), MARGIN = 1, FUN = paste0,collapse = '')
            specCov <- c(specCov, paste0(notSpec, '_Fixed__Effects__Comps'))
            covInPres <- covInPres[!covInPres %in% specCov]
            
          }
          
          if (!is.null(data$.__enclos_env__$private$speciesIntercepts)) {
          
            if (!data$.__enclos_env__$private$speciesIntercepts)  covInPres <- covInPres[!covInPres %in% paste0(notSpec,'_intercept')]
          
          }
          
          if (!is.null(data$.__enclos_env__$private$speciesSpatial)) {
            
            if (data$.__enclos_env__$private$speciesSpatial == 'copy') {
              
              rmSpat <-  grepl(paste0(notSpec, collapse = '|'), covInPres) & grepl('_spatial', covInPres) 
              covInPres <- covInPres[!rmSpat]
              
            }
            
          }
          
        }
          #if bias
          
        #covInPres <- intersect(oldIn[[pd]], formula_terms)
        ##Get all old vars in test like and thin with formula terms
        predForm <- formula(paste0('~(', paste(covInPres, collapse = ' + '), ')'))
          
        #Change this to testLike[[pd]]$data
        testPredicts <- suppressWarnings(predict(trainedModel,testLike[[pd]]$data, formula = predForm)) # testData[[1]][[pd]]
        #IF cp add log(1) or NA to ipoints?
        #if (testLike[[pd]]$family == 'cp') {
          
        #  nIPS <- nrow(data$.__enclos_env__$private$IPS)
        #  testLike[[pd]]$data$olikhoodvar <- c(testPredicts$mean, rep(0, nIPS))
          
        #} else 
        testLike[[pd]]$data$olikhoodvar <- testPredicts$mean
        
        }
        
        #Need to look at formula, and get the correct comps here
        foldOptions <- data$.__enclos_env__$private$optionsINLA
        foldOptions$control.family <- foldOptions$control.family[sourcePred]
        
        optionsTest <- append(options, foldOptions)
        compsIntercepts <- paste0('testIntercept', 1:length(testData[[1]]),'(1)')
        compPreds <- formula(paste0('~ - 1 + olikhoodvar(main = olikhoodvar, model = "linear") + ', paste0(compsIntercepts, collapse = ' + ')))

        testModel <- try(inlabru::bru(components = compPreds,
                                   testLike, options = optionsTest))
        
        if (inherits(testModel, 'try-error')) results[[paste(dataToUse, collapse = ' and ')]][[paste0('fold',fold)]] <- NA
        else results[[paste(dataToUse, collapse = ' and ')]][[paste0('fold',fold)]] <- testModel$mlik[[1]] #mean(testModel$residuals$deviance.residuals)##trainedModel$mlik[[1]] - 
        
        
      }
    
    }
  }
    
  }
  
  comps <- formula(paste0(' ~ ', paste0(gsub('\\(.*$', '', data$.__enclos_env__$private$Components), collapse = ' + ')))
  
  if (method == 'DIC') {
  
  results <- append(results, list(Formula = comps))
  class(results) <- c('blockedCV', 'list')
  
  }
  else {
    
    class(results) <- c('blockedCVpred', 'list')
    
  }
  results

}


#' Export class blockedCV
#' 
#' @export

setClass('blockedCV')

#' Print for blockedCV
#' 
#' @export print.blockedCV

#' Export print.blockedCV
#' @title Print function for \code{blockedCV}.
#' @param x A blockedCV object.
#' @param ... Unused argument.
#' 
#' @exportS3Method 

print.blockedCV <- function(x, ...) {
  
  cat('Spatial block cross-validation score:')
  cat('\n\n')
  cat('Formula: ')

  cat(deparse1(x$Formula))
  cat('\n\n')
  x$Formula <- NULL
  mean.deviance <- sapply(x, function(y) y$mean.deviance)
  p.eff <- sapply(x, function(y) y$p.eff)
  dic <- sapply(x, function(y) y$dic)
  
  dataobj <- data.frame(mean.deviance = mean.deviance,
                        p.eff = p.eff,
                        dic = dic)
  
  row.names(dataobj) <- paste0('fold ', 1:nrow(dataobj))
  
  print.data.frame(dataobj)
  
  cat('\nmean DIC score: ')
  cat(mean(dataobj$dic, na.rm = TRUE))


}

#' Export class blockedCVpred
#' 
#' @export

setClass('blockedCVpred')

#' Print for blockedCVpred
#' 
#' @export print.blockedCVpred

#' Export print.blockedCVpred
#' @title Print function for \code{blockedCV}.
#' @param x A blockedCV object.
#' @param ... Unused argument.
#' 
#' @exportS3Method 

print.blockedCVpred <- function(x, ...) {
  
  cat('Spatial block cross-validation score:')
  cat('\n\n')
  cat('Log marginal likelihoods obtained by using the predictions from the model as an offset:\n\n')

  foldFrame <- do.call(rbind.data.frame, x)
  foldFrame <- data.frame(foldFrame, mean = rowMeans(foldFrame))
  colnames(foldFrame) <- c(paste0('fold ', 1:(ncol(foldFrame) - 1)), 'mean')
  
  #foldFrame <- data.frame(Data.included = row.names(foldFrame),
  #                        foldFrame)
  #names(foldFrame)[1] <- 'Data included'
  #row.names(foldFrame) <- NULL
  print.data.frame(foldFrame)
  
}
