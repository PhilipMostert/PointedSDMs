#' Function to run marked species distribution models with data from presence only and presence absence data. (and then later account for over dispersion).
#'
#' @param ... Point process datasets with coordinates of species, and optionally marks and covariates explaining the coordinates.
#' @param spatialcovariates Data frame of the spatial covariates accompanied by their associated coordinates.
#' @param marks Should the model be a marked point process. Defaults to \code{FALSE}.
#' @param markfamily Assumed distribution of the marks. Defaults to \code{"gaussian"}.
#' @param inclmarks. A vector of which marks should be included in the model. Defaults to \code{NULL}.
#' @param coords Vector of the names of the coordinates used in datasets. Defaults to \code{c('X','Y')} (For now should be standardized)
#' @param presname Names of presences column in observs. Defaults to \code{"NPres"}. Note that this column can also be logical.
#' @param trialname Names of column of number of columns in observs. Defaults to \code{NULL}.
#' @param inclcoords Should coordinates be used in data. Defaults to \code{FALSE}.
#' @param mesh An inla.mesh object. Defaults to \code{NULL}.
#' @param meshpars List of mesh parameters. Requires the following items: "cut.off", "max.edge" and "offset". Defaults to \code{NULL}.
#' @param ips Integration points. Defaults to \code{NULL}.
#' @param bdry Polygon of boundary for region, of class Polygon. If \code{NULL}, draws a boundary around the points.
#' @param proj Projection to use if data is not a projection. Defaults to utm (hopefully).
#' @param residuals Should residuals for each dataset be calculated. Defaults to \code{TRUE}.
#' @param predictions Boolean: should predictions (on the linear scale) be made? Defaults to \code{FALSE}.
#' @param intercept Include joint intercept in the model. Defaults to \code{NULL}.
#' @param options A bru_options options object or a list of options passed on to bru_options()
#' @param namespde Name for spatial term used for the points in the model. Defaults to \code{"MySPDE"}. Set to \code{NULL} if no spatial term wanted.
#' @param markspde Name for spatial term used for the marks in the model. Defaults to \code{"MarkSPDE"}. Set to \code{NULL} if no spatial term wanted.
#' @param poformula Formula given to the presence only datasets. Defaults to \code{NULL}.
#' @param paformula Formula given to the presence/absence datasets. Defaults to \code{NULL}.
#' @param tol Tolerance parameter for SpatialPixelsDataFrame. Defaults to \code{NULL}.
#' 
#' @import sp
#' @import INLA
#' @import inlabru
#' @import rgeos ## Need for mesh creation??

bru_sdm = function(..., spatialcovariates, marks = FALSE, markfamily = 'gaussian',
                   inclmarks = NULL, coords = c('X','Y'), presname = 'NPres', trialname= NULL,
                   inclcoords = FALSE, mesh = NULL, meshpars = NULL ,ips = NULL, 
                   bdry = NULL, proj = CRS("+proj=longlat +ellps=WGS84"), predictions = FALSE,
                   residuals = TRUE, intercept = TRUE, namespde = 'MySPDE', markspde = 'MarkSPDE',
                   options = list(), poformula = NULL, paformula = NULL, tol = NULL) {
  
  
  if (is.null(spatialcovariates)) stop("Spatial covariates not provided.")
  
  if (presname == 'NPres') warning('Default presence names is being used for presence/absence data.\nPlease change if this is incorrect.')
  
  if (length(presname) > 1) stop("More than one name given for presences column.")
  
  if (length(trialname) > 1) stop("More than one name given for number of trials column.")
  
  if (length(coords) != 2) stop("Coordinates must have two components.")
  
  if (as.character(proj@projargs) == "+proj=longlat +ellps=WGS84 +no_defs") warning("Default CRS is being used. Please change if incorrect.")
  
  if (!is.null(poformula)) {
    
    if (as.character(poformula[2]) != 'resp') stop("Response variable for presence only datasets should be 'coordinates'.")
    
  }
  
  if(!is.null(paformula)) {
    
    if (as.character(paformula[2]) != 'resp') stop("Response variable for presence absence datasets should be 'resp'.")
    
  }
  
  if (class(spatialcovariates) == 'data.frame' & is.null(tol)) stop('Please provide a tolerance parameter to convert the spatial covariates to a SpatialPixelsDataFrame.')
  
  if (class(proj) != 'CRS') stop("Proj needs to be a CRS object.")
  
  if (is.null(bdry) & is.null(ips)) stop("Either boundary or intepration points need to be provided.")
  
  if (is.null(mesh) & is.null(meshpars)) stop("Either a mesh or mesh parameters need to be provided.")
  
  if (is.null(mesh)) {
    
    if (sum(names(meshpars)%in%c( "cut.off", "max.edge", "offset")) < 3)
      
      stop("Meshpars requires three items in the list: cut.off, max.edge and offset.")
    
  }
  
  if (!marks) names_marks <- NULL; data_marks <- NULL
  
  if(!marks & !is.null(inclmarks)) {
    
    warning('Marks to include is non null but include marks is set to FALSE.\nMarks are now being set to TRUE.')
    marks <- TRUE
    
  } 
  
  
  datasets = list(...)
  
  coords_in = unlist(lapply(datasets, function(dat) {
    
    if (class(dat) == 'data.frame') coords%in%names(dat)
    
  }))
  
  if (!all(coords_in)) stop("At least one dataset does not have coordinates in it.\nEither check your datasets or change your coordinates argument.")
  
  
  data_names <- setdiff(as.character(match.call(expand.dots=TRUE)), 
                        as.character(match.call(expand.dots=FALSE)))
  
  #Separate PO and PA data by inclusion/exclusion of 'trialname'.
  
  data_attributes <- lapply(datasets,function(dat) {
    if (inherits(dat,"Spatial")) {
      if (class(dat) == "SpatialPoints") {
        
        dat = sp::SpatialPointsDataFrame(coords = sp::coordinates(dat),
                                         data = data.frame(resp = rep(1,nrow(coordinates(dat)))),
                                         proj4string = proj)
        attr(dat,'family') <- 'poisson'
        attr(dat,'data_type') <- 'Present only'
        dat
        
      } 
      else
        if (presname%in%colnames(dat@data)){
          
          dat <- sp::SpatialPointsDataFrame(coords = dat@coords,
                                            data = as.data.frame(dat@data),
                                            proj4string = proj)
          names(dat)[names(dat) == presname] <- 'resp' 
          dat@data[,'resp'] <- as.numeric(dat@data[,'resp'])
          if (!is.null(trialname)) {
            if (trialname%in%colnames(dat@data)) attr(dat,'Ntrials') <- dat@data[,trialname]
            else attr(dat,'Ntrials') <- 1
          }
          attr(dat,'family') <- 'binomial'
          attr(dat,'data_type') <- 'Present absent'
          dat
          
        }
      else { 
        
        dat <- sp::SpatialPointsDataFrame(coords = dat@coords,
                                          data = as.data.frame(dat@data),
                                          proj4string = proj)
        dat$resp = 1
        attr(dat,'family') <- 'poisson'
        attr(dat,'data_type') <- 'Present only'
        dat
        
      }
      
    }
    else
      if (class(dat) == 'data.frame') {
        if (presname%in%colnames(dat)) {
          if (ncol(dat) == 1) stop("Only one column provided in data frame.\nNeed two columns for coordinates and one for presence name\nfor presence/absence data.")
          else
            if (ncol(dat) == 2) stop("Only two columns provided for presence/absence data. Either coordinates is of length one or presence name not given.")
          else {
            
            names <- names(dat)[!names(dat)%in%c(coords)]
            dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                              data = as.data.frame(dat[,!names(dat)%in%coords]),
                                              proj4string = proj) #dat[,!names(x)%in%coords]
            colnames(dat@data) <- names
            names(dat)[names(dat) == presname] <- 'resp' 
            dat@data[,'resp'] <- as.numeric(dat@data[,'resp'])
            if (!is.null(trialname)) {
              if (trialname%in%colnames(dat@data)) attr(dat,'Ntrials') <- dat@data[,trialname]
              else attr(dat,'Ntrials') <- 1
            }
            attr(dat,'family') <- 'binomial'
            attr(dat,'data_type') <- 'Present absent'
            dat
            
          }
        }
        else {
          if (ncol(dat) == 1) stop('Data Frame provided only has one column. Coordinates require two columns.')
          else
            if (ncol(dat) == 2){
              
              dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                                data = data.frame(resp = rep(1,nrow(dat))),
                                                proj4string = proj)
              attr(dat,'family') <- 'poisson'
              attr(dat,'data_type') <- 'Present only'
              dat
              
            }
          else {
            
            names <- names(dat)[!names(dat)%in%c(coords)]
            dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                              data = as.data.frame(dat[,!names(dat)%in%coords]),
                                              proj4string = proj)
            colnames(dat@data) <- names
            dat$resp = 1
            attr(dat,'family') <- 'poisson'
            attr(dat,'data_type') <- 'Present only'
            dat
            
          }
          
        }
      }
  })
  
  names(data_attributes) <- data_names
  
  if (inclcoords) {
    ##Should I include?
    for (i in 1:length(data_attributes)) {
      if (class(data_attributes[[i]]) == 'SpatialPoints') {
        
        data_attributes[[i]] <- sp::SpatialPointsDataFrame(coords = data_attributes[[i]]@coords,
                                                           data = as.data.frame(data_attributes[[i]]@coords),
                                                           proj4string = proj)
        
      }
      else {
        
        data_attributes[[i]]@data[,coords] <- data_attributes[[i]]@coords
        
      }
    }
    
  }
  
  if (marks) {
    
    data_marks = list()
    #Make a unique SpatialPointsDataframe for each mark (to be run on spatial covariates).
    #Incorporate standardized list of marks/covariates for data.
    for (i in 1:length(data_attributes)) {
      
      ind <- 1 + length(data_marks)
      
      if (class(data_attributes[[i]]) == 'SpatialPoints') data_marks[[ind]] <- FALSE
      else {
        
        names = names(data_attributes[[i]])[!names(data_attributes[[i]])%in%c('resp',coords,trialname)]
        
        if (!is.null(inclmarks)) names <- names[names%in%inclmarks]
        
        if (is.null(names) | identical(names,character(0))) data_marks[[ind]] <- FALSE
        
        else
          if(length(names) == 1) {
            
            mark <- sp::SpatialPointsDataFrame(coords = coordinates(data_attributes[[i]]),
                                               data = as.data.frame(data_attributes[[i]]@data[,names]),
                                               proj4string = proj)
            colnames(mark@data) <- paste0(names(data_attributes)[i],'_',names)
            attr(mark,'family') <- markfamily
            capital_markfamily <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2", 
                                       markfamily, perl = TRUE)
            attr(mark,'data_type') <- paste0(capital_markfamily,' mark')
            data_marks[[ind]] <- mark
            names(data_marks)[[ind]] <- paste0(names(data_attributes)[i],'_',names)
            
          }
        
        else
          if(length(names) > 1) {
            for(j in 1:length(names)) {
              
              ind <- ind + j - 1
              mark <- sp::SpatialPointsDataFrame(coords = coordinates(data_attributes[[i]]),
                                                 data = as.data.frame(data_attributes[[i]]@data[,names[j]]),
                                                 proj4string = proj)
              colnames(mark@data) <- paste0(names(data_attributes)[i],'_',names[j])
              attr(mark,'family') <- markfamily
              capital_markfamily <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2", 
                                         markfamily, perl = TRUE)
              attr(mark,'data_type') <- paste0(capital_markfamily,' mark')
              data_marks[[ind]] <- mark
              names(data_marks)[[ind]] <- paste0(names(data_attributes)[i],'_',names[j])
              
            }
            
          }
        
      }
    }
    
    data_marks[sapply(data_marks,is.logical)] <- NULL
    
    if (length(data_marks) == 0) stop("Either marks have been set to TRUE and no datasets contain marks, or marks to include only contains marks not present in any dataset.")
    
    names_marks <- sapply(data_marks, function(mark) colnames(mark@data))
    
  }
  
  if (is.null(mesh)) {
    
    warning("Mesh not provided. Will try to create own mesh.")
    
    #Make mesh same way as PointedSDMs
    if (is.null(bdry)) {
      
      if (inherits(spatialcovariates, "Spatial")) spatcoords <- sp::SpatialPoints(coords = spatialcovariates@coords,
                                                                                  proj4string = proj)
      else spatcoords <- sp::SpatialPoints(coords = spatialcovariates[,coords],
                                           proj4string  = proj)
      
      bstart <- min(c(diff(sort(unique(spatcoords@coords[,1]))), diff(sort(unique(spatcoords@coords[,2])))))
      
      poly.tmp <- rgeos::gBuffer(spatcoords, width=bstart, byid=TRUE)
      
      bdry <- rgeos::gBuffer(rgeos::gUnaryUnion(poly.tmp), width=bstart)
      
    }
    
    else {
      if (class(bdry)!="SpatialPolygons") {
        
        bdry <- sp::SpatialPolygons(Srl=list(Polygons(srl=list(bdry), ID="eek")))
        
      } 
      else {
        if (!is.projected(bdry)) bdry <- spTransform(bdry, CRSobj = proj)
      }
      
    }
    
    region.bdry <- inla.sp2segment(bdry)
    
    #Is there a nice way to add other parameters to inla.mesh.2d?
    mesh <- inla.mesh.2d(boundary=region.bdry, 
                         cutoff=meshpars$cutoff,
                         max.edge=meshpars$max.edge, 
                         offset=meshpars$offset)
    
  }
  
  if (is.null(ips)) {
    
    warning('Integration points not provided. Will try to create own points')
    
    ips <- ipoints(samplers = bdry,
                   domain = mesh)
    
  }
  
  #Should I do the same for raster data??
  #Does this do the same thing as 'GetNearestCovariate'?
  #When inlabru update comes, change SpatialPointsDataFrame part
  #SpatialGridDataFrame?
  
  if (inherits(spatialcovariates,'Spatial')) {
    if (ncol(spatialcovariates) == 1) {
      if(class(spatialcovariates) == 'SpatialPixelsDataFrame') {
        
        proj4string(spatialcovariates) <- proj
        spatnames <- names(spatialcovariates@data)
        assign(spatnames, spatialcovariates)}
      
      else {
        
        spatnames <- names(spatialcovariates)
        spatcoords <- spatialcovariates@coords
        spatdata <- spatialcovariates@data[,!colnames(spatialcovariates@data)%in%coords]
        spatpix <- sp::SpatialPixelsDataFrame(points = spatcoords,
                                              data = data.frame(spatdata), 
                                              tolerance = tol,
                                              proj4string = proj)
        names(spatpix@data) <- spatnames
        assign(names(spatpix@data),spatpix)
      }
    }
    else {
      
      spatcoords <- spatialcovariates@coords
      spatdata <-  spatialcovariates@data[,!colnames(spatialcovariates@data)%in%coords]
      spatnames <- names(spatdata)
      
      for (i in 1:ncol(spatdata)) {
        
        spatpix <- sp::SpatialPixelsDataFrame(points = spatcoords,
                                              data = data.frame(spatdata[,i]), 
                                              tolerance = tol,
                                              proj4string = proj)
        colnames(spatpix@data) = colnames(spatdata)[i]
        assign(names(spatpix@data),spatpix)
        
      }
    }
  }
  
  else
    if (class(spatialcovariates) == 'data.frame') {
      
      warning("Spatialcovariates is of class 'data.frame'.\nWill convert it to a SpatialPixelsDataFrame.")
      spatcoords <- spatialcovariates[,coords]
      spatdata <- as.data.frame(spatialcovariates[,!colnames(spatialcovariates)%in%coords])
      spatnames <- names(spatialcovariates)[!colnames(spatialcovariates)%in%coords]
      
      for (i in 1:ncol(spatdata)) {
        spatpix <- sp::SpatialPixelsDataFrame(points = spatcoords,
                                              data = data.frame(spatdata[,i]), 
                                              tolerance = tol,
                                              proj4string = proj)
        colnames(spatpix@data) = spatnames[i]
        
        assign(names(spatpix@data),spatpix)
        
      }
      
    }
  
  spde2 = inla.spde2.matern(mesh)
  
  ##Construct joint components for the likelihoods.
  ##Will need to change with inclusion of separate covariates.
  
  if (is.null(poformula) | is.null(paformula)){
    
    components_joint <- formula(paste(c('~ 0',paste0(spatnames,'(main = ',spatnames,', model = "linear")')), collapse = '+'))
    
    if (inclcoords) {
      
      components_joint <- update(components_joint, paste0('~ . +',coords, collapse = '+'))
    }
    
    if (!is.null(namespde)) {
      
      components_joint <- update(components_joint, paste0(' ~ . +',namespde,"(main = coordinates, model = spde2)"))
      
    }
    
    if (intercept) {
      
      #Find way to add individual intercepts nicely
      components_joint <- update(components_joint, ~ . + Intercept(1))
    }
    
  }
  
  likelihoods = list()
  
  family <- sapply(data_attributes, function(x) attributes(x)$family)
  
  trials <- sapply(data_attributes, function(x){
    
    if (!is.null(attributes(x)$Ntrials)) attributes(x)$Ntrials
    else 1
    
  })
  
  ##Take out any brackets from 'components_joint'.
  ##I.e (for now) run coordinates only on spatial covariates (and optional others).
  form_elements <- gsub(" *\\(.*?\\) *", "",components_joint)
  
  formula <- sapply(family, function(fam) {
    if (!is.null(poformula) & fam == 'poisson') {
      
      formula <- poformula
    }
    else
      if (is.null(poformula) & fam == 'poisson') {
        
        formula <- formula(paste0(c("resp ~", form_elements[2]),collapse = " ")) 
        
      }
    
    else 
      if(!is.null(paformula) & fam == 'binomial') {
        
        formula <- paformula
        
      }
    else{
      
      formula <- formula(paste0(c("resp ~", form_elements[2]),collapse = " "))
      
    }
    
  })
  
  for (i in 1:1) {
    
    lhoods <- like(formula = formula[[i]], ##Add tag to this likelihood somehow?
                   family = family[i],
                   data = data_attributes[[i]],
                   mesh = mesh,
                   ips = ips,
                   Ntrials = trials[i])
    likelihoods <- like_list(lhoods)
    
    if (length(family) > 1) {
      for (j in 2:length(family)) {
        
        lhoods <- like(formula = formula[[j]],
                       family = family[j],
                       data = data_attributes[[j]],
                       mesh = mesh,
                       ips = ips,
                       Ntrials = trials[j])
        likelihoods[[j]] <- lhoods
        
      }
    }
    
    likelihoods
  }
  
  if (marks) {
    
    family_marks <- sapply(data_marks, function(x) attributes(x)$family)
    formula_marks <- list()
    likelihoods_marks <- list()
    
    if (!is.null(markspde)) {
      
      #Somehow add copy feauture to marks SPDE? Not working
      #Add grouping? Not working
      components_joint <- update(components_joint, paste0(' ~ . +',markspde,"(main = coordinates, model = spde2)"))
      
    }    
    
    for (i in 1:length(family_marks)) {
      
      formula_marks[[i]] <- formula(paste0(c(names_marks[i],'~',form_elements[2]),collapse = " "))
      #Can I just remove the coordinates spde from the marks formula?
      formula_marks[[i]] <- update(formula_marks[[i]], paste0(' . ~ .  -', namespde))
      
    }
    
    for (k in 1:length(family_marks)) {
      
      lhoods <- like(formula = formula_marks[[k]],
                     family = family_marks[k],
                     data = data_marks[[k]],
                     mesh = mesh,
                     ips = ips)
      likelihoods_marks[[k]] <- lhoods
      
      
    }
    for (l in 1:length(likelihoods_marks)) {
      
      #Better way to do this?
      likelihoods[[l+length(family)]] <- likelihoods_marks[[l]]
      
    }
    
  }
  
  ips$resp = 0
  proj4string(ips) = proj
  
  #Run integration points only on spatialcovariates?
  like_ip = like(formula = formula(paste0(c("resp ~", form_elements[2]),collapse = " ")),
                 family = 'poisson',
                 mesh = mesh,
                 E = ips$weight,
                 data = ips)
  
  likelihoods[[length(likelihoods) + 1]] = like_ip
  
  names(likelihoods) <- c(data_names,names_marks, 'like_ip')
  
  
  
  model_joint <- bru(components = components_joint,
                     likelihoods, options = list(control.compute = list(cpo = FALSE, dic = FALSE, waic = FALSE), control.inla = list(int.strategy = "eb")))
  #return(model_joint)
  
  if (residuals) {
    
    name_resp <- c()
    
    for (i in 1:(length(likelihoods) - 1)) {
      
      name_resp[i] <- gsub("\\(|\\)","",likelihoods[[i]]$formula[2])
      
    }
    
    fitted_residuals = list()
    for (i in 1:(length(likelihoods) - 1)) {
      res <- c()
      for (j in (model_joint$bru_iinla$inla_stack$data$index[[i]][1]):tail(model_joint$bru_iinla$inla_stack$data$index[[i]], n = 1)) {
        
        res[j] <- INLA::inla.emarginal(function(x) x, model_joint$marginals.fitted.values[[j]])
        
        
      }
      
      fitted_residuals[[i]] <- as.vector(na.omit(res))
      
    }
    
    
    residuals = list()
    
    for (k in 1:length(fitted_residuals)) {
      
      residuals[[k]] <- likelihoods[[k]]$data@data[,name_resp[k]] - fitted_residuals[[k]]
      
    }
    
    names(residuals) <- c(data_names,names_marks) #,'like_ip')
    #residuals[['like_ip']] <- NULL
    model_joint[['model_residuals']] = residuals
    
    
  }
  
  data_type <- sapply(c(data_attributes,data_marks), function(x) attributes(x)[['data_type']])
  names(data_type) <- c(data_names,names_marks)
  model_joint[['data_type']] <- data_type
  
  class(model_joint) <- c('bru_sdm',class(model_joint))
  return(model_joint)
  
  
}


##To do:

#Predictions of model.
#Can't just use predict(model, pixels, ~ exp(predictors))?
#Multiple likelihood problem.
#Make predictions for each distribution family?
#Would that mean running each distribution in a separate model?

#Modeling multiple species in a dataset.
#Find nice way to do this.
#Would this be done with a categorical variable?
#Add new parameter to model a single (or subset vector of species).
#Do this via attributes?

#Stack data
#Should I stack each dataset by family.
#Reason: for poisson each dataset has its own integration points.
#A lot slower than just stacking.
