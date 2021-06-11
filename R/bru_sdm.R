#' Function to run marked species distribution models with data from presence only and presence absence data. (and then later account for over dispersion).
#'
#' @param ... Point process datasets with coordinates of species, and optionally marks and covariates explaining the coordinates.
#' @param spatialcovariates Data frame of the spatial covariates accompanied by their associated coordinates.
#' @param marks Should the model be a marked point process. Defaults to \code{FALSE}.
#' @param markfamily Assumed distribution of the marks. Defaults to \code{"gaussian"}.
#' @param inclmarks. A vector of which marks should be included in the model. Defaults to \code{NULL}.
#' @param speciespresence Name of column used if multiple species are implmented in the model. Defaults to \code{NULL}.
#' @param inclspecies A vector of which species should be included in the model. Defaults to \code{NULL}.
#' @param coords Vector of the names of the coordinates used in datasets. Defaults to \code{c('X','Y')} (For now should be standardized).
#' @param poresp Name for the response variable for the presence only datasets. Defaults to \code{NULL}. If no presence only response is found in dataset, a vector of 1's will be used.  
#' @param paresp Name for the response variable for the presence absence datasets. Defaults to \code{NULL}. Note that this column may also be logical.
#' @param trialname Names of column of number of columns in observs. Defaults to \code{NULL}.
#' @param inclcoords Should coordinates be used in data. Defaults to \code{FALSE}.
#' @param mesh An inla.mesh object. Defaults to \code{NULL}.
#' @param meshpars List of mesh parameters. Requires the following items: "cut.off", "max.edge" and "offset". Defaults to \code{NULL}.
#' @param ips Integration points. Defaults to \code{NULL}.
#' @param bdry Polygon of boundary for region, of class Polygon. If \code{NULL}, draws a boundary around the points.
#' @param proj Projection to use if data is not a projection. Defaults to utm (hopefully).
#' @param residuals Should residuals for each dataset be calculated. Defaults to \code{TRUE}.
#' @param predictions Boolean: should predictions (on the linear scale) be made? Defaults to \code{FALSE}.
#' @param intercept Include joint intercept in the model. Defaults to \code{FALSE}.
#' @param indvintercepts Include individual intercepts for each dataset in the model. Defaults to \code{TRUE}.
#' @param options A bru_options options object or a list of options passed on to bru_options()
#' @param pointsspatial Should spatial effects be used for the points in the model. Defaults to \code{TRUE}.
#' @param marksspatial Should spatial effects be used for the marks in the model. Defaults to \code{TRUE}.
#' @param poformula Formula given to the presence only datasets. Defaults to \code{NULL}.
#' @param paformula Formula given to the presence/absence datasets. Defaults to \code{NULL}.
#' @param tol Tolerance parameter for SpatialPixelsDataFrame. Defaults to \code{NULL}.
#' 
#' @import sp
#' @import INLA
#' @import inlabru
#' @import rgeos

bru_sdm = function(..., spatialcovariates, marks = FALSE, markfamily = 'gaussian',
                   inclmarks = NULL, speciespresence = NULL, inclspecies = NULL,
                   coords = c('X','Y'), poresp = NULL, paresp = NULL,trialname= NULL,
                   inclcoords = FALSE, mesh = NULL, meshpars = NULL, 
                   ips = NULL, bdry = NULL, proj = CRS("+proj=longlat +ellps=WGS84"),
                   predictions = FALSE, residuals = TRUE, intercept = FALSE,
                   indivintercepts = TRUE, pointsspatial = TRUE, marksspatial = TRUE, 
                   options = list(), poformula = NULL, paformula = NULL, tol = NULL) {
  
  
  if (is.null(spatialcovariates)) stop("Spatial covariates not provided.")
  
  if (is.null(poresp) | is.null(paresp)) stop("Either the precense only or the precense absence response is null.")
  
  if (poresp == paresp) stop("Please provide different names for presence only and presence absence datasets.")
  
  if (length(paresp) > 1 | length(poresp) > 1) stop("More than one name given for presences column.")
  
  if (length(trialname) > 1) stop("More than one name given for number of trials column.")
  
  if (length(coords) != 2) stop("Coordinates must have two components.")
  
  if(ncol(spatialcovariates) > 1 & is.null(tol)) stop("Tolerance parameter not provided.")
  
  if (as.character(proj@projargs) == "+proj=longlat +ellps=WGS84 +no_defs") warning("Default CRS is being used. Please change if incorrect.")
  
  if (!is.null(poformula)) {
    
    if (as.character(poformula[2]) != poresp) stop(paste("Response variable for presence only datasets should be: ", poresp,'.', sep = ""))
    
  }
  
  if(!is.null(paformula)) {
    
    if (as.character(paformula[2]) != paresp) stop(paste("Response variable for presence absence datasets should be: ", paresp,'.', sep = ""))    
  }
  
  if (class(spatialcovariates) == 'data.frame' & is.null(tol)) stop('Please provide a tolerance parameter to convert the spatial covariates to a SpatialPixelsDataFrame.')
  
  if (class(proj) != 'CRS') stop("Proj needs to be a CRS object.")
  
  if (is.null(mesh) & is.null(meshpars)) stop("Either a mesh or mesh parameters need to be provided.")
  
  if (is.null(mesh)) {
    
    if (sum(names(meshpars)%in%c( "cutoff", "max.edge", "offset")) < 3)
      
      stop("Meshpars requires three items in the list: cut.off, max.edge and offset.")
    
  }

  if (!marks) {
    
    names_marks <- NULL 
    data_marks <- NULL
    
    if (!is.null(speciespresence) & marksspatial)
    marksspatial <- TRUE
    } else markspatial <- FALSE

  if (!marks & !is.null(inclmarks)) {
    
    warning('Marks to include is non null but include marks is set to FALSE.\nMarks are now being set to TRUE.')
    marks <- TRUE
    
  }
  
  if (is.null(speciespresence) & !is.null(inclspecies)) {
    
    warning('Speciespresence is set to FALSE but a list of species to include is provided. Setting speciespresence to TRUE')
    speciespresence <- TRUE
    
  }
  
  if (is.null(speciespresence)) {
    
    names_all_species <- NULL 
    species_names <- NULL
    
  }
  
  datasets = list(...)
  
  coords_in = unlist(lapply(datasets, function(dat) {
    
    if (class(dat) == 'data.frame') coords%in%names(dat)
    else
    if (inherits(dat, 'Spatial')) {
      x_coord <- colnames(dat@coords)[1]
      y_coord <- colnames(dat@coords)[2]
      coords%in%c(x_coord,y_coord)
    }
    
  }))
  
  if (!all(coords_in)) stop("At least one dataset does not have coordinates in it.\nEither check your datasets or change your coordinates argument.")
  
  
  data_names <- setdiff(as.character(match.call(expand.dots=TRUE)), 
                        as.character(match.call(expand.dots=FALSE)))
  
  #Separate PO and PA data by inclusion/exclusion of 'trialname'.

  data_attributes <- lapply(datasets,function(dat) {
    if (inherits(dat,"Spatial")) {
      if (class(dat) == "SpatialPoints") {
        
        dat <- sp::SpatialPointsDataFrame(coords = sp::coordinates(dat),
                                         data = data.frame(resp = rep(1,nrow(coordinates(dat)))),
                                         proj4string = proj)
        names(dat) <- poresp
        attr(dat,'family') <- 'poisson'
        attr(dat,'data_type') <- 'Present only'
        dat
        
      } 
      else #if class == SpatialPointsDataFrame
        if (paresp%in%colnames(dat@data)){
          
          dat <- sp::SpatialPointsDataFrame(coords = dat@coords,
                                            data = as.data.frame(dat@data),
                                            proj4string = proj)
          dat@data[,paresp] <- as.numeric(dat@data[,paresp])
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
        if (!poresp%in%colnames(dat@data)) {
          dat@data[,poresp] <- 1
          
        }
        attr(dat,'family') <- 'poisson'
        attr(dat,'data_type') <- 'Present only'
        dat
        
      }
      
    }
    else #is not a spatial object
      if (class(dat) == 'data.frame') {
        if (paresp%in%colnames(dat)) {
          if (ncol(dat) == 1) stop("Only one column provided in data frame.\nNeed two columns for coordinates and one for presence name\nfor presence/absence data.")
          else
            if (ncol(dat) == 2) stop("Only two columns provided for presence/absence data. Either coordinates is of length one or presence name not given.")
          else {
            
            names <- names(dat)[!names(dat)%in%c(coords)]
            dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                              data = as.data.frame(dat[,!names(dat)%in%coords]),
                                              proj4string = proj) #dat[,!names(x)%in%coords]
            colnames(dat@data) <- names
            dat@data[,paresp] <- as.numeric(dat@data[,paresp])
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
              names(dat) <- poresp
              attr(dat,'family') <- 'poisson'
              attr(dat,'data_type') <- 'Present only'
              dat
              
            }
          else  {
            
            names <- names(dat)[!names(dat)%in%c(coords)]
            dat <- sp::SpatialPointsDataFrame(coords = dat[,coords],
                                              data = as.data.frame(dat[,!names(dat)%in%coords]),
                                              proj4string = proj)
            colnames(dat@data) <- names
            if (!poresp%in%colnames(dat@data)) {
              dat[,poresp] <- 1
            }
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

    data_attributes[[i]]@data[,coords] <- data_attributes[[i]]@coords
        
    }
    
  }
  
  #Issues with running species and marks separately
   #If marks = TRUE and speciespresence = FALSE and species is given as an available mark
   #Then species are run incorrectly
   #Could fix marks to run multinomial
   #But then I would need to get rid of the inclspecies variable.
  if (!is.null(speciespresence)) {
    
    data_species <- lapply(data_attributes, function(dat) {
      
      if (speciespresence%in%names(dat)) {
        
        if (is.null(inclspecies)) {
          
          if (attributes(dat)$data_type == 'Present only') {
            
          dat@data[,speciespresence] <- factor(dat@data[,speciespresence])
          dat@data[,'species_response'] <- dat@data[,poresp]
          
          if (length(unique(dat@data[,speciespresence])) > 1) {
            dat@data[,'phi'] <- rep(1,nrow(dat))
            attr(dat,'family') <- 'poisson'
          }
          #Seems to work with only 1 species
          #Make all poisson for now
          #Is this the correct way to do multinomial?
          else attr(dat, 'family') <- 'poisson'#'binomial'
          attr(dat,'multispecies') <- TRUE
          attr(dat,'species_included') <- dat@data[,speciespresence]
          dat 
          
          }
          
          else {
           #if data_type == Present absence  
            dat@data[,speciespresence] <- factor(dat@data[,speciespresence])
            dat@data[,'species_response'] <- dat@data[,paresp]
            dat@data[,'phi'] <- rep(1,nrow(dat))
            attr(dat, 'family') <- 'poisson'#'binomial'            
            attr(dat,'multispecies') <- TRUE
            attr(dat,'species_included') <- dat@data[,speciespresence]
            dat
            
          }
          
        } 
        else {
          #if include data is not null
          dat <- dat[dat@data[,speciespresence]%in%inclspecies,]
          if (length(dat) == 0) NULL
          
          else
            if (length(unique(dat@data[,speciespresence])) == 1) {
            #Multinomial seems to work with only one species, might remove later  
            dat@data[,speciespresence] <- factor(dat@data[,speciespresence])
           
            if (attributes(dat)$data_type == 'Present only') {
            dat@data[,'species_response'] <- dat@data[,poresp]
            }
            else dat@data[,'species_response'] <- dat@data[,paresp]
            dat@data[,'phi'] <- rep(1,nrow(dat))
            attr(dat, 'family') <- 'poisson'#'binomial'
            attr(dat,'multispecies') <- TRUE
            attr(dat,'species_included') <- dat@data[,speciespresence]
            dat
              
            }
          
          else {
            
            dat@data[,speciespresence] <- factor(dat@data[,speciespresence])
            dat@data[,'species_response'] <- dat@data[,poresp]
            dat@data[,'phi'] <- rep(1,nrow(dat))
            attr(dat, 'family') <- 'poisson'
            attr(dat,'multispecies') <- TRUE
            attr(dat,'species_included') <- dat@data[,speciespresence]
            dat 
            
            }
          }
      }
      else NULL
      
    })
    
    ##Throw error if all attributes(dat)$species_included is NULL
    ##I.e !is.null(speciespresence) but column name not found in any dataset
    
    data_species[sapply(data_species,is.null)] <- NULL
    
    if (length(data_species) == 0) stop("Species to include is given but species presence name column not given in any datasets.")
    
    names(data_species) <- paste0(names(data_species),'_',speciespresence)
    species_names <- names(data_species)
    
    if (!is.null(inclspecies)) {
      all_species_in <- lapply(data_species, function(dat) {
        
        if (attributes(dat)$multispecies) {
          if (length(dat) == 0) FALSE
          else "Multispecies included"
        }
        else FALSE
      })
      
      all_species_in[sapply(all_species_in,is.logical)] <- NULL
      
      if (length(all_species_in) == 0) stop('None of the species specified are found in any datasets.')
    }
    
    names_all_species <- unique(unlist(lapply(data_species, function(dat) {
      
      names <- as.character(attributes(dat)$species_included)
      names
      
    })))
    
    ##Remove 'speciespresence' from data_attributes??
    ##Not sure why but fixes issue
    data_attributes <- lapply(data_attributes, function(dat){
      
      if (speciespresence%in%names(dat)) {
      
        dat@data[,speciespresence] <- NULL
        dat
        
        }
      else dat
      
    })
    
  }
  
  if (marks) {
    
    data_marks = list()
    #Make a unique SpatialPointsDataframe for each mark (to be run on spatial covariates).
    #Incorporate standardized list of marks/covariates for data.
    #Include non numeric marks as well
    for (i in 1:length(data_attributes)) {
      
      ind <- 1 + length(data_marks)
      
      if (class(data_attributes[[i]]) == 'SpatialPoints') data_marks[[ind]] <- FALSE
      else {
        
        names = names(data_attributes[[i]])[!names(data_attributes[[i]])%in%c(poresp,paresp,coords,trialname,speciespresence)]
        
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
          if (length(names) > 1) { #I.e. more than one mark in a dataset
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
    names(data_marks) <- names_marks
  }
  
  #Do I need this?
  #Seems to work fine if I add it to the components

  #if (indivintercepts) {
  #  for (i in 1:length(data_attributes)) {
  #    
  #    data_attributes[[i]]@data[,paste0(names(data_attributes)[i],'_intercept')] <- 1
  #    
  #  }
  #  if (marks) {
  #    for (i in 1:length(data_marks)) {
  #    
  #    data_marks[[i]]@data[,paste0(names(data_marks)[i],'_intercept')] <- 1
  #    
  #  }
  #  }
  #  
  #}

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
  #Easiest way to fix is by spatdata <- as(rasterdata, 'SpatialPixelsDataFrame')
  #Is there loss in data this way??
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
  
  spde2 <- inla.spde2.matern(mesh)
  
  ##Construct joint components for the likelihoods.
  ##Will need to change with inclusion of separate covariates.
  
  if (is.null(poformula) | is.null(paformula)){
    
    components_joint <- formula(paste(c('~ 0',paste0(spatnames,'(main = ',spatnames,', model = "linear")')), collapse = '+'))
    
    if (inclcoords) {
      
      components_joint <- update(components_joint, paste0('~ . +',coords, collapse = '+'))
    }
    
    if (intercept) {
      
      components_joint <- update(components_joint, ~ . + Intercept(1))
    
      }
    
    
  }
  
  likelihoods = list()
  
  family <- unlist(sapply(data_attributes, function(x) attributes(x)$family))

  trials <- sapply(data_attributes, function(x){
    
    if (!is.null(attributes(x)$Ntrials)) attributes(x)$Ntrials
    else 1
    
  }) 

  ##Take out any brackets from 'components_joint'.
  ##I.e (for now) run coordinates only on spatial covariates (and optional others).
  form_elements <- gsub(" *\\(.*?\\) *", "",components_joint)

  formula <- mapply(function(fam,ind) {
    if (!is.null(poformula) & fam == 'poisson') {
      
      formula <- poformula
    }
    else
      if (is.null(poformula) & fam == 'poisson') {
        
        formula <- formula(paste0(c(poresp,'~', form_elements[2]),collapse = " ")) 
        
      }
    
    else 
      if(!is.null(paformula) & fam == 'binomial') {
        
        formula <- paformula
        
      }
    else
       if(is.null(paformula) & fam == 'binomial'){
      
      formula <- formula(paste0(c(paresp,"~", form_elements[2]),collapse = " "))
      
    }
    
    if (indivintercepts) {
       
        formula <- update(formula,paste0(' ~ . +', paste0(data_names[ind],'_intercept'), collapse = ' + '))
       
    }
    else formula
    
    if (pointsspatial) {
      
    formula <- update(formula, paste0('~ . +',data_names[[ind]],'_spde'))

    }
    else formula
    
  }, fam = family, ind = 1:length(family))
  
  for (i in 1:1) {

    lhoods <- like(formula = formula[[i]], ##Add tag to this likelihood somehow?
                   family = family[i],
                   data = data_attributes[[i]],
                   mesh = mesh,
                   ips = ips,
                   Ntrials = trials[i])
    likelihoods <- like_list(lhoods)
    
    
    if (length(family) > 1) { #Better way of doing this??
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
    
    for (i in 1:length(family_marks)) {
      
      formula_marks[[i]] <- formula(paste0(c(names_marks[i],'~',form_elements[2]),collapse = " "))
      
      if (marksspatial) {
        
        formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +", names_marks[i],'_spde'))
      
      }
      
      if (indivintercepts) {
        
        formula_marks[[i]] <- update(formula_marks[[i]],paste0(' ~ . +', paste0(names_marks[i],'_intercept'), collapse = ' + '))
        
      }
      
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
  
  if(!is.null(speciespresence)) {
    
    family_species <- unlist(sapply(data_species, function(x) attributes(x)$family))
    formula_species <- list()
    
    for (i in 1:length(family_species)){
    
    formula_species[[i]] <- formula(paste0(c('species_response','~',form_elements[2]),collapse = " "))
    
    if (family_species[[i]] == 'poisson') {
      
      formula_species[[i]] <- update(formula_species[[i]], '. ~ .  + phi') 
      formula_species[[i]] <- update(formula_species[[i]], paste('. ~ . +', speciespresence))
      
    }
    
    if (family_species[[i]] == 'binomial') { #Keep this for now
     
      formula_species[[i]] <- update(formula_species[[i]], paste('. ~ . +', unique(data_species[[i]]@data[,speciespresence])))
      
    }
    
    if (marksspatial) {
      
      formula_species[[i]] <- update(formula_species[[i]],paste0('. ~ . +',species_names[i],'_spde'))
    
      }
    
    }

    likelihoods_species <- list()

    for (k in 1:length(data_species)) {
      #Add NTrials above
      lhoods <- like(formula = formula_species[[k]],
                     family = family_species[k],
                     data = data_species[[k]],
                     mesh = mesh,
                     ips = ips)
      likelihoods_species[[k]] <- lhoods
      
  }
    n <- length(likelihoods)
    for (l in 1:length(likelihoods_species)) {
      
      #Better way to do this?
      likelihoods[[l + n]] <- likelihoods_species[[l]]
      
    }
    
  }
  ips$int_resp <- 0
  #ips <- spTransform(ips, proj) <- doesn't work if ips is not projected
  proj4string(ips) <- proj # <- is this fine?
  #Run integration points only on spatialcovariates?
  like_ip = like(formula = formula(paste0(c('int_resp ~ 0', c(spatnames)) ,collapse = '+')), #formula(paste0(c("resp ~", form_elements[2]),collapse = " ")),
                 family = 'poisson',
                 mesh = mesh,
                 E = ips$weight,
                 data = ips)
  
  likelihoods[[length(likelihoods) + 1]] = like_ip

  names(likelihoods) <- c(data_names,names_marks, species_names, 'like_ip')
  
  if (indivintercepts) {
    
    components_joint <- update(components_joint, paste0(' ~ . +', paste0(c(data_names,names_marks),'_intercept(1)'), collapse = ' + '))
    
  }
  
  if (pointsspatial) {
    
    components_joint <- update(components_joint, paste('. ~ . +',paste0(data_names,'_spde(main = coordinates, model = spde2)',collapse = ' + ')))
    
  }
  
  if (marksspatial) {
    if (!is.null(names_marks)) { #I.e. if marks is null but species is non null. Should I add a seperate random for species?
    components_joint <- update(components_joint, paste('. ~ . +',paste0(names_marks,'_spde(main = coordinates, model = spde2)',collapse = ' + ')))
    }
    
  }
  
  if (!is.null(speciespresence)) {
    if ('poisson'%in%family_species) {
      
   #components_joint <- update(components_joint, paste0('. ~ . +' ,speciespresence,'(main =',speciespresence,' ,model = "factor_full")'))
   components_joint <- update(components_joint, paste0('. ~ . +' ,speciespresence,'(main =',speciespresence,',model = "iid",constr = FALSE, fixed=TRUE)'))
   components_joint <- update(components_joint, ' . ~ .  + phi(main = phi, model = "iid", initial = -10, fixed = TRUE)')
   
    }
    
    if ('binomial'%in%family_species) { #Keep for now
      for (i in 1:length(data_species)) {
        if (attributes(data_species[[i]])$family == 'binomial') {
          
          components_joint <- update(components_joint, paste0(' . ~ . +',unique(attributes(data_species[[i]])$species_included),'(1)'))

        }
      }
  
    }
    
    if (marksspatial) {
      
      components_joint <- update(components_joint, paste(' . ~ . +',paste0(species_names,'_spde(main = coordinates, model = spde2)', collapse = ' + ')))
    }
  
  }
  
  #Add options parameter here
  model_joint <- bru(components = components_joint,
                     likelihoods) #, options = list(control.compute = list(cpo = FALSE, dic = FALSE, waic = FALSE), control.inla = list(int.strategy = "eb")))

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
  
  model_joint[['species_in_model']] <- names_all_species
  
  class(model_joint) <- c('bru_sdm',class(model_joint))
  return(model_joint)
  
  
}
