#' Function to organize data to input into bru_sdm
#' 
#' @param ... Point process datasets with coordinates of species, and optionally marks and covariates explaining the coordinates. May be the dataset objects as \code{SpatialPoints*} or \code{data.frame} objects, or may be a list of \code{SpatialPoints} or \code{data.frame} objects.
#' @param poresp Name for the response variable for the presence only datasets. Defaults to \code{NULL}. If no presence only response is found in dataset, a vector of 1's will be used.  
#' @param paresp Name for the response variable for the presence absence datasets. Defaults to \code{NULL}. Note that this column may also be logical.
#' @param trialname Name of column in data denoting the number of trials used in a binomial process for the points. Defaults to \code{NULL}.
#' @param marktrialname Name of column in data denoting the number of trials used in a binomial process for the points. Defaults to \code{NULL}.
#' @param coords Vector of the names of the coordinates used in datasets. Defaults to \code{c('X','Y')} (For now should be standardized).
#' @param proj Projection to use if data is not a projection. Defaults to utm (hopefully).
#' @param marks Should the model be a marked point process. Defaults to \code{FALSE}.
#' @param inclmarks. A vector of which marks should be included in the model. Defaults to \code{NULL}.
#' @param markfamily Assumed distribution of the marks. May be either a single character string or a named list/vector of each mark's distribution in the form: <mark name> = <distribution family>. Defaults to \code{"gaussian"}.
#' @param pointcovariates The names of the non-spatial covariates in the model that are attached to the dataset. Defaults to \code{NULL}.
#' @param speciesname Name of the species name variable used in the model. Defaults to \code{NULL}.
#' @param ips Integration points. Defaults to \code{NULL}.
#' @param mesh An inla.mesh object. Defaults to \code{NULL}.
#' @param meshpars List of mesh parameters. Requires the following items: "cut.off", "max.edge" and "offset". Defaults to \code{NULL}.
#' @param boundary Polygon of boundary for region, of class Polygon. If \code{NULL}, draws a boundary around the points.
#' 
#' @export


organize_data <- function(..., poresp = NULL, paresp = NULL,
                          trialname = NULL, marktrialname = NULL,
                          coords = NULL, proj = NULL,
                          marks = FALSE, inclmarks = NULL,
                          markfamily = 'gaussian', 
                          pointcovariates = NULL, speciesname = NULL,
                          ips = NULL, mesh = NULL, meshpars = NULL,
                          boundary = NULL) {
  
  if (is.null(poresp) | is.null(paresp)) stop('Both poresp and paresp must be non-null.')
  
  if (poresp == paresp) stop('Presence only response cannot be the same as present absence response. Please change one of them.')
  
  if (length(coords) != 2) stop('Coordinates needs to have two components')
  
  if (!is.null(inclmarks) & !marks) {
    
  warning('Marks is set to false but inclmarks is non-null. Setting marks to TRUE.')
  marks <- TRUE
    
  }
  
  if (is.null(mesh) & is.null(meshpars)) stop("Either a mesh or mesh parameters need to be provided.")
  
  if (is.null(mesh)) {
    
  if (sum(names(meshpars)%in%c( "cut.off", "max.edge", "offset")) < 3)
      
  stop("Meshpars requires three items in a list: cut.off, max.edge and offset.")
    
  }
  
  datasets <- list(...)
  
  datasets_class <- lapply(datasets, class)
  
  if (length(datasets_class) == 1 & class(datasets_class) == 'list') {

  datasets <- unlist(datasets)
  datasets_class <- lapply(datasets, class)
  data_list <- TRUE
  
  } else data_list <- FALSE

  if (any(!datasets_class%in%c('SpatialPointsDataFrame','SpatialPoints', 'data.frame'))) {
    
  stop('Datasets need to be either a SpatialPoints* object or a data frame.')
    
  }
  
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
  
  if (data_list) {
  
  data_names <-  setdiff(gsub('list[(]|[)]','',as.character(match.call(expand.dots=TRUE))), 
                        gsub('list[(]|[)]','',as.character(match.call(expand.dots=FALSE))))

  data_names <- unlist(strsplit(x = data_names, split = ', '))
  
  if (length(data_names) != length(datasets)) {
    
  warning('Issues with naming the datasets from a list. Will create generic dataset names.')
  data_names <- paste0('dataset_',seq_len(length(datasets))) ## Fix this  
    
  }
  
  }
  
  else {
    
  data_names <- setdiff(as.character(match.call(expand.dots=TRUE)), 
                        as.character(match.call(expand.dots=FALSE)))
  }
  
  data_points <- lapply(datasets, function(data) {
    
  data <- as.data.frame(data)
  data_vars <- names(data)
    
  if (paresp%in%data_vars) {
      
  dat <- sp::SpatialPointsDataFrame(coords = data[,coords],
                                    data = data.frame(data[,!data_vars%in%coords]),
                                    proj4string = proj)
      
  if (ncol(dat@data) == 1) names(dat@data) <- paresp
      
  if (!is.null(trialname)) {
        
  if (trialname%in%data_vars) attr(dat, 'Ntrials') <- data.frame(dat@data[,trialname])[,1]
  else attr(dat, 'Ntrials') <- NULL
        
  }
      
  attr(dat, 'family') <- 'binomial'
  attr(dat, 'data_type') <- 'Present absence'
  attr(dat, 'Ntrials') <- NULL
      
  dat
      
  }
  else
  if (poresp%in%data_vars) {
        
  dat <- sp::SpatialPointsDataFrame(coords = data[,coords],
                                    data = data.frame(data[,!data_vars%in%coords]),
                                    proj4string = proj)
        
  if (ncol(dat@data) == 1) names(dat@data) = poresp
        
  if (any(dat@data[,poresp] > 1)) {
          
  attr(dat, 'Ntrials') <- NULL
  attr(dat, 'family') <- 'poisson'
  attr(dat, 'data_type') <- 'Present only abundance'
          
  dat
          
  }
  else {
          
  attr(dat, 'Ntrials') <- NULL
  attr(dat, 'family') <- 'cp'
  attr(dat, 'data_type') <- 'Present only' 
          
  dat
          
  }
        
  }
  else {
      
  data[,poresp] <- rep(1,nrow(data))
  data_vars <- c(data_vars,poresp)
      
  dat <- sp::SpatialPointsDataFrame(coords = data[,coords],
                                    data = data.frame(data[,!data_vars%in%coords]),
                                    proj4string = proj)
      
  if (ncol(dat@data) == 1) names(dat@data) <- poresp
      
  attr(dat, 'Ntrials') <- NULL
  attr(dat, 'family') <- 'cp'
  attr(dat, 'data_type') <- 'Present only'
      
  dat
      
  }
    
  })

  if (!is.null(speciesname)) {
    
  all_species <- sapply(data_points, function(dat) {
      
  speciesname%in%names(dat@data)
      
  })
    
  if (!all(all_species)) stop('All datasets are required to have the species name variable included.')
    
  data_points <- lapply(data_points, function(dat) {
      
  dat@data[,speciesname] <- factor(dat@data[,speciesname])
  dat
      
  })
    
  }
  
  if (!is.null(pointcovariates)) {
    
  all_covariates <- sapply(data_points, function(dat) {
      
  pointcovariates%in%names(dat@data)
      
  })
   
  if (!all(all_covariates)) stop('All datasets are required to have all the points covariates included.')
    
  }
  
  names(data_points) <- data_names
  
  if (marks) {
    
  data_marks = list()
    
  for (i in 1:length(data_points)) {
      
  ind <- 1 + length(data_marks)
      
  if (class(data_points[[i]]) == 'SpatialPoints') data_marks[[ind]] <- FALSE
      
  else {
        
  names = names(data_points[[i]])[!names(data_points[[i]])%in%c(poresp, paresp, coords, trialname,
                                                                speciesname, marktrialname, pointcovariates)]
        
  if (!is.null(inclmarks)) names <- names[names%in%inclmarks]      
        
  class_marks <- sapply(data_points[[i]]@data[names], class)
        
  if (is.null(names) | identical(names,character(0))) data_marks[[ind]] <- FALSE
        
  else
          
  for(j in 1:length(names)) {
            
  index <- ind + j - 1
            
  if (class_marks[j] == 'character'| class_marks[j] == 'factor') {
              
  if (attributes(data_points[[i]])$family == 'cp' | attributes(data_points[[i]])$family == 'poisson') mark_response <- data_points[[i]]@data[,poresp]
              
  else mark_response <- data_points[[i]]@data[,paresp]
              
  mark <- sp::SpatialPointsDataFrame(coords = coordinates(data_points[[i]]),
                                     data = data.frame(factor((data_points[[i]]@data[,names[j]]))),
                                     proj4string = proj)
  colnames(mark@data) <- names[j]
              
  mark@data[,paste0(names[j],'_phi')] <- rep(1,nrow(mark@coords))
  mark@data[,paste0(names[j],'_response')] <- mark_response
  mark@data[,'mark_response_weights'] <- mark_response
  mark@data[,'placeholder_for_factor'] <- mark@data[,names[j]]
              
  weights <- mark@data %>% 
             dplyr::group_by(placeholder_for_factor) %>% 
             dplyr::summarize(n_factor = sum(mark_response_weights)) %>%
             dplyr::mutate(weight = sum(n_factor)/n_factor) %>%
             dplyr::select(weight) %>%
             data.frame()
              
  attr(mark,'weights') <- weights$weight[as.numeric(mark@data[,names[j]])]
  mark@data[,'placeholder_for_factor'] <- NULL
  mark@data[,'mark_response_weights'] <- NULL
  attr(mark, 'Ntrials') <- rep(1,nrow(mark@coords))
  attr(mark,'response') <- paste0(names[j],'_response')
  attr(mark,'family') <- 'poisson'
  attr(mark,'data_type') <- 'Multinomial mark'
  attr(mark,'mark_name') <- names[j]
  attr(mark, 'phi') <- paste0(names[j],'_phi')
  attr(mark,'dataset') <- names(data_points)[i]
              
  data_marks[[index]] <- mark
  names(data_marks)[[index]] <- paste0(names(data_points)[i],'_',names[j])
              
  }
  else
              
  if (class_marks[j] == 'numeric' | class_marks[j] == 'integer') {
    
  mark <- sp::SpatialPointsDataFrame(coords = coordinates(data_points[[i]]),
                                     data = as.data.frame(data_points[[i]]@data[,names[j]]),
                                     proj4string = proj)
  colnames(mark@data) <- names[j]
                
  if (!is.null(names(markfamily))) {
                  
  if (names[j]%in%names(markfamily)) {
                    
  attr(mark,'family') <- markfamily[names(markfamily) == names[j]] 
  capital_markfamily <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2", 
                             markfamily[names(markfamily) == names[j]], perl = TRUE)
  attr(mark,'data_type') <- paste0(capital_markfamily,' mark')
                    
  }
  else {  
                    
  warning(names[j], ' has not been assigned a family. Will assign it "gaussian"')  
  attr(mark,'family') <- 'gaussian'
  attr(mark,'data_type') <- 'Gaussian'    
                    
  }
                  
  }
  else {
                  
  attr(mark,'family') <- markfamily
                  
  capital_markfamily <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2", 
                             markfamily, perl = TRUE)
  attr(mark,'data_type') <- paste0(capital_markfamily,' mark')
                  
  }
                
  attr(mark,'response') <- names[j]
  
  if (!is.null(marktrialname)) {
      
  if (attributes(mark)$family == 'binomial') {
  
  if (marktrialname%in%names(data_points[[i]]@data)) {
    
  attr(mark, 'Ntrials') <- data.frame(data_points[[i]]@data[,marktrialname])[,1]
    
  }
  else attr(mark, 'Ntrials') <- rep(1, nrow(mark@data))
    
  }
  else attr(mark, 'Ntrials') <- rep(1, nrow(mark@data))
    
  }
  else attr(mark, 'Ntrials') <- rep(1, nrow(mark@data))
    
  
  attr(mark,'mark_name') <- names[j]
  attr(mark,'phi') <- NULL
  attr(mark, 'weights') <- rep(1,nrow(mark@coords))
  attr(mark,'dataset') <- names(data_points)[i]
  data_marks[[index]] <- mark
  names(data_marks)[[index]] <- paste0(names(data_points)[i],'_',names[j])
               
  }
  else FALSE
            
  }
        
  }
      
  }
  data_marks[sapply(data_marks,is.logical)] <- NULL
    
  if (length(data_marks) == 0) stop("Either marks have been set to TRUE and no datasets contain marks, or marks to include only contains marks not present in any dataset.")
    
  names_marks <- sapply(data_marks, function(mark) attributes(mark)$mark_name)
    
  response_marks <- sapply(data_marks, function(dat) attributes(dat)$response)
    
  datasets_numeric_marks <- unlist(sapply(data_marks, function(mark) {
      
  if (attributes(mark)$data_type != 'Multinomial mark') attributes(mark)$dataset
      
  }))
    
  multinom_incl <- sapply(data_marks, function(mark) attributes(mark)$data_type == 'Multinomial mark')
    
  phi_vars <- sapply(data_marks, function(mark) attributes(mark)$phi)
    
  if (any(multinom_incl)) {
      
  multinom_vars <- unique(unlist(sapply(data_marks, function(mark) {
        
  if(attributes(mark)$data_type == 'Multinomial mark') attributes(mark)$mark_name
        
  })))
      
  datasets_multinom_marks <- unlist(sapply(data_marks, function(mark) {
        
  if (attributes(mark)$data_type == 'Multinomial mark') attributes(mark)$dataset
        
  }))
      
  data_points <- lapply(data_points, function(dat){
        
  if (any(multinom_vars%in%names(dat))) {
          
  dat@data[,multinom_vars] <- NULL
  dat
          
  }
  else dat
        
  })
      
  }
    
  else {
      
  multinom_vars <- NULL
  datasets_multinom_marks <- NULL
      
  }
    
  if (all(multinom_incl)) {
      
  datasets_numeric_marks <- NULL
  marksintercept <- FALSE
      
  }
    
  }
  else {
    
  data_marks <- NULL
  multinom_incl <- NULL
  multinom_vars <- NULL
  response_marks <- NULL
    
  }
  
  if (is.null(mesh)) {
    
  warning("Mesh not provided. Will try to create own mesh.")
    
  #Make mesh same way as PointedSDMs
  if (is.null(boundary)) {
      
  if (inherits(spatialcovariates, "Spatial")) spatcoords <- sp::SpatialPoints(coords = spatialcovariates@coords,
                                                                              proj4string = proj)
  else spatcoords <- sp::SpatialPoints(coords = spatialcovariates[,coords],
                                       proj4string  = proj)
      
  bstart <- min(c(diff(sort(unique(spatcoords@coords[,1]))), diff(sort(unique(spatcoords@coords[,2])))))
      
  poly.tmp <- rgeos::gBuffer(spatcoords, width=bstart, byid=TRUE)
      
  boundary <- rgeos::gBuffer(rgeos::gUnaryUnion(poly.tmp), width=bstart)
      
  }
    
  else {
  if (class(boundary)!="SpatialPolygons") {
        
  boundary <- sp::SpatialPolygons(Srl=list(Polygons(srl=list(boundary), ID="eek")))
        
  } 
  else {
        
  if (!is.projected(boundary)) boundary <- spTransform(boundary, CRSobj = proj)
        
  }
      
  }
    
  region.bdry <- inla.sp2segment(boundary)
    
  mesh <- inla.mesh.2d(boundary=region.bdry, 
                       cutoff=meshpars$cutoff,
                       max.edge=meshpars$max.edge, 
                       offset=meshpars$offset)
  
  mesh$crs <- proj
    
  }
  
  if (is.null(ips)) {
    
  warning('Integration points not provided. Will try to create own points')
    
  ips <- inlabru::ipoints(samplers = boundary,
                          domain = mesh)
    
  }
  
  colnames(ips@coords) <- coords
  proj4string(ips) <- proj
  
  
  family_index <- sapply(data_points, function(x) attributes(x)$family == 'poisson' | attributes(x)$family == 'cp')
  
  if (any(family_index) & !all(family_index)) {
    
  PO_data <- data_points[family_index]
  PA_data <- data_points[!family_index]
    
  }
  else
    
  if (all(family_index)) {
      
  PO_data <- data_points
  PA_data <- NULL
      
  }
  else {
    
  PO_data <- NULL
  PA_data <- data_points
  
  }
  
  object <- new('bru_sdm_data',
                PO_data = PO_data,
                PA_data = PA_data,
                Mark_data = data_marks,
                ips = ips,
                mesh = mesh)
  
  attr(object,'Points_response') <- c(poresp, paresp)
  attr(object,'Points_family') <- sapply(data_points, function(dat) attributes(dat)$family)
  attr(object,'Points_trials') <- lapply(data_points, function(dat) attributes(dat)$Ntrials)
  attr(object,'Points_data_type') <- sapply(data_points, function(dat) attributes(dat)$data_type)
  
  attr(object,'Marks') <- marks
  attr(object,'Mark_family') <- sapply(data_marks, function(dat) attributes(dat)$family)
  attr(object,'Mark_data_type') <- sapply(data_marks, function(dat) attributes(dat)$data_type)
  attr(object,'Mark_name') <- sapply(data_marks, function(dat) attributes(dat)$mark_name)
  attr(object,'Mark_phi') <- sapply(data_marks, function(dat) attributes(dat)$phi)
  attr(object,'Mark_dataset') <- sapply(data_marks, function(dat) attributes(dat)$dataset)
  attr(object,'Mark_weight') <- lapply(data_marks, function(dat) attributes(dat)$weight)
  attr(object,'Mark_response') <- response_marks
  attr(object, 'Mark_trials') <- lapply(data_marks, function(dat) attributes(dat)$Ntrials)
  attr(object, 'Multinom_incl') <- multinom_incl
  attr(object, 'Multinom_vars') <- multinom_vars
  attr(object, 'Sources_of_information') <- unname(c(data_names, sapply(data_marks, function(dat) attributes(dat)$dataset)))
  
  attr(object, 'Pointcovariates') <- pointcovariates
  
  attr(object, 'Species') <- speciesname
  
  return(object)
  
}
