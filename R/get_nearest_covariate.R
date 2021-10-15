#' Function to get covariate values nearest the data/integration point.
#' 
#' @param points A SpatialPoints* containing the locations where the covariate needs to be found.
#' @param spatialcovariates A SpatialPoints* containing the covariate values.
#' @param covariatestokeep A vector of covariate names to keep in the model. Defaults to \code{NULL}.
#' @param componentstokeep A vector of components to keep. Defaults to \code{c(poresp,paresp,'weight')}.
#' @param coords Vector of the names of the coordinates used in datasets.
#' @param proj Projection to use if data is not a projection.
#' @param attributestokeep A vector of attributes to keep from the points object. Defaults to \code{c('Ntrials','family','data_type')}.
#' 

get_nearest_covariate <- function(points, spatialcovariates, covariatestokeep = NULL,
                                  componentstokeep = c(poresp,paresp,'weight'),
                                  coords, proj,
                                  attributestokeep = c('Ntrials','family','data_type')) {
  
  if (!inherits(spatialcovariates,'Spatial')) stop("Spatialcovariates are required to be spatial.")
  
  if (!inherits(points,'Spatial')) stop("Points are required to be spatial.")
  
  atts <- attributes(points)[attributestokeep]

  coords_data <- as.data.frame(points)
  spatnames <- names(spatialcovariates@data)
  
  if (!is.null(covariatestokeep)) spatnames <- covariatestokeep
  
  colnames(spatialcovariates@coords) <- coords
  spatcovs <- cbind(spatialcovariates@coords,spatialcovariates@data)
  
  spatcovs$ID <- seq(nrow(spatcovs))
  
  closest <- RANN::nn2(spatcovs[coords], coords_data[coords], k = 1)
  
  closest <- as.data.frame(closest) %>%
  dplyr::rename(ID = nn.idx)
  
  joined <- dplyr::inner_join(closest,spatcovs, by = 'ID') %>%
  dplyr::select(-c(nn.dists,coords[1],coords[2],ID))
  
  allpts <- dplyr::bind_cols(coords_data, joined)
  
  points.df <- sp::SpatialPointsDataFrame(coords = allpts[coords],
                                          data = allpts[,names(allpts)%in%c(spatnames,componentstokeep)],
                                          proj4string = proj)
  
  missing_vals <- sapply(points.df@data, is.na)
  
  if (any(missing_vals)) {
    
  if (class(spatialcovariates) == 'SpatialPointsDataFrame') {
      
  spatialcovariates <- as(spatialcovariates, 'SpatialPixelsDataFrame')
  proj4string(spatialcovariates) <- proj
    
  }
  
  for (name in spatnames) { 
    
  points.df@data[,name] <- inlabru::bru_fill_missing(data = spatialcovariates,
                                                     where = points,
                                                     values = points.df@data[,name])
  
  }
  
  }
 
  attributes(points.df) <- append(attributes(points.df), atts)
  
  return(points.df)
  
  }