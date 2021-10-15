#' Function to create model matrix to run covariates by species.
#' 
#' @param datasets A list of datasets with a species variable
#' @param species Species variable name in the datasets.
#' @param covariates Spatial covariates used in the model.
#' @param covariatesbydataset A named list of covariates by dataset. Defaults to \code{NULL}.
#' @param proj Projection to use if data is not a projection.
#' @param attributestokeep A vector of attributes to keep from the points object. Defaults to \code{NULL}.
#' @param covariatestokeep A vector of covariate names to keep in the model. Defaults to \code{NULL}.
#' @param componentstokeep A vector of components to keep. Defaults to \code{NULL}.
#' @param coords Vector of the names of the coordinates used in datasets.
#' @param proj Projection to use if data is not a projection.
#' @param attributestokeep A vector of attributes to keep from the points object. Defaults to \code{NULL}.
#'

model_matrix_maker <- function(datasets, species, covariates,
                               covariatesbydataset = NULL,
                               covariatestokeep = NULL,
                               componentstokeep = NULL,
                               coords = NULL,
                               proj = NULL,
                               attributestokeep = NULL) {
  
  ##Will probably also need to run this for marks?
  
  for (dataset in names(datasets)) {
  
  species_in_dataset <- as.character(unique(datasets[[dataset]]@data[,species]))
  
  if (!is.null(covariates)) {  
      
  colnames(covariates@coords) <- coords  
    
  if (dataset%in%names(covariatesbydataset)) covariatestokeep <- covariatesbydataset[[dataset]]
  
  else covariatestokeep <- names(covariates)
  
  datasets[[dataset]] <- get_nearest_covariate(points = datasets[[dataset]],
                                                         spatialcovariates = covariates,
                                                         covariatestokeep = covariatestokeep,
                                                         componentstokeep = componentstokeep,
                                                         coords = coords,
                                                         proj = proj,
                                                         attributestokeep = attributestokeep)
  
  matrix_formula <- formula(paste(' ~ - 1', species,  paste0(species, ':', names(covariates), collapse = ' + '), sep = ' + '))
  
  }
    
  else {
    
  matrix_formula <- formula(paste(' ~ - 1', species, sep = ' + '))
  
  }  
  
  if (length(unique(datasets[[dataset]]@data[,species])) != 1) {
  
  species_matrix <- data.frame(model.matrix(matrix_formula, datasets[[dataset]]@data))

  unique_species <- as.character(unique(datasets[[dataset]]@data[,species]))

  names(species_matrix)[names(species_matrix)%in%paste0(species, unique_species)] <- paste0(unique_species, '_intercept')
  
  names(species_matrix)[names(species_matrix)%in%as.vector(outer(paste0(species,unique_species,'.'),names(covariates), FUN = 'paste0'))] <- as.vector(outer(paste0(unique_species,'_'),names(covariates), FUN = 'paste0'))
  
  species_matrix[, names(species_matrix)%in%names(covariates)] <- NULL

  datasets[[dataset]]@data <- dplyr::bind_cols(datasets[[dataset]]@data, species_matrix)
  
  }
  else {
  
  unique_species <- as.character(unique(datasets[[dataset]]@data[,species]))
  
  datasets[[dataset]]@data[, paste0(unique_species,'_intercept')] <- 1  
  
  names(datasets[[dataset]]@data)[names(datasets[[dataset]]@data) == covariatestokeep] <- paste0(unique_species,'_',covariatestokeep)
    
  }

  }
  
  return(datasets)

}

## After this:
 # Need to add all these created interactions to the relevant formulas/include components
  # probably not difficult: just need to add dataset:covariate name
  # Add all of these interactions to the components
  # Then go to the leave one out function and remove all dodge components

#' Function to create model matrix for the integration points.
#' @param ips Is the matrix being constructed for the integration points. Defaults to \code{False}.
#' @param covariates Spatial covariates used in the model.
#' @param species Name of the species variable name used in the model.
#' @param all_species If integration points are considered, what are the names of all the species. Defaults to \code{NULL}.
#' @param covariatestokeep A vector of covariate names to keep in the model. Defaults to \code{NULL}.
#' @param componentstokeep vector of components to keep. Defaults to \code{NULL}.
#' @param coords Vector of the names of the coordinates used in datasets.
#' @param proj Projection to use if data is not a projection.
#' @param attributestokeep A vector of attributes to keep from the points object. Defaults to \code{NULL}.
#' 
ips_model_matrix_maker <- function(ips, covariates, species,
                                   all_species,
                                   covariatestokeep = NULL,
                                   componentstokeep = NULL,
                                   coords,
                                   proj,
                                   attributestokeep = NULL) {
  
  ips_coords <- do.call(rbind, replicate(length(all_species), ips@coords, simplify = FALSE))
  
  ips_weight <- do.call(rbind, replicate(length(all_species), ips@data, simplify = FALSE))
  
  coords_rows <- nrow(ips@coords)
  ips_fact <- rep(all_species, each = coords_rows)
  
  ips_data <- data.frame(ips_weight, ips_fact)
  
  ips <- sp::SpatialPointsDataFrame(coords = ips_coords,
                                    data = ips_data,
                                    proj = proj,
                                    match.ID = FALSE)
  
  names(ips@data) <- c(names(ips_weight), species)
  
  if (!is.null(covariates)) {
    
  colnames(covariates@coords) <- coords  
  

  ips <- get_nearest_covariate(points = ips,
                               spatialcovariates = covariates,
                               covariatestokeep = covariatestokeep,
                               componentstokeep = componentstokeep,
                               coords = coords,
                               proj = proj,
                               attributestokeep = attributestokeep)
  
  for (cov in names(covariates)) {
  
  ips@data[,paste0(all_species,'_',cov)] <- ips@data[,cov]    
    
  }
  #matrix_formula <- formula(paste(' ~ - 1', species, paste0(species, ':', names(covariates), collapse = ' + '), sep = ' + '))
  
  }
  
  else {
    
  #matrix_formula <- formula(paste(' ~ - 1', species, sep = ' + '))
    
  }
  
  if (length(unique(ips@data[,species])) != 1) {
    
  #ips_matrix <- data.frame(model.matrix(matrix_formula, ips@data))
  
  unique_species <- as.character(unique(all_species))
  
  #names(ips_matrix)[names(ips_matrix)%in%paste0(species, unique_species)] <- unique_species
  
  ips@data[,paste0(unique_species,'_intercept')] <- 0
  
  ips@data[, names(ips@data)%in%names(covariates)] <- NULL
  
  #ips@data <- dplyr::bind_cols(ips@data, ips_matrix)
  
  } else {
    
  unique_species <- as.character(unique(all_species))
  
  ips@data[, paste0(unique_species,'_intercept')] <- 0  

  names(ips@data)[names(ips@data) == covariatestokeep] <- paste0(species, unique_species,'_',covariatestokeep)
    
  }
  
  ips@data[,species] <- rep(seq_len(length(all_species)), each = coords_rows)
  
  return(ips)
}