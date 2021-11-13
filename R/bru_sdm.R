#' Function to run marked species distribution models with data from marked presence only and presence absence data.
#' 
#' @param data A bru_sdm_data object created with \code{organize_data}.
#' @param spatialcovariates Data frame of the spatial covariates accompanied by their associated coordinates. Defaults to \code{NULL}.
#' @param covariatestoinclude A vector of spatial covariate names to include in the model. Defaults to \code{NULL}.
#' @param covariatesbydataset A named list which includes which covariates are modeled to each dataset. Defaults to \code{NULL}.
#' @param specieseffects Calculate effects for the species. Defaults to \code{FALSE}.
#' @param pointsintercept Include individual intercepts for each point process in the model. Defaults to \code{TRUE}.
#' @param marksintercept Include individual intercepts for each mark process in the model. Defaults to \code{TRUE}.
#' @param spatialdatasets A vector of which datasets have spatial effects. Defaults to \code{NULL} which implies all datasets have spatial effects.
#' @param spdemodel inla.spde model used in the model. May be a named list where the name of the spde object is the name of the associated dataset. Default \code{NULL} uses \code{inla.spde2.matern}.
#' @param pointsspatial Should spatial effects be used for the points in the model. Defaults to \code{TRUE}.
#' @param marksspatial Should spatial effects be used for the marks in the model. Defaults to \code{TRUE}.
#' @param sharedspatial Should a spatial effect be shared across datasets. Defaults to \code{FALSE}.
#' @param speciesmodel INLA \code{control.group} model to use. Defaults to \code{list(model = "exchangeable")}.
#' @param tolerance Tolerance parameter to convert covariates to SpatialPixelsDataFrame. Defaults to \code{NULL}.
#' @param options INLA or inlabru options to be used in the model.
#' 
#' @export

bru_sdm <- function(data, spatialcovariates = NULL, covariatestoinclude = NULL,
                    covariatesbydataset = NULL, specieseffects = FALSE, pointsintercept = TRUE,
                    marksintercept = TRUE, sharedspatial = FALSE, spdemodel = NULL, 
                    pointsspatial = TRUE, marksspatial = TRUE,
                    spatialdatasets = NULL, speciesmodel = list(model = "exchangeable"),
                    tolerance = NULL, options = list()) {
  
  ##Add another argument called ** EXTRA COVARIATES **
  ##ie non-spatial; non-marks values attached to the datasets
  ##That act as additional covariates within the model
  ##e.g. effort covariate and coordinates
  ##Then make go to bru_sdm and add 0's onto the 
  ##ips such that they are modeled as well.
  ##Maybe just call it pointcovariates??
  
  #if (!is.null(attributes(data)$Pointcovariates)) point_covs_incl <- TRUE
  #else point_covs_incl <- FALSE
  
  ##Go to model matrix maker and say if point_covs_incl add those extra cols??
  ##Or just do it completely sep; prob easier to just add columns of 0's afterwards
  
  ##Then add to formulas
  ##Then add to components_joint

  if (class(data)[1] != 'bru_sdm_data') stop('Please supply data formed by the "organize_data" function.')

  proj <- data@ips@proj4string
  ##Change covariate coords name prior to running GNC
  coords <- colnames(data@ips@coords)
  
  data_points <- append(data@PO_data, data@PA_data)
  data_names <- names(data_points)
  points_family <- sapply(data_points, function(data) attributes(data)$family)
  points_response <- attributes(data)$Points_response
  points_response_variables <- sapply(data_points, function(data) {
    
    if(attributes(data)$family == 'cp') return('coordinates')
    else
    if(attributes(data)$family == 'binomial') return(points_response[2])
    else return(points_response[1])
    
  })

  if (!is.null(spatialdatasets)) {
    
  if (!all(spatialdatasets%in%data_names)) stop('At least one of the datasets speciefied to include spatial effects were not included in the model.')  
    
  }
  
  if (attributes(data)$Marks) {
    
  data_marks <- data@Mark_data
  names_marks <- names(data_marks)
  family_marks <- attributes(data)$Mark_family
  mark_weights <- attributes(data)$Mark_weight
  response_marks <- attributes(data)$Mark_response
  multinom_incl <- attributes(data)$Multinom_incl
  multinom_vars <- attributes(data)$Multinom_vars
  marks_name <- sapply(data_marks, function(x) attributes(x)$mark_name)
  mark_type <- sapply(data_marks, function(x) attributes(x)$data_type)
  mark_phi <- sapply(data_marks, function(x) attributes(x)$phi)
    
  }
  else {
    
  names_marks <- NULL
  data_marks <- NULL
  multinom_incl <- NULL
  multinom_vars <- NULL
  response_marks <- NULL
    
  }
  
  pointcovariates_incl <- attributes(data)$Pointcovariates
  
  if (!is.null(covariatesbydataset)) {
    
  if (any(!names(covariatesbydataset)%in%data_names)) stop('covariatesbydataset includes a dataset not available')  
  
  covs_for_all_datasets <- unique(unlist(covariatesbydataset))
  
  if (!all(covs_for_all_datasets%in%names(spatialcovariates))) stop('Covariates supplied to dataset are not included in the spatialcovariates object.')
    
  }
  
  if (!is.null(spatialcovariates)) {
    
  if (class(spatialcovariates) == 'RasterLayer' | class(spatialcovariates) == 'RasterBrick' | class(spatialcovariates) == 'RasterStack') {
      
  spatialcovariates <- as(spatialcovariates, 'SpatialPixelsDataFrame')
      
  }
    
  if (class(spatialcovariates) == 'data.frame') {
      
  spatialcovariates <- sp::SpatialPointsDataFrame(coords = spatialcovariates[,coords],
                                                  data = spatialcovariates[,!names(spatialcovariates)%in%coords],
                                                  proj4string = proj)
  
  if (is.null(tolerance)) {
      
  spatialcovariates <- as(spatialcovariates, 'SpatialPixelsDataFrame')
  
  }
  else {
    
  spatialcovariates <- sp::SpatialPixelsDataFrame(points = spatialcovariates@coords,
                                                  data = spatialcovariates@data,
                                                  proj4string = proj, tolerance = tolerance)  
    
  }
      
  }
  
  spatnames <- names(spatialcovariates@data)
  spatdata_class <- sapply(spatialcovariates@data, class)
 
  if (!is.null(covariatestoinclude)) {
      
  spatdata_class <- spatdata_class[spatnames%in%covariatestoinclude] 
  spatnames <- spatnames[(spatnames%in%covariatestoinclude)]
      
  if (is.null(spatnames) | identical(spatnames,character(0))) stop('covariatestoinclude contains covariate names not found in spatialcovariates')
      
  }
    
  for (name in spatnames) {
  
  if (is.null(tolerance)) {
        
  pixels_df <- sp::SpatialPixelsDataFrame(points = spatialcovariates@coords,
                                          data = data.frame(spatialcovariates@data[,name]),
                                          proj4string = proj)
  
  } 
  else {
    
  pixels_df <- sp::SpatialPixelsDataFrame(points = spatialcovariates@coords,
                                          data = data.frame(spatialcovariates@data[,name]),
                                          proj4string = proj, tolerance = tolerance)  
    
  }
    
  names(pixels_df) <- name
  assign(name,pixels_df)
      
  }
    
  }
  else spatnames <- NULL
  
  species <- attributes(data)$Species
  
  if (!is.null(species)) {
  
  if (specieseffects) {  
    
  species_dataset <- lapply(data_points, function(data) {
    
  data@data[,species]  
    
  }) 
  
  all_species <- unlist(species_dataset)
  
  numeric_species <- as.numeric(all_species)
  
  data_points <- model_matrix_maker(datasets = data_points, species = species,
                                    covariates = spatialcovariates,
                                    allspecies = as.character(unique(all_species)),
                                    componentstokeep = c(points_response, species, 'weight', pointcovariates_incl),
                                    coords = coords,
                                    attributestokeep = c('Ntrials', 'data_type'),
                                    covariatesbydataset = covariatesbydataset,
                                    proj =  proj)

  
  for (k in 1:length(data_points)) {
      
  if (k == 1) { 
        
  length_var <- (1:length(data_points[[1]]))
        
  }
  else {
        
  length_var <- (length(data_points[[k-1]]) + 1):(length(data_points[[k-1]]) + length(data_points[[k]]))
        
  }
    
  data_points[[k]]@data[,species] <- numeric_species[length_var]
  
  }
  
  data@ips <- ips_model_matrix_maker(ips = data@ips, covariates = spatialcovariates, allspecies = as.character(unique(all_species)),
                                     coords = coords, proj =  proj,
                                     species = species, componentstokeep = c(points_response, species, 'weight'))
  }
    
  } else specieseffects <- FALSE
  
  if (!is.null(pointcovariates_incl)) {
  
  pointcovariates_dataframe <-  data.frame(matrix(data=0, nrow = nrow(data@ips@data), ncol = length(pointcovariates_incl)))
  names(pointcovariates_dataframe) <- pointcovariates_incl
  
  data@ips@data <- cbind(data@ips@data, pointcovariates_dataframe)
    
  }
  
 
  if (is.null(spdemodel)) {
    
  spdemodel <- inla.spde2.matern(data@mesh)
    
  }
  else
    
  if (is.list(spdemodel[[1]])) { 
      
  if (length(names(spdemodel)) > 1) {
        
  if (is.null(names(spdemodel))) stop('Please provide a named list of spatial objects where the name of the object is the associated datasets name.') 
        
  if (length(names(spdemodel)) < length(names(data_points))) {
          
  names_out <- names(data_points)[!names(data_points)%in%names(spdemodel)]
          
  for (name in names_out)
            
  spdemodel[[name]] <- inla.spde2.matern(data@mesh)   
          
  }
        
  if (length(names(spdemodel)) == length(names(data_points))) {
          
  if (!all(names(spdemodel)%in%names(data_points))) stop('Names provided in spdemodel are not the same as the dataset names.')  
          
  }    
        
  }
      
  for (name in names(spdemodel)) {
        
  assign(paste0(name,'_spde'),spdemodel[[name]])  
        
  }  
      
  }  
  
  components_joint <- formula( ~ - 1)
  
  if (!is.null(spatialcovariates)) {
    
  for (cov in 1:length(spatdata_class)) {
      
  if (spatdata_class[cov] == 'numeric' | spatdata_class[cov] == 'integer') {
        
  components_joint <- update(components_joint, paste(c(' ~ . +', paste0(spatnames[cov],'(main = ', spatnames[cov], ', model = "linear")'))))
        
  }
  else
        
  if (pointsintercept | marksintercept) {
          
  components_joint <- update(components_joint, paste(c(' ~ . +', paste0(spatnames[cov],'(main = ', spatnames[cov], ', model = "factor_contrast")'))))
          
  } 
  else {
        
  components_joint <- update(components_joint, paste(c(' ~ . +', paste0(spatnames[cov],'(main = ', spatnames[cov], ', model = "factor_full")'))))
        
  }
      
  }
    
  }
  
  #form_elements <- gsub(" *\\(.*?\\) *", "",components_joint)
  formula <- formula_maker(response = points_response_variables, dataset = data_names,
                           covariates = spatnames, pointsspatial = pointsspatial,
                           pointsintercept = pointsintercept,
                           spatialdatasets = spatialdatasets, covariatesbydataset = covariatesbydataset,
                           species = species, speciesdataset = species_dataset,
                           pointcovariates = pointcovariates_incl)
  
  include <- list()
 
  for (i in 1:length(formula)) {
    
  variables <- all.vars(formula[[i]])
  include[[i]] <- variables[!variables%in%c(points_response,'coordinates')]
  formula[[i]] <- as.formula(paste(variables[!variables%in%include[[i]]], '~ .'))
    
  }
  
  for (i in 1:1) {
    
  lhoods <- inlabru::like(formula = formula[[i]],
                          family = points_family[i],
                          data = data_points[[i]],
                          mesh = data@mesh,
                          ips = data@ips,
                          Ntrials = attributes(data_points[[i]])$Ntrials,
                          include = include[[i]])
    
  likelihoods <- inlabru::like_list(lhoods)
    
  if (length(points_family) > 1) { #Better way of doing this??
      
  for (j in 2:length(points_family)) {
        
  lhoods <- inlabru::like(formula = formula[[j]],
                          family = points_family[j],
                          data = data_points[[j]],
                          mesh = data@mesh,
                          ips = data@ips,
                          Ntrials = attributes(data_points[[j]])$Ntrials,
                          include = include[[j]])
        
  likelihoods[[j]] <- lhoods
        
  }
      
  }
    
  likelihoods
    
  }
  
  if (attributes(data)$Marks) {
  
  likelihoods_marks <- list()
  
  formula_marks <- mark_formula_maker(responsemarks = response_marks, namesdatasetmarks = names_marks,
                                      marksname = marks_name, covariates = spatnames,
                                      covariatesbydataset = covariatesbydataset, marksspatial = marksspatial,
                                      sharedspatial = sharedspatial, spatialdatasets = spatialdatasets,
                                      marksintercept = marksintercept, marktype = mark_type, markphi = mark_phi)
  
  include_marks <- list()
    
  for (i in 1:length(formula_marks)) {
      
  variables <- all.vars(formula_marks[[i]])
  include_marks[[i]] <- variables[!variables%in%c(response_marks,coords)]
  formula_marks[[i]] <- as.formula(paste(response_marks[i], '~ .'))
      
  }
    
  for (k in 1:length(family_marks)) {
      
  lhoods <- inlabru::like(formula = formula_marks[[k]],
                          family = family_marks[[k]],
                          data = data_marks[[k]],
                          mesh = data@mesh,
                          Ntrials = attributes(data_marks[[k]])$Ntrials,
                          ips = data@ips,
                          E = mark_weights[[k]],
                          include = include_marks[[k]])
  likelihoods_marks[[k]] <- lhoods
      
  }
    
  n <- length(likelihoods)
    
  for (l in 1:length(likelihoods_marks)) {
      
  likelihoods[[l + n]] <- likelihoods_marks[[l]]
      
  }
    
  }
  
  if (pointsspatial) {
    
  if (sharedspatial) {
    
  if (is.list(spdemodel[[1]])) stop('Shared spatial model should only have one spde model.')  
    
  components_joint <- update(components_joint, paste(' ~ . +','shared_spatial(main = coordinates, model = spdemodel)'))  
    
  }  
  else  
  if (is.list(spdemodel[[1]])) { 
      
  for (name in names(data_points)) {
  
  if (!is.null(spatialdatasets)) {
      
  if (name%in%spatialdatasets) {
           
  components_joint <- update(components_joint, paste(' ~ . +',paste0(name,'_spde(main = coordinates, model =',name,'_spde)')))
        
  }
    
  }
  else {
    
  components_joint <- update(components_joint, paste(' ~ . +',paste0(name,'_spde(main = coordinates, model =',name,'_spde)')))
    
  }
  
  }  
      
  } 
  else {  
  
  if (!is.null(spatialdatasets)) {
    
  components_joint <- update(components_joint, paste(' ~ . +',paste0(spatialdatasets,'_spde(main = coordinates, model = spdemodel)',collapse = ' + ')))
    
  }
  else {  
  
  components_joint <- update(components_joint, paste(' ~ . +',paste0(data_names,'_spde(main = coordinates, model = spdemodel)',collapse = ' + ')))
      
  }
  
  }    
    
  }  
  
  if (pointsintercept) {
    
  if (specieseffects) {
    
  if (length(unique(all_species)) > 1) {  
    
  components_joint <- update(components_joint, paste0(' ~ . +', paste(paste0(unique(all_species),'_intercept'), collapse = ' + ')))
  
  }
  else components_joint <- update(components_joint, paste0(' ~ . + intercept(1)'))
    
  }
  else {  
  components_joint <- update(components_joint, paste0(' ~ . +', paste0(data_names,'_intercept(1)'), collapse = ' + '))
  
  }
    
  }
  
  if (specieseffects) {
    
  if (length(unique(all_species)) > 1) {  

  if (!is.null(spatnames)) {  
    
  for (name in data_names) {
    
  if (!is.null(covariatesbydataset)) {
      
  if (name%in%names(covariatesbydataset)) {
        
  incl_cov <- covariatesbydataset[[name]]
        
  } else {
        
  incl_cov <- spatnames
        
  }
      
  } 
  else incl_cov <- spatnames
  
  if (specieseffects) {
    
  cov_list <- paste(as.vector(outer(paste0(as.character(unique(species_dataset[[name]])),'_'),incl_cov, 'paste0')), collapse = ' + ')

  components_joint <- update(components_joint, paste('~ . +', cov_list))
  
  }  
  } 
  
  }
    
  }  
    
  components_joint <- update(components_joint, paste0('~ . +', species,'_spde(main = coordinates, model = spdemodel, group = ', species,', ngroup = ', max(numeric_species),', control.group = ', list(speciesmodel), ')'))
    
  }
  
  if (!is.null(pointcovariates_incl)) {
    
  components_joint <- update(components_joint, paste(' ~ . +', paste(pointcovariates_incl, collapse = ' + ')))    
   
  }  
  
  if (attributes(data)$Marks) {
    
  if (marksspatial) {
    
  if (sharedspatial) {
    
  if (!pointsspatial) {
    
  components_joint <- update(components_joint, paste(' ~ . +', 'shared_spatial(main = coordinates, model = spdemodel)'))  
    
  }  
    
  }  
  else {  
  if (!is.null(spatialdatasets)) {
    
  marks_spatial_names <- names(data_marks)[gsub('_.*$',"",names(data_marks))%in%spatialdatasets]  
  marks_spatial_datasets <- gsub('_.*$', '', marks_spatial_names)
    
  }
  else  {
    
  marks_spatial_names <- names(data_marks)
  marks_spatial_datasets <- attributes(data)$Mark_dataset
  
  }
     
  if (pointsspatial) {
  
  components_joint <- update(components_joint, paste(' ~ . +',paste0(marks_spatial_names,'_spde(main = coordinates, copy = ', paste0("\"", marks_spatial_datasets,'_spde',"\""),',  hyper = list(beta = list(fixed = FALSE)))'),collapse = ' + '))
        
  }
  else {
        
  components_joint <- update(components_joint, paste(' ~ . +',paste0(marks_spatial_names,'_spde(main = coordinates, model = spdemodel)',collapse = ' + ')))
        
  }
      
  }
    
  }  
    
  if (marksintercept) {
      
  for (i in 1:length(data_marks)) {
        
  if (attributes(data_marks[[i]])$data_type != "Multinomial mark") {
          
  components_joint <- update(components_joint, paste0(' ~ . +', paste0(c(names(data_marks)[[i]]),'_intercept(1)'), collapse = ' + '))
          
  }
        
  }
      
  }
    
  if (any(multinom_incl)) {
      
  factor_vars <- sapply(data_marks, function(name) attributes(name)$mark_name)
  factor_vars <- unique(factor_vars[multinom_incl])
  components_joint <- update(components_joint, paste('  ~ . + ', paste0(factor_vars,'(main = ', factor_vars, ', model = "iid",constr = FALSE, fixed=TRUE)', collapse = ' + ')))
      
  phi_vars <- sapply(data_marks, function(name) attributes(name)$phi)
  phi_vars <- unique(phi_vars[multinom_incl])
  components_joint <- update(components_joint, paste('  ~ . +', paste0(phi_vars, '(main = ',phi_vars, ', model = "iid", initial = -10, fixed = TRUE)', collapse = ' + ')))
      
  }
    
  }
  
  for (i in 1:(length(likelihoods))) {
    
  if (likelihoods[[i]]$response == points_response[[2]]) options[['control.family']][[i]] <- list(link = 'cloglog')
    
  else options[['control.family']][[i]] <- list(link = 'default')
    
  }
  
  names(likelihoods) <- c(data_names,names_marks)

  model_joint <- inlabru::bru(components = components_joint,
                              likelihoods, options = options)

  data_type <- sapply(c(data_points,data_marks), function(x) attributes(x)$data_type)
  names(data_type) <- c(data_names,attributes(data)$Mark_name)
  model_joint[['data_type']] <- data_type
  model_joint[['dataset_names']] <- data_names
  
  if (any(multinom_incl)) { 
    
  model_joint[['multinom_vars']] <- multinom_vars
    
  }
  
  model_joint[['sources_of_information']] <- unname(c(data_names, sapply(data_marks, function(dat) attributes(dat)$dataset)))
  
  model_joint[['components']] <- components_joint
  
  model_joint[['bru_sdm_options']] <- options
  
  model_joint[['spatial_covariates_used']] <- spatnames
  
  if (!is.null(data_marks)) {
    
  model_joint[['marks_used']] <- sapply(data_marks, function(x) attributes(x)$mark_name)
  names(model_joint[['marks_used']]) <- sapply(data_marks, function(x) attributes(x)$dataset)
    
  }
  else model_joint[['marks_used']] <- NULL
  
  if (pointsspatial | marksspatial) {
    
  if (is.null(spatialdatasets)) {
    
  model_joint[['spatial_datasets']] <- data_names   
    
  }
  else model_joint[['spatial_datasets']] <- spatialdatasets
    
  }
  else model_joint[['spatial_datasets']] <- NULL
  
  if (specieseffects) {
    
  model_joint[['species_in']] <- species_dataset  
  attr(model_joint, 'Species')  <- species
  attr(model_joint, 'Speciesmodel') <- speciesmodel
  
  }
  else model_joint[['species_in']] <- NULL
  
  class(model_joint) <- c('bru_sdm',class(model_joint))
  
  return(model_joint)
  
}
