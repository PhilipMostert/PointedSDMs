#' Function to run marked species distribution models with data from marked presence only and presence absence data.
#' 
#' @param data A bru_sdm_data object created with /code{organize_data}.
#' @param spatialcovariates Data frame of the spatial covariates accompanied by their associated coordinates. Defaults to \code{NULL}.
#' @param covariatestoinclude A vector of spatial covariate names to include in the model. Defaults to \code{NULL}.
#' @param pointsintercept Include individual intercepts for each point process in the model. Defaults to \code{TRUE}.
#' @param marksintercept Include individual intercepts for each mark process in the model. Defaults to \code{TRUE}.
#' @param spdemodel inla.spde model used in the model. May be a named list where the name of the spde object is the name of the associated dataset. Default \code{NULL} uses \code{inla.spde2.matern}.
#' @param pointsspatial Should spatial effects be used for the points in the model. Defaults to \code{TRUE}.
#' @param marksspatial Should spatial effects be used for the marks in the model. Defaults to \code{TRUE}.
#' @param timemodel Time-series model to use. Defaults to \code{list(model = 'ar1')}.
#' @param options INLA or inlabru options to be used in the model.

bru_sdm <- function(data, spatialcovariates = NULL, covariatestoinclude = NULL,
                    pointsintercept = TRUE, marksintercept = TRUE,
                    spdemodel = NULL, pointsspatial = TRUE, marksspatial = TRUE,
                    timemodel = list(model = 'ar1'), options = list()) {
  
  if (class(data)[1] != 'bru_sdm_data') stop('Please supply data formed by the "organize_data" function.')
  
  proj <- data@ips@proj4string
  coords <- colnames(data@ips@coords)
  
  data_points <- append(data@PO_data, data@PA_data)
  data_names <- names(data_points)
  points_family <- sapply(data_points, function(data) attributes(data)$family)
  points_response <- attributes(data)$Points_response
  
  if (attributes(data)$Marks) {
    
  data_marks <- data@Mark_data
  names_marks <- names(data_marks)
  family_marks <- attributes(data)$Mark_family
  mark_weights <- attributes(data)$Mark_weight
  response_marks <- attributes(data)$Mark_response
  multinom_incl <- attributes(data)$Multinom_incl
  multinom_vars <- attributes(data)$Multinom_vars
    
  }
  else {
    
  names_marks <- NULL
  data_marks <- NULL
  multinom_incl <- NULL
  multinom_vars <- NULL
  response_marks <- NULL
    
  }
  
  time_var <- attributes(data)$Timevariable
  
  if (!is.null(time_var)) {
    
  ##leave a var which translates factor var to numeric  
  numeric_time <- as.numeric(unlist(sapply(data_points, function(data) {
      
  data@data[,time_var]  
      
  })))
    
    
  for (k in 1:length(data_points)) {
      
  if (k == 1) { 
        
  length_var <- (1:length(data_points[[1]]))
        
  }
  else {
        
  length_var <- (length(data_points[[k-1]]) + 1):(length(data_points[[k-1]]) + length(data_points[[k]]))
        
  }
      
  data_points[[k]]@data[,time_var] <- numeric_time[length_var]   
      
  }
    
  ips_coords <- do.call(rbind, replicate(max(numeric_time), data@ips@coords, simplify = FALSE))
    
  ips_weight <- do.call(rbind, replicate(max(numeric_time), data@ips@data, simplify = FALSE))
    
  ips_fact <- rep(1:max(numeric_time), each = nrow(data@ips@coords))
    
  ips_data <- data.frame(ips_weight, ips_fact)
    
  data@ips <- sp::SpatialPointsDataFrame(coords = ips_coords,
                                         data = ips_data,
                                         proj = proj,
                                         match.ID = FALSE)
    
  names(data@ips@data) <- c(names(ips_weight), paste0(time_var))
    
    
  }
  
  if (!is.null(spatialcovariates)) {
    
  if (class(spatialcovariates) == 'RasterLayer') {
      
  spatialcovariates <- as(spatialcovariates, 'SpatialPixelsDataFrame')
      
  }
    
  if (class(spatialcovariates) == 'data.frame') {
      
  spatialcovariates <- sp::SpatialPointsDataFrame(coords = spatialcovariates[,coords],
                                                  data = spatialcovariates[,!names(spatialcovariates)%in%coords],
                                                  proj4string = proj)
      
  spatialcovariates <- as(spatialcovariates, 'SpatialPixelsDataFrame')
      
  }
    
  spatnames <- names(spatialcovariates@data)
  spatdata_class <- sapply(spatialcovariates@data, class)
    
  if (!is.null(covariatestoinclude)) {
      
  spatdata_class <- spatdata_class[spatnames%in%covariatestoinclude] 
  spatnames <- spatnames[(spatnames%in%covariatestoinclude)]
      
  if (is.null(spatnames) | identical(spatnames,character(0))) stop('covariatestoinclude contains covariate names not found in spatialcovariates')
      
  }
    
  for (name in spatnames) {
      
  pixels_df <- sp::SpatialPixelsDataFrame(points = spatialcovariates@coords,
                                          data = data.frame(spatialcovariates@data[,name]),
                                          proj4string = proj)
  names(pixels_df) <- name
  assign(name,pixels_df)
      
  }
    
  }
  else spatnames <- NULL
  
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
  
  form_elements <- gsub(" *\\(.*?\\) *", "",components_joint)
  
  formula <- mapply(function(fam,index) {
    
  if (fam == 'cp') {
      
  formula <- formula(paste0(c('coordinates','~', form_elements[2]),collapse = " ")) 
      
  }
  else
  if (fam == 'poisson') {
        
  formula <- formula(paste0(c(points_response[1],'~', form_elements[2]),collapse = " ")) 
        
  }
  else
  if (fam == 'binomial') {
        
  formula <- formula(paste0(c(points_response[2],"~", form_elements[2]),collapse = " "))
        
  }
    
  if (pointsintercept) {
      
  formula <- update(formula,paste0(' ~ . +', paste0(data_names[index],'_intercept'), collapse = ' + '))
      
  }
  else formula
    
  if (pointsspatial) {
      
  formula <- update(formula, paste0('~ . +',data_names[[index]],'_spde'))
      
  }
  else formula
    
  if (!is.null(time_var)) {
      
  formula <- update(formula, paste0('~ . +',time_var,'_spde'))
      
      
  }
  else formula
    
  }, fam = points_family, index = 1:length(points_family))
  
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
    
  formula_marks <- list()
  likelihoods_marks <- list()
    
  for (i in 1:length(response_marks)) {
      
  formula_marks[[i]] <- formula(paste0(c(response_marks[i],'~',form_elements[2]),collapse = " "))
      
  if (marksspatial) {
        
  formula_marks[[i]] <- update(formula_marks[[i]], paste0(" . ~ . +",names(data_marks)[[i]],'_spde'))
        
  }
      
  if (marksintercept) {
        
  if (attributes(data_marks[[i]])$data_type != 'Multinomial mark') {
          
  formula_marks[[i]] <- update(formula_marks[[i]],paste0('. ~ . +', paste0(names(data_marks)[[i]],'_intercept'), collapse = ' + '))
          
  }
        
  }
      
  if (attributes(data_marks[[i]])$data_type == 'Multinomial mark') {
        
  formula_marks[[i]] <- update(formula_marks[[i]], paste0(paste0(names_marks[i],'_response'), ' ~ . + ', paste(attributes(data_marks[[i]])$mark_name, attributes(data_marks[[i]])$phi, sep = ' + ')))
        
  }
      
  }
    
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
                          mesh = mesh,
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
    
  if (is.list(spdemodel[[1]])) { 
      
  for (name in names(data_points)) {
        
  components_joint <- update(components_joint, paste(' ~ . +',paste0(name,'_spde(main = coordinates, model =',name,'_spde)')))
        
  }    
      
  } 
  else {  
      
  components_joint <- update(components_joint, paste(' ~ . +',paste0(data_names,'_spde(main = coordinates, model = spdemodel)',collapse = ' + ')))
      
  }
    
  }  
  
  if (pointsintercept) {
    
  components_joint <- update(components_joint, paste0(' ~ . +', paste0(data_names,'_intercept(1)'), collapse = ' + '))
    
  }
  
  if (!is.null(time_var)) {
    
  components_joint <- update(components_joint, paste0('~ . +', time_var,'_spde(main = coordinates, model = spdemodel, group = ', time_var,', ngroup = ', max(numeric_time),', control.group = ', timemodel,')'))
  #if bottom one remove ips things
  #components_joint <- update(components_joint, paste0('~ . +', time_var,'_spde(main = ', time_var,', model = ','\"', timemodel,'\"',')'))
    
  }
  
  if (attributes(data)$Marks) {
    
  if (marksspatial) {
      
  if (pointsspatial) {
        
  components_joint <- update(components_joint, paste(' ~ . +',paste0(names(data_marks),'_spde(main = coordinates, copy = ', paste0("\"",attributes(data)$Mark_dataset,'_spde',"\""),',  hyper = list(beta = list(fixed = FALSE)))'),collapse = ' + '))
        
  }
  else {
        
  components_joint <- update(components_joint, paste(' ~ . +',paste0(names(data_marks),'_spde(main = coordinates, model = spdemodel)',collapse = ' + ')))
        
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
    
  model_joint[['marks_used']] <- response_marks
  names(model_joint[['marks_used']]) <- sapply(data_marks, function(x) attributes(x)$dataset)
    
  }
  else model_joint[['marks_used']] <- NULL
  
  class(model_joint) <- c('bru_sdm',class(model_joint))
  
  return(model_joint)
  
}
