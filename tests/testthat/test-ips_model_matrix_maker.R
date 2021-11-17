testthat::test_that('Test ips_model_matrix_maker makes a matrix for the integration points', {
  
  ##Make integration points
  library(inlabru)
  data("gorillas")
  ips <- ipoints(samplers = mesh)
  covariates <- gorillas$gcov$waterdist
  
  #species variable name
  species_var_name <- 'species'
  
  #make three arbitrary species
  allspecies <- c('Species_a','Species_b','Species_c')
  
  ips_matrix <- ips_model_matrix_maker(ips = ips, covariates = covariates, componentstokeep = c(species_var_name,'weight'),
                                       species = species_var_name, allspecies = allspecies,
                                       coords = colnames(ips@coords), proj = ips@proj4string)
  
  ##ie weight is retained in the data of the ips matrix
  expect_equal(names(ips_matrix@data)[names(ips_matrix@data) == 'weight'], 'weight')
  
  #ie species variable name is retained in the ips matrix
  expect_equal(names(ips_matrix@data)[names(ips_matrix@data) == species_var_name], species_var_name)
  
  expect_setequal(names(ips_matrix@data), c("weight", "species", "Species_a_waterdist", "Species_b_waterdist", "Species_c_waterdist", "Species_a_intercept",
                                            "Species_b_intercept", "Species_c_intercept"))
  
  })