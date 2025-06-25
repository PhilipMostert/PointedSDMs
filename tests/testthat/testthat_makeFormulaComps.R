test_that('makeFormulaComps makes components based on the formula given.', {
  
  form <- ~ cov + I(cov^2) + sqrt(cov)
  speciesnames <- c('bird', 'fish', 'dog')
  
  biasNospecies <- PointedSDMs:::makeFormulaComps(form = form,species = FALSE, speciesnames = NULL, type = 'Bias')
  expect_true(biasNospecies == "Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")")
  
  covNospecies <- PointedSDMs:::makeFormulaComps(form = form,species = FALSE, speciesnames = NULL, type = 'Cov')
  expect_true(covNospecies == "Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")")
  
  biasSpecies <- PointedSDMs:::makeFormulaComps(form = form,species = 'community', speciesnames = speciesnames, type = 'Bias')
  expect_setequal(biasSpecies, c("bird_Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                 "fish_Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                 "dog_Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")"))
  
  covSpecies <- PointedSDMs:::makeFormulaComps(form = form,species = 'community', speciesnames = speciesnames, type = 'Cov')
  expect_setequal(covSpecies, c("bird_Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                 "fish_Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                 "dog_Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")"))
  
  biasSpecies <- PointedSDMs:::makeFormulaComps(form = form,species = 'stack', speciesnames = speciesnames, type = 'Bias')
  expect_setequal(biasSpecies, c("bird_Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                 "fish_Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                 "dog_Bias__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")"))
  
  covSpecies <- PointedSDMs:::makeFormulaComps(form = form,species = 'stack', speciesnames = speciesnames, type = 'Cov')
  expect_setequal(covSpecies, c("bird_Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                "fish_Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")",
                                "dog_Fixed__Effects__Comps(main = ~cov + I(cov^2) + sqrt(cov) - 1, model = \"fixed\")"))
}
)