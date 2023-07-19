test_that('checkVar can identify a variable name in the data.', {
  
  ##Set up a model
  ##Set up arbitrary data
  projection <- '+proj=tmerc'
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  xy <- cbind(x, y)
  SpatialPoly <- st_sfc(st_polygon(list(xy)), crs = projection)
  
  
  ##Old coordinate names

  #Make random points
  #Random presence only dataset
  PO <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PO) <- 'geometry'
  #Generate random variable
  PO$var <- runif(n = nrow(PO))
  
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PA) <- 'geometry'
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  #Generate random variable
  PA$var <- runif(n = nrow(PA))
  
  ##Make PA a data.frame object
  PA$long <- st_coordinates(PA)[,1]
  PA$lat <- st_coordinates(PA)[,2]
  st_geometry(PA) <- NULL
  PA <- data.frame(PA)
  
  spData <- list(PO,PA)  
  
  #Set new coordinate names
  
  mod <- checkVar(data = spData, var = 'var')
  expect_true(mod)
  
  ##Change var argument
  
  mod2 <- checkVar(data = spData, var = 'notvar')
  expect_false(mod2)
  
  ##Change variable name of PO
  names(PO) <- 'notvar'
  spData <- list(PO,PA)
  
  mod3 <- checkVar(data = spData, var = 'var')
  expect_false(mod3)
  
  ##Change coords of PO back to coordNames; change coord names of PA
  names(PO) <- 'var'
  names(PA)[names(PA) == 'var'] <- 'notvar'
  
  spData <- list(PO,PA)
  mod4 <- checkVar(data = spData, var = 'var')
  expect_false(mod4)
  
})