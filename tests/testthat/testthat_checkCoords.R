test_that('checkCoords can identify the coordinates in the data.', {
  
  ##Set up a model
  ##Set up arbitrary data
  projection <- '+proj=tmerc'
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  xy <- cbind(x, y)
  SpatialPoly <- st_sfc(st_polygon(list(xy)), crs = projection)
  
  ##Old coordinate names
  coordNames <- c('x','y')
  #Make random points
  #Random presence only dataset
  PO <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PO) <- 'geometry'
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, 100, crs = projection))
  st_geometry(PA) <- 'geometry'
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  
  
  ##Make PA a data.frame object
  PA$x <- st_coordinates(PA)[,1]
  PA$y <- st_coordinates(PA)[,2]
  st_geometry(PA) <- NULL
  PA <- data.frame(PA)
  
  spData <- list(PO,PA)  
  
  #Set new coordinate names
  
  mod <- checkCoords(data = spData, coords = coordNames)
  expect_true(mod)
  
  ##Change coords argument
  
  mod2 <- checkCoords(data = spData, coords = c('long', 'lat'))
  expect_false(mod2)
  
  ##Change coords of PO
  colnames(PA)[names(PA) %in% c('x', 'y')] <- c('long','lat')
  spData <- list(PO,PA)
  
  mod3 <- checkCoords(data = spData, coords = coordNames)
  expect_false(mod3)
  
  
})