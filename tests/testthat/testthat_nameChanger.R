testthat::test_that('nameChanger changes the variable name of a list of data.', {
  
  ##Generate arbitrary data:
  #Arbitrary projection
  projection <- '+proj=tmerc'
  
  #Make random shape to generate points on
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  SpatialPoly <- st_sfc(st_polygon(list(xy)), crs = projection)
  ##Old coordinate names
  #Make random points
  #Random presence only dataset
  
  PO <- st_as_sf(st_sample(x = SpatialPoly, 100))
  st_geometry(PO) <- 'geometry'
  ##Add random variable
  PO$oldvar <- runif(n = nrow(PO))
  
  #Random presence absence dataset
  PA <- st_as_sf(st_sample(SpatialPoly, 100, type = 'random'))
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA), replace = TRUE)
  PA$oldvar <- runif(n = nrow(PA))
  PA$long <- st_coordinates(PA)[,1]
  PA$lat <- st_coordinates(PA)[,2]
  st_geometry(PA) <- NULL
  ##Make PA a data.frame object
  PA <- data.frame(PA)
  
  spData <- list(PO,PA)  
  
  #Set new coordinate names
  newvar <- 'newvar'
  
  mod <- nameChanger(data = spData, oldName = 'oldvar', newName =  'newvar')
  
 expect_true('newvar'%in%names(mod[[1]]))
 expect_false('oldvar'%in%names(mod[[1]]))
 expect_true('newvar'%in%names(mod[[2]]))
 expect_false('olv_var'%in%names(mod[[2]]))
  
})