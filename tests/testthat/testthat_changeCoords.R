testthat::test_that('changeCoords changes the coordinate variables of a list of data.', {
  
  ##Generate arbitrary data:
  #Arbitrary projection
  projection <- CRS('+proj=tmerc')
  
  #Make random shape to generate points on
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)
  
  Poly = Polygon(xy)
  Poly = Polygons(list(Poly),1)
  SpatialPoly = SpatialPolygons(list(Poly), proj4string = projection)
  
  ##Old coordinate names
  oldcoords <- c('x','y')
  #Make random points
  #Random presence only dataset
  PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  colnames(PO@coords) <- oldcoords
  
  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
  colnames(PA@coords) <- oldcoords
  
  ##Make PA a data.frame object
  PA <- data.frame(PA)
  
  spData <- list(PO,PA)  
  
  #Set new coordinate names
  newcoords <- c('long', 'lat')
  
  mod <- changeCoords(data = spData, oldcoords = oldcoords, newcoords = newcoords)
  
  expect_equal(colnames(mod[[1]]@coords), newcoords)
  expect_true(all(newcoords%in%names(mod[[2]])))
  expect_false(any(oldcoords%in%names(mod[[2]])))
  
})