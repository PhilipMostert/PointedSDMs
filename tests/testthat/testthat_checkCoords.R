test_that('checkCoords can identify the coordinates in the data.', {
  
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
  coordNames <- c('x','y')
  #Make random points
  #Random presence only dataset
  PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  colnames(PO@coords) <- coordNames
  
  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
  colnames(PA@coords) <- coordNames
  
  ##Make PA a data.frame object
  PA <- data.frame(PA)
  
  spData <- list(PO,PA)  
  
  #Set new coordinate names
  
  mod <- checkCoords(data = spData, coords = coordNames)
  expect_true(mod)
  
  ##Change coords argument
  
  mod2 <- checkCoords(data = spData, coords = c('long', 'lat'))
  expect_false(mod2)
  
  ##Change coords of PO
  colnames(PO@coords) <- c('long','lat')
  spData <- list(PO,PA)
  
  mod3 <- checkCoords(data = spData, coords = coordNames)
  expect_false(mod3)
  
  ##Change coords of PO back to coordNames; change coord names of PA
  colnames(PO@coords) <- coordNames
  names(PA)[names(PA) == coordNames] <- c('long', 'lat')
  
  spData <- list(PO,PA)
  mod4 <- checkCoords(data = spData, coords = coordNames)
  expect_false(mod4)
  
})