testthat::test_that('nameChanger changes the variable name of a list of data.', {
  
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
  #Make random points
  #Random presence only dataset
  PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  ##Add random variable
  PO$oldvar <- runif(n = nrow(PO@coords))
  
  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
  PA$oldvar <- runif(n = nrow(PA@coords))

  ##Make PA a data.frame object
  PA <- data.frame(PA)
  
  spData <- list(PO,PA)  
  
  #Set new coordinate names
  newvar <- 'newvar'
  
  mod <- nameChanger(data = spData, oldName = 'oldvar', newName =  'newvar')
  
 expect_true('newvar'%in%names(mod[[1]]@data))
 expect_false('oldvar'%in%names(mod[[1]]@data))
 expect_true('newvar'%in%names(mod[[2]]))
 expect_false('olv_var'%in%names(mod[[2]]))
  
})