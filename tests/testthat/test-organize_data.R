

testthat::test_that('organize_data gives the correct output given arguments', {
  
 #Arbitrary projection
projection <- CRS('+proj=tmerc')

#Make random shape to generate points on
x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
xy <- cbind(x, y)

Poly = Polygon(xy)
Poly = Polygons(list(Poly),1)
SpatialPoly = SpatialPolygons(list(Poly), proj4string = projection)

#Make random points
 #Random presence only dataset
PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
 #Add species name
PO$species <- sample(x = c('a','b'), size = nrow(PO@coords), replace = TRUE)
 #Add mark
PO$mark1 <- rgamma(n = nrow(PO@coords),2,3)

 #Random presence absence dataset
PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
 #Add species name
PA$species <- sample(x = c('c','d'), size = nrow(PA@coords), replace = TRUE)
 #Add mark
PA$mark2 <- rnorm(n = nrow(PA@coords))
 #Add response name
PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
 #Add trial name
PA$trial <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)
 #Add multionomial mark
PA$multi <- sample(c('x','y'), size = nrow(PA@coords), replace = TRUE)

obj <- organize_data(PO, PA, poresp = 'POresp', paresp = 'PAresp',
                     trialname = 'trial', coords = colnames(PO@coords), proj = projection,
                     marks = TRUE, markfamily = c(mark1 = 'gamma', mark2 = 'gaussian'),
                     speciesname = 'species',
                     meshpars = list(cutoff=0.8, max.edge=c(1, 3), offset=c(1,1)),
                     boundary = SpatialPoly)

expect_s4_class(obj,'bru_sdm_data')
expect_equal(class(obj@ips)[1], 'SpatialPointsDataFrame')
expect_equal(class(obj@mesh), 'inla.mesh')
expect_equal(class(obj@PO_data), 'list')
expect_equal(class(obj@PA_data), 'list')
expect_equal(class(obj@Mark_data), 'list')

expect_equal(sapply(obj@PO_data, class), c(PO = 'SpatialPointsDataFrame'))
expect_equal(sapply(obj@PA_data, class), c(PA = 'SpatialPointsDataFrame'))
expect_equal(sapply(obj@Mark_data, class), c(PO_mark1 = 'SpatialPointsDataFrame',
                                             PA_mark2 = 'SpatialPointsDataFrame',
                                             PA_multi = 'SpatialPointsDataFrame'))

expect_equal(names(obj@PO_data), 'PO')
expect_equal(names(obj@PA_data), 'PA')
expect_equal(names(obj@Mark_data), c('PO_mark1','PA_mark2', 'PA_multi'))

expect_equal(attributes(obj)$Mark_data_type, c(PO_mark1 = 'Gamma mark',
                                               PA_mark2 = 'Gaussian mark',
                                               PA_multi = 'Multinomial mark'))
expect_equal(attributes(obj)$Multinom_incl, c(PO_mark1 = FALSE,
                                              PA_mark2 = FALSE,
                                              PA_multi = TRUE))

expect_equal(attributes(obj)$Points_family, c(PO = 'cp', PA = 'binomial'))

expect_null(attributes(obj)$Points_trials$PO)
expect_identical(attributes(obj)$Points_trials$PA, PA$trial)

})
