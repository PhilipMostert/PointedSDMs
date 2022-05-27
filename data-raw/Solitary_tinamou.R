
library(usethis)
devtools::install_github('oharar/PointedSDMs')
library(PointedSDMs)


Projection <- CRS("+proj=longlat +ellps=WGS84")
data('SolTin_covariates')
#Using script from the PointedSDMs vignette to make data

Forest <- SpatialPixelsDataFrame(points = PointedSDMs::SolTin_covariates[,1:2],
                                 data = data.frame(Forest = PointedSDMs::SolTin_covariates[,'Forest']), tol = 0.1,
                                 proj4string = Projection)

Forest <- raster(Forest)
#Forest <- disaggregate(Forest, 5, method='bilinear')

NPP <- SpatialPixelsDataFrame(points = PointedSDMs::SolTin_covariates[,1:2],
                              data = data.frame(NPP = PointedSDMs::SolTin_covariates[,'NPP']), tol = 0.1,
                              proj4string = Projection)

NPP <- raster(NPP)
#NPP <- disaggregate(NPP, 5, method='bilinear')

Altitude <- SpatialPixelsDataFrame(points = PointedSDMs::SolTin_covariates[,1:2],
                                   data = data.frame(Altitude = PointedSDMs::SolTin_covariates[,'Altitude']), tol = 0.1,
                                   proj4string = Projection)

Altitude <- raster(Altitude)
#Altitude <- disaggregate(Altitude, 5, method='bilinear')

data('SolTin_ebird')
data('SolTin_gbif')
data('SolTin_parks')
SolTin_parks$Present <- as.numeric(SolTin_parks$Present)


region.mask <- as.owin(cbind(SolTin_covariates[,c("X","Y")], In=rep(TRUE,nrow(SolTin_covariates))),
                    step=c(0.25, 0.3))
region.mask$m[is.na(region.mask$m)] <- FALSE
Region.poly <- simplify.owin(as.polygonal(region.mask), dmin=0.5)
PolyPoints <- cbind(Region.poly$bdry[[1]]$x[c(1,length(Region.poly$bdry[[1]]$x):1)],
                    Region.poly$bdry[[1]]$y[c(1,length(Region.poly$bdry[[1]]$y):1)])
Pgon <- Polygons(list(region=Polygon(coords=PolyPoints)), ID="region")
Region <- SpatialPolygons(list(Pgon), proj4string = Projection)


Meshpars <- list(cutoff=0.8, max.edge=c(1, 3), offset=c(1,1))
Mesh <- MakeSpatialRegion(data=NULL, bdry=Region, meshpars=Meshpars,
                          proj = Projection)


SolitaryTinamou <- list(datasets = list(eBird = SolTin_ebird, Parks = SolTin_parks, Gbif = SolTin_gbif),
                        covariates = list(Forest = Forest, NPP = NPP, Altitude = Altitude),
                        region = Region,
                        mesh = Mesh$mesh)

usethis::use_data(SolitaryTinamou, overwrite = TRUE, compress = 'xz')

