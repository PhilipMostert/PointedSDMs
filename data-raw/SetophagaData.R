library(raster)
library(sf)
library(USAboundaries)
library(elevatr)
library(FedData)
library(dplyr)
library(censusapi)
library(spocc)

# coordinate reference system to use throughout
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Get outline of PA
PA <- USAboundaries::us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, "Spatial")

if (!file.exists("Data/BBA.csv")) {
  
  # Local Location of file with data from Miller et al. (2019)
  # Data downloaded from here:
  Miller.url <- "https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13110&file=mee313110-sup-0001-supplementA.zip"
  Miller.file <- "Data/mee313110-sup-0001-supplementa.zip"
  
  # download the Miller file if needed
  if (!file.exists(Miller.file)) {
    download.file(Miller.url, Miller.file)
  }
  
  load(unzip(Miller.file, files = "DI_Data.Rdata"))
  
  # sum counts (counts undertaken in 5 time intervals at a single spatial point)
  BBA_Wren <- bba %>%
    mutate(total = dplyr::select(., v1:v5) %>% rowSums(na.rm = TRUE)) %>%
    dplyr::select(-c(v1:v5)) %>%
    mutate(present = if_else(total == 0, FALSE, TRUE)) %>%
    dplyr::rename(X = Longitude, Y = Latitude)
  
  write.csv(BBA_Wren, file = "Data/BBA.csv")
} else {
  BBA_Wren <- read.csv(file = "Data/BBA.csv")
}

BBA_sp <- SpatialPointsDataFrame(
  coords = BBA_Wren[, c("X", "Y")],
  data = BBA_Wren[, c("present", "point")],
  proj4string = crs(proj))

if (!file.exists("Data/BBS.csv")) {
 
  routes.url <- 'https://www.sciencebase.gov/catalog/file/get/625f151ed34e85fa62b7f926?f=__disk__95%2F4e%2F1a%2F954e1a5ec1112c34bba91ad1fce6bd2d3c9df90e'
  routes.file <- 'Data/routes.zip'
  
 routes <- download.file(routes.url, routes.file)
   
 routes <- read.csv(unzip(routes.file)) %>% 
    filter(StateNum == 72) %>% select(Longitude, Latitude, Route)
  
 BBS.url <- 'https://www.sciencebase.gov/catalog/file/get/625f151ed34e85fa62b7f926?f=__disk__b5%2Fde%2F6b%2Fb5de6b3d221373125e90a85b5fcba4b7059c7f94'
 BBS.file <- 'Data/50_stops.zip'
 
 stops <- download.file(BBS.url, BBS.file) 
 
 stops <- unzip(BBS.file)
 
 BBS_Wren <- read.csv(unzip(stops[grep('Fifty8.zip', stops)]))  %>% 
   filter(AOU == 06540, Year >= 2005, Year <= 2009, StateNum == 72)
  
  BBS_Wren <- BBS_Wren %>% mutate(NPres = rowSums(dplyr::select(., starts_with("stop")) > 0)) %>%
    mutate(Ntrials = rowSums(!is.na(dplyr::select(., starts_with("stop")))))
  
  BBS_Wren <- BBS_Wren %>% group_by(Route) %>%
    summarise(
      Ntrials = sum(Ntrials),
      NPres = sum(NPres)) 
  
  routes_in <- intersect(BBS_Wren$Route, routes$Route)
  routes <- routes[routes$Route %in% routes_in,]
  
  BBS_Wren <- BBS_Wren %>% left_join(routes)
  
  BBS_Wren <- sp::SpatialPointsDataFrame(coords = data.frame(BBS_Wren[, c('Longitude', 'Latitude')]),
                                         data = data.frame(Ntrials = BBS_Wren$Ntrials, NPres = BBS_Wren$NPres),
                                         proj4string = proj)
  
  BBS_Wren$Ntrials <- NULL; names(BBS_Wren) <- 'Counts'
  
  write.csv(BBS_Wren, file = "Data/BBS.csv")

  } else {
  
    BBS_Wren <- read.csv(file = "Data/BBS.csv")

  }

elev_raster <- elevatr::get_elev_raster(PA, z = 6, clip = "locations")

PA@proj4string <- proj
NLCD_canopy <- get_nlcd(
  template = PA,
  year = 2011,
  dataset = "canopy",
  label = "PA_lc")

NLCD_canopy <- projectRaster(from = NLCD_canopy, to = elev_raster)
NLCD_canopy_raster <- mask(NLCD_canopy, PA)

SetophagaData <- list(BBS = BBS, BBA = BBA, 
                      elev_raster = elev_raster, 
                      NLCD_canopy_raster = NLCD_canopy_raster)

usethis::use_data(SetophagaData, overwrite = TRUE, compress = 'xz')
