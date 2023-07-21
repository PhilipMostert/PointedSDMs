library(raster)
library(sf)
library(sp)
library(USAboundaries)
library(elevatr)
library(FedData)
library(dplyr)
library(censusapi)
library(spocc)

# coordinate reference system to use throughout
#proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj <- '+proj=longlat +datum=WGS84 +no_defs'
# Get outline of PA
PA <- USAboundaries::us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
st_crs(PA) <- proj

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
    mutate(NPres = if_else(total == 0, 0, 1)) %>%
    mutate(Species_name = 'caerulescens') %>%
    dplyr::rename(X = Longitude, Y = Latitude)
  
  write.csv(BBA_Wren, file = "Data/BBA.csv")
} else {
  BBA_Wren <- read.csv(file = "Data/BBA.csv")
}

BBA <- st_as_sf(
  x = BBA_Wren[, c("NPres", "Species_name", 'X', 'Y')],
  coords = c("X", "Y"),
  crs = proj)

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
      #Ntrials = sum(Ntrials),
      Counts = sum(NPres)) %>%
    mutate(Species_name = 'caerulescens')
    
  routes_in <- intersect(BBS_Wren$Route, routes$Route)
  routes <- routes[routes$Route %in% routes_in,]
  
  BBS_Wren <- BBS_Wren %>% left_join(routes)
  BBS_Wren$Route <- NULL
  
  BBS <- st_as_sf(x = BBS_Wren,
                  coords = c('Longitude', 'Latitude'),
                  crs = proj)
  
  write.csv(BBS_Wren, file = "Data/BBS.csv")

  } else {
  
    BBS_Wren <- read.csv(file = "Data/BBS.csv")

  }

elev_raster <- terra::rast(elevatr::get_elev_raster(PA, z = 6, clip = "locations"))
elev_raster <- terra::project(elev_raster, proj)
#PA@proj4string <- sp::CRS(proj)
NLCD_canopy <- get_nlcd(
  template = PA,
  year = 2011,
  dataset = "landcover",
  label = "PA_lc")

NLCD_canopy <- raster::projectRaster(from = NLCD_canopy, to = elev_raster)
NLCD_canopy_raster <- raster::mask(NLCD_canopy, PA)

NLCD_canopy_raster <- terra::rast(NLCD_canopy_raster)
NLCD_canopy_raster <- terra::resample(NLCD_canopy_raster, elev_raster)

terra::writeRaster(
  c(elev_raster, NLCD_canopy_raster),
  filename = 'inst/extdata/SetophagaCovariates.tif',
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW")
)

SetophagaData <- list(BBS = BBS, BBA = BBA)

usethis::use_data(SetophagaData, overwrite = TRUE, compress = 'xz')
