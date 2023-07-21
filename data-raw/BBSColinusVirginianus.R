library(sp)
library(sf)
library(dplyr)

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
AL <- USAboundaries::us_states(states = "Alabama")
AL <- as(AL, "sf")

if (!file.exists('data/BBS_Colinus_virginianus')) {
  
  routes.url <- 'https://www.sciencebase.gov/catalog/file/get/625f151ed34e85fa62b7f926?f=__disk__95%2F4e%2F1a%2F954e1a5ec1112c34bba91ad1fce6bd2d3c9df90e'
  routes.file <- 'Data/routes.zip'
  
  routes <- download.file(routes.url, routes.file)
  
  routes <- read.csv(unzip(routes.file)) %>% 
    filter(Active == 1, StateNum == 2)
  
  
  BBS.url <- 'https://www.sciencebase.gov/catalog/file/get/625f151ed34e85fa62b7f926?f=__disk__b5%2Fde%2F6b%2Fb5de6b3d221373125e90a85b5fcba4b7059c7f94'
  BBS.file <- 'Data/50_stops.zip'
  
  stops <- download.file(BBS.url, BBS.file) 
  
  stops <- unzip(BBS.file)
  
  BBS_Colinus_virginianus <- read.csv(unzip(stops[grep('Fifty1.zip', stops)]))  %>% 
    filter(AOU == 2890, Year >= 2015, Year <= 2017, StateNum == 2) %>%
    mutate(NPres = rowSums(dplyr::select(., starts_with("stop")) > 0)) %>%
    mutate(Ntrials = rowSums(!is.na(dplyr::select(., starts_with("stop")))))  %>% 
    group_by(Route) %>%
    summarise(
      Year = Year,
      Ntrials = sum(Ntrials),
      NPres = sum(NPres))
  
  routes_in <- intersect(BBS_Colinus_virginianus$Route, routes$Route)
  routes <- routes[routes$Route %in% routes_in,]
  
  BBS_Colinus_virginianus <- BBS_Colinus_virginianus %>% left_join(routes) %>%
    filter(!is.na(Latitude), !is.na(Longitude))
  
  BBS_Colinus_virginianus <- st_as_sf(x = BBS_Colinus_virginianus[, c('Longitude', 'Latitude', 
                                                                      'Year', 'Ntrials', 'NPres')],
                                      coords = c('Longitude', 'Latitude'),
                                      crs = proj)
  
  BBSColinusVirginianus <- BBS_Colinus_virginianus[unlist(st_intersects(AL, BBS_Colinus_virginianus)),]
  
  write.csv(BBS_Colinus_virginianus, file = "Data/BBSColinusVirginianus.csv")
  
} else {
  
  BBS_Colinus_virginianus <- read.csv('Data/BBSColinusVirginianus.csv')
  
}

usethis::use_data(BBSColinusVirginianus, overwrite = TRUE, compress = 'xz')
