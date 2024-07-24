##############################
## CONNECTIVITY - TELEMETRY
#############################
# Created: 05-14-24
# By: Sophia Horigan
# Contact: shorigan@uchicago.edu
# Last Updated: 05-14-24

# This script takes in MoveBank telemetry data for Eidolon dupreanum and computes home range size
# to use for intermingling in the metapopulation model

rm(list=ls())

library(lubridate)
library(dplyr)
library(sp)
library(sf)
library(adehabitatHR)
library(ggplot2)
library(gganimate)
library(ggmap)
library(geosphere)


#--------------------------------------------------------------------------------------------
## SET WD
homedir <- "/Users/shorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/WIP/"
setwd(homedir)

#--------------------------------------------------------------------------------------------
## LOAD DATA
points.df <- read.csv(paste0(homedir,"/Input/Movebank-AllTags-AllSensorTypes-06182024.csv"), header = TRUE)

#--------------------------------------------------------------------------------------------
## CLEAN DATA

# cleanup
# convert timestamp column 
points.df$timestamp <- ymd_hms(points.df$timestamp)
# change time zone to Madagascar
points.df$timestamp <- with_tz(points.df$timestamp, "Africa/Addis_Ababa")
# add site column
points.df <- points.df %>% 
  mutate(site = case_when(grepl("MAN", individual.local.identifier) ~ 'Ambositra', # pteropus
                          grepl("ANA", individual.local.identifier) ~ 'Analambotaka', # pteropus
                          grepl("TSI", individual.local.identifier) ~ 'Marotsipohy', # pteropus
                          grepl("MARO", individual.local.identifier) ~ 'Marovitsika', # pteropus
                          grepl("HAR", individual.local.identifier) ~ 'Nosy Hara', # eidolon
                          grepl("LOR", individual.local.identifier) ~ 'Ambositra', # eidolon
                          grepl("NAT", individual.local.identifier) ~ 'Mangroves', # pteropus
                          grepl("VHL", individual.local.identifier) ~ 'Vahialava', # pteropus
                          grepl("KEL", individual.local.identifier) ~ 'Angavokely', # eidolon
                          grepl("WAY", individual.local.identifier) ~ 'Ankarana')) # eidolon

table(points.df$individual.local.identifier)
# drop any that have less than 5 data points (right now just 1)
points.df <- points.df[!(points.df$individual.local.identifier %in% "LOR002"),]


# separate Pteropus and Eidolon
pteropus.points <- points.df %>%
  filter(individual.taxon.canonical.name == "Pteropus rufus")
  
eidolon.points <- points.df %>%
  filter(individual.taxon.canonical.name == "Eidolon dupreanum")


#--------------------------------------------------------------------------------------------
## FORMAT FOR SPATIAL ANALYSIS

# To format for spatial analyses, we need to remove NA values from the coordinate columns
eidolon.points <- eidolon.points[!is.na(eidolon.points$location.lat) & !is.na(eidolon.points$location.long),]

# DAY VS NIGHT 
dusk_hour = 17 # made up 
dawn_hour = 6

eidolon.points <- eidolon.points %>%
  mutate(daytime = ifelse(hour(timestamp) < dusk_hour & hour(timestamp) > dawn_hour, "day", "night"))

tmp <- eidolon.points[, c("site", "location.long", "location.lat", "daytime")]

eidolon.only.day <- tmp %>% 
  filter(daytime == "day")

eidolon.only.night <- tmp %>% 
  filter(daytime == "night")


# The dataframe should only have 3 columns (x, y, and an identifier) for home ranges
eidolon.sp <- eidolon.points[, c("individual.local.identifier", "location.long", "location.lat")]
eidolon.site <- eidolon.points[, c("site", "location.long", "location.lat")]
eidolon.only.day.sp <- eidolon.only.day[, c("site", "location.long", "location.lat")]
eidolon.only.night.sp <- eidolon.only.night[, c("site", "location.long", "location.lat")]

# Turn into a spatial points dataframe (class: SpatialPointsDataFrame)
coordinates(eidolon.sp) <- c("location.long", "location.lat")
coordinates(eidolon.site) <- c("location.long", "location.lat")
coordinates(eidolon.only.day.sp) <- c("location.long", "location.lat")
coordinates(eidolon.only.night.sp) <- c("location.long", "location.lat")

# Examine the structure of our SpatialPointsDataFrame
str(eidolon.sp)
str(eidolon.site)
str(eidolon.only.day.sp)
str(eidolon.only.night.sp)

# convert to UTM
locations_sf <- points.df %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)

locations_projected <- locations_sf %>% 
  st_transform(crs = "EPSG:32738") # see https://epsg.io/32604 for more info...


# Set coordinate system & projection
crs_wgs84 <- CRS("+proj=longlat +datum=WGS84")
class(crs_wgs84)
slot(eidolon.sp, "proj4string") <- crs_wgs84
slot(eidolon.site, "proj4string") <- crs_wgs84
slot(eidolon.only.day.sp, "proj4string") <- crs_wgs84
slot(eidolon.only.night.sp, "proj4string") <- crs_wgs84

# convert to UTM
eidolon.only.night.sp <- spTransform(eidolon.only.night.sp, CRS("+proj=utm +zone=38 +south ellps=WGS84"))
eidolon.only.night.sp

#--------------------------------------------------------------------------------------------
## INTERMINGLING
# CALCULATE KERNAL DENSITIY AREA

# INDIVIDUAL
# calculate UD
eidolon.kernel <- kernelUD(eidolon.only.night.sp, h="href")

# get contour
list <- vector(mode = "list", length = length(eidolon.kernel))

for (i in 1:length(eidolon.kernel)){
  list[[i]] <- as.image.SpatialGridDataFrame(eidolon.kernel[[i]])
  image(list[[i]])
  contour(list[[i]], add=TRUE)
}

# calculate area
# grid size issues
UD.area <- kernel.area(eidolon.kernel, percent = seq(10, 90, by=5), unout = "m2")
plot(UD.area)

UD.area

## ----------------------------
# make contours
UD.rad <- sqrt(UD.area/pi)

# --------
# BY SITE # NEED TO UPDATE
# calculate kernel densities
eidolon.kernel.site <- kernelUD(eidolon.site, h = "href")
image(eidolon.kernel.site)

# get volume
eidolon.kernel.site.vud <- getvolumeUD(eidolon.kernel.site)

# get contour
levels.site <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
list.site <- vector(mode = "list", length = length(eidolon.kernel.site.vud))

for (i in 1:length(eidolon.kernel.site.vud)){
  list.site[[i]] <- as.image.SpatialGridDataFrame(eidolon.kernel.site.vud[[i]])
}

# plot
for (i in 1:length(eidolon.kernel.site.vud)){
  plot(eidolon.kernel.site.vud[[i]])
  contour(list.site[[i]], add = TRUE, levels = levels.site)
}

# calculate home range area
# grid size issues
homerange <- kernel.area(eidolon.kernel.site, percent = seq(10, 90, by = 5))
plot(homerange)

homerange
# average the area of each homerange %? or take random draws from each bat as a representative
# random draw - keeps correlations between bats.
# or develop range of area for 50% occupation..?

#--------------------------------------------------------------------------------------------
## DISPSERSAL

# 1. Count the number of putative roosts used across the total number of days
ggplot(data = eidolon.daytime[eidolon.daytime], aes(location.long, location.lat)) + 
  geom_point(aes(colour = daytime)) +
  facet_wrap(individual.local.identifier ~ ., scales = 'free')

# Make sure points are in order (in case they weren't before)
eidolon.only.day <- tmp %>%
  arrange(individual.local.identifier, timestamp) %>% # arrange by animal, and ascending date
  filter(!is.na(location.long), !is.na(location.lat)) # remove NA's

# subset by ID
har002 <- eidolon.only.day %>%
  filter(individual.local.identifier == 'HAR002')
har003 <- eidolon.only.day %>%
  filter(individual.local.identifier == 'HAR003')
kel815 <- eidolon.only.day %>%
  filter(individual.local.identifier == 'KEL815')
kel818 <- eidolon.only.day %>%
  filter(individual.local.identifier == 'KEL818')
way155 <- eidolon.only.day %>%
  filter(individual.local.identifier == 'WAY155')
way179 <- eidolon.only.day %>%
  filter(individual.local.identifier == 'WAY179')

# COUNT : NUMBER OF ROOST SWITCHES
# har002 
plot(har002$timestamp, har002$location.long)
plot(har002$timestamp, har002$location.lat)
for (i in 1:length(har002$location.long)){ # points > 1000m?
  print(i)
  print(i+1)
  print(distm(c(har002$location.long[i], har002$location.lat[i]), c(har002$location.long[i+1], har002$location.lat[i+1]), fun = distGeo))
}
num_switch_har002 = 1
# duration
x = interval(har002$timestamp[1], har002$timestamp[length(har002$timestamp)])
y = int_length(x) # in seconds
d = y/86400 # days
d # 10 days
daily_prob_disp_har002 = num_switch_har002/d
daily_prob_disp_har002 

# har003 - 0 roost switches
plot(har003$timestamp, har003$location.long)
plot(har003$timestamp, har003$location.lat)
num_switch_har003 = 1
# duration
x = interval(har003$timestamp[1], har003$timestamp[length(har003$timestamp)])
y = int_length(x) # in seconds
d = y/86400
d # 1
daily_prob_disp_har003 = num_switch_har003/d
daily_prob_disp_har003 

# kel815 
plot(kel815$timestamp, kel815$location.long)
plot(kel815$timestamp, kel815$location.lat)
for (i in 1:length(kel815$location.long)){ # points > 1000m?
  print(i)
  print(i+1)
  print(distm(c(kel815$location.long[i], kel815$location.lat[i]), c(kel815$location.long[i+1], kel815$location.lat[i+1]), fun = distGeo))
}
num_roost_switch_kel815 = 1
# duration
x = interval(kel815$timestamp[1], kel815$timestamp[length(kel815$timestamp)])
y = int_length(x) # in seconds
d = y/86400
d # 60
daily_prob_disp_kel815 = num_roost_switch_kel815/d
daily_prob_disp_kel815 # 0.1011129

# kel818 - 2 roost switches
plot(kel818$timestamp, kel818$location.long)
plot(kel818$timestamp, kel818$location.lat)
num_roost_switch_kel815 = 2
# duration
x = interval(kel818$timestamp[1], kel818$timestamp[length(kel818$timestamp)])
y = int_length(x) # in seconds
d = y/86400
d # 48
daily_prob_disp_kel818 = num_roost_switch_kel815/d
daily_prob_disp_kel818 # 0.07221715

# way155 - 2 roost switch
plot(way155$timestamp, way155$location.long)
plot(way155$timestamp, way155$location.lat)
# duration
x = interval(way155$timestamp[1], way155$timestamp[length(way155$timestamp)])
y = int_length(x) # in seconds
d = y/86400
d # 24
daily_prob_disp_way155 = 1/d
daily_prob_disp_way155 # 0.04166667

# way179 - 2 roost switches
plot(way179$timestamp, way179$location.long)
plot(way179$timestamp, way179$location.lat)
# duration
x = interval(way179$timestamp[1], way179$timestamp[length(way179$timestamp)])
y = int_length(x) # in seconds
d = y/86400
d # 36
daily_prob_disp_way179 = 2/d
daily_prob_disp_way179 # 0.05405622

# average daily prob of dispersal all 5 bats
a <- c(daily_prob_disp_har002, daily_prob_disp_har003, daily_prob_disp_kel815, daily_prob_disp_kel818)
ave_prob_disp = mean(a) # 0.02292944
ave_prob_disp
# biweekly dispersal probability
# what is the probability that each bat disperses at least once in a two week period?
biweek_prob = 1 - (1 - ave_prob_disp)^14
# is this nearest neighbor? should I make it scaled by distance between the sites?
biweek_prob
#--------------------------------------
# MAKE TRACK PLOTS AND ANIMATIONS
# Make spatial
har002.sp <- har002
coordinates(har002.sp) <- c("location.long", "location.lat")
har003.sp <- har003
coordinates(har003.sp) <- c("location.long", "location.lat")
kel815.sp <- kel815
coordinates(kel815.sp) <- c("location.long", "location.lat")
kel818.sp <- kel818
coordinates(kel818.sp) <- c("location.long", "location.lat")
way155.sp <- way155
coordinates(way155.sp) <- c("location.long", "location.lat")
way179.sp <- way179
coordinates(way179.sp) <- c("location.long", "location.lat")

# Set coordinate system & projection
crs_wgs84 <- CRS(SRS_string = "EPSG:4326") # WGS 84 has EPSG code 4326
class(crs_wgs84)

slot(har002.sp, "proj4string") <- crs_wgs84
slot(har003.sp, "proj4string") <- crs_wgs84
slot(kel815.sp, "proj4string") <- crs_wgs84
slot(kel818.sp, "proj4string") <- crs_wgs84
slot(way155.sp, "proj4string") <- crs_wgs84
slot(way179.sp, "proj4string") <- crs_wgs84


# Make back into dataframe (but include date for our animation)
# ggmap and gganimate use dataframes for plotting
points.geo.har002 <- as.data.frame(har002.sp@coords)
points.geo.har002$id <- har002.sp@data$individual.local.identifier # add individual identifier
points.geo.har002$date <- as.POSIXct(har002.sp@data$timestamp, format = "%m/%d/%y %H:%M") # Important! the variable for revealing in the animation must be

points.geo.har003 <- as.data.frame(har003.sp@coords)
points.geo.har003$id <- har003.sp@data$individual.local.identifier # add individual identifier
points.geo.har003$date <- as.POSIXct(har003.sp@data$timestamp, format = "%m/%d/%y %H:%M") # Important! the variable for revealing in the animation must be

points.geo.kel815 <- as.data.frame(kel815.sp@coords)
points.geo.kel815$id <- kel815.sp@data$individual.local.identifier # add individual identifier
points.geo.kel815$date <- as.POSIXct(kel815.sp@data$timestamp, format = "%m/%d/%y %H:%M") # Important! the variable for revealing in the animation must be

points.geo.kel818 <- as.data.frame(kel818.sp@coords)
points.geo.kel818$id <- kel818.sp@data$individual.local.identifier # add individual identifier
points.geo.kel818$date <- as.POSIXct(kel818.sp@data$timestamp, format = "%m/%d/%y %H:%M") # Important! the variable for revealing in the animation must be

points.geo.way155 <- as.data.frame(way155.sp@coords)
points.geo.way155$id <- way155.sp@data$individual.local.identifier # add individual identifier
points.geo.way155$date <- as.POSIXct(way155.sp@data$timestamp, format = "%m/%d/%y %H:%M") # Important! the variable for revealing in the animation must be

points.geo.way179 <- as.data.frame(way179.sp@coords)
points.geo.way179$id <- way179.sp@data$individual.local.identifier # add individual identifier
points.geo.way179$date <- as.POSIXct(way179.sp@data$timestamp, format = "%m/%d/%y %H:%M") # Important! the variable for revealing in the animation must be



register_google(key = "AIzaSyCwa0OOmg7nRXgOrBZgBxgvmdn1h_bIO7g")


har002.map <- get_map(location = c(lon = mean(har002.sp@coords[,1]), 
                                  lat = mean(har002.sp@coords[,2])), 
                     source = "google", 
                     zoom = 12,
                     maptype = 'satellite')
ggmap(har002.map)

har003.map <- get_map(location = c(lon = mean(har003.sp@coords[,1]), 
                                   lat = mean(har003.sp@coords[,2])), 
                      source = "google", 
                      zoom = 12,
                      maptype = 'satellite')
ggmap(har003.map)

kel815.map <- get_map(location = c(lon = mean(kel815.sp@coords[,1]), 
                                   lat = mean(kel815.sp@coords[,2])), 
                      source = "google", 
                      zoom = 12,
                      maptype = 'satellite')
ggmap(kel815.map)

kel818.map <- get_map(location = c(lon = mean(kel818.sp@coords[,1]), 
                                   lat = mean(kel818.sp@coords[,2])), 
                      source = "google", 
                      zoom = 12,
                      maptype = 'satellite')
ggmap(kel818.map)

way155.map <- get_map(location = c(lon = mean(way155.sp@coords[,1]), 
                                   lat = mean(way155.sp@coords[,2])), 
                      source = "google", 
                      zoom = 12,
                      maptype = 'satellite')
ggmap(way155.map)

way179.map <- get_map(location = c(lon = mean(way179.sp@coords[,1]), 
                                   lat = mean(way179.sp@coords[,2])), 
                      source = "google", 
                      zoom = 12,
                      maptype = 'satellite')
ggmap(way179.map)

#---------------
# HAR002
har002.paths <- ggmap(har002.map) + 
  geom_point(data = points.geo.har002, aes(x = location.long, y = location.lat, colour = as.character(id))) +
  geom_path(data = points.geo.har002, aes(x = location.long, y = location.lat, colour = as.character(id), group = as.character(id))) +
  #theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") 
  
har002.paths

path.animate.plot.har002 <- har002.paths +
  transition_reveal(along = date) +
  labs(title = 'Date: {frame_along}')

anananimate(path.animate.plot.har002,
        fps = 5, renderer = gifski_renderer())

# HAR003
har003.paths <- ggmap(har003.map) + 
  geom_point(data = points.geo.har003, aes(x = location.long, y = location.lat, colour = as.character(id))) +
  geom_path(data = points.geo.har003, aes(x = location.long, y = location.lat, colour = as.character(id), group = as.character(id))) +
  #theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") 

har003.paths

path.animate.plot.har003 <- har003.paths +
  transition_reveal(along = date) +
  labs(title = 'Date: {frame_along}')  

animate(path.animate.plot.har003,
        fps = 5, renderer = gifski_renderer())

# KEL815
kel815.paths <- ggmap(kel815.map) + 
  geom_point(data = points.geo.kel815, aes(x = location.long, y = location.lat, colour = as.character(id))) +
  geom_path(data = points.geo.kel815, aes(x = location.long, y = location.lat, colour = as.character(id), group = as.character(id))) +
  #theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") 

kel815.paths

path.animate.plot.kel815 <- kel815.paths +
  transition_reveal(along = date) +
  labs(title = 'Date: {frame_along}')  

animate(path.animate.plot.kel815,
        fps = 5, renderer = gifski_renderer())


# KEL818
kel818.paths <- ggmap(kel818.map) + 
  geom_point(data = points.geo.kel818, aes(x = location.long, y = location.lat, colour = as.character(id))) +
  geom_path(data = points.geo.kel818, aes(x = location.long, y = location.lat, colour = as.character(id), group = as.character(id))) +
  #theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") 

kel818.paths

path.animate.plot.kel818 <- kel818.paths +
  transition_reveal(along = date) +
  labs(title = 'Date: {frame_along}')

animate(path.animate.plot.kel818,
        fps = 5, renderer = gifski_renderer())

# WAY155
way155.paths <- ggmap(way155.map) + 
  geom_point(data = points.geo.way155, aes(x = location.long, y = location.lat, colour = as.character(id))) +
  geom_path(data = points.geo.way155, aes(x = location.long, y = location.lat, colour = as.character(id), group = as.character(id))) +
  #theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") 

way155.paths

path.animate.plot.way155 <- way155.paths +
  transition_reveal(along = date) +
  labs(title = 'Date: {frame_along}')

animate(path.animate.plot.way155,
        fps = 5, renderer = gifski_renderer())

# WAY179
way179.paths <- ggmap(way179.map) + 
  geom_point(data = points.geo.way179, aes(x = location.long, y = location.lat, colour = as.character(id))) +
  geom_path(data = points.geo.way179, aes(x = location.long, y = location.lat, colour = as.character(id), group = as.character(id))) +
  #theme(legend.position = c(0.15, 0.80)) +
  labs(x = "Longitude", y = "Latitude") 

way179.paths

path.animate.plot.way179 <- way179.paths +
  transition_reveal(along = date) +
  labs(title = 'Date: {frame_along}')

animate(path.animate.plot.way179,
        fps = 5, renderer = gifski_renderer())


# SAVE
anim_save(path.animate.plot,
          file = "animatedpaths.gif")










