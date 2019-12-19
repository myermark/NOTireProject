####---
# Title: TireProjectModeling.R
#
# Author: Mark Myer
#
# Date: 12/18/2019
#
# Purpose: To create models for mosquito larvae modeling in tires
#
# R version 3.6.1 Action of the Toes
####---

#Load packages
library(INLA)
library(rgdal)
library(sf)
library(ggmap)
library(dplyr)
library(raster)

#Import data 
dat.selected <- readRDS(file = "Mosquito_Variables_Selected.rds")

#Plot the sampling locations
border <- readOGR("/Volumes/Mark Drive/NOLA/GIS/New Orleans/CityLimits.shp")
border <- spTransform(border, crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) #reproject to decimal degrees in WGS84

#Visualize using ggmaps
nola_stamen <- get_stamenmap(bbox = border@bbox, zoom = 11, maptype = "terrain")
map <- ggmap(nola_stamen)
loc <- unique(data.frame(long = dat.selected[[1]]$LatY, lat = dat.selected[[1]]$LongX)) #Note that lat and long are reversed in the dataset. 
map + geom_point(data = loc, pch=21, stroke = 1, aes(x=long, y= lat))  + 
  guides(size = F) +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Tire Sample Locations") +
  theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
        plot.title = element_text(size=20, face = "bold"))
