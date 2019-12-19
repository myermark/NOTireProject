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
