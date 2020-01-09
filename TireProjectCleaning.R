####---
# Title: TireProjectCleaning.R
#
# Author: Mark Myer
#
# Date: 12/11/2019
#
# Purpose: To import and clean tire project data for analysis. 
#
# R version 3.6.1 Action of the Toes
####---

#Load packages
library(raster)
library(dplyr)
library(forcats)
library(ggplot2)
library(readr)
library(readxl)
library(tidyr)
library(sp)
library(rgdal)

#Import data
tire_rawdat <- read_excel("MMEdit_Tire_Data_ELC_Updated_Nov20_2019.xlsx") %>% mutate(MosqCount = as.numeric(MosqCount))
tiresite_vars <- read_excel("MMEdit_Tire_Data_ELC_Updated_Nov20_2019.xlsx", sheet = 3)

#The names of latitude and longitude are reversed. Need to fix. 
names(tire_rawdat)[which(names(tire_rawdat) == "LongX")] = "Placeholder"
names(tire_rawdat)[which(names(tire_rawdat) == "LatY")] = "LongX"
names(tire_rawdat)[which(names(tire_rawdat) == "Placeholder")] = "LatY"

names(tire_rawdat)[which(names(tire_rawdat) == "Adj_X")] = "Placeholder"
names(tire_rawdat)[which(names(tire_rawdat) == "Adj_Y")] = "Adj_X"
names(tire_rawdat)[which(names(tire_rawdat) == "Placeholder")] = "Adj_Y"

#Match site variables to observations
tire_fulldat <- merge(tire_rawdat, tiresite_vars, by="WayPt_ID")

#Change site variable to just a number ("2" rather than "WP2")
tire_fulldat <- mutate(tire_fulldat, WayPt_ID = substring(WayPt_ID, 3))

#Standardize the number of mosquitoes found by the volume of water sampled (mosquitoes per liter rounded to the nearest mosquito) and put it in front
tire_fulldat <- mutate(tire_fulldat, MosqPerL = round(MosqCount/water_L, 0)) %>% dplyr::select(MosqPerL, everything()) 

#Convert the locations to a shapefile then reproject to kilometers to make INLA happy
coordinates(tire_fulldat) <- ~Adj_X + Adj_Y
proj4string(tire_fulldat) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
tire_fulldat <- spTransform(tire_fulldat, crs("+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))

tire_fulldat <- data.frame(tire_fulldat)

#Make epiweek start at 1 for INLA
tire_fulldat$INLAWeek <- tire_fulldat$EpiWeek - (min(tire_fulldat$EpiWeek) -1)

#Save the modified data
write.csv(tire_fulldat, file="TireData121219.csv")

