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
library(forcats)
library(ggplot2)
library(readr)
library(readxl)
library(tidyr)
library(sp)
library(rgdal)
library(dplyr)
<<<<<<< HEAD
source("/Volumes/Mark Drive/Papers and Textbooks/Highstat Guide to INLA/HighstatLibV11.R") #Remember to cite these helper functions from Highstat
=======
>>>>>>> e7d7cf6dc2ed7bb1703f2a1455d8b3eade273db1

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

#Convert the land cover variables to fractions of the total buffer area
<<<<<<< HEAD
radius = 2 #2km
buffer_area = pi * radius^2
tire_fulldat <- tire_fulldat %>% mutate(Open.Water = Open.Water / buffer_area,
                                        Developed.Open.Space = Developed.Open.Space / buffer_area,
                                        Developed.Low.Intensity = Developed.Low.Intensity / buffer_area,
                                        Developed.Medium.Intensity = Developed.Medium.Intensity / buffer_area,
                                        Developed.High.Intensity = Developed.High.Intensity / buffer_area,
                                        Barren.Land = Barren.Land / buffer_area,
                                        Deciduous.Forest = Deciduous.Forest / buffer_area, 
                                        Shrub.Scrub = Shrub.Scrub / buffer_area, 
                                        Grassland.Herbaceous = Grassland.Herbaceous / buffer_area, 
                                        Pasture.Hay = Pasture.Hay/ buffer_area, 
                                        Cultivated.Crops = Cultivated.Crops/ buffer_area, 
                                        Woody.Wetlands = Woody.Wetlands/ buffer_area, 
                                        Emergent.Herbaceous.Wetlands = Emergent.Herbaceous.Wetlands/ buffer_area)
=======
tire_fulldat <- tire_fulldat %>% mutate(Open.Water = Open.Water / Area_km2,
                                        Developed.Open.Space = Developed.Open.Space / Area_km2,
                                        Developed.Low.Intensity = Developed.Low.Intensity / Area_km2,
                                        Developed.Medium.Intensity = Developed.Medium.Intensity / Area_km2,
                                        Developed.High.Intensity = Developed.High.Intensity / Area_km2,
                                        Barren.Land = Barren.Land / Area_km2,
                                        Deciduous.Forest = Deciduous.Forest / Area_km2, 
                                        Shrub.Scrub = Shrub.Scrub / Area_km2, 
                                        Grassland.Herbaceous = Grassland.Herbaceous / Area_km2, 
                                        Pasture.Hay = Pasture.Hay/ Area_km2, 
                                        Cultivated.Crops = Cultivated.Crops/ Area_km2, 
                                        Woody.Wetlands = Woody.Wetlands/ Area_km2, 
                                        Emergent.Herbaceous.Wetlands = Emergent.Herbaceous.Wetlands/ Area_km2)
>>>>>>> e7d7cf6dc2ed7bb1703f2a1455d8b3eade273db1

#Combine the detritus types into organic/nonorganic
tire_fulldat <- tire_fulldat %>% mutate(org_debris = seeds + sticks + leaf + algae_moss + shells + grass, 
                                        inorg_debris = rubber + man.made + rocks)

#Change the factor variables to 0/1 rather than y/n
tire_fulldat <- tire_fulldat %>% mutate(veg_in_tire = ifelse(veg_in_tire == "y", 1, 0),
                                        cop = ifelse(cop == "y", 1, 0), 
                                        ostra = ifelse(ostra == "y", 1, 0),
                                        daph = ifelse(daph == "y", 1, 0))  %>% 
                                 mutate(veg_in_tire = factor(veg_in_tire ),
                                        cop = factor(cop), 
                                        ostra = factor(ostra),
                                        daph = factor(daph)) 
  

#Change the name of the veg cover variable to eliminate the trailing period
tire_fulldat <- tire_fulldat %>% mutate(veg_cover = veg_cover.) %>% dplyr::select(-c(veg_cover.))

#Change the order of variables so dependent vars are at the end of the dataset
tire_fulldat <- tire_fulldat %>% dplyr::select(MosqPerL, WayPt_ID, LongX, LatY, Adj_X, Adj_Y, Day, Month, EpiWeek, INLAWeek, MosqCount, everything())

#Remove variables from consideration that we are entirely uninterested in, or are derivatives of others, or researcher aid info that is not intended as covariate
tire_fulldat <- tire_fulldat %>% dplyr::select(-c(Address, 
                                           Sampledtire, 
                                           water_L, 
                                           no_cover, 
                                           open., 
                                           w_cover, 
                                           seeds, 
                                           sticks, 
                                           leaf, 
                                           algae_moss, 
                                           shells, 
                                           grass, 
                                           rubber, 
                                           man.made, 
                                           rocks, 
                                           Buffer, 
                                           Area_km2, 
                                           Barren.Land,
                                           Deciduous.Forest, 
                                           Shrub.Scrub, 
                                           Grassland.Herbaceous, 
                                           Pasture.Hay, 
                                           Cultivated.Crops, 
                                           PovPast12M, 
                                           Unemp, 
                                           Employed, 
                                           NoSchool, 
                                           SomeColleg, 
                                           AssDegree, 
                                           Bachelors, 
                                           Masters, 
                                           PhD, 
                                           ProfDegree, 
                                           Rented_Unocupied, 
                                           SoldUnocup, 
                                           OtherVacan))

<<<<<<< HEAD

=======
>>>>>>> e7d7cf6dc2ed7bb1703f2a1455d8b3eade273db1
#Save the modified data
write.csv(tire_fulldat, file="TireData121219.csv")


