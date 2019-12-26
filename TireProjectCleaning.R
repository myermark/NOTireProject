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
library(dplyr)
library(forcats)
library(ggplot2)
library(readr)
library(readxl)
library(tidyr)

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

#Standardize the number of mosquitoes found by the volume of water sampled (mosquitoes per liter) and put it in front
tire_fulldat <- mutate(tire_fulldat, MosqPerL = MosqCount/water_L) %>% dplyr::select(MosqPerL, everything()) %>% dplyr::select(-c(MosqCount))

#Convert the locations to a shapefile then reproject to kilometers to make INLA happy
coordinates(tire_fulldat) <- ~Adj_X + Adj_Y
proj4string(tire_fulldat) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
tire_fulldat <- spTransform(tire_fulldat, crs("+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))

tire_fulldat <- data.frame(tire_fulldat)

#Make epiweek start at 1 for INLA
tire_fulldat$INLAWeek <- tire_fulldat$EpiWeek - (min(tire_fulldat$EpiWeek) -1)

#Save the modified data
write.csv(tire_fulldat, file="TireData121219.csv")

#Plot each species' counts by week 
for(n in 1:length(splist)) {
  temp <- filter(tire_fulldat, MosqSpp == splist[n])
  tiff(filename=paste0("./Figures/", temp$MosqSpp[1], "_CountByWeek.tiff"), width = 5, height = 4, units = "in", compression = "lzw", type = "cairo", res = 300)
  plot(temp$EpiWeek, temp$MosqPerL, ylim =c(0, max(tire_fulldat$MosqPerL)),
        pch=16, xlab = "Week", ylab = "Larvae Per Liter", main = temp$MosqSpp[1])
  dev.off()
}
rm(temp)

#Determine the relative count of each species captured 
#In total
spp_totals <- tire_fulldat %>% group_by(MosqSpp) %>% summarise(total = sum(MosqCount))
tiff(filename="./Figures/YearlyTotals.tiff", width = 10, height = 8, units = "in", compression = "lzw", type = "cairo", res = 300)
barplot(total ~ fct_reorder(MosqSpp, total, .desc=TRUE), data = spp_totals, 
          ylim = c(0, 5000), ylab = "Mosquito Larvae", xlab = "", main = "Yearly Total")
dev.off()

#By week
tiff(filename="./Figures/WeeklyTotals.tiff", width = 5, height = 4, units = "in", compression = "lzw", type = "cairo", res = 300)
ggplot(tire_fulldat %>% group_by(MosqSpp, EpiWeek) %>% summarise(total = sum(MosqPerL)), aes(x = EpiWeek, y = total, colour=MosqSpp)) + 
  geom_point(pch = 16, size = 2) + 
  labs(title = "Total by Week", x = "Week", y = "Larvae Per Liter", colour = "Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
