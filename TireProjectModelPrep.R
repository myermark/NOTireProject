####---
# Title: TireProjectModelPrep.R
#
# Author: Mark Myer
#
# Date: 1/7/2020
#
# Purpose: To create models for mosquito larvae modeling in tires
#
# R version 3.6.1 Action of the Toes
####---


#Load packages----
library(INLA)
library(fields)
library(rgdal)
library(sf)
library(sp)
library(ggmap)
library(raster)
library(gstat)
library(tigris)
library(formattable)
library(dplyr)
source("/Volumes/Mark Drive/Papers and Textbooks/Highstat Guide to INLA/HighstatLibV11.R") #Remember to cite these helper functions from Highstat

#Import data---- 
#Import data
tire_fulldat <- read.csv("TireData121219.csv") %>% dplyr::select(., -c(X, comments, DispX, DispY, optional)) %>% na.omit()
dat.selected <- readRDS(file = "Mosquito_Variables_Selected.rds")
splist <- unique(tire_fulldat$MosqSpp)

#Visualize data and prepare for modeling----
#Plot the sampling locations
border <- readOGR("/Volumes/Mark Drive/NOLA/GIS/New Orleans/CityLimits.shp")
border <- spTransform(border, crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) #reproject to decimal degrees in WGS84

#Visualize using ggmaps
nola_stamen <- get_stamenmap(bbox = border@bbox, zoom = 11, maptype = "terrain")
map <- ggmap(nola_stamen)
loc <- unique(data.frame(long = dat.selected[[1]]$LongX, lat = dat.selected[[1]]$LatY)) 

# Distance between samples
D <- dist(unique(data.frame(long = dat.selected[[1]]$Adj_X, lat = dat.selected[[1]]$Adj_Y)))
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sampling locations (km)",
     ylab = "Frequency")

map + geom_point(data = loc, pch=21, stroke = 1, aes(x=long, y= lat))  + 
  guides(size = F) +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Tire Sample Locations") +
  theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
        plot.title = element_text(size=20, face = "bold"))

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
#Plot the relative mosquito counts per tire sample site
for(i in 1:length(dat.selected)) {
  temp  <- dflist[[i]] %>% group_by(WayPt_ID) %>% summarize(LatY = mean(LatY, na.rm=T), LongX = mean(LongX, na.rm=T), MosqCount = sum(MosqCount, na.rm=T))
  #tiff(filename = paste0("./Figures/Counts/", names(dat.selected)[i], " MappedCounts.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
  print(map + geom_point(data = arrange(temp, MosqCount), pch=19, stroke = 1, aes(x=LongX, y= LatY, size = MosqCount, col = MosqCount))  + 
          guides(size = F) +
          labs(x = "Longitude", y = "Latitude") + 
          ggtitle(paste0(names(dflist)[i], " Total Counts")) +
          scale_color_gradient(high = "red", low = "blue") +
          theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
                axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
                plot.title = element_text(size=20, face = "bold"))
  )
  #dev.off()
}

#Check the distribution of responses to determine link function for GLMM (if any)
#All responses
for(i in 1:length(dat.selected)) {
  tiff(filename = paste0("./Figures/Response Distributions/", names(dat.selected)[i], "WithZeroes.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
  plot(table(dat.selected[[i]]$MosqPerL), main=paste(names(dat.selected)[i], "All Responses"), xlab = "Mosquito Larvae Per Liter", ylab = "Times Occurring")
  dev.off()
}
#Zero-excluded
for(i in 1:length(dat.selected)) {
  x = filter(dat.selected[[i]], MosqPerL != 0)
  tiff(filename = paste0("./Figures/Response Distributions/", names(dat.selected)[i], "ZeroesExcluded.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
  plot(table(x$MosqPerL), main=paste(names(dat.selected)[i], "Zeroes Excluded"), xlab = "Mosquito Larvae Per Liter", ylab = "Times Occurring")
  dev.off()  
}

#Determine the percentage of zero responses by species
zero.counts <- data.frame(species = splist, zeroes = rep(NA, times = length(splist)), percent.zero =  rep(NA, times = length(splist)))
for(i in 1:length(dat.selected)) {
  zero.counts$species[i] = splist[i]
  zero.counts$zeroes[i] = sum(dat.selected[[i]]$MosqPerL == 0)
  zero.counts$percent.zero[i] = round(sum(dat.selected[[i]]$MosqPerL == 0) / nrow(dat.selected[[i]]), 2)  
}

formattable(zero.counts, align =c("l", rep("r", ncol(zero.counts) - 1)))

#Plot the response vs. each predictor variable by species
for(i in 1:length(dat.selected)) {
len <- length(dat.selected[[i]])
covars <- c(names(dat.selected[[i]])[len-4], 
               names(dat.selected[[i]])[len-3], 
               names(dat.selected[[i]])[len-2], 
               names(dat.selected[[i]])[len-1], 
               names(dat.selected[[i]])[len])
tiff(filename = paste0("./Figures/Variable Plots/", names(dat.selected)[i], "VarScatterplot.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
MyMultipanel.ggp2(Z = dat.selected[[i]], 
                  varx = covars, 
                  vary = "MosqPerL", 
                  ylab = "Mosquito Larvae Per Liter",
                  main = names(dat.selected)[i],
                  addSmoother = FALSE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)
dev.off()
}
#Define the nonspatial GLM model formulas
nonspatial.formulas <- lapply(1:length(dat.selected), function(i) {
  formula <- list()
  len <- length(dat.selected[[i]])
  formula[[i]] <- paste("MosqPerL ~ ",
                        names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                        "+",
                        names(dat.selected[[i]])[len-3], 
                        "+",
                        names(dat.selected[[i]])[len-2], 
                        "+",
                        names(dat.selected[[i]])[len-1], 
                        "+",
                        names(dat.selected[[i]])[len]
  )
  return(formula[[i]])
}
)

#Fit nonspatial GLM models
nonspatial.results <- lapply(1:length(dat.selected), function (i) {
  nonspatial.result <- list()
  nonspatial.result[[i]] <- glm(as.formula(nonspatial.formulas[[i]]), data = dat.selected[[i]], family = poisson)
  return(nonspatial.result[[i]])
}
)
names(nonspatial.results) <- splist

nonspatial.resids <- lapply(1:length(nonspatial.results), function (i) {
  resid <- list()
  resid[[i]] <- dat.selected[[i]]$MosqPerL - predict(nonspatial.results[[i]]) 
}
)
names(nonspatial.resids) <- splist

#Assess spatial autocorrelation using a variogram
plot.vario <- function (residuals, Xkm, Ykm, cutoff) {
  MyData1 <- data.frame(E = residuals, 
                        Xkm = Xkm, 
                        Ykm = Ykm)
  coordinates(MyData1) <- c("Xkm", "Ykm")
  V = variogram(E ~ 1, data = MyData1, cressie = F, cutoff)
  return(V)
}

variograms = list()
variograms <- lapply(1:length(dat.selected), function(i) {
  plot.vario(
    residuals = nonspatial.resids[[i]],
    Xkm = dat.selected[[i]]$Adj_X,
    Ykm = dat.selected[[i]]$Adj_Y,
    cutoff = 10
  )
}
)

for (i in 1:length(dat.selected)) {
  tiff(filename = paste0("./Figures/Variograms/", names(dat.selected)[i], " Vario.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
  lo <- loess(variograms[[i]]$gamma ~ variograms[[i]]$dist)
  xl <- seq(min(variograms[[i]]$dist),max(variograms[[i]]$dist), (max(variograms[[i]]$dist) - min(variograms[[i]]$dist))/1000)
  plot(x = variograms[[i]]$dist, y = variograms[[i]]$gamma, main = names(dat.selected)[[i]], xlab = "Distance", ylab = "Semivariance")
  lines(xl, predict(lo, xl), lwd = 1.5, col = "blue")
  dev.off()
}

#Save environment for later loading by the modeling script
rm(D, i, variograms, plot.vario, lo, covars, len, xl, loc, map, nola_stamen)
save.image()
