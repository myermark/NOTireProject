####---
# Title: TireProjectModelEval.R
#
# Author: Mark Myer
#
# Date: 1/7/2020
#
# Purpose: To validate and evaluate mosquito larvae models in tires
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
library(dplyr)
library(raster)
library(gstat)
library(tigris)

#Load data----
load(".Rdata")

#Perform model evaluation and validation
#Plot observed vs. fits for the gamma portion of the model 
lapply(1:length(dat.selected), function (i) {
  lapply(1:length(modlist), function (j) {
    tryCatch(
      do.call(function (i,j) { 
        N <- nrow(dat.selected[[i]])
        mu <- results.gam[[i]][[j]]$summary.fitted.values[1:N, "mean"]
        r <- results.gam[[i]][[j]]$summary.hyperpar[1, "mean"]
        varY <- mu^2 / r
        E <- (dat.selected[[i]]$MosqPerL - mu) / (sqrt(varY))
        
        dat.temp <- dat.selected[[i]]
        dat.temp$mu <- mu
        dat.temp$E <- E
        dat.pos <- filter(dat.temp, MosqPerL > 0)
        
        #tiff(filename = paste0("./Figures/Validation/", names(dat.selected)[i], modlist[j], " ObservedVsFit.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
        
        plot(y = dat.pos$MosqPerL, 
             x = dat.pos$mu, 
             xlab = "Fitted values",
             ylab = "Observed data", 
             main = paste(names(dat.selected)[i], modlist[j], "Observed vs. Fits"),
             xlim = range(dat.temp$MosqPerL), 
             ylim = range(dat.temp$MosqPerL))
        #dev.off()
      }, args = list(i,j))
    , error = function (e) print("Error"))
  })
})

#Combine the models to get a joint mean and variance
#Plot the residuals vs the fitted values for each model 
lapply(1:length(dat.selected), function (i) {
  lapply(1:length(modlist), function (j) {
    tryCatch({
      pi = results.bin[[i]][[j]]$summary.fitted.values[1:N, "mean"]
      mu = results.gam[[i]][[j]]$summary.fitted.values[1:N, "mean"]
      exp = pi * mu
      r <- results.gam[[i]][[j]]$summary.hyperpar[1, "mean"]
      var= (mu ^ 2) * (pi * r + pi - pi^2 * r) / r
      E = as.numeric((dat.selected[[i]]$MosqPerL - exp) / sqrt(var))
      #tiff(filename = paste0("./Figures/Validation/", names(dat.selected)[i], modlist[j], " ResidVsFit.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
      plot(y = E, 
           x = exp, 
           xlab = "Fitted values",
           ylab = "Residuals", 
           xlim = range(dat.selected[[i]]$MosqPerL),
           main = paste(names(dat.selected)[i], modlist[j], "Residuals vs. Fits") 
      )
      #dev.off()
    }, error = function(e) print("Error"))
  })
})

#Plot the spatial field 
PlotField2 <- function(field, mesh, ContourMap, xlim, ylim, Add=FALSE, MyMain, ...){
  stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
  if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
  
  # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, 
                              xlim = xlim, 
                              ylim = ylim, 
                              dims = c(300, 300))
  # The function inla.mesh.project can then 
  # be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  # And plot the whole thing
  image.plot(list(x = proj$x, 
                  y = proj$y,
                  z = field.proj), 
             xlim = xlim, 
             ylim = ylim,
             asp = 1,
             add = Add,
             main = MyMain,
             axes = FALSE,
             ...)  
}

#Visualize using road shapefile from TIGERLINE
nola_roads <- roads(state="Louisiana", county = "Orleans")
nola_roads <- spTransform(nola_roads, crs("+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
border_km <- spTransform(border, crs("+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
par(oma=c( 0,0,0,0), mar = c(4,4,1,1)) # margin of 4 spaces width at right hand side

for(i in 1:length(dat.selected)) {
  km_loc <- unique(data.frame(long = dat.selected[[i]]$Adj_X, lat = dat.selected[[i]]$Adj_Y))
  #Plot binomial models
  #tiff(filename = paste0("./Figures/Spatial/", names(dat.selected)[i], " Binomial.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
  w.bin <- results.bin[[i]]$spatial$summary.random$spatial$mean
  PlotField2(field = w.bin, 
             mesh = meshes[[i]], 
             xlim = range(meshes[[i]]$loc[,1]), 
             ylim = range(meshes[[i]]$loc[,2]),
             MyMain = paste(names(dat.selected)[i], "Occurrence Model (Binomial)")
  )
  axis(1); axis(2)
  points(x = km_loc[,1],
         y = km_loc[,2], 
         cex = 0.5, 
         col = "black", 
         pch = 16)
  plot(nola_roads, add=T, lwd = 0.75)
  plot(border_km, add=T)
  #dev.off()
  
  #Plot gamma models
  #tiff(filename = paste0("./Figures/Spatial/", names(dat.selected)[i], " Gamma.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw", type = "cairo")
  w.gam <- results.gam[[i]]$spatial$summary.random$spatial$mean
  PlotField2(field = w.gam, 
             mesh = meshes[[i]], 
             xlim = range(meshes[[i]]$loc[,1]), 
             ylim = range(meshes[[i]]$loc[,2]),
             MyMain = paste(names(dat.selected)[i], "Abundance Model (Gamma)")
  )
  axis(1); axis(2)
  points(x = km_loc[,1],
         y = km_loc[,2], 
         cex = 0.5, 
         col = "black", 
         pch = 16)
  plot(nola_roads, add=T, lwd = 0.75)
  plot(border_km, add=T)
  #dev.off()
}