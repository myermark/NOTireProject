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

#Load packages----
library(INLA)
library(rgdal)
library(sf)
library(ggmap)
library(dplyr)
library(raster)
library(gstat)

#Import data---- 
dat.selected <- readRDS(file = "Mosquito_Variables_Selected.rds")

#Visualize data and prepare for modeling----
#Plot the sampling locations
border <- readOGR("/Volumes/Mark Drive/NOLA/GIS/New Orleans/CityLimits.shp")
border <- spTransform(border, crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) #reproject to decimal degrees in WGS84

#Visualize using ggmaps
nola_stamen <- get_stamenmap(bbox = border@bbox, zoom = 11, maptype = "terrain")
map <- ggmap(nola_stamen)
loc <- unique(data.frame(long = dat.selected[[1]]$LongX, lat = dat.selected[[1]]$LatY)) 
map + geom_point(data = loc, pch=21, stroke = 1, aes(x=long, y= lat))  + 
  guides(size = F) +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Tire Sample Locations") +
  theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
        plot.title = element_text(size=20, face = "bold"))

#Check the distribution of responses to determine link function for GLMM (if any)
#All responses
lapply(1:length(dat.selected), function(i){
  hist(dat.selected[[i]]$MosqPerL, main=paste(names(dat.selected)[i], "All Responses"), xlab = "Mosquito Larvae Per Liter")
})
#Zero-excluded
lapply(1:length(dat.selected), function(i){
  x = filter(dat.selected[[i]], MosqPerL != 0)
  hist(x$MosqPerL, main=paste(names(dat.selected)[i], "Zeroes Excluded"), xlab = "Mosquito Larvae Per Liter")
})

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
  nonspatial.result[[i]] <- glm(as.formula(nonspatial.formulas[[i]]), data = dat.selected[[i]], family = gaussian)
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
  print(plot(variograms[[i]], main = names(dat.selected)[[i]], col = 1,))
}

#Define the INLA model ----
#Create mesh
meshes <- lapply(1:length(dat.selected), function (i) {
  loc <- list()
  bnd <- list()
  mesh <- list()
  loc[[i]] <- cbind(dat.selected[[i]]$Adj_X, dat.selected[[i]]$Adj_Y) 
  bnd[[i]] <-inla.nonconvex.hull(loc[[i]]) #Makes a nonconvex hull around the points
  mesh[[i]] <-inla.mesh.2d(boundary=bnd[[i]], cutoff=0.5, max.edge=c(1,4)) #Max-edge is based on estimated correlation distance from semivariograms
  return(mesh[[i]])
  }
)

#Make the SPDE
spdes <- lapply(1:length(meshes), function (i) {
  spde <- list()
  spde[[i]] <- inla.spde2.pcmatern(
  mesh=meshes[[i]], alpha=2, ### mesh and smoothness parameter
  prior.range=c(1, 0.05), ### P(range<1km)=0.05
  prior.sigma=c(0.5, 0.05) ### P(sigma>0.5)=0.05
    )
  return(spde[[i]])
  }
)

#Create the projector matrix 
As <- lapply(1:length(meshes), function (i) {
  loc <- list()
  A <- list()
  loc[[i]] <- cbind(dat.selected[[i]]$Adj_X, dat.selected[[i]]$Adj_Y)
  A[[i]] <- inla.spde.make.A(meshes[[i]], loc[[i]])
  return(A[[i]])
  }
)

#Create the stacks 
#Binomial stack, where anything that is not zero is coded as 1. Presence/absence.
stacks.bin <- lapply(1:length(dat.selected), function (i) {
  stack <- list()
  stack[[i]] <- inla.stack(
    tag="Fit",
    data=list(y=ifelse(dat.selected[[i]]$MosqPerL == 0, 0, 1)),
    A=list(As[[i]],1),
    effects=list(list(spatial = 1:spdes[[i]]$n.spde), data.frame(dat.selected[[i]])
                 # cbind(dplyr::select(data.frame(dat.selected[[i]]), 1:(ncol(dat.selected[[i]]) - 5)), 
                 #       (dplyr::select(data.frame(dat.selected[[i]]), (ncol(dat.selected[[i]]) - 4):ncol(dat.selected[[i]])) %>%
                 #         mutate_if(is.numeric, scale))
                 #       )
      )
    )
  return(stack[[i]])
  }
)

stacks.gam <- lapply(1:length(dat.selected), function (i) {
  stack <- list()
  stack[[i]] <- inla.stack(
    tag="Fit",
    data=list(y=ifelse(dat.selected[[i]]$MosqPerL == 0, NA, dat.selected[[i]]$MosqPerL)),
    A=list(As[[i]],1),
    effects=list(list(spatial = 1:spdes[[i]]$n.spde), data.frame(dat.selected[[i]])
                 # cbind(dplyr::select(data.frame(dat.selected[[i]]), 1:(ncol(dat.selected[[i]]) - 5)), 
                 #       (dplyr::select(data.frame(dat.selected[[i]]), (ncol(dat.selected[[i]]) - 4):ncol(dat.selected[[i]])) %>%
                 #          mutate_if(is.numeric, scale))
                 # )
    )
  )
  return(stack[[i]])
}
)

#Define the INLA formulas
formulas <- lapply(1:length(dat.selected), function(i) {
  formula <- list()
  len <- length(dat.selected[[i]])
  formula[[i]] <- paste("y ~ -1 + ",
                        names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                        "+",
                        names(dat.selected[[i]])[len-3], 
                        "+",
                        names(dat.selected[[i]])[len-2], 
                        "+",
                        names(dat.selected[[i]])[len-1], 
                        "+",
                        names(dat.selected[[i]])[len], 
                        "+ f(spatial, model=spdes[[i]])"#, #Spatial model component
                        #"+ f(WayPt_ID, model = 'iid')"#, #Random intercept by waypoint site
                        #"+ f(INLAWeek, model = 'ar1', hyper = list(theta1=list(prior='pc.prec', param=c(0.5,0.5)), theta2=list(prior='pc.cor1', param=c(0.9,0.9))))" #Temporal model component
  )
  return(formula[[i]])
}
)

#Define prior for precision
prec.prior <- list(prior='pc.prec', param=c(0.5, 0.5))

#Run the binomial part of the models
results.bin <- lapply(1:length(dat.selected), function (i) {
  inlas.bin <- list()
  inlas.bin[[i]] <- inla(as.formula(formulas[[i]]),
                      family="binomial",
                      Ntrials = 1,
                      data=inla.stack.data(stacks.bin[[i]]),
                      control.predictor=list(compute=TRUE, A=inla.stack.A(stacks.bin[[i]]), link = 1), 
                      control.inla=list(int.strategy='auto', correct = TRUE, correct.factor = 10),
                      control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE,po=TRUE,config=TRUE),
                      control.fixed=list(expand.factor.strategy ='inla'),
                      verbose = F)
  return(inlas.bin[[i]])
  }
) 
names(results.bin) <- splist

#Run the gamma part of the models
results.gam <- lapply(1:length(dat.selected), function (i) {
  inlas.gam <- list()
  inlas.gam[[i]] <- inla(as.formula(formulas[[i]]),
                         family="gamma",
                         data=inla.stack.data(stacks.gam[[i]]),
                         control.predictor=list(compute=TRUE, A=inla.stack.A(stacks.gam[[i]]), link = 1), 
                         control.inla=list(int.strategy='auto', correct = TRUE, correct.factor = 10),
                         control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE,po=TRUE,config=TRUE),
                         control.fixed=list(expand.factor.strategy ='inla'),
                         verbose = F)
  return(inlas.gam[[i]])
}
) 
names(results.gam) <- splist

#Check the results
summary(results.bin$`A. aeg`)
summary(results.bin$`A. albo`)
summary(results.bin$`Cx. salinarius`)
summary(results.bin$`Cx. quinq`)
summary(results.bin$`A. crucians`)
summary(results.bin$`Cx. restuans`)
summary(results.bin$`Cx. nigripalpus`)

summary(results.gam$`A. aeg`)
summary(results.gam$`A. albo`)
summary(results.gam$`Cx. salinarius`)
summary(results.gam$`Cx. quinq`)
summary(results.gam$`A. crucians`)
summary(results.gam$`Cx. restuans`)
summary(results.gam$`Cx. nigripalpus`)
