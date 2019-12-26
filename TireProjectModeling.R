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
loc <- unique(data.frame(long = dat.selected[[1]]$LongX, lat = dat.selected[[1]]$LatY)) 
map + geom_point(data = loc, pch=21, stroke = 1, aes(x=long, y= lat))  + 
  guides(size = F) +
  labs(x = "Longitude", y = "Latitude") + 
  ggtitle("Tire Sample Locations") +
  theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
        plot.title = element_text(size=20, face = "bold"))

#Check the distribution of responses to determine link function for GLMM (if any)
lapply(1:length(dat.selected), function(i){
  hist(dat.selected[[i]]$MosqPerL, main=names(dat.selected)[i], xlab = "Mosquito Larvae Per Liter")
})

#Define the INLA model 
meshes <- lapply(1:length(dat.selected), function (i) {
  loc <- list()
  bnd <- list()
  mesh <- list()
  loc[[i]] <- cbind(dat.selected[[i]]$Adj_X, dat.selected[[i]]$Adj_Y) 
  bnd[[i]] <-inla.nonconvex.hull(loc[[i]]) #Makes a nonconvex hull around the points
  mesh[[i]] <-inla.mesh.2d(boundary=bnd[[i]], cutoff=1, max.edge=c(2,4))
  return(mesh[[i]])
  }
)

#Make the SPDE
spdes <- lapply(1:length(meshes), function (i) {
  spde <- list()
  spde[[i]] <- inla.spde2.pcmatern(
  mesh=meshes[[i]], alpha=2, ### mesh and smoothness parameter
  prior.range=c(10, 0.9), ### P(range<10km)=0.9
  prior.sigma=c(1, 0.5) ### P(sigma>1)=0.5
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

#Create the stack 
stacks <- lapply(1:length(dat.selected), function (i) {
  stack <- list()
  stack[[i]] <- inla.stack(
    tag="Fit",
    data=list(y=dat.selected[[i]]$MosqPerL),
    A=list(As[[i]],1),
    effects=list(list(spatial = 1:spdes[[i]]$n.spde), data.frame(dat.selected[[i]]))
    )
  return(stack[[i]])
  }
)

#Define the formulas
formulas <- lapply(1:length(dat.selected), function(i) {
  formula <- list()
  len <- length(dat.selected[[i]])
  formula[[i]] <- paste("y ~ 0 +",
                         names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                         "+",
                         names(dat.selected[[i]])[len-3], 
                         "+",
                         names(dat.selected[[i]])[len-2], 
                         "+",
                         names(dat.selected[[i]])[len-1], 
                         "+",
                         names(dat.selected[[i]])[len], 
                         "+ f(spatial, model=spdes[[i]])"#, 
                         #"+ f(INLAWeek, model = 'ar1', hyper = list(theta1=list(prior='pc.prec', param=c(0.5,0.5)), theta2=list(prior='pc.cor1', param=c(0.9,0.9))))"
                         )
  return(formula[[i]])
  }
)

#Define prior for precision
prec.prior <- list(prior='pc.prec', param=c(0.5, 0.5))

#Run the models
results <- list()
#Ae. aegypti
results[[1]]<- inla(as.formula(formulas[[1]]),
                  family="nbinomial",
                  data=inla.stack.data(stacks[[1]]),
                  control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[1]]), link = 1), #link = 1 scales the fitted values with the logit link
                  control.family=list(link='log', hyper=list(theta=prec.prior)),
                  control.inla=list(int.strategy='auto'),
                  control.compute=list(dic=TRUE,cpo=TRUE),
                  control.fixed=list(expand.factor.strategy ='inla'),
                  verbose = F)
#Ae. albopictus
results[[2]]<- inla(as.formula(formulas[[2]]),
                    family="nbinomial",
                    data=inla.stack.data(stacks[[2]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[2]]), link = 1), #link = 1 scales the fitted values with the logit link
                    control.family=list(link='log', hyper=list(theta=prec.prior)),
                    control.inla=list(int.strategy='auto'),
                    control.compute=list(dic=TRUE,cpo=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)
#Cx. salinarius
results[[3]]<- inla(as.formula(formulas[[3]]),
                    family="poisson",
                    data=inla.stack.data(stacks[[3]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[3]]), link = 1), #link = 1 scales the fitted values with the logit link
                    control.family=list(link='log'),#, hyper=list(theta=prec.prior)),
                    control.inla=list(int.strategy='auto'),
                    control.compute=list(dic=TRUE,cpo=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)
#Cx. quinquefasciatus
results[[4]]<- inla(as.formula(formulas[[4]]),
                    family="nbinomial",
                    data=inla.stack.data(stacks[[4]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[4]]), link = 1), #link = 1 scales the fitted values with the logit link
                    control.family=list(link='log', hyper=list(theta=prec.prior)),
                    control.inla=list(int.strategy='auto'),
                    control.compute=list(dic=TRUE,cpo=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)
#A. crucians
results[[5]]<- inla(as.formula(formulas[[5]]),
                    family="nbinomial",
                    data=inla.stack.data(stacks[[5]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[5]]), link = 1), #link = 1 scales the fitted values with the logit link
                    control.family=list(link='log', hyper=list(theta=prec.prior)),
                    control.inla=list(int.strategy='auto'),
                    control.compute=list(dic=TRUE,cpo=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)
#Cx. restuans
results[[6]]<- inla(as.formula(formulas[[6]]),
                    family="nbinomial",
                    data=inla.stack.data(stacks[[6]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[6]]), link = 1), #link = 1 scales the fitted values with the logit link
                    control.family=list(link='log', hyper=list(theta=prec.prior)),
                    control.inla=list(int.strategy='auto'),
                    control.compute=list(dic=TRUE,cpo=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)
#Cx. nigripalpus
results[[7]]<- inla(as.formula(formulas[[7]]),
                    family="poisson",
                    data=inla.stack.data(stacks[[7]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks[[7]]), link = 1), #link = 1 scales the fitted values with the logit link
                    control.family=list(link='log'),# hyper=list(theta=prec.prior)),
                    control.inla=list(int.strategy='auto'),
                    control.compute=list(dic=TRUE,cpo=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)

names(results) <- splist

#Check the results
summary(results$`A. aeg`)
summary(results$`A. albo`)
summary(results$`Cx. salinarius`)
summary(results$`Cx. quinq`)
summary(results$`A. crucians`)
summary(results$`Cx. restuans`)
summary(results$`Cx. nigripalpus`)
