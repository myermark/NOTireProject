####---
# Title: TireProjectModelRun.R
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
library(dplyr)
library(raster)
library(gstat)
library(tigris)


#Import data---- 
load(".Rdata")

#Define the INLA model ----
#Create mesh
meshes <- lapply(1:length(dat.selected), function (i) {
  loc <- list()
  bnd <- list()
  mesh <- list()
  loc[[i]] <- cbind(dat.selected[[i]]$Adj_X, dat.selected[[i]]$Adj_Y) 
  bnd[[i]] <-inla.nonconvex.hull(loc[[i]]) #Makes a nonconvex hull around the points
  mesh[[i]] <-inla.mesh.2d(boundary=bnd[[i]], cutoff=0.25, max.edge=c(1,4)) #Max-edge is based on estimated correlation distance from semivariograms
  return(mesh[[i]])
  }
)

#Make the SPDE
spdes <- lapply(1:length(meshes), function (i) {
  spde <- list()
  spde[[i]] <- inla.spde2.pcmatern(
  mesh=meshes[[i]], alpha=2, ### mesh and smoothness parameter
  prior.range=c(1, 0.05), ### P(range<1km)=0.05
  prior.sigma=c(sqrt(log(max(dat.selected[[i]]$MosqPerL+1))/3), 0.05) 
  # P(sigma > ) = 0.05
  # MosqPerL = exp(u_i) 
  # some u_i have to be as large as sigma_u to cover max(dat.selected$MosqPerL)
  # If u_i ~ N(0, sigma_u^2) then it is unlikely that sigma_u > log(max(dat.selected$MosqPerL))
  # So this takes the 99.7th quantile (3 standard deviations) of N(0, sigma_u^2) as the max value for the prior
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
    effects=list(list(spatial = 1:spdes[[i]]$n.spde), # data.frame(dat.selected[[i]])
                 cbind(dplyr::select(data.frame(dat.selected[[i]]), 1:(ncol(dat.selected[[i]]) - 5)),
                       (dplyr::select(data.frame(dat.selected[[i]]), (ncol(dat.selected[[i]]) - 4):ncol(dat.selected[[i]])) %>%
                         mutate_if(is.numeric, scale))
                       )
      )
    )
  return(stack[[i]])
  }
)

#Negative binomial stack, where we model the counts of positive outcomes
stacks.nb <- lapply(1:length(dat.selected), function (i) {
  stack <- list()
  stack[[i]] <- inla.stack(
    tag="Fit",
    data=list(y=ifelse(dat.selected[[i]]$MosqPerL == 0, NA, dat.selected[[i]]$MosqPerL)),
    A=list(As[[i]],1),
    effects=list(list(spatial = 1:spdes[[i]]$n.spde), #data.frame(dat.selected[[i]])
                 cbind(dplyr::select(data.frame(dat.selected[[i]]), 1:(ncol(dat.selected[[i]]) - 5)),
                       (dplyr::select(data.frame(dat.selected[[i]]), (ncol(dat.selected[[i]]) - 4):ncol(dat.selected[[i]])) %>%
                          mutate_if(is.numeric, scale))
                 )
    )
  )
  return(stack[[i]])
}
)

#Define the INLA formulas
formulas <- data.frame(nonspatial = rep(NA, times = 7), randint = rep(NA, times = 7), spatial = rep(NA, times = 7), temporal = rep(NA, times = 7), spatiotemporal = rep(NA, times = 7))
for (i in 1:length(dat.selected)) {
  len <- length(dat.selected[[i]])
  formulas$nonspatial[i] <- paste("y ~ -1 + ",
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
  formulas$randint[i]<- paste("y ~ -1 + ",
                                   names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                                   "+",
                                   names(dat.selected[[i]])[len-3], 
                                   "+",
                                   names(dat.selected[[i]])[len-2], 
                                   "+",
                                   names(dat.selected[[i]])[len-1], 
                                   "+",
                                   names(dat.selected[[i]])[len], 
                                   "+ f(WayPt_ID, model = 'iid')"#, #Random intercept by waypoint site
  )
  formulas$spatial[i] <- paste("y ~ -1 + ",
                                   names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                                   "+",
                                   names(dat.selected[[i]])[len-3], 
                                   "+",
                                   names(dat.selected[[i]])[len-2], 
                                   "+",
                                   names(dat.selected[[i]])[len-1], 
                                   "+",
                                   names(dat.selected[[i]])[len], 
                                   "+ f(spatial, model=spdes[[", i, "]])"#, #Spatial model component
  )
  formulas$temporal[i]<- paste("y ~ -1 + ",
                                   names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                                   "+",
                                   names(dat.selected[[i]])[len-3], 
                                   "+",
                                   names(dat.selected[[i]])[len-2], 
                                   "+",
                                   names(dat.selected[[i]])[len-1], 
                                   "+",
                                   names(dat.selected[[i]])[len], 
                                   "+ f(INLAWeek, model = 'ar1', hyper = list(theta1=list(prior='pc.prec', param=c(0.5,0.5)), theta2=list(prior='pc.cor1', param=c(0.9,0.9))))" #Temporal model component
  )
  formulas$spatiotemporal[i] <- paste("y ~ -1 + ",
                                   names(dat.selected[[i]])[len-4], #This pastes the last 5 variable names together into the formula
                                   "+",
                                   names(dat.selected[[i]])[len-3], 
                                   "+",
                                   names(dat.selected[[i]])[len-2], 
                                   "+",
                                   names(dat.selected[[i]])[len-1], 
                                   "+",
                                   names(dat.selected[[i]])[len], 
                                   "+ f(spatial, model=spdes[[", i, "]])", #Spatial model component
                                   "+ f(INLAWeek, model = 'ar1', hyper = list(theta1=list(prior='pc.prec', param=c(0.5,0.5)), theta2=list(prior='pc.cor1', param=c(0.9,0.9))))" #Temporal model component
  )
}
#Name the models
modlist <- c("Nonspatial", "Randint", "Spatial", "Temporal", "Spatiotemporal")

#Define prior for precision
prec.prior <- list(prior='pc.prec', param=c(0.5, 0.5))

#Run the binomial part of the models
results.bin<- list()
for(i in 1:length(dat.selected)) {
  mod.templist <- list()
  for(j in 1:length(modlist)) {
     temp <- try(
              inla(as.formula(formulas[i,j]),
                    family="binomial",
                    Ntrials = 1,
                    data=inla.stack.data(stacks.bin[[i]]),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(stacks.bin[[i]]), link = 1), 
                    control.inla=list(int.strategy='auto', correct = TRUE, correct.factor = 10),
                    control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE,po=TRUE,config=TRUE),
                    control.fixed=list(expand.factor.strategy ='inla'),
                    verbose = F)
              )
              if(inherits(temp, "try-error"))
              {
                mod.templist[[j]] <- NA
                next
              }
     mod.templist[[j]] <- temp
  }
  names(mod.templist) = modlist
  results.bin[[i]] <- mod.templist
}
names(results.bin) <- splist

#Run the negative binomial part of the models
results.nb<- list()
for(i in 1:length(dat.selected)) {
  mod.templist <- list()
  for(j in 1:length(modlist)) {
    temp <- try(
              inla(as.formula(formulas[i,j]),
                           family="zeroinflatednbinomial0",
                           data=inla.stack.data(stacks.nb[[i]]),
                           control.predictor=list(compute=TRUE, A=inla.stack.A(stacks.nb[[i]]), link = 1), 
                           control.inla=list(int.strategy='auto', correct = TRUE, correct.factor = 10),
                           control.family=list(hyper = list(prec.prior)), 
                           control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE,po=TRUE,config=TRUE),
                           control.fixed=list(expand.factor.strategy ='inla'),
                           verbose = F)
            )
            if(inherits(temp, "try-error"))
            {
              mod.templist[[j]] <- NA
              next
            }
            mod.templist[[j]] <- temp
  }
  names(mod.templist) = modlist
  results.nb[[i]] <- mod.templist
}
names(results.nb) <- splist

#Clean up environment
rm(i,j,len,mod,temp,mod.templist)
#Save image for model evaluation
save.image()

#TESTING AREA-----
i.test <- 2
#Create mesh
  loc.test <- cbind(dat.selected[[i.test]]$Adj_X, dat.selected[[i.test]]$Adj_Y) 
  bnd.test <-inla.nonconvex.hull(loc.test) #Makes a nonconvex hull around the points
  mesh.test <-inla.mesh.2d(boundary=bnd.test, cutoff=0.25, max.edge=c(1,4)) 
  
#Make the SPDE
count.test = 0
test.range = list()
for(i in seq(0.5, 2.5, by=0.25)) {
  for(j in seq(0.5, 2.5, by = 0.25)) {
count.test = count.test + 1
  spde.test <- inla.spde2.pcmatern(
    mesh=mesh.test, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.5, 0.05), ### P(range<0.5km)=0.05
    prior.sigma=c(7, 0.05) ### P(sigma>0.5)=0.05
  )

#Create the projector matrix 
A.test <- inla.spde.make.A(mesh.test, loc.test)

stack.test <- inla.stack(
  tag="Fit",
  data=list(y=ifelse(dat.selected[[i.test]]$MosqPerL == 0, NA, dat.selected[[i.test]]$MosqPerL)),
  A=list(A.test,1),
  effects=list(list(spatial = 1:spde.test$n.spde), #data.frame(dat.selected[[i.test]])
               cbind(dplyr::select(data.frame(dat.selected[[i.test]]), 1:(ncol(dat.selected[[i.test]]) - 5)),
                     (dplyr::select(data.frame(dat.selected[[i.test]]), (ncol(dat.selected[[i.test]]) - 4):ncol(dat.selected[[i.test]])) %>%
                        mutate_if(is.numeric, scale))
               )
  )
)

len.test <- length(dat.selected[[i.test]])
formula.test <- paste("y ~ -1 + ",
                             names(dat.selected[[i.test]])[len.test-4], #This pastes the last 5 variable names together into the formula
                             "+",
                             names(dat.selected[[i.test]])[len.test-3], 
                             "+",
                             names(dat.selected[[i.test]])[len.test-2], 
                             "+",
                             names(dat.selected[[i.test]])[len.test-1], 
                             "+",
                             names(dat.selected[[i.test]])[len.test], 
                             "+ f(spatial, model=spde.test)"#, #Spatial model component
                             #"+ f(WayPt_ID, model = 'iid')"#, #Random intercept by waypoint site
                             #"+ f(INLAWeek, model = 'ar1', hyper = list(theta1=list(prior='pc.prec', param=c(0.5,0.5)), theta2=list(prior='pc.cor1', param=c(0.9,0.9))))" #Temporal model component
)


inla.test <- inla(as.formula(formula.test),
                  family="zeroinflated.nbinomial0",
                  data=inla.stack.data(stack.test),
                  control.predictor=list(compute=TRUE, A=inla.stack.A(stack.test), link = 1), 
                  control.inla=list(int.strategy='auto', correct = TRUE, correct.factor = 10),
                  control.family=list(hyper = list(theta = prec.prior)), 
                  control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE,po=TRUE,config=TRUE),
                  control.fixed=list(expand.factor.strategy ='inla'),
                  verbose = F)

summary(inla.test)

test.range[[count.test]]<- c(inla.test$summary.hyperpar$mean[2], inla.test$waic$waic)
  }
}

w.test <- inla.test$summary.random$spatial$mean
PlotField2(field = w.test, 
           mesh = mesh.test, 
           xlim = range(mesh.test$loc[,1]), 
           ylim = range(mesh.test$loc[,2]),
           MyMain = paste(names(dat.selected)[i.test], "Test Model (Negative Binomial)")
)
axis(1); axis(2)
points(x = km_loc[,1],
       y = km_loc[,2], 
       cex = 0.5, 
       col = "black", 
       pch = 16)
plot(nola_roads, add=T, lwd = 0.75)
plot(border_km, add=T)

df.test<- data.frame(nonspatial = seq(1:7), randint = seq(8:14), spatial = seq(1:7), temporal = seq(1:7), spatiotemporal = seq(1:7))

rm(list = ls(pattern = ".test$"))

