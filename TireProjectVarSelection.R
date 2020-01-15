####---
# Title: TireProjectVarSelection.R
#
# Author: Mark Myer
#
# Date: 12/12/2019
#
# Purpose: To select predictor variables for mosquito larvae modeling in tires
#
# R version 3.6.1 Action of the Toes
####---

#Load packages
library(BayesVarSel)
library(caret)
library(corrgram)
library(corrplot)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
source("/Volumes/Mark Drive/Papers and Textbooks/Highstat Guide to INLA/HighstatLibV11.R") #Remember to cite these helper functions from Highstat

#Import data
tire_fulldat <- read.csv("TireData121219.csv") %>% dplyr::select(., -c(X, DispX, DispY, comments, optional)) %>% na.omit()

#Separate dataset into dataframes by species (maybe not necessary)
splist <- unique(tire_fulldat$MosqSpp)
aeg <- filter(tire_fulldat, MosqSpp == "A. aeg")
albo <- filter(tire_fulldat, MosqSpp == "A. albo")
sal <- filter(tire_fulldat, MosqSpp == "Cx. salinarius")
quin <- filter(tire_fulldat, MosqSpp == "Cx. quinq")
cru <- filter(tire_fulldat, MosqSpp == "A. crucians")
res <- filter(tire_fulldat, MosqSpp == "Cx. restuans")
nigr <- filter(tire_fulldat, MosqSpp == "Cx. nigripalpus")
dflist <- list(aeg, albo, sal, quin, cru, res, nigr) #List of dataframes for easier use of "apply" functions
names(dflist) = splist
saveRDS(dflist, file = "dflist.rds")
rm(aeg, albo, sal, quin, cru, res, nigr)

#Get correlations between every IV and the dependent variable for each species separately
dflist.num = list()
dflist.cor = list()
for (n in 1:length(splist)) {
  dflist.num[[n]] = select_if(dflist[[n]], is.numeric) %>% select_if(~ length(unique(.)) > 1) #Remove columns where all values are the same. These cannot affect the model due to zero variance
  dflist.cor[[n]] = cor(dflist.num[[n]])
}
names(dflist.num) = splist
names(dflist.cor) = splist

#Create a giant correlogram for each species
# for (n in 1:length(dflist.num)) {
#   tiff(filename=paste0("Figures/Correlograms/",names(dflist.num)[n],"corrplot.tiff"), width = 15, height = 15, units = "in", pointsize = 8, res = 96, compression = "lzw", type = "cairo")
#   corrplot(dflist.cor[[n]], method = "number", tl.cex = 1.75, number.cex = 0.75)
#   dev.off()
# }

#Obtain correlated variables (>0.5) to consider for removal 
dflist.remove <- list()
for (n in 1:length(dflist.cor)) {
dflist.remove[[n]] <- findCorrelation(cor(dplyr::select(dflist.num[[n]], -c(MosqPerL, WayPt_ID, LongX, LatY, Adj_X, Adj_Y, Day, Month, EpiWeek, INLAWeek, MosqCount))), cutoff = 0.5, names = TRUE)
}
names(dflist.remove) = splist

#Remove correlated variables, leaving one from each pair that is the least correlated with other variables in the dataset
dflist.reduced <- list()
for (n in 1:length(dflist.remove)) {
  dflist.reduced[[n]] <-dplyr::select(dflist.num[[n]], -c(dflist.remove[[n]])) #%>%
    #cbind(dplyr::select(dflist[[n]], c(veg_in_tire, cop, ostra, daph))) #Rebind factor variables
}
names(dflist.reduced) = splist

#Eliminate more variables using Variance Inflation Factors 
for(i in 1:length(dflist.reduced)) {
dplyr::select_if(dflist.reduced[[i]], is.numeric) %>% dplyr::select(-c(MosqPerL, WayPt_ID, LongX, LatY, Adj_X, Adj_Y, Day, Month, EpiWeek, INLAWeek, MosqCount)) %>% corvif()
}

#Eliminate variables that are "0" for more than 95% of outcomes - these are effectively outliers
colSums(tire_fulldat == 0) / nrow(tire_fulldat) 
#Oligo, tris, noto, dycitid are all 95% zero or more. 

for(i in 1:length(dflist.reduced)) {
  dflist.reduced[[i]] <- dplyr::select(dflist.reduced[[i]], -c(Developed.Medium.Intensity, Emergent.Herbaceous.Wetlands, oligo, noto, tris, dyticid))
}


#Conduct Bayesian variable selection to identify the probabilities of variable inclusion in the best models excluding spatial and temporal
bvs.list <- list()
bvs.vars <- list()
for (n in 1:length(dflist.reduced)) {
bvs.list[[n]] <- Bvs(MosqPerL ~ . , data = dplyr::select(dflist.reduced[[n]], -c(WayPt_ID, LongX, LatY, Adj_X, Adj_Y, Day, Month, EpiWeek, INLAWeek, MosqCount)), time.test=F)
bvs.vars[[n]] <- bvs.list[[n]]$inclprob[order(bvs.list[[n]]$inclprob, decreasing=T)]
}
names(bvs.list) = splist
names(bvs.vars) = splist

# #Remove the 'y' from the end of dummy-coded factor variables
# for (n in 1:length(bvs.vars)) {
#   names(bvs.vars[[n]]) <- sub("y$", "", names(bvs.vars[[n]]))
#   names(bvs.vars[[n]]) <- sub("Developed.Medium.Intensit", "Developed.Medium.Intensity", names(bvs.vars[[n]])) #Put the 'y' back on this one. Stupid, I know.
# }

#Reduce the datasets to the top 5 variables from each list and bind with the spatial/temporal ones held out
dat.selected <- list()
for (n in 1:length(bvs.vars)) {
  dat.selected[[n]] <- cbind(dplyr::select(dflist[[n]], c(MosqPerL, WayPt_ID, LongX, LatY, Adj_X, Adj_Y, Day, Month, EpiWeek, INLAWeek, MosqCount)), dplyr::select(dflist[[n]], names(bvs.vars[[n]][1:5])))
}
names(dat.selected) = splist

saveRDS(dat.selected, file = "Mosquito_Variables_Selected.rds")
