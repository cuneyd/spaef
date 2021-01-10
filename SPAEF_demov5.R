
#This is a R script to calculate spatial efficiency (SPAEF) which can be used to compare spatial patterns in two raster maps 
#Note that the "no values (here -9999)" in 2D maps are cleared first. Observed and Simulated two arrays (2D) are provided to the function.

#Detailed explanation goes here

#The newly proposed spatial efficiency metric (SPAEF) is proven to be robust 
#and easy to interpret due to its three distinct and complementary components of correlation, variance and histogram matching.

#Created on Wed Feb 6 2019
#@ authors:                 Mehmet Cüneyd Demirel
#@ author's website:        http://www.space.geus.dk/
#@ author's webpage:        http://akademi.itu.edu.tr/demirelmc/
#@ author's email id:       demirelmc@itu.edu.tr
#
#A libray with R functions for calculation of spatial efficiency (SPAEF) metric.
#
#Literature:
#
# [1] Demirel, M. C., Mai, J., Mendiguren, G., Koch, J., Samaniego, L., & Stisen, S. (2018). Combining satellite data and appropriate objective functions for improved spatial pattern performance of a distributed hydrologic model. Hydrology and Earth System Sciences, 22(2), 1299-1315. https://doi.org/10.5194/hess-22-1299-2018
# [2] Koch, J., Demirel, M. C., & Stisen, S. (2018). The SPAtial EFficiency metric (SPAEF): multiple-component evaluation of spatial patterns for optimization of hydrological models. Geoscientific Model Development, 11(5), 1873-1886. https://doi.org/10.5194/gmd-11-1873-2018

#clean Global Environment
rm(list=ls()) 

# WD: working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #For Rstudio users

###############automated setwd for all R users
script.dir <- getSrcDirectory(function(x) {x})
print(script.dir)
setwd(script.dir)
##############

# setwd("D:/SPAE_4public")


#install.packages("pracma") #histc is in pracma, uncomment and install first
library(pracma)

# read ascii
mask <- read.delim('./map_files/mask_1km.asc', header = FALSE, sep = ",", dec = ".")
dimens=dim(mask)
mask=array(as.numeric(unlist(mask)), dim=c(dimens[1], dimens[2]))

map1 <- read.delim('./map_files/obs.asc', header = FALSE, sep = "\t", dec = ".")  #observed data comes here
map1=array(as.numeric(unlist(map1)), dim=c(dimens[1], dimens[2]))
map1=map1[as.logical(mask)]

map2 <- read.delim('./map_files/sim_1.asc', header = FALSE, sep = "\t", dec = ".") #You can test other asc files here
map2=array(as.numeric(unlist(map2)), dim=c(dimens[1], dimens[2]))
map2=map2[as.logical(mask)]

obs <- map1[ map1 != -9999 ]
sim <- map2[ map2 != -9999 ]


#CORR
alpha=cor(obs,sim)

#coefficient of variation
cv_obs=sd(obs)/mean(obs);
cv_sim=sd(sim)/mean(sim);

beta=cv_sim/cv_obs;

#HISTOmatch
obs=(obs-mean(obs))/sd(obs)
sim=(sim-mean(sim))/sd(sim)

bins=floor(sqrt(length(obs)))

h1 <- hist(obs, breaks=bins, freq=TRUE, plot=TRUE)
h2 <- hist(sim, breaks=bins, freq=TRUE, plot=TRUE) #False makes Density instead of frequency, try it

a=histc(obs, h1$breaks)
b=histc(sim, h1$breaks)
c=cbind(a$cnt, b$cnt)
d <- pmin(c[,1],c[,2])
overlap=sum(d)
histogram_match=overlap/sum(a$cnt)
gamma=histogram_match

spaef = 1- sqrt( (alpha-1)^2 + (beta-1)^2 + (gamma-1)^2 )

#print(paste0("SPAEF: ", round(spaef, digits = 6)))
cat("SPAEF: ", spaef)