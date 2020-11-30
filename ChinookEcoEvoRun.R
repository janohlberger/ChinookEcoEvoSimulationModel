##########################################################################
##                                                                      ##
##   Run model to simulate trends in Chinook salmon age-size structure  ##
##                                                                      ##
##########################################################################
pkgs<-c("plyr","msm")
if(length(setdiff(pkgs,rownames(installed.packages())))>0) { install.packages(setdiff(pkgs,rownames(installed.packages()))) }
lapply(pkgs,library,character.only=T)
source("ChinookEcoEvoParameters.R")
source("ChinookEcoEvoCode.R")
ChinookEcoEvoModel(parameters) 
# save(list=ls(),file="Rworkspace.RData")
##########################################################################