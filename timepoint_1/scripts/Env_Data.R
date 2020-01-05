#Title: Environmental Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191028
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 

pH.1 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_pH_20697595.csv", header=F, skip=14)
pH.2 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_pH_20697596.csv", header=F, skip=14)
pH.3 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_pH_20697597.csv", header=F, skip=14)

tmp.1 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_TEMP_20444027.csv", header=F, skip=3)
tmp.2 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_TEMP_20444038.csv", header=F, skip=3)
tmp.3 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_TEMP_20444041.csv", header=F, skip=3)

sal.1 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_Cond_20715881.csv", header=T, skip=1)
sal.2 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_Cond_20715882.csv", header=F, skip=3)
sal.3 <- read.csv("RAnalysis/Data/Logger_Calibrations/Cross_Calibration/20191028_Cond_20715883.csv", header=F, skip=3)

