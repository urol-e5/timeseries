#Title: Protein Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191104
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(tidyr)
library(dplyr)
library(plotrix)
library(ggplot2)

Licor.Data <- read.csv("data/1_Env_data/Cross_Calibration/20200111/1_20200111_Licor_Light_Calib.csv", header=FALSE)
Site1 <- read.csv("data/1_Env_data/Cross_Calibration/20200111/Logger1_FCDC4664E473_1578805794987.csv", header=FALSE, skip=81)
Site2 <- read.csv("data/1_Env_data/Cross_Calibration/20200111/Logger2_EA62F175DD93_1578805769597.csv", header=FALSE, skip=83)
Site3 <- read.csv("data/1_Env_data/Cross_Calibration/20200111/Logger3_ECB51D4FB2D6_1578805825453.csv", header=FALSE, skip=83)
Site1 <- Site1[1:49,]
Site2 <- Site2[1:49,]
Site3 <- Site3[1:49,]
Licor.Data$um.m2.s1 <- (Licor.Data$V2*10^4)/(5*60)

Data <- cbind(Licor.Data, Site1$V3,Site2$V3,Site3$V3)
colnames(Data) <- c("Date.Time", "Intg.Licor", "Inst.Licor", "Logger1", "Logger2", "Logger3")
  
S1 <- coef(lm(Logger1 ~ Inst.Licor, data=Data))
S2 <- coef(lm(Logger2 ~ Inst.Licor, data=Data))
S3 <- coef(lm(Logger3 ~ Inst.Licor, data=Data))

S1.Data <- read.csv("data/1_Env_data/1_FCDC4664E473_1578713594885_Site1_Light.csv")
S2.Data <- read.csv("data/1_Env_data/1_EA62F175DD93_1578713736764_Site2_Light.csv")
S3.Data <- read.csv("data/1_Env_data/1_ECB51D4FB2D6_1578713762658_Site3_Light.csv")

S1.Data$light <-(S1.Data$data1*S1[2])+S1[1]