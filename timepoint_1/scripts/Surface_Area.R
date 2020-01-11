#Title: Surface Area By Wax Dipping 
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191027
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 


#Surface area by wax dipping Veal et al 

library(tidyverse)

wax.data <- read.csv("data/1_Wax_dipping.csv", header=TRUE)

wax.data$delta.mass.g <- wax.data$weight2.g-wax.data$weight1.g
stnds <- subset(wax.data, Sample=="Standard")
stnds <- stnds[-1,]
stnds$rad <- stnds$Diameter/2
stnds$surface.area.cm2 <- 4*pi*(stnds$rad)^2
stnd.curve <- lm(surface.area.cm2~delta.mass.g, data=stnds)
plot(surface.area.cm2~delta.mass.g, data=stnds)


stnd.curve$coefficients
summary(stnd.curve)$r.squared

#Calculate surface area
smpls <- subset(wax.data, Sample=="Coral")
smpls$surface.area.cm2 <- stnd.curve$coefficients[2] * smpls$delta.mass.g + stnd.curve$coefficients[1]

smpls <- smpls %>%
  select(-Sample, -Diameter)


write.csv(smpls, "output/coral_surface_area.csv")

