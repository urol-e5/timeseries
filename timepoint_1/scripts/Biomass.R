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

Data <- read.csv("data/1_Biomass.csv")
Meta <- read.csv("data/1_Sample_Info.csv")
Data <-merge(Data, Meta, by="colony_id")

Data$delta.mass.g <- Data$dry.pan.mass.g - Data$initial.mass.g

#different volumes for sym and host
#sym = 5ml
#host = 4ml
Data <- Data %>%
  mutate(delta.mass.g.ml = case_when(partner=="Host" ~ delta.mass.g/4, partner=="Sym" ~ delta.mass.g/5))

#Load homogenate volume
homog.vol <- read.csv("data/1_homogenate_vols.csv", header=TRUE)
Data <- merge(Data, homog.vol, by="colony_id")


Data$mass.g <- Data$delta.mass.g.ml * Data$homog_vol_ml
Data$AFDW.g <- Data$dry.pan.mass.g - Data$burnt.pan.mass.g

#calculate surface area standard curve
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
plot(smpls$surface.area.cm2)
range(smpls$surface.area.cm2)
range(stnds$surface.area.cm2)

Data <- merge(Data, smpls, by="colony_id")

Data$drymass.g.cm2 <- Data$mass.g / Data$surface.area.cm2
Data$drymass.mg.cm2 <- Data$drymass.g.cm2 * 1000

Data$AFDW.g.cm2 <- Data$AFDW.g / Data$surface.area.cm2
Data$AFDW.mg.cm2 <- Data$AFDW.g.cm2 * 1000

Data%>%
  group_by(Species, Site)%>%
  summarise(mean.value = mean(drymass.mg.cm2), se = std.error(drymass.mg.cm2)) %>%
  ggplot(aes(x = Site, y = mean.value, group = Species, color = Species))+
  ylab("Dry Biomass mg cm-2")+
  geom_point(size = 3)+
  geom_errorbar(aes(x = Site, ymin = mean.value-se, ymax = mean.value+se), width = 0.5)+
  facet_grid(~Species, scales = "free_y")
