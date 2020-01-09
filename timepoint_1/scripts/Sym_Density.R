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

Data <- read.csv("data/1_Sym_Counts_Data.csv")
Meta <- read.csv("data/1_Sample_Info.csv")
Data <-merge(Data, Meta, by="Plug.Number")

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

range(smpls$surface.area.cm2)
range(stnds$surface.area.cm2)

Data <- merge(Data, smpls, by="Plug.Number")

#Load homogenate volume
homog.vol <- read.csv("data/1_homogenate_vols.csv", header=TRUE)

Data <- merge(Data, homog.vol, by="Plug.Number")

Data$avg.cells <- rowMeans(Data[,3:8])
#Data$avg.cells[2] <- rowMeans(Data[2,4:7])
#Concentration = Number of cells x 10,000 / Number of squares
Data$cells.ml <- Data$avg.cells * 10000 /Data$Squares.Counted
Data$cells <- Data$cells.ml * Data$homog_vol_ml
Data$cells.cm2 <- Data$cells/Data$surface.area.cm2
Data$cells.106.cm2 <- Data$cells.cm2/10^6
Data <- subset(Data, Species!="BK")
plot(Data$cells.106.cm2 ~Data$Plug.Number, las=2)


Data%>%
  group_by(Species, Site)%>%
  summarise(mean.value = mean(cells.106.cm2), se = std.error(cells.106.cm2)) %>%
  ggplot(aes(x = Site, y = mean.value, group = Species, color = Species))+
  geom_point(size = 3)+
  geom_errorbar(aes(x = Site, ymin = mean.value-se, ymax = mean.value+se), width = 0.5)+
  facet_grid(~Species, scales = "free_y")


