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

Data <- read.csv("data/1_surface_area.csv", header=TRUE)

Stnds <- read.csv("RAnalysis/Data/wax.standards_1.csv", header=TRUE)
Stnds$delta.mass.g <- Stnds$weight2.g-Stnds$weight1.g
#Stnds$area.cm2 <- ((pi*(Stnds$diameter.cm)^2) + (2*pi*(Stnds$diameter.cm)))
stnd.curve <- lm(surface.area.cm2~delta.mass.g, data=Stnds)
plot(surface.area.cm2~delta.mass.g, data=Stnds)

stnd.curve$coefficients

fecundity <- subset(Data, Sample.Type=="Fecundity")
fecundity$delta.mass.g <- fecundity$waxedmass.g-fecundity$mass.g
fecundity$SA.cm2 <- stnd.curve$coefficients[2] * fecundity$delta.mass.g + stnd.curve$coefficients[1]


E5.Data <- subset(Data, Sample.Type=="E5")
E5.Data$delta.mass.g <- E5.Data$waxedmass.g-E5.Data$mass.g

E5.Data$SA.cm2 <- stnd.curve$coefficients[2] * E5.Data$delta.mass.g + stnd.curve$coefficients[1]
write.csv(E5.Data, "RAnalysis/Data/E5surface_area.csv")

