#Title: Biomass Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191031
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 

Data <- read.csv("RAnalysis/Data/Biomass_Data.csv")
Meta <- read.csv("RAnalysis/Data/Sample_Metadata.csv")
Surface.Area <- read.csv("RAnalysis/Data/E5surface_area.csv") 
Data <-merge(Data, Meta, by="Tube.ID")
Data <- merge(Data, Surface.Area,  by="Sample.ID")


Data$delta.mass.g <- Data$dry.pan.mass.g - Data$pan.mass.g
Data$mass.g.ml <- Data$delta.mass.g / Data$Vol.added.ml
Data$mass.g <- Data$mass.g.ml * Data$Homog.Vol.ml.x
Data$AFDW.g <- Data$dry.pan.mass.g - Data$burnt.pan.mass.g

Data$drymass.g.cm2 <- Data$mass.g / Data$SA.cm2
Data$drymass.mg.cm2 <- Data$drymass.g.cm2 * 1000

Data$AFDW.g.cm2 <- Data$AFDW.g / Data$SA.cm2
Data$AFDW.mg.cm2 <- Data$AFDW.g.cm2 * 1000


Data$group <- paste0(Data$Species,"_",Data$Site)

myColors <- c("grey", "orange", "orange", "cyan", "cyan", "cyan", "green", "green", "green" )

pdf("RAnalysis/Output/Biomass.pdf")
par(mfrow=c(1,2))
boxplot(Data$drymass.mg.cm2 ~ Data$group, las=2, ylab=expression(paste("Dry Biomass mg cm-"^"2")),col=myColors, xaxt='n', xlab="") 
legend("topleft", legend = c("Garden A.pulchra","A.pulchra", "Poc.meandrina", "Por.lutea" ) , 
       col = c("grey","orange", "cyan","green" ) , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))
boxplot(Data$AFDW.mg.cm2 ~ Data$group, las=2, ylab=expression(paste("AFDW mg cm-"^"2")),col=myColors, xaxt='n', xlab="") 
dev.off()

jpeg("RAnalysis/Output/Biomass.jpg")
par(mfrow=c(1,2))
boxplot(Data$drymass.mg.cm2 ~ Data$group, las=2, ylab=expression(paste("Dry Biomass mg cm-"^"2")),col=myColors, xaxt='n', xlab="") 
legend("topleft", legend = c("Garden A.pulchra","A.pulchra", "Poc.meandrina", "Por.lutea" ) , 
       col = c("grey","orange", "cyan","green" ) , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))
boxplot(Data$AFDW.mg.cm2 ~ Data$group, las=2, ylab=expression(paste("AFDW mg cm-"^"2")),col=myColors, xaxt='n', xlab="") 
dev.off()

