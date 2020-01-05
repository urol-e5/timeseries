#Title: Protein Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191104
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 

Data <- read.csv("RAnalysis/Data/Sym_Counts_Data.csv")
Meta <- read.csv("RAnalysis/Data/Sample_Metadata.csv")
Surface.Area <- read.csv("RAnalysis/Data/E5surface_area.csv") 
Data <-merge(Data, Meta, by="Tube.ID")
Data <- merge(Data, Surface.Area,  by="Sample.ID")

Data$avg.cells <- rowMeans(Data[,4:9])
Data$avg.cells[2] <- rowMeans(Data[2,4:7])
#Concentration = Number of cells x 10,000 / Number of squares
Data$cells.ml <- Data$avg.cells * 10000 /Data$Squares.Counted
Data$cells <- Data$cells.ml * Data$Homog.Vol.ml.x
Data$cells.cm2 <- Data$cells/Data$SA.cm2
Data$cells.106.cm2 <- Data$cells.cm2/10^6

Data$group <- paste0(Data$Species,"_",Data$Site)

myColors <- c("grey", "orange", "orange", "cyan", "cyan", "cyan", "green", "green", "green" )

pdf("RAnalysis/Output/Sym_Density.pdf")
#par(mfrow=c(1,2))
boxplot(Data$cells.106.cm2 ~ Data$group, las=2, ylab=expression(paste("Cells x 10^6 cm-"^"2")),col=myColors, xaxt='n', xlab="") 
legend("topleft", legend = c("Garden A.pulchra","A.pulchra", "Poc.meandrina", "Por.lutea" ) , 
       col = c("grey","orange", "cyan","green" ) , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))
dev.off()

jpeg("RAnalysis/Output/Sym_Density.jpg")
#par(mfrow=c(1,2))
boxplot(Data$cells.106.cm2 ~ Data$group, las=2, ylab=expression(paste("Cells x 10^6 cm-"^"2")),col=myColors, xaxt='n', xlab="") 
legend("topleft", legend = c("Garden A.pulchra","A.pulchra", "Poc.meandrina", "Por.lutea" ) , 
       col = c("grey","orange", "cyan","green" ) , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))
dev.off()
