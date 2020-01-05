#Title: Protein Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191031
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 

Host.Data <- read.csv("RAnalysis/Data/20191030_host_total_protein.csv")
Host.Map <- read.csv("RAnalysis/Data/20191030_host_total_protein_WELLS.csv")
Meta <- read.csv("RAnalysis/Data/Sample_Metadata.csv")
Surface.Area <- read.csv("RAnalysis/Data/E5surface_area.csv") 

Host<-merge(Host.Data, Host.Map, by="Well")
Stnds <- subset(Host, Sample.Type=="Standard")
Host <-merge(Host, Meta, by="Tube.ID")
Host <- merge(Host, Surface.Area,  by="Sample.ID")

Stnds.eq <- lm(A562~Concentration, data=Stnds) 
plot(A562~Concentration, data=Stnds)

coefficients(Stnds.eq)

Host$protein.ug.ml <- coefficients(Stnds.eq)[2]*Host$A562 + coefficients(Stnds.eq)[1]
Host$protein.ug <- Host$protein.ug.ml * Host$Homog.Vol.ml.x
#added 0.725ml of base and acid to 1ml of homogenate for solubilization for 1.37931 dilution factor
Host$protein.ug <- Host$protein.ug * 1/0.725
Host$protein.ug.cm2 <- Host$protein.ug/Host$SA.cm2
Host$group <- paste0(Host$Species,"_",Host$Site)
Host.avgs <- aggregate(protein.ug.cm2~ Sample.ID, data=Host, FUN=mean)
Host.avgs <- merge(Host.avgs, Meta, by="Sample.ID")
Host.avgs$group <- paste0(Host.avgs$Species,"_",Host.avgs$Site)
boxplot(Host.avgs$protein.ug.cm2 ~ Host.avgs$group, las=2)  

# HOLOBIONT DATA
Holo.Data <- read.csv("RAnalysis/Data/20191031_holobiont_total_protein.csv")
Holo.Map <- read.csv("RAnalysis/Data/20191031_holobiont_total_protein_WELLS.csv")

Holo<-merge(Holo.Data, Holo.Map, by="Well")
Stnds <- subset(Holo, Sample.Type=="Standard")
Holo <-merge(Holo, Meta, by="Tube.ID")
Holo <- merge(Holo, Surface.Area,  by="Sample.ID")

Stnds.eq <- lm(A562~Concentration, data=Stnds) 
plot(A562~Concentration, data=Stnds)

coefficients(Stnds.eq)

Holo$protein.ug.ml <- coefficients(Stnds.eq)[2]*Holo$A562 + coefficients(Stnds.eq)[1]
Holo$protein.ug <- Holo$protein.ug.ml * Holo$Homog.Vol.ml.x
#added 0.420ml of base and acid to 0.5ml of homogenate for solubilization for 1.190476 dilution factor
Holo$protein.ug <- Holo$protein.ug * 0.5/0.420
Holo$protein.ug.cm2 <- Holo$protein.ug/Holo$SA.cm2
Holo$group <- paste0(Holo$Species,"_",Holo$Site)
Holo.avgs <- aggregate(protein.ug.cm2~ Sample.ID, data=Holo, FUN=mean)
Holo.avgs <- merge(Holo.avgs, Meta, by="Sample.ID")
Holo.avgs$group <- paste0(Holo.avgs$Species,"_",Holo.avgs$Site)
boxplot(Holo.avgs$protein.ug.cm2 ~ Holo.avgs$group, las=2)  


myColors <- c("grey", "orange", "orange", "cyan", "cyan", "cyan", "green", "green", "green" )

pdf("RAnalysis/Output/Protein.pdf")
par(mfrow=c(1,2))
boxplot(Host.avgs$protein.ug.cm2 ~ Host.avgs$group, las=2, ylab=expression(paste("Host Total Protein µg cm-"^"2")), ylim=c(0,0.65),col=myColors, xaxt='n', xlab="") 
legend("topleft", legend = c("Garden A.pulchra","A.pulchra", "Poc.meandrina", "Por.lutea" ) , 
       col = c("grey","orange", "cyan","green" ) , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))

boxplot(Holo.avgs$protein.ug.cm2 ~ Holo.avgs$group, las=2, ylab=expression(paste("Holobiont Total Protein µg cm-"^"2")), ylim=c(0,0.65),col=myColors, xaxt='n', xlab="") 
dev.off()

jpeg("RAnalysis/Output/Protein.jpg")
par(mfrow=c(1,2))
boxplot(Host.avgs$protein.ug.cm2 ~ Host.avgs$group, las=2, ylab=expression(paste("Host Total Protein µg cm-"^"2")), ylim=c(0,0.65),col=myColors, xaxt='n', xlab="") 
legend("topleft", legend = c("Garden A.pulchra","A.pulchra", "Poc.meandrina", "Por.lutea" ) , 
       col = c("grey","orange", "cyan","green" ) , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))

boxplot(Holo.avgs$protein.ug.cm2 ~ Holo.avgs$group, las=2, ylab=expression(paste("Holobiont Total Protein µg cm-"^"2")), ylim=c(0,0.65),col=myColors, xaxt='n', xlab="") 
dev.off()
