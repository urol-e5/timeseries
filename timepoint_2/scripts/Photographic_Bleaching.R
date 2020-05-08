#Title: Bleaching estimates from pictures
#Project: NSF BSF
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20181206
#See Readme file for details

rm(list=ls()) #clears workspace 

if ("vegan" %in% rownames(installed.packages()) == 'FALSE') install.packages('vegan') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggpubr') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("multcompView" %in% rownames(installed.packages()) == 'FALSE') install.packages('multcompView') 


#Read in required libraries
##### Include Versions of libraries
library("vegan")
library("ggpubr")
library("gridExtra")
library("plyr") 
library("lsmeans")
library("multcompView")


#Required Data files

# Set Working Directory:
setwd("~/MyProjects/Holobiont_Integration/RAnalysis/") #set working

Data <- read.csv("Data/Bleaching ImageJ_20200211.csv", header=T, sep=",", na.string="NA") #read in data file
Sample.Info <- read.csv("Data/Master_Fragment_Sheet.csv", header=T, sep=",", na.string="NA") #read in data file 
Sample.Info <- Sample.Info[,c(7,5,9)]
Tank.Info <- read.csv("Data/Tank_to_Treatment.csv", header=T, sep=",", na.string="NA") #read in data file 


A <- Sample.Info$PLUG.ID
B <- Data$PLUG.ID
bad.id.tanks <- which(B %in% A ==FALSE)
bad.id.tanks


#Data <- merge(Data, Sample.Info, by="PLUG.ID")
Data <-na.omit(Data)
#Data <- merge(Data, Tank.Info, by="Tank")

Data$Red.Norm.Coral <- Data$Red.Coral/Data$Red.Standard #normalize to color standard
Data$Green.Norm.Coral <- Data$Green.Coral/Data$Green.Standard #normalize to color standard
Data$Blue.Norm.Coral <- Data$Blue.Coral/Data$Blue.Standard #normalize to color standard

# xxx <- subset(Data, Timepoint=="Time5")
# xxx <- subset(xxx, Species.x=="Pocillopora")
# xxx <- subset(xxx, Treatment.x=="ATHC")
# par(mfrow=c(1,3))
# plot(xxx$Red.Norm.Coral ~ xxx$Tank)
# plot(xxx$Green.Norm.Coral ~ xxx$Tank)
# plot(xxx$Blue.Norm.Coral ~ xxx$Tank)


#Red.Norm.Coral <- Data$Red.Coral/Data$Red.Standard #normalize to color standard
#Green.Norm.Coral <- Data$Green.Coral/Data$Green.Standard #normalize to color standard
#Blue.Norm.Coral <- Data$Blue.Coral/Data$Blue.Standard #normalize to color standard

blch.scor <- as.matrix(cbind(Data$Red.Norm.Coral,Data$Green.Norm.Coral,Data$Blue.Norm.Coral)) #create matrix
rownames(blch.scor) <- Data$PLUG.ID #name columns in dataframe

dist <- vegdist(blch.scor, method="euclidean") #calculate distance matrix of color scores

PCA.color <- princomp(dist) #run principal components Analysis
summary(PCA.color) # view variance explained by PCs

Blch <- as.data.frame(PCA.color$scores[,1]) #extract PC1
Blch$PLUG.ID <- rownames(blch.scor)
#Blch <- merge(Blch, Data, by="PLUG.ID")
Blch  <- cbind(Blch, Data$Timepoint, Data$Treatment, Data$Species) #make a dataframe of PC1 and experiment factors
colnames(Blch) <- c("Bleaching.Score", "PLUG.ID", "Timepoint", "Treatment", "Species")
Blch$Group <- paste(Blch$Timepoint, Blch$Treatment, Blch$Species)
Blch$SpGroup <- paste(Blch$Treatment, Blch$Species)

# write.table(Blch,"~/MyProjects/Holobiont_Integration/RAnalysis/Output/Bleaching_Score.csv",sep=",", row.names=FALSE)

write.table(Blch, file = "~/MyProjects/Holobiont_Integration/RAnalysis/Output/Bleaching_Score.csv", append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

x<- Blch #assign values to test data frame to look for outliers
par(mar=c(10,4,2,2)) #bottom, left, top and right margins respectively
boxplot(Bleaching.Score ~ Group, data = x, lwd = 1, ylab = 'PC1Color', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint

## Sample Sizes 
Bleaching.Score <- read.csv("Output/Bleaching_Score.csv", header=T, sep=",", na.string="NA") #read in data file


##### Repeat investigation until min >-25
#min <- which(x$Bleaching.Score==min(x$Bleaching.Score))
#x[min,]
#x<- x[-min,]
#####

pdf("~/MyProjects/Holobiont_Integration/RAnalysis/Output/Photographic_Bleaching.pdf")
par(mar=c(10,4,2,2)) #bottom, left, top and right margins respectively
boxplot(Bleaching.Score ~ Group, data = Blch, lwd = 1, ylab = 'PC1Color', las=2) #plot boxplot of PC1 color score by Genotype and timepoint
stripchart(Bleaching.Score ~ Group, vertical = TRUE, data = Blch, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue', cex=0.2) #include all datapoints in blue overlaid on boxplots
text(x= 0.5, y= min(Blch$Bleaching.Score), labels= "pale") #add text to indicate dark and pale on graphic
text(x= 0.5, y= max(Blch$Bleaching.Score)+2, labels= "dark") #add text to indicate dark and pale on graphic
dev.off()

mod1 <- aov(sqrt(Bleaching.Score+200) ~ Treatment*Timepoint*Species, data=Blch) #run an ANOVA by Genotype
hist(residuals(mod1)) #look at normality of data
boxplot(residuals(mod1)) #look at normality of data
summary(mod1)


marginal = lsmeans(mod1, ~ Treatment*Timepoint*Species)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")

CLD <- CLD[order(CLD$Timepoint, CLD$Treatment, CLD$Species),]
CLD

All.Means <- ddply(Blch, c('Timepoint','Species', 'Treatment'), summarize,
                    mean= mean(Bleaching.Score, na.rm=T), #mean 
                    N = sum(!is.na(Bleaching.Score)), # sample size
                    se = sd(Bleaching.Score, na.rm=T)/sqrt(N)) #SE
All.Means
All.Means$se[is.na(All.Means$se)] <- 0
All.Means$Group <- paste(All.Means$Timepoint, All.Means$Treatment, All.Means$Species)
All.Means$SpGroup <- paste(All.Means$Treatment, All.Means$Species)

cols <- c("lightblue", "blue", "salmon", "red3")
#All.Means$Timepoint <- factor(All.Means$Timepoint, levels = c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7", "Week8", "Week9", "Week10", "Week11", "Week12", "Week13", "Week14", "Week15", "Week16"))
All.Means$Timepoint <- factor(All.Means$Timepoint, levels = c("Time1", "Time2", "Time3", "Time4", "Time5", "Time6", "Time7", "Time8", "Time9", "Time10", "Time11", "Time12", "Time13", "Time14", "Time15", "Time16"))

rec <- data.frame(x1=c(1,3,1,5,4))


Fig.All <- ggplot(All.Means, aes(x=Timepoint, y=mean, group=SpGroup)) + 
  geom_line(aes(linetype= Species, colour=Treatment, group=SpGroup), position = position_dodge(width = 0.1), alpha=0.5) + # colour, group both depend on cond2
  geom_errorbar(aes(ymin=All.Means$mean-All.Means$se, ymax=All.Means$mean+All.Means$se), colour="black", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Species), size = 1, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Timepoint") +
  ylab(expression(paste("Bleaching Score"))) +
  ylim(-110,30) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position= "right") + #remove legend background
  ggtitle("Photographic Bleaching") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 10, 
                                  hjust = 0))
Fig.All



# paling <- ggline(Blch, x = "Timepoint", y = "Bleaching.Score", color = "Group", group="Species",
#                   las=2,
#                   add = c("mean_se"),
#                   palette = c("lightblue", "lightblue",  "blue",  "blue", "pink", "pink", "red","red"))
# paling 


Mcap.Blch <- subset(Blch, Species=="Montipora")
Pact.Blch <- subset(Blch, Species=="Pocillopora")

Mcap.Means <- ddply(Mcap.Blch, c('Timepoint','Species', 'Treatment'), summarize,
                  mean= mean(Bleaching.Score, na.rm=T), #mean 
                  N = sum(!is.na(Bleaching.Score)), # sample size
                  se = sd(Bleaching.Score, na.rm=T)/sqrt(N)) #SE
Mcap.Means

Bleaching_mont.aov <- aov(mean ~ Treatment, data = Mcap.Means)
summary(Bleaching_mont.aov)
TukeyHSD(Bleaching_mont.aov)

Mcap.Means$Timepoint <- factor(Mcap.Means$Timepoint, levels = c("Time1", "Time2", "Time3", "Time4", "Time5", "Time6", "Time7", "Time8", "Time9", "Time10", "Time11", "Time12", "Time13", "Time14", "Time15", "Time16", "Time17"))

pdf("~/MyProjects/Holobiont_Integration/RAnalysis/Output/Montipora_Bleaching.pdf")
Fig.MC <- ggplot(Mcap.Means, aes(x=Timepoint, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=Mcap.Means$mean-Mcap.Means$se, ymax=Mcap.Means$mean+Mcap.Means$se), colour="black", width=.1, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  geom_line(aes(colour=Treatment, group=Treatment), position = position_dodge(width = 0.2), alpha=0.4) + # colour, group both depend on cond2
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Weeks") +
  ylab(expression(paste("Bleaching Score"))) +
  ylim(-20,31) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("M. capitata") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  geom_vline(xintercept = c(8.45), colour="red", linetype="dotted") + 
  geom_vline(xintercept = c(8.55), colour="blue", linetype="dotted") + 
  geom_vline(xintercept = c(16.55), colour="red", linetype="dotted") + 
  geom_vline(xintercept = c(16.45), colour="blue", linetype="dotted") +
  geom_vline(xintercept = c(0.5), colour="red", linetype="dotted") 
  # geom_text(aes(x=1, y=-30, label="2-month Stress"), size=2.5, colour="red", angle=90, vjust=-1, hjust=0) +
  # geom_text(aes(x=9, y=-30, label="2-month Recovery"), size=2.5, colour="blue", angle=90, vjust=-1, hjust=0) +
  # geom_text(aes(x=17, y=-30, label="Recurrent Stress"), size=2.5, colour="red", angle=90, vjust=-1, hjust=0) 
dev.off()

Fig.MC


Pact.Means <- ddply(Pact.Blch, c('Timepoint','Species', 'Treatment'), summarize,
                    mean= mean(Bleaching.Score, na.rm=T), #mean pnet
                    N = sum(!is.na(Bleaching.Score)), # sample size
                    se = sd(Bleaching.Score, na.rm=T)/sqrt(N)) #SE
Pact.Means

Bleaching_poc.aov <- aov(mean ~ Treatment, data = Pact.Means)
summary(Bleaching_poc.aov)
TukeyHSD(Bleaching_poc.aov)

Pact.Means$Timepoint <- factor(Pact.Means$Timepoint, levels = c("Time1", "Time2", "Time3", "Time4", "Time5", "Time6", "Time7", "Time8", "Time9", "Time10", "Time11", "Time12", "Time13", "Time14", "Time15", "Time16"))

Fig.PA <- ggplot(Pact.Means, aes(x=Timepoint, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=Pact.Means$mean-Pact.Means$se, ymax=Pact.Means$mean+Pact.Means$se), colour="black", width=.1, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  geom_line(aes(colour=Treatment, group=Treatment), position = position_dodge(width = 0.2), alpha=0.4) + # colour, group both depend on cond2
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Weeks") +
  ylab(expression(paste("Bleaching Score"))) +
  ylim(-20,100) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("P. acuta") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) + 
  geom_vline(xintercept = c(8.45), colour="red", linetype="dotted") + 
  geom_vline(xintercept = c(8.55), colour="blue", linetype="dotted") + 
  geom_vline(xintercept = c(16.45), colour="blue", linetype="dotted") +
  geom_vline(xintercept = c(0.5), colour="red", linetype="dotted") 
Fig.PA



Mcap.bch <- ggline(Mcap.Blch, x = "Timepoint", y = "Bleaching.Score", color = "Treatment",
                  ylim=c(-60,10),
                  add = c("mean_se", "jitter"),
                  title = "M. capitata",
                  palette = c("lightblue", "darkblue","pink", "red"))

Pact.bch <- ggline(Pact.Blch, x = "Timepoint", y = "Bleaching.Score", color = "Treatment",
                   ylim=c(-60,10),
                   add = c("mean_se", "jitter"),
                   title = "P. acuta",
                   palette = c("lightblue", "darkblue","pink", "red"))

Bch.Figs <- arrangeGrob(Fig.All, ncol=1)
ggsave(file="Output/Photo_Bch.pdf", Bch.Figs, width = 4, height = 3, units = c("in"))


## Color Distance: https://peerj.com/articles/6398/
if ("colordistance" %in% rownames(installed.packages()) == 'FALSE') install.packages('colordistance') 


## Visualizing Data (ES)
if ("hrbrthemes" %in% rownames(installed.packages()) == 'FALSE') install.packages('hrbrthemes') 
if ("GGally" %in% rownames(installed.packages()) == 'FALSE') install.packages('GGally') 
if ("viridis" %in% rownames(installed.packages()) == 'FALSE') install.packages('viridis') 
if ("reshape" %in% rownames(installed.packages()) == 'FALSE') install.packages('reshape') 

library(hrbrthemes)
library(GGally)
library(viridis)
library(dplyr)
library(tidyr)
library(reshape)

iris

Montipora.Blch <- Blch %>% dplyr::select("Timepoint", "Bleaching.Score", "PLUG.ID", "Treatment", "Species") %>% filter(Species %in% "Montipora") 
Montipora.NA <- reshape(Montipora.Blch, direction="wide", idvar=c("PLUG.ID", "Treatment", "Species"), timevar = "Timepoint")
Montipora <- na.omit(Montipora.NA) # ommit NAs

# produced error: couldn't find pivot_wider function Montipora.Blch.wide <- pivot_wider(Montipora.Blch, names_from=Timepoint, values_from=Bleaching.Score)
# this goes from wide to long: Montipora <- melt(Montipora.Blch, id=c("Timepoint", "Bleaching.Score"))

ggparcoord(Montipora,
           columns = 4:19, groupColumn = 2, order = "anyClass",
           showPoints = TRUE, scale = "globalminmax",
           alphaLines = 0.3
          ) + 
  scale_color_viridis(discrete=TRUE) + theme_classic()

# using Mcap.means and Pact.means data frames to plot means 


