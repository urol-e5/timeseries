#Title: Photosynthesis Irradiance Curves
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191024
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
#if ("phytotools" %in% rownames(installed.packages()) == 'FALSE') install.packages('phytotools') 

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
#library('phytotools')


##### PI Curve Rate Calculation #####
path.p<-"RAnalysis/Data/All_Resp/" #the location of all your respirometry files 

#bring in the files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*12, ncol=4)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "µmol.L.sec", "Temp")

#Load Sample meta info Info
Sample.Info <- read.csv(file="RAnalysis/Data/All_PI_Curve_Sample_Info.csv", header=T) #read sample.info data
Sample.Info$Fragment.ID <- paste(Sample.Info$Fragment.ID,"_",Sample.Info$Light_Level, sep = "")

#subset the data by light step using time breaks in the data
for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]),  header=T, sep=",", skip=1, na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- as.data.frame(cbind(Photo.Data1$Time, Photo.Data1$Value, Photo.Data1$Temp)) #subset columns of interest
  colnames(Photo.Data1) <- c("Time", "Value", "Temp")
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", "HST") -7*60*60 #convert time from character to time
  Photo.Data1$Time <-format(strptime(Photo.Data1$Time, format="%Y-%m-%d %H:%M:%S"), "%H:%M:%S")
  levs <- which(Sample.Info$Sample.ID ==sub("_.*", "", file.names[i]))
  brk <- unique(as.POSIXct(Sample.Info$Start.time[levs],format="%H:%M","HST"))
  brk <-format(strptime(brk, format="%Y-%m-%d %H:%M:%S"), "%H:%M:%S")
  Step1 <- subset(Photo.Data1, Photo.Data1$Time > brk[1] & Photo.Data1$Time < brk[2])
  Step2 <- subset(Photo.Data1, Photo.Data1$Time > brk[2] & Photo.Data1$Time < brk[3])
  Step3 <- subset(Photo.Data1, Photo.Data1$Time > brk[3] & Photo.Data1$Time < brk[4])
  Step4 <- subset(Photo.Data1, Photo.Data1$Time > brk[4] & Photo.Data1$Time < brk[5])
  Step5 <- subset(Photo.Data1, Photo.Data1$Time > brk[5] & Photo.Data1$Time < brk[6])
  Step6 <- subset(Photo.Data1, Photo.Data1$Time > brk[6] & Photo.Data1$Time < brk[7])
  Step7 <- subset(Photo.Data1, Photo.Data1$Time > brk[7] & Photo.Data1$Time < brk[8])
  Step8 <- subset(Photo.Data1, Photo.Data1$Time > brk[8] & Photo.Data1$Time < brk[9])
  Step9 <- subset(Photo.Data1, Photo.Data1$Time > brk[9] & Photo.Data1$Time < brk[10])
  Step10 <- subset(Photo.Data1, Photo.Data1$Time > brk[10] & Photo.Data1$Time < brk[11])
  Step11 <- subset(Photo.Data1, Photo.Data1$Time > brk[11] & Photo.Data1$Time < brk[12])
  #Step12 <- subset(Photo.Data1, Photo.Data1$Time > brk[12])
  lt.levs <- list(Step1,Step2,Step3,Step4,Step5,Step6,Step7,Step8,Step9,Step10,Step11) #list levels of segmentation
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data)[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data)[1] #list length of trimmed data
    Photo.Data$sec <- as.numeric(1:n) #set seconds by one from start to finish of run
    Photo.Data$Value <- as.numeric(as.character(Photo.Data$Value)) #save O2 data as numeric
    Photo.Data$Temp <- as.numeric(as.character(Photo.Data$Temp)) #save O2 data as numeric
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("~/MyProjects/E5_ROL/RAnalysis/Output/PI_Curves/",rename,"_",j,"thinning.pdf"))
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'),  axes=FALSE) #plot data as a function of time
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data$Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    
    #set data reduction thinning parameter
    thin <- 20
    #save original unthinned data
    Photo.Data.orig<-Photo.Data
    
    Photo.Data   <-  thinData(Photo.Data.orig , by=thin)$newData1 #thin data 
    Photo.Data $sec <- as.numeric(rownames(Photo.Data)) #maintain numeric values for time
    Photo.Data$Temp <- NA # add a new column to fill with the thinned data
    Photo.Data$Temp <-  thinData(Photo.Data.orig, xy=c(2,3),by=thin)$newData1[,2] #thin data  for the temp values
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    
    Regs  <-  rankLocReg(xall=Photo.Data$sec, yall=Photo.Data$Value, alpha=0.2, 
                         method="pc", verbose=TRUE) 
    pdf(paste0("~/MyProjects/E5_ROL/RAnalysis/Output/PI_Curves/",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
  }
}


#view rates
Photo.R
write.csv(Photo.R,"~/MyProjects/E5_ROL/RAnalysis/Output/All_PI_Curve_rates.csv")

#Merge rates with sample info
Data <- merge(Photo.R, Sample.Info, by="Fragment.ID")

#correct for the size of the chamber
Data$micromol.s <- Data$µmol.L.sec * Data$Chamber.Vol.L

#calculate the average of the blanks from each time step
blnks <- subset(Data, Sample.Type=="Blank")
blnks <-mean.blnks <- aggregate(micromol.s ~ Light_Level, data=blnks, FUN=mean)
Data <- merge(Data, blnks, by="Light_Level")
colnames(Data)[colnames(Data) == 'micromol.s.x'] <- 'sample.micromol.s'
colnames(Data)[colnames(Data) == 'micromol.s.y'] <- 'blank.micromol.s'

#subtract the average of the blanks from each time step
Data$corr.micromol.s <- Data$sample.micromol.s - Data$blank.micromol.s
#Data$corr.micromol.s <- Data$micromol.s
Data$micromol.cm2.s <- Data$corr.micromol.s/Data$Surf.Area.cm2
Data$micromol.cm2.h <- Data$micromol.cm2.s*3600

Data <- subset(Data, Species!="Blank")
write.csv(Data,"~/MyProjects/E5_ROL/RAnalysis/Data/All_PI_Curve_rates.csv")

AP.geno <- subset(Data, Species=="Apul")
AP <- subset(Data, Species=="Apulchra")
PL <- subset(Data, Species=="Porites")
PM <- subset(Data, Species=="Pocillopora")
#write.csv(PA,"~/MyProjects/Holobiont_Integration/RAnalysis/Output/20180914_PI_Curve_rates_Pacuta.csv")
#MC <- subset(Data, Species=="Mcapitata")
#write.csv(MC,"~/MyProjects/Holobiont_Integration/RAnalysis/Output/20180914_PI_Curve_rates_Mcapitata.csv")


##### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980) #####

pdf("~/MyProjects/E5_ROL/RAnalysis/Output/All_NLLS_PICurves.pdf", width=10, height=5)
par(mfrow=c(2,2))

# ACROPORA Sites
#Plot curves
PAR <- as.numeric(AP.geno$Light_Value)
Pc <- as.numeric(AP.geno$micromol.cm2.h)

#pdf("~/MyProjects/E5_ROL/RAnalysis/Output/20191024_NLLS_Acropora_PICurves.pdf")
#par(mfrow=c(1,1))
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="A) A. pulchra Garden", adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PA = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                      start=list(Am=(max(Pc)-min(Pc)),AQY=0.05,Rd=-min(Pc),theta=.001)) 

my.fit.PA <- summary(curve.nlslrc.PA ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PA <- curve((1/(2*summary(curve.nlslrc.PA)$coef[4,1]))*(summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1]-sqrt((summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1])^2-4*summary(curve.nlslrc.PA)$coef[2,1]*summary(curve.nlslrc.PA)$coef[4,1]*summary(curve.nlslrc.PA)$coef[1,1]*x))-summary(curve.nlslrc.PA)$coef[3,1],lwd=2,col="blue",add=T)
#dev.off()

#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit.PA$parameters[1]

#AQY (AP.genoparent quantum yield) alpha  
AQY <- my.fit.PA$parameters[2]

#Rd (dark respiration)
Rd <- my.fit.PA$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
PI.Output <- as.data.frame(rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic))
row.names(PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
colnames(PI.Output) <- c("AP.geno")

# ACROPORA Sites
#Plot curves
PAR <- as.numeric(AP$Light_Value)
Pc <- as.numeric(AP$micromol.cm2.h)

#pdf("~/MyProjects/E5_ROL/RAnalysis/Output/20191024_NLLS_Acropora_PICurves.pdf")
#par(mfrow=c(1,1))
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="B) A. pulchra Sites", adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PA = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                      start=list(Am=(max(Pc)-min(Pc)),AQY=0.05,Rd=-min(Pc),theta=0.01)) 

my.fit.PA <- summary(curve.nlslrc.PA ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PA <- curve((1/(2*summary(curve.nlslrc.PA)$coef[4,1]))*(summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1]-sqrt((summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1])^2-4*summary(curve.nlslrc.PA)$coef[2,1]*summary(curve.nlslrc.PA)$coef[4,1]*summary(curve.nlslrc.PA)$coef[1,1]*x))-summary(curve.nlslrc.PA)$coef[3,1],lwd=2,col="blue",add=T)
#dev.off()

#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit.PA$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY <- my.fit.PA$parameters[2]

#Rd (dark respiration)
Rd <- my.fit.PA$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
PI.Output <- as.data.frame(rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic))
row.names(PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
colnames(PI.Output) <- c("Apulchra")



# PORITES
#Plot curves
PAR <- as.numeric(PL$Light_Value)
Pc <- as.numeric(PL$micromol.cm2.h)
#pdf("~/MyProjects/E5_ROL/RAnalysis/Output/20191024_NLLS_Porites_PICurves.pdf")
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="C) Porites", adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PA = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                      start=list(Am=(max(Pc)-min(Pc)),AQY=0.05,Rd=-min(Pc),theta=0.9)) 

my.fit.PA <- summary(curve.nlslrc.PA ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PA <- curve((1/(2*summary(curve.nlslrc.PA)$coef[4,1]))*(summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1]-sqrt((summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1])^2-4*summary(curve.nlslrc.PA)$coef[2,1]*summary(curve.nlslrc.PA)$coef[4,1]*summary(curve.nlslrc.PA)$coef[1,1]*x))-summary(curve.nlslrc.PA)$coef[3,1],lwd=2,col="blue",add=T)
#dev.off()
#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit.PA$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY <- my.fit.PA$parameters[2]

#Rd (dark respiration)
Rd <- my.fit.PA$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
PI.Output.PL <- as.data.frame(rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic))
row.names(PI.Output.PL) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
colnames(PI.Output.PL) <- c("Porites")
PI.Output <- cbind(PI.Output, PI.Output.PL)

# POCILLIPORA
#Plot curves

PAR <- as.numeric(PM$Light_Value)
Pc <- as.numeric(PM$micromol.cm2.h)
#pdf("~/MyProjects/E5_ROL/RAnalysis/Output/20191024_NLLS_Pocillopora_PICurves.pdf")
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="D) Pocillopora", adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PA = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                      start=list(Am=(max(Pc)-min(Pc)),AQY=0.05,Rd=-min(Pc),theta=0.1)) 

my.fit.PA <- summary(curve.nlslrc.PA ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PA <- curve((1/(2*summary(curve.nlslrc.PA)$coef[4,1]))*(summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1]-sqrt((summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1])^2-4*summary(curve.nlslrc.PA)$coef[2,1]*summary(curve.nlslrc.PA)$coef[4,1]*summary(curve.nlslrc.PA)$coef[1,1]*x))-summary(curve.nlslrc.PA)$coef[3,1],lwd=2,col="blue",add=T)
#dev.off()
#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit.PA$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY <- my.fit.PA$parameters[2]

#Rd (dark respiration)
Rd <- my.fit.PA$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
PI.Output.PM <- as.data.frame(rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic))
row.names(PI.Output.PM) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
colnames(PI.Output.PM) <- c("Pocillopora")
PI.Output <- cbind(PI.Output, PI.Output.PM)
PI.Output <- round(PI.Output, digits = 3)


# legend('topright', ncol = 4L, cex=0.9,
#        legend = c('parameter', rownames(PI.Output), 
#                   colnames(PI.Output)[1], PI.Output[1:6,1], 
#                   colnames(PI.Output)[2], PI.Output[1:6,2], 
#                   colnames(PI.Output)[3], PI.Output[1:6,3]))
dev.off()

# pdf("~/MyProjects/E5_ROL/RAnalysis/Output/20191024_NLLSR_PICurves.pdf")
# par(mfrow=c(3,1))
# #Plot input data and model fit
# plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-0.4,1.0),cex.lab=0.8,cex.axis=0.8,cex=1, main="A) A. pulchra", adj = 0.05) #set plot info
# mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
# mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels
# lines(curve.fitting.PA ,lwd=2,col="blue") #add fit line
# lines(cbind(curve.fitting.PA$x, PI.Output[2]),lwd=2,col="yellow") #add Pnet max line
# abline(v=Ik,lwd=2,col="gray", lty=3)
# #abline(PI.Output[3], PI.Output[4],lwd=2,col="green") #add alpha line
# #legend( max(PAR)-300, min(Pc)+5, legend=c("Curve Fit", "Pnet max", "Ik", "alpha"), col=c("blue", "yellow", "red", "green"), lty=1, cex=1, box.lty = 0 ) #add a legend
# dev.off()
# 
# PI.Output.PA <- as.data.frame(PI.Output)
# colnames(PI.Output.PA) <- c("Apulchra")
# #row.names(PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
# PI.Output.PA
# 
# write.csv(PI.Output.PA,"~/MyProjects/E5_ROL/RAnalysis/Output/20191024_PI_Curve_Output.csv")
