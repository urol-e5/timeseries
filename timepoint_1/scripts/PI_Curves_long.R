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
path.p<-"data/1_PICurves/" #the location of all your respirometry files 

#bring in the files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*12, ncol=4)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "µmol.L.sec", "Temp")

#Load Sample meta info Info
Sample.Info <- read.csv(file="data/1_PI_Curve_Sample_Info.csv", header=T) #read sample.info data
Sample.Info$Fragment.ID <- paste(Sample.Info$Fragment.ID,"_",Sample.Info$Light_Level, sep = "")

#subset the data by light step using time breaks in the data
for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]),  header=T, sep=",", skip=1, na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- as.data.frame(cbind(Photo.Data1$Time, Photo.Data1$Value, Photo.Data1$Temp)) #subset columns of interest
  colnames(Photo.Data1) <- c("Time", "Value", "Temp")
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", "HST") #-7*60*60 #convert time from character to time
  Photo.Data1$Time <-format(strptime(Photo.Data1$Time, format="%Y-%m-%d %H:%M:%S"), "%H:%M:%S")
  levs <- which(Sample.Info$Plug.Number ==sub("_.*", "", file.names[i]))
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
  Step10 <- subset(Photo.Data1, Photo.Data1$Time > brk[10]) #& Photo.Data1$Time < brk[11])
  #Step11 <- subset(Photo.Data1, Photo.Data1$Time > brk[11] & Photo.Data1$Time < brk[12])
  #Step12 <- subset(Photo.Data1, Photo.Data1$Time > brk[12])
  lt.levs <- list(Step1,Step2,Step3,Step4,Step5,Step6,Step7,Step8,Step9,Step10) #list levels of segmentation
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    #n<-dim(Photo.Data)[1] #identify length of data
    #Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data)[1] #list length of trimmed data
    Photo.Data$sec <- as.numeric(1:n) #set seconds by one from start to finish of run
    Photo.Data$Value <- as.numeric(as.character(Photo.Data$Value)) #save O2 data as numeric
    Photo.Data$Temp <- as.numeric(as.character(Photo.Data$Temp)) #save O2 data as numeric
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/reg_graphs/",rename,"_",j,"thinning.pdf"))
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
    pdf(paste0("output/reg_graphs/",rename,"_",j,"regression.pdf"))
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
write.csv(Photo.R,"output/All_PI_Curve_rates_not_normalized.csv")

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

#Correct for surface area
Data$micromol.cm2.s <- Data$corr.micromol.s/Data$surface.area.cm2
Data$micromol.cm2.h <- Data$micromol.cm2.s*3600

Data <- subset(Data, Species!="Blank")
write.csv(Data,"output/All_PI_Curve_rates.csv")


Data <- read.csv("output/All_PI_Curve_rates.csv", sep=",", header=TRUE)


# vvv RC vvv

# Define function for PI curve
nls_PI <- function(Light_Value, micromol.cm2.h, data, start) {
  nls(micromol.cm2.h ~ (1/(2*theta))*(AQY*Light_Value+Am-sqrt((AQY*Light_Value+Am)^2-4*AQY*theta*Am*Light_Value))-Rd, data = data, start = start)
}

# Make 'safe' function that will return NA if nls gets an error
nls_PI_safe <- safely(nls_PI, otherwise = NA_real_)

# Fit the PI curve for each coral, and pull out model parameters
Data2 <- Data %>%
  group_by(Plug.Number) %>%
  nest() %>%
  mutate(
    startpars = map(data, ~ list(
      Am = max(.$micromol.cm2.h) - min(.$micromol.cm2.h),
      Rd = -min(.$micromol.cm2.h),
      AQY = 0.005, theta = 0.001)),
    nls_fit = map2(data, startpars, ~ nls_PI_safe(data = .x, start = .y)),
    nls_res = map(nls_fit, ~ .$result),
  nls_class = map_chr(nls_res, ~ class(.)),
   nls_pars = map(nls_fit, ~ broom::tidy(.$result)))

# Unnest nls parameters and join with sample metadata
params <- Data2 %>% 
  unnest(nls_pars) %>%
  left_join(distinct(select(Data, Plug.Number, Species, Site)))

# Plot parameters for each species at each site
params %>%
  ggplot(aes(x = Site, y = estimate, color = Site)) +
  geom_boxplot() +
  facet_grid(term ~ Species, scales = "free_y")

# Plot individual PI curves and nls fits
Data2 <- Data2 %>%
  filter(nls_class == "nls") %>%
  mutate(Light_Value = list(Light_Value = seq(0:max(Data$Light_Value))),
         micromol.cm2.h = map2(nls_res, newdat, ~ predict(.x, newdata = data.frame(Light_Value = seq(0:max(Data$Light_Value))))))

preds <- Data2 %>% unnest(Light_Value, micromol.cm2.h)

Data2 <- Data2 %>%
  mutate(plot = pmap(list(dat = data, pars = nls_pars), function(dat, pars) {
    ggplot(dat, aes(x = Light_Value, y = micromol.cm2.h)) +
      geom_point() +
      geom_line(data = preds)
    }))



Data2 %>% pull(plot)



# ^^^ RC ^^^




AP.S2 <- subset(Data, Species=="Acropora" & Site =="Site2")
AP.S3 <- subset(Data, Species=="Acropora" & Site =="Site3")
PL.S2 <- subset(Data, Species=="Porites" & Site =="Site2")
PL.S3 <- subset(Data, Species=="Porites" & Site =="Site3")
PM.S2 <- subset(Data, Species=="Pocillopora" & Site =="Site2")
PM.S3 <- subset(Data, Species=="Pocillopora" & Site =="Site3")

##### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980) #####

#pdf("output/All_NLLS_PICurves.pdf", width=10, height=5)
#par(mfrow=c(1,3))

# ACROPORA 
#Plot curves
PAR.S2 <- as.numeric(AP.S2$Light_Value)
Pc.S2 <- as.numeric(AP.S2$micromol.cm2.h)

PAR.S3 <- as.numeric(AP.S3$Light_Value)
Pc.S3 <- as.numeric(AP.S3$micromol.cm2.h)


pdf("output/NLLS_Acropora_PICurves.pdf")
#par(mfrow=c(1,1))
plot(PAR.S2,Pc.S2, col="red", xlab="", ylab="", xlim=c(0,max(PAR.S2)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="A) Acropora", adj = 0.05) #set plot info
points(PAR.S3,Pc.S3, col="blue", xlab="", ylab="", xlim=c(0,max(PAR.S3)), ylim=c(-1, 2), adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.AP.S2 = nls(Pc.S2 ~ (1/(2*theta))*(AQY*PAR.S2+Am-sqrt((AQY*PAR.S2+Am)^2-4*AQY*theta*Am*PAR.S2))-Rd,
                      start=list(Am=(max(Pc.S2)-min(Pc.S2)),AQY=0.05,Rd=-min(Pc.S2),theta=0.0001)) 

curve.nlslrc.AP.S3 = nls(Pc.S3 ~ (1/(2*theta))*(AQY*PAR.S3+Am-sqrt((AQY*PAR.S3+Am)^2-4*AQY*theta*Am*PAR.S3))-Rd,
                         start=list(Am=(max(Pc.S3)-min(Pc.S3)),AQY=0.05,Rd=-min(Pc.S3),theta=0.001)) 

my.fit.AP.S2 <- summary(curve.nlslrc.AP.S2 ) #summary of model fit
my.fit.AP.S3 <- summary(curve.nlslrc.AP.S3 ) #summary of model fit

#draw the curve using the model fit
curve.fitting.AP.S2 <- curve((1/(2*summary(curve.nlslrc.AP.S2)$coef[4,1]))*(summary(curve.nlslrc.AP.S2)$coef[2,1]*x+summary(curve.nlslrc.AP.S2)$coef[1,1]-sqrt((summary(curve.nlslrc.AP.S2)$coef[2,1]*x+summary(curve.nlslrc.AP.S2)$coef[1,1])^2-4*summary(curve.nlslrc.AP.S2)$coef[2,1]*summary(curve.nlslrc.AP.S2)$coef[4,1]*summary(curve.nlslrc.AP.S2)$coef[1,1]*x))-summary(curve.nlslrc.AP.S2)$coef[3,1],lwd=2,col="red",add=T)
curve.fitting.AP.S3 <- curve((1/(2*summary(curve.nlslrc.AP.S3)$coef[4,1]))*(summary(curve.nlslrc.AP.S3)$coef[2,1]*x+summary(curve.nlslrc.AP.S3)$coef[1,1]-sqrt((summary(curve.nlslrc.AP.S3)$coef[2,1]*x+summary(curve.nlslrc.AP.S3)$coef[1,1])^2-4*summary(curve.nlslrc.AP.S3)$coef[2,1]*summary(curve.nlslrc.AP.S3)$coef[4,1]*summary(curve.nlslrc.AP.S3)$coef[1,1]*x))-summary(curve.nlslrc.AP.S3)$coef[3,1],lwd=2,col="blue",add=T)

dev.off()

#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross.S2 <- my.fit.AP.S2$parameters[1]
Pmax.gross.S3 <- my.fit.AP.S3$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY.S2 <- my.fit.AP.S2$parameters[2]
AQY.S3 <- my.fit.AP.S3$parameters[2]

#Rd (dark respiration)
Rd.S2 <- my.fit.AP.S2$parameters[3]
Rd.S3 <- my.fit.AP.S3$parameters[3]

# Ik light saturation point
Ik.S2 <- Pmax.gross.S2/AQY.S2
Ik.S3 <- Pmax.gross.S3/AQY.S3

# Ic light compensation point
Ic.S2 <- Rd.S2/AQY.S2
Ic.S3 <- Rd.S3/AQY.S3

# Net photosynthetic rates
Pmax.net.S2 <- Pmax.gross.S2 - Rd.S2
Pmax.net.S3 <- Pmax.gross.S3 - Rd.S3

# #output parameters into a table
# PI.Output <- as.data.frame(rbind(Pmax.gross.S2, Pmax.net, -Rd, AQY,Ik,Ic))
# colnames(PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
# rownames(PI.Output) <- c("Site1", "Site2")



# PORITES
#Plot curves
PAR.S2 <- as.numeric(PL.S2$Light_Value)
Pc.S2 <- as.numeric(PL.S2$micromol.cm2.h)

PAR.S3 <- as.numeric(PL.S3$Light_Value)
Pc.S3 <- as.numeric(PL.S3$micromol.cm2.h)


pdf("output/NLLS_Porites_PICurves.pdf")
#par(mfrow=c(1,1))
plot(PAR.S2,Pc.S2, col="red", xlab="", ylab="", xlim=c(0,max(PAR.S2)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="B) Porites", adj = 0.05) #set plot info
points(PAR.S3,Pc.S3, col="blue", xlab="", ylab="", xlim=c(0,max(PAR.S3)), ylim=c(-1, 2), adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PL.S2 = nls(Pc.S2 ~ (1/(2*theta))*(AQY*PAR.S2+Am-sqrt((AQY*PAR.S2+Am)^2-4*AQY*theta*Am*PAR.S2))-Rd,
                         start=list(Am=(max(Pc.S2)-min(Pc.S2)),AQY=0.05,Rd=-min(Pc.S2),theta=0.01)) 

curve.nlslrc.PL.S3 = nls(Pc.S3 ~ (1/(2*theta))*(AQY*PAR.S3+Am-sqrt((AQY*PAR.S3+Am)^2-4*AQY*theta*Am*PAR.S3))-Rd,
                         start=list(Am=(max(Pc.S3)-min(Pc.S3)),AQY=0.05,Rd=-min(Pc.S3),theta=0.001)) 

my.fit.PL.S2 <- summary(curve.nlslrc.PL.S2 ) #summary of model fit
my.fit.PL.S3 <- summary(curve.nlslrc.PL.S3 ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PL.S2 <- curve((1/(2*summary(curve.nlslrc.PL.S2)$coef[4,1]))*(summary(curve.nlslrc.PL.S2)$coef[2,1]*x+summary(curve.nlslrc.PL.S2)$coef[1,1]-sqrt((summary(curve.nlslrc.PL.S2)$coef[2,1]*x+summary(curve.nlslrc.PL.S2)$coef[1,1])^2-4*summary(curve.nlslrc.PL.S2)$coef[2,1]*summary(curve.nlslrc.PL.S2)$coef[4,1]*summary(curve.nlslrc.PL.S2)$coef[1,1]*x))-summary(curve.nlslrc.PL.S2)$coef[3,1],lwd=2,col="red",add=T)
curve.fitting.PL.S3 <- curve((1/(2*summary(curve.nlslrc.PL.S3)$coef[4,1]))*(summary(curve.nlslrc.PL.S3)$coef[2,1]*x+summary(curve.nlslrc.PL.S3)$coef[1,1]-sqrt((summary(curve.nlslrc.PL.S3)$coef[2,1]*x+summary(curve.nlslrc.PL.S3)$coef[1,1])^2-4*summary(curve.nlslrc.PL.S3)$coef[2,1]*summary(curve.nlslrc.PL.S3)$coef[4,1]*summary(curve.nlslrc.PL.S3)$coef[1,1]*x))-summary(curve.nlslrc.PL.S3)$coef[3,1],lwd=2,col="blue",add=T)

dev.off()

#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross.S2 <- my.fit.PL.S2$parameters[1]
Pmax.gross.S3 <- my.fit.PL.S3$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY.S2 <- my.fit.PL.S2$parameters[2]
AQY.S3 <- my.fit.PL.S3$parameters[2]

#Rd (dark respiration)
Rd.S2 <- my.fit.PL.S2$parameters[3]
Rd.S3 <- my.fit.PL.S3$parameters[3]

# Ik light saturation point
Ik.S2 <- Pmax.gross.S2/AQY.S2
Ik.S3 <- Pmax.gross.S3/AQY.S3

# Ic light compensation point
Ic.S2 <- Rd.S2/AQY.S2
Ic.S3 <- Rd.S3/AQY.S3

# Net photosynthetic rates
Pmax.net.S2 <- Pmax.gross.S2 - Rd.S2
Pmax.net.S3 <- Pmax.gross.S3 - Rd.S3

# #output parameters into a table
# PI.Output.PL <- as.data.frame(rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic))
# row.names(PI.Output.PL) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
# colnames(PI.Output.PL) <- c("Porites")
# PI.Output <- cbind(PI.Output, PI.Output.PL)

# POCILLIPORA
#Plot curves
PAR.S2 <- as.numeric(PM.S2$Light_Value)
Pc.S2 <- as.numeric(PM.S2$micromol.cm2.h)

PAR.S3 <- as.numeric(PM.S3$Light_Value)
Pc.S3 <- as.numeric(PM.S3$micromol.cm2.h)


pdf("output/NLLS_Pocillopora_PICurves.pdf")
#par(mfrow=c(1,1))
plot(PAR.S2,Pc.S2, col="red", xlab="", ylab="", xlim=c(0,max(PAR.S2)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="C) Pocillopora", adj = 0.05) #set plot info
points(PAR.S3,Pc.S3, col="blue", xlab="", ylab="", xlim=c(0,max(PAR.S3)), ylim=c(-1, 2), adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PM.S2 = nls(Pc.S2 ~ (1/(2*theta))*(AQY*PAR.S2+Am-sqrt((AQY*PAR.S2+Am)^2-4*AQY*theta*Am*PAR.S2))-Rd,
                         start=list(Am=(max(Pc.S2)-min(Pc.S2)),AQY=0.05,Rd=-min(Pc.S2),theta=0.00001)) 

curve.nlslrc.PM.S3 = nls(Pc.S3 ~ (1/(2*theta))*(AQY*PAR.S3+Am-sqrt((AQY*PAR.S3+Am)^2-4*AQY*theta*Am*PAR.S3))-Rd,
                         start=list(Am=(max(Pc.S3)-min(Pc.S3)),AQY=0.05,Rd=-min(Pc.S3),theta=0.0001)) 

my.fit.PM.S2 <- summary(curve.nlslrc.PM.S2 ) #summary of model fit
my.fit.PM.S3 <- summary(curve.nlslrc.PM.S3 ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PM.S2 <- curve((1/(2*summary(curve.nlslrc.PM.S2)$coef[4,1]))*(summary(curve.nlslrc.PM.S2)$coef[2,1]*x+summary(curve.nlslrc.PM.S2)$coef[1,1]-sqrt((summary(curve.nlslrc.PM.S2)$coef[2,1]*x+summary(curve.nlslrc.PM.S2)$coef[1,1])^2-4*summary(curve.nlslrc.PM.S2)$coef[2,1]*summary(curve.nlslrc.PM.S2)$coef[4,1]*summary(curve.nlslrc.PM.S2)$coef[1,1]*x))-summary(curve.nlslrc.PM.S2)$coef[3,1],lwd=2,col="red",add=T)
curve.fitting.PM.S3 <- curve((1/(2*summary(curve.nlslrc.PM.S3)$coef[4,1]))*(summary(curve.nlslrc.PM.S3)$coef[2,1]*x+summary(curve.nlslrc.PM.S3)$coef[1,1]-sqrt((summary(curve.nlslrc.PM.S3)$coef[2,1]*x+summary(curve.nlslrc.PM.S3)$coef[1,1])^2-4*summary(curve.nlslrc.PM.S3)$coef[2,1]*summary(curve.nlslrc.PM.S3)$coef[4,1]*summary(curve.nlslrc.PM.S3)$coef[1,1]*x))-summary(curve.nlslrc.PM.S3)$coef[3,1],lwd=2,col="blue",add=T)

dev.off()

#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross.S2 <- my.fit.PM.S2$parameters[1]
Pmax.gross.S3 <- my.fit.PM.S3$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY.S2 <- my.fit.PM.S2$parameters[2]
AQY.S3 <- my.fit.PM.S3$parameters[2]

#Rd (dark respiration)
Rd.S2 <- my.fit.PM.S2$parameters[3]
Rd.S3 <- my.fit.PM.S3$parameters[3]

# Ik light saturation point
Ik.S2 <- Pmax.gross.S2/AQY.S2
Ik.S3 <- Pmax.gross.S3/AQY.S3

# Ic light compensation point
Ic.S2 <- Rd.S2/AQY.S2
Ic.S3 <- Rd.S3/AQY.S3

# Net photosynthetic rates
Pmax.net.S2 <- Pmax.gross.S2 - Rd.S2
Pmax.net.S3 <- Pmax.gross.S3 - Rd.S3

