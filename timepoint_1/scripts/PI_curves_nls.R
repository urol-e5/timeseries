#Title: Photosynthesis Irradiance Curves
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191024
#See Readme file for details

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
  group_by(colony)id %>%
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
  left_join(distinct(select(Data, colony_id, Species, Site)))

# Plot parameters for each species at each site
params %>%
  ggplot(aes(x = Site, y = estimate, color = Site)) +
  geom_boxplot() +
  facet_grid(term ~ Species, scales = "free_y")

# Plot individual PI curves and nls fits
dd <- Data2 %>%
  filter(nls_class == "nls") %>%
  mutate(Light_Value = list(Light_Value = seq(0:max(Data$Light_Value))),
         micromol.cm2.h = map(nls_res, ~ predict(.x, newdata = data.frame(Light_Value = seq(0:max(Data$Light_Value))))))

preds <- dd %>% unnest(Light_Value, micromol.cm2.h) %>%
  nest(-Plug.Number, .key = "preds") %>%
  full_join(select(dd, Plug.Number, data, nls_pars))

preds <- preds %>%
  mutate(plot = pmap(list(dat = data, pars = nls_pars, preds = preds, id = as.character(Plug.Number)), function(dat, pars, preds, id) {
    ggplot(dat, aes(x = Light_Value, y = micromol.cm2.h)) +
      geom_point() +
      geom_line(data = preds) +
      geom_hline(yintercept = c(0, pull(filter(pars, term == "Am"), estimate)), lty = c(3,2)) +
      labs(data = id, title = as.character(id))
  }))


preds %>% pull(plot)
###

Data %>%
  group_by(Plug.Number) %>%
  do(nls_PI_safe(data = ., start = list(Am = 1, Rd = 0.2, AQY = 0.005, theta = 0.001)))



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
                         start=list(Am=(max(Pc.S3)-min(Pc.S3)),AQY=0.05,Rd=-min(Pc.S3),theta=0.0001)) 

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
                         start=list(Am=(max(Pc.S2)-min(Pc.S2)),AQY=0.05,Rd=-min(Pc.S2),theta=0.001)) 

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
                         start=list(Am=(max(Pc.S3)-min(Pc.S3)),AQY=0.05,Rd=-min(Pc.S3),theta=0.001)) 

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