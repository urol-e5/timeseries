#Title: Calculating calcification rates following the total alkalinity anomaly method
#Project: NSF E5 UROL
#Author: NJ Silbiger
#Edited by: DM Becker, AS Huffmyer, E5 team members
#Date Last Modified: 20210504

#clear workspace
rm(list=ls())

##Install packages
#load packages
library(broom)
library(purrr)
library(lubridate)
library(tidyverse)
library(nlstools)
library(here)
library(stringr)

#set wd
here()

#bring in calcification data file with TA and chamber pH, temp, salinity measurements

calc.data <- read.csv("data/2_calcification/2_TA_data.csv")

#create new data frame of just the initial data (initial bottles taken before each run for initial TA measurements), pull out initial data from sample.type

rows.initial <- which(calc.data$sample.type == "Initial") #tells you all the rows that you have with initial sample type

initial <- calc.data[rows.initial,] #shows you the rows with the initial sample type

calc.data <- calc.data[-rows.initial,] #to remove the rows with initial data

#remove and create new data frame with just blanks

rows.blanks <- which(calc.data$sample.type == "Blank") #tells you all the rows that you have wwith blanks

blanks <-calc.data[rows.blanks,] #shows you the rows with the blank sample type

calc.data <- calc.data[-rows.blanks,] #to remove the rows with blank data

#need to join sample data frame, join initial with your calc.data, only by temperature, pull out the columns we need

initial <- initial[, c("date","colony_ID", "salinity.chamber", "salinity.lab", "TA", "mV.chamber")]
names(initial)[3:6] <- paste0(names(initial)[3:6], "_initial") #use this to rename all of our columns

#join blanks and carb chem data frame

calc.data <- left_join(calc.data, initial) #joining the initials to the data frame for carb chem
blanks <- left_join(blanks,initial)

#figure out delta TA, initial-final

blanks$delta.TA.blank <- blanks$TA_initial - blanks$TA

#getting the averages of blanks for each temperature and each date
mean.blanks <- blanks %>% 
  group_by(date) %>%
  summarise(mean.blanks=mean(delta.TA.blank))

calc.data <- left_join(calc.data, mean.blanks) #bring in mean blanks to calc.data

#need to join in SA, time data by colony ID, before calculating NEC

sample.data <- read.csv("../timepoint_2/output/2_surface_area.csv") #bring in SA and volume data sheet

SA <- sample.data[, c("colony_id", "surf.area.cm2")] #pull out the necessary columns and treatment 

calc.data2 <- left_join(calc.data, SA) # join carb chem and SA data

#bring in the time data from resp.data sheet
resp.data <- read.csv("data/2_calcification/2_DeltaTA_metadata.csv")

#pull out columns that we want to use for our bind to the calc.data2 sheet

time.data <- resp.data[, c("colony_id","TA.Start.Time", "TA.Stop.Time")] 

full.calc.data <- left_join(calc.data2, time.data)

#adjust the time information and format

#convert time from character to time
full.calc.data$start.time <- strptime(as.character(full.calc.data$TA.Start.Time), "%I:%M:%S %p")

#convert time from character to time
full.calc.data$stop.time <- strptime(as.character(full.calc.data$TA.Stop.Time), "%I:%M:%S %p")

#calculate the net ecosystem calcification rate

full.calc.data$deltaTA<- (full.calc.data$TA_initial - full.calc.data$TA) - full.calc.data$mean.blanks
full.calc.data$timediff <- as.numeric((full.calc.data$stop.time - full.calc.data$start.time)) 

#convert volume from cm3 to m3
full.calc.data$volume <- full.calc.data$volume * 0.000001

#equation to calculate NEC rates 

full.calc.data$umol.cm2.hr <- (full.calc.data$deltaTA/2)*(1.023)*(full.calc.data$volume/full.calc.data$surf.area.cm2)*(1/full.calc.data$timediff)*(1/1000)


#anything that it <0 make it zero
calcification$umol.cm2.hr[calcification$umol.cm2.hr<0]<-0

#log x+ 1 you will have to log x + whatever the difference is here, a line for dissolution and calcification
calcification$umol.cm2.hr <- log(calcification$umol.cm2.hr+1) 

write.csv(calcification, 'data/2_calcification/2_calcification_rates.csv') 
