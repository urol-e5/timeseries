##Practicing calcification rates 

rm(list=ls())

##Install packages
# load packages
library(broom)
library(purrr)
library(lubridate)
library(tidyverse)
library(nlstools)
library(here)
library(stringr)

#set wd
here()

#bring in files 

calc.data <- read.csv("TPC_curves/Data/CarbChem2020-06-27.csv")

#create new data frame of just the initial data, pulll out initiail data from sample.type

rows.initial <- which(calc.data$Sample.Type == "Initial") #tells you all the rows that you have wwith initial

initial <- calc.data[rows.initial,]

calc.data <- calc.data[-rows.initial,] #to remove the rows with initial data

#remove and create new data frame with just blanks

rows.blanks <- which(calc.data$Sample.Type == "Blank") #tells you all the rows that you have wwith blanks

blanks <-calc.data[rows.blanks,]

calc.data <- calc.data[-rows.blanks,] #to remove the rows with blank data

#need to join sample data frame, join initial with your calc.data, only by temperature, pull out the columns we need

initial <- initial[, c("date","temp.Cat", "salinity.insitu", "salinity.lab", "TA", "pH")]
names(initial)[3:6] <- paste0(names(initial)[3:6], "_initial") #use this to rename all of our columns

#join blanks and carb chem data frame

calc.data <- left_join(calc.data, initial) #joing the initials to my data frame for carb chem
blanks <- left_join(blanks,initial)

#figure out delta TA, initial-final

blanks$delta.TA.blank <- blanks$TA_initial - blanks$TA

#getting the averages of blanks for each temperature and each date
mean.blanks <- blanks %>% 
  group_by(date, temp.Cat) %>%
  summarise(mean.blanks=mean(delta.TA.blank))

calc.data <- left_join(calc.data, mean.blanks) #bring in mean blanks to calc.data

#need to join in SA,volume data and time data by fragment ID, before calculating NEC

sample.data <- read.csv("Respirometry/Data/sample_info_TT.csv") #bring in SA and volume data sheet

SA <- sample.data[, c("fragment.ID", "surf.area.cm2", "volume")] #pull out the necessary columns and treatment 

SA2 <- SA %>%
  dplyr::mutate(light_dark = str_extract(fragment.ID, "[^_]+$")) #split the fragment.ID at the underscore to delete the Dark values

#steps to organize and delete D values
rows.D <- which(SA2$light_dark == "D") #tells you all the rows that you have with blanks
D <-SA2[rows.D ,] 
SA <- SA2[-rows.D,] #to remove the rows with dark data, now have data sheet with only light, SA and volume values
SA[,4] = NULL #had to delete the light-dark column to allow to just join by fragment.ID

calc.data2 <- left_join(calc.data, SA)

#bring in the time data from resp.data sheet
resp.data <- read.csv("Respirometry/Data/resp_data_TT.csv")

#need to organize and first delete the dark information
rows.dark <- which(resp.data$light_dark == "dark") #tells you all the rows that you have with blanks
dark <-resp.data[rows.dark ,] 
resp.data <- resp.data[-rows.dark,]

#pull out columns that we want to use for our bind to the calc.data2 sheet

time.data <- resp.data[, c("fragment.ID","treatment", "temp.Cat", "start.time", "stop.time")] 

full.calc.data <- left_join(calc.data2, time.data)

#adjust the time information and format
#full.calc.data$start.time <- hms(full.calc.data$start.time)#convert time from character to time
full.calc.data$start.time <- strptime(as.character(full.calc.data$start.time), "%I:%M:%S %p")
#full.calc.data$stop.time <- hms(full.calc.data$stop.time) #convert time from character to time
full.calc.data$stop.time <- strptime(as.character(full.calc.data$stop.time), "%I:%M:%S %p")

#calculate the NEC rate

full.calc.data$deltaTA<- (full.calc.data$TA_initial - full.calc.data$TA) - full.calc.data$mean.blanks
full.calc.data$timediff <- as.numeric((full.calc.data$stop.time - full.calc.data$start.time)) 

#make a row that has rate.type for "C" assigned to light 

full.calc.data$rate.type <-ifelse(full.calc.data$light_dark=='light', "C")
view(full.calc.data)

#convert volume from cm3 to m3
#full.calc.data$volume <- full.calc.data$volume * 0.000001

#equation to calculate NEC rates 

full.calc.data$umol.cm2.hr <- (full.calc.data$deltaTA/2)*(1.023)*(full.calc.data$volume/full.calc.data$surf.area.cm2)*(1/full.calc.data$timediff)*(1/1000)

#reformat calc data sheet here to match the exact format for my TPC photo.T data sheet
#use split function to create new column that lists frag.ID without the "_L"

full.calc.data$individual.ID <- strsplit(full.calc.data$fragment.ID, '_L')

#need to bring in the Temp.C linked to NP for calcification, edit data frame by removing all GP and R data so I can use the NP Temp.C values in my new data frame

adjust.calc.data <- read.csv("TPC_curves/Data/Photo.T.csv") #bring in photo.T sheet with Temp.C values
rows.GP <- which(adjust.calc.data$rate.type == "GP") #tells you all the rows that you have wwith GP
GP <- adjust.calc.data[-rows.GP,] #to remove the rows with GP data

rows.R <- which(GP$rate.type == "R") #tells you all the rows that you have wwith R
NP.only.data <- GP[-rows.R,] #to remove the rows with R data

#pull out only columns I need to link Temp.C values with fragment ID
Temp.C.data <- NP.only.data[, c("fragment.ID", "temp.Cat", "Temp.C")] #pull out the necessary columns and treatment 
full.calc.data2 <- full.calc.data[, c("individual.ID", "temp.Cat", "treatment", "umol.cm2.hr", "rate.type", "light_dark", "fragment.ID")] #pull out the necessary columns and treatment 
calc.temp.c <- left_join(Temp.C.data, full.calc.data2)

#remove fragments that where not used for calcicifation



calcification <- calc.temp.c[, c("individual.ID", "temp.Cat", "treatment", "Temp.C", "umol.cm2.hr", "rate.type", "light_dark", "fragment.ID")] #pull out the necessary columns and treatment 

#replace the _L with and _C for all fragment IDs

calcification$fragment.ID <- str_replace_all(calcification$fragment.ID, fixed("_L"), "_C")
calcification$individual.ID <- as.character(calcification$individual.ID)

#anything that it <0 make it zero
calcification$umol.cm2.hr[calcification$umol.cm2.hr<0]<-0

#log x+ 1 you will have to log x + whatever the difference is here, a line for dissolution and calcification
calcification$umol.cm2.hr <- log(calcification$umol.cm2.hr+1) 

write.csv(calcification, 'TPC_curves/Data/calcification1.csv') 
