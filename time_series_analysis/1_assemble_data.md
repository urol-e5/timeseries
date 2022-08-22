Assemble master data frame for time series analysis
================
Ariana S Huffmyer, E5 RoL Team
04/14/2022

# Set Up

``` r
## install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("lme4")) install.packages("lme4")
if (!require("lmerTest")) install.packages("lmerTest")
if (!require("car")) install.packages("car")
if (!require("effects")) install.packages("effects")
if (!require("ggfortify")) install.packages("ggfortify")
if (!require("cowplot")) install.packages("cowplot")
if (!require("vegan")) install.packages("vegan")
if (!require("corrr")) install.packages("corrr")
if (!require("ggcorrplot")) install.packages("ggcorrplot")
if (!require("GGally")) install.packages("GGally")
if (!require("broom")) install.packages("broom")
if (!require("cowplot")) install.packages("cowplot")

# load packages
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(car)
library(effects)
library(ggfortify)
library(cowplot)
library(vegan)
library(corrr)
library(ggcorrplot)
library(GGally)
library(broom)
library(cowplot)
```

# Load and manipulate data

## Loading data files

Load all .csv files from output of all timepoints for each biological
response

Note: I had to remake some of these files in each timepoint because some
files were saving with capitol letters (?) and therefore not loading
into the master dataframe.

``` r
biomass_files<-list.files("../", pattern = "biomass.csv", recursive=T, full.names=T)
antioxidant_files<-list.files("../", pattern = "antioxidant_capacity.csv", recursive=T, full.names=T)
pi_curve_pars_files<-list.files("../", pattern = "pi_curve_pars_nls.csv", recursive=T, full.names=T)
surface_area_files<-list.files("../", pattern = "surface_area.csv", recursive=T, full.names=T)
protein_files<-list.files("../", pattern = "protein.csv", recursive=T, full.names=T)
symb_densities_files<-list.files("../", pattern = "symb_densities.csv", recursive=T, full.names=T)
chlorophyll_files<-list.files("../", pattern = "chlorophyll.csv", recursive=T, full.names=T)
calcification_files<-list.files("../", pattern = "calcification_output.csv", recursive=T, full.names=T)
#need to add calcification once complete
```

## Read data files

Load all data frames.

``` r
#biomass 
biomass_dataset <- data.frame()

for (i in 1:length(biomass_files)){
  biomass_df <- read.csv(biomass_files[i]) #each file will be read in
  biomass_dataset <- rbind(biomass_dataset, biomass_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#antioxidant capacity 
antioxidant_dataset <- data.frame()

for (i in 1:length(antioxidant_files)){
  antioxidant_df <- read.csv(antioxidant_files[i]) #each file will be read in
  antioxidant_dataset <- rbind(antioxidant_dataset, antioxidant_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#PI curve pars
pi_curve_pars_dataset <- data.frame()

for (i in 1:length(pi_curve_pars_files)){
  pi_curve_pars_df <- read.csv(pi_curve_pars_files[i]) #each file will be read in
  pi_curve_pars_dataset <- rbind(pi_curve_pars_dataset, pi_curve_pars_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#surface area 
surface_area_dataset <- data.frame()

for (i in 1:length(surface_area_files)){
  surface_area_df <- read.csv(surface_area_files[i]) #each file will be read in
  surface_area_dataset <- rbind(surface_area_dataset, surface_area_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#protein 
protein_dataset <- data.frame()

for (i in 1:length(protein_files)){
  protein_df <- read.csv(protein_files[i]) #each file will be read in
  protein_dataset <- rbind(protein_dataset, protein_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#symb densities 
symb_densities_dataset <- data.frame()

for (i in 1:length(symb_densities_files)){
  symb_densities_df <- read.csv(symb_densities_files[i]) #each file will be read in
  symb_densities_dataset <- rbind(symb_densities_dataset, symb_densities_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#chlorophyll files 
chlorophyll_dataset <- data.frame()

for (i in 1:length(chlorophyll_files)){
  chlorophyll_df <- read.csv(chlorophyll_files[i]) #each file will be read in
  chlorophyll_dataset <- rbind(chlorophyll_dataset, chlorophyll_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#calcification files 
calcification_dataset <- data.frame()

for (i in 1:length(calcification_files)){
  calcification_df <- read.csv(calcification_files[i]) #each file will be read in
  calcification_dataset <- rbind(calcification_dataset, calcification_df) #for each iteration, bind the new data to the building dataset
}
```

``` r
#remove files that are not needed from loops
rm(list = ls(pattern = "*_df"))
```

## Generate master data frame

Read in tag metadata frame.

``` r
#Load tag metadata sheet 
tags<-read.csv("../metadata/coral_id_metadata.csv")
```

Prepare datasets for merging by renaming columns and spreading
dataframes. Each file needs to have one line per colony per response
variable in order to merge and should be without site/species columns as
these will be added from metadata sheet.

``` r
#antioxidant data 
antioxidant_dataset<-antioxidant_dataset%>%
  select(colony_id, cre.umol.mgprot, timepoint)

#biomass data 
biomass_dataset<- biomass_dataset %>%
  filter(!colony_id=="1xPBS")%>%
  nest(value_col = c(AFDW.mg.cm2, DW.mg.cm2)) %>%
  spread(key = partner, value = value_col) %>%
  unnest(Host, Sym, .sep = '_')%>%
  select(colony_id, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Host_DW.mg.cm2, Sym_DW.mg.cm2, timepoint)

#protein data
protein_dataset<-protein_dataset%>%
  select(colony_id, prot_ug.cm2, timepoint)
  
#symbiont densitites
symb_densities_dataset<-symb_densities_dataset%>%
  select(colony_id, cells.cm2, timepoint)
```

Join all data frames together.

``` r
#can also use "join_all"
master <- full_join(full_join(full_join(full_join(full_join(full_join(full_join(
  antioxidant_dataset,
  biomass_dataset, by=c("colony_id", "timepoint")),
  chlorophyll_dataset, by=c("colony_id", "timepoint")),
  protein_dataset, by=c("colony_id", "timepoint")),
  surface_area_dataset, by=c("colony_id", "timepoint")), 
  symb_densities_dataset, by=c("colony_id", "timepoint")), 
  pi_curve_pars_dataset, by=c("colony_id", "timepoint")), 
  calcification_dataset, by=c("colony_id", "timepoint"))

head(master)
```

    ##   colony_id cre.umol.mgprot  timepoint Host_AFDW.mg.cm2 Sym_AFDW.mg.cm2
    ## 1   ACR-139      0.07678208 timepoint1        1.4072814       0.9942749
    ## 2   ACR-140      0.16359919 timepoint1        0.8512491       0.8185088
    ## 3   ACR-145      0.11473069 timepoint1        0.8806505       0.3906111
    ## 4   ACR-150      0.08877581 timepoint1        1.0759520       1.2938664
    ## 5   ACR-165      0.10763769 timepoint1        1.6758688       0.8437534
    ## 6   ACR-173      0.12499656 timepoint1        1.2823276       0.9486607
    ##   Host_DW.mg.cm2 Sym_DW.mg.cm2 chla.ug.cm2 chlc2.ug.cm2 prot_ug.cm2
    ## 1      293.15578      375.3623    2.788699    2.5416344    581.4959
    ## 2      147.35613      191.1504    1.352630    0.8816772    391.1590
    ## 3       78.49784      113.6264    1.790804    1.2134689    308.9662
    ## 4      128.25125      161.1390    2.801127    1.4677382    557.2875
    ## 5      197.68680      258.6335    2.770324    1.0237056    377.5935
    ## 6      260.60473      326.5465    2.991138    3.2424896    476.5595
    ##   surface.area.cm2 cells.cm2        Am         AQY        Rd  umol.cm2.hr
    ## 1         16.99731  628007.8 1.0423651 0.002247138 0.4110558 0.0069629397
    ## 2         24.43468  403797.7 0.7395432 0.002678070 0.4154000 0.0074498259
    ## 3         12.67245  684950.4 0.7092019 0.002476605 0.2058895 0.0000000000
    ## 4         32.30627  909111.4 2.0342383 0.003363935 0.8398342 0.0021854881
    ## 5         14.60735  756468.6 1.0504973 0.004994669 0.3808200 0.0002372078
    ## 6         22.92706  800364.3 0.7474726 0.001818345 0.4458121 0.0000000000

Identify rows that have no data for any metric (typos, blanks, missing
colonies). This helps with cleaning the data set.

``` r
master<-master%>%
  relocate(timepoint, .after = colony_id)%>%
  mutate(master, "RowSum" = rowSums(is.na(master[ , 3:16])))%>%
  filter(!RowSum>=12) %>% #remove rows with all na's (rowSums=11)
  select(!RowSum)
```

Correct for colony ID’s that changed over time.

Control for colony id, site name, and timepoint names.

``` r
#Add column for month of year for each timepoint
master$month[master$timepoint == "timepoint1"] <- "January 2020"
master$month[master$timepoint == "timepoint2"] <- "March 2020"
master$month[master$timepoint == "timepoint3"] <- "September 2020"
master$month[master$timepoint == "timepoint4"] <- "November 2020"
```

``` r
#replace colony id with the original tag number if the tag was changed during the timeseries, as shown by a new tag id at the end of the time series
master$colony_id_corr<-tags$Original_Tag[match(master$colony_id, tags$Tag_ID_Nov)] 

#if there is an na (because the colony tag id did not change over the timesereis), populate with the colonyid number from colony_id column
master$colony_id_corr[is.na(master$colony_id_corr)] <- master$colony_id[is.na(master$colony_id_corr)]

#reorder columns
master<-master%>%
    relocate(month, .after = timepoint)%>%
   relocate(colony_id_corr, .after = colony_id)

#change the one colony id manually that was changed in the middle of the series
master<-master%>%
  mutate(colony_id_corr=as.character(colony_id_corr))%>%
  mutate(colony_id_corr=if_else(colony_id=="POR-206", "POR-251", colony_id_corr))%>%
  mutate(colony_id_corr=as.factor(colony_id_corr))

#match site name by site number and species from tag sheet 
master$site<-tags$Site[match(master$colony_id_corr, tags$Original_Tag)]
master$species<-tags$Coral_Species[match(master$colony_id_corr, tags$Original_Tag)]

#add a nutrient level for each site  
master<-master%>%
  mutate(nutrient = if_else(site == "Hilton", "Medium", 
                            if_else(site == "Mahana", "Low", 
                                    if_else(site == "Manava", "High", "NA"))))%>%
  mutate(site_code = paste(site, nutrient))

master$site_code <- factor(master$site_code, levels = c("Mahana Low", "Hilton Medium", "Manava High"))

#rearrange columns
master<-master%>%
  relocate(site, .after = month)%>%
  relocate(species, .after = colony_id_corr)%>%
  relocate(nutrient, .after = site)%>%
  relocate(site_code, .after = nutrient)
```

140 colonies in this dataset. This number is expected.

``` r
length(levels(master$colony_id_corr))
```

    ## [1] 140

Check for duplicates, will return TRUE if there are no duplicates.

``` r
master<-master%>%
  mutate(code=as.factor(paste(colony_id_corr, "-", timepoint)))

length(unique(master$code)) == nrow(master)
```

    ## [1] FALSE

``` r
length(unique(master$code))
```

    ## [1] 448

``` r
nrow(master)
```

    ## [1] 525

``` r
master<-master%>%
  unique()

length(unique(master$code)) == nrow(master)
```

    ## [1] TRUE

``` r
length(unique(master$code))
```

    ## [1] 448

``` r
nrow(master)
```

    ## [1] 448

Export list of duplicate values to examine.

``` r
duplicates<-master[duplicated(master$code),]
duplicates<-duplicates%>%
   droplevels()
list<-levels(duplicates$code)
 
 master%>%
   filter(code %in% list)%>%
   write_csv(file="output/timeseries_duplicates_qc.csv")
```

There are no duplicates.

# Normalization and write master file

In order to examine the effect of normalizer on all response variables,
normalize each response to: (1) Host protein (ug/cm2)  
(2) Host AFDW (mg/cm2)  
(3) Host surface area (cm2)

To do this, back-normalize variables that are already normalized and
divide by each new normalizer.

``` r
master<-master%>%
  mutate(mgprot = (prot_ug.cm2*surface.area.cm2)/1000)%>% #get absolute values of normaliers for each coral sample, change from ug to mg
  mutate(mgafdw = Host_AFDW.mg.cm2*surface.area.cm2)%>% #get absolute values of normaliers for each coral sample
  mutate(cells = cells.cm2 * surface.area.cm2)%>% #get absolute number of symbiont cells
  mutate(cells.mgprot = (cells.cm2 * surface.area.cm2)/mgprot)%>% #convert cells per cm2 to cells per mg protein
  mutate(cells.mgAFDW = (cells.cm2 * surface.area.cm2)/mgafdw)%>% #converts cells per cm2 to cells per mg afdw
  mutate(chla.ug.mgprot = (chla.ug.cm2 * surface.area.cm2)/mgprot)%>% #normalize chla to protein
  mutate(chla.ug.mgAFDW = (chla.ug.cm2 * surface.area.cm2)/mgafdw)%>% #normalize chla to afdw
  mutate(chlc2.ug.mgprot = (chlc2.ug.cm2 * surface.area.cm2)/mgprot)%>% #normalize chlc2 to protein
  mutate(chlc2.ug.mgAFDW = (chlc2.ug.cm2 * surface.area.cm2)/mgafdw)%>% #normalize chlc2 to afdw
  mutate(prot_mg.cm2=prot_ug.cm2/1000)%>% #change protein to mg
  mutate(prot_mg.mgafdw=(prot_mg.cm2*surface.area.cm2)/mgafdw)%>% #normalize protein to afdw
  mutate(cre.umol.mgafdw= (cre.umol.mgprot*mgprot)/mgafdw)%>% #normalize tac to afdw
  mutate(cre.umol.cm2= (cre.umol.mgprot*mgprot)/surface.area.cm2)%>% #normalize tac to surface area
  mutate(chla.ug.cell = (chla.ug.cm2 * surface.area.cm2)/cells) %>% #normalize to chl per cell
  mutate(chlc2.ug.cell = (chlc2.ug.cm2 * surface.area.cm2)/cells) %>% #normalize to chl per cell
  rename(calc.umol.cm2.hr = umol.cm2.hr)%>%
  mutate(calc.umol.mgprot.hr = (calc.umol.cm2.hr * surface.area.cm2)/mgprot)%>%
  mutate(calc.umol.mgAFDW.hr = (calc.umol.cm2.hr * surface.area.cm2)/mgafdw)
```

Calculate additional metrics of host:symbiont biomass ratio and total
cholorophyll (a + c2).

``` r
master <- master %>% 
  mutate(Ratio_AFDW.mg.cm2=Sym_AFDW.mg.cm2/(Sym_AFDW.mg.cm2+Host_AFDW.mg.cm2))%>%
  mutate(Total_Chl=chla.ug.mgAFDW+chlc2.ug.mgAFDW)%>%
  mutate(Total_Chl_cell=chla.ug.cell+chlc2.ug.cell)

head(master)
```

    ##   colony_id colony_id_corr  species  timepoint        month   site nutrient
    ## 1   ACR-139        ACR-139 Acropora timepoint1 January 2020 Manava     High
    ## 2   ACR-140        ACR-140 Acropora timepoint1 January 2020 Manava     High
    ## 3   ACR-145        ACR-145 Acropora timepoint1 January 2020 Manava     High
    ## 4   ACR-150        ACR-150 Acropora timepoint1 January 2020 Manava     High
    ## 5   ACR-165        ACR-165 Acropora timepoint1 January 2020 Manava     High
    ## 6   ACR-173        ACR-173 Acropora timepoint1 January 2020 Manava     High
    ##     site_code cre.umol.mgprot Host_AFDW.mg.cm2 Sym_AFDW.mg.cm2 Host_DW.mg.cm2
    ## 1 Manava High      0.07678208        1.4072814       0.9942749      293.15578
    ## 2 Manava High      0.16359919        0.8512491       0.8185088      147.35613
    ## 3 Manava High      0.11473069        0.8806505       0.3906111       78.49784
    ## 4 Manava High      0.08877581        1.0759520       1.2938664      128.25125
    ## 5 Manava High      0.10763769        1.6758688       0.8437534      197.68680
    ## 6 Manava High      0.12499656        1.2823276       0.9486607      260.60473
    ##   Sym_DW.mg.cm2 chla.ug.cm2 chlc2.ug.cm2 prot_ug.cm2 surface.area.cm2 cells.cm2
    ## 1      375.3623    2.788699    2.5416344    581.4959         16.99731  628007.8
    ## 2      191.1504    1.352630    0.8816772    391.1590         24.43468  403797.7
    ## 3      113.6264    1.790804    1.2134689    308.9662         12.67245  684950.4
    ## 4      161.1390    2.801127    1.4677382    557.2875         32.30627  909111.4
    ## 5      258.6335    2.770324    1.0237056    377.5935         14.60735  756468.6
    ## 6      326.5465    2.991138    3.2424896    476.5595         22.92706  800364.3
    ##          Am         AQY        Rd calc.umol.cm2.hr                 code
    ## 1 1.0423651 0.002247138 0.4110558     0.0069629397 ACR-139 - timepoint1
    ## 2 0.7395432 0.002678070 0.4154000     0.0074498259 ACR-140 - timepoint1
    ## 3 0.7092019 0.002476605 0.2058895     0.0000000000 ACR-145 - timepoint1
    ## 4 2.0342383 0.003363935 0.8398342     0.0021854881 ACR-150 - timepoint1
    ## 5 1.0504973 0.004994669 0.3808200     0.0002372078 ACR-165 - timepoint1
    ## 6 0.7474726 0.001818345 0.4458121     0.0000000000 ACR-173 - timepoint1
    ##      mgprot mgafdw    cells cells.mgprot cells.mgAFDW chla.ug.mgprot
    ## 1  9.883867  23.92 10674444      1079987     446256.0       4.795732
    ## 2  9.557846  20.80  9866667      1032311     474359.0       3.458005
    ## 3  3.915359  11.16  8680000      2216910     777777.8       5.796116
    ## 4 18.003882  34.76 29370000      1631315     844936.7       5.026359
    ## 5  5.515639  24.48 11050000      2003394     451388.9       7.336789
    ## 6 10.926109  29.40 18350000      1679463     624149.7       6.276526
    ##   chla.ug.mgAFDW chlc2.ug.mgprot chlc2.ug.mgAFDW prot_mg.cm2 prot_mg.mgafdw
    ## 1       1.981621        4.370855       1.8060598   0.5814959      0.4132052
    ## 2       1.588994        2.254012       1.0357452   0.3911590      0.4595118
    ## 3       2.033501        3.927514       1.3779234   0.3089662      0.3508386
    ## 4       2.603394        2.633718       1.3641297   0.5572875      0.5179483
    ## 5       1.653067        2.711132       0.6108507   0.3775935      0.2253121
    ## 6       2.332585        6.803955       2.5285969   0.4765595      0.3716363
    ##   cre.umol.mgafdw cre.umol.cm2 chla.ug.cell chlc2.ug.cell calc.umol.mgprot.hr
    ## 1      0.03172675   0.04464847 4.440548e-06  4.047138e-06        0.0119741844
    ## 2      0.07517576   0.06399330 3.349772e-06  2.183463e-06        0.0190455160
    ## 3      0.04025196   0.03544790 2.614502e-06  1.771616e-06        0.0000000000
    ## 4      0.04598128   0.04947365 3.081170e-06  1.614476e-06        0.0039216528
    ## 5      0.02425207   0.04064329 3.662179e-06  1.353269e-06        0.0006282095
    ## 6      0.04645326   0.05956830 3.737221e-06  4.051267e-06        0.0000000000
    ##   calc.umol.mgAFDW.hr Ratio_AFDW.mg.cm2 Total_Chl Total_Chl_cell
    ## 1        0.0049477948         0.4140127  3.787681   8.487686e-06
    ## 2        0.0087516402         0.4901961  2.624740   5.533235e-06
    ## 3        0.0000000000         0.3072626  3.411425   4.386118e-06
    ## 4        0.0020312133         0.5459770  3.967524   4.695646e-06
    ## 5        0.0001415432         0.3348730  2.263918   5.015449e-06
    ## 6        0.0000000000         0.4252199  4.861182   7.788488e-06

Subset data for D. Conetta Acropora colonies for separate analysis.

``` r
conetta_colonies<-c("ACR-140", "ACR-150", "ACR-173", "ACR-185", "ACR-187", "ACR-225", "ACR-228", "ACR-234", "ACR-237", "ACR-244", "ACR-343", "ACR-368", "ACR-379", "ACR-389", "ACR-390", "ACR-193", "ACR-180", "ACR-178", "ACR-243", "ACR-256", "ACR-374", "ACR-350", "ACR-351", "ACR-363")
conetta<-master%>%
  filter(colony_id_corr %in% conetta_colonies)%>%
  arrange(colony_id_corr)%>%
  write.csv(., "Output/conetta_data.csv")
```

Write master file to csv.

``` r
write_csv(master, "Output/master_timeseries.csv")
```

# Generating summary by groups

Summarize all response values by group (species + timepoint + site).

``` r
summary_colony<-master%>%
  select(., -c(colony_id, timepoint, site_code))%>%
  group_by(colony_id_corr, month, site, nutrient, species)%>%
  summarise(across(.cols=everything(), ~mean(.x, na.rm = TRUE)))%>%
  drop_na("site")%>%
  write_csv(., "Output/Summary_Colony_Responses.csv")

summary_group<-master%>%
  select(., -c(colony_id, timepoint, site_code, colony_id_corr))%>%
  group_by(month, site, nutrient, species)%>%
  summarise(across(.cols=everything(), ~mean(.x, na.rm = TRUE)))%>%
  drop_na("site")%>%
  write_csv(., "Output/Summary_Responses.csv")
```
