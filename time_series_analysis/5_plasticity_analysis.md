Plasticity analysis of E5 time series biological data
================
Ariana S Huffmyer, E5 RoL Team
08/22/2022

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

# Load dataframe

Load in master dataframe generated from 1_assemble_data.Rmd.

``` r
master<-read.csv("Output/master_timeseries.csv")

#reorder site levels 
master$side_code<-as.factor(master$site_code)
master$site_code<-fct_relevel(master$site_code, "Mahana Low", "Hilton Medium", "Manava High")
```

# Population - level analysis

Calculate the plasticity of physiological profiles (ratio of distance
between centroids to the average spread of points in the PCA) at the
population level produced in 4_multivariate_analysis.Rmd).

## All responses

### Acropora

Generate PCA using scaled (scaled_acr_afdw) from dataframe
(acr_data_afdw).

``` r
acr_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

acr_data_afdw<-acr_data_afdw[complete.cases(acr_data_afdw), ]

scaled_acr_afdw<-prcomp(acr_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 

acr_info<-acr_data_afdw[c(2,4)]

acr_data<-scaled_acr_afdw%>%
  augment(acr_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

acr.centroids<-acr_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate sd of mean distance between all points and the centroid of
all points (spread)

``` r
#calculate mean centroid location
mean.centroid.acr <- acr.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points
spread.acr<-acr_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.acr[1])^2+(PC2.mean-mean.centroid.acr[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #calculate the sd of distances
  summarise(distance=mean(distance)) #summarize across groups

spread.acr$distance
```

    ## [1] 0.6409587

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.acr<-acr.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.acr.all<-distance.acr%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.acr$distance)%>%
  mutate(fraction="All Responses")%>%
  mutate(species="Acropora")
plasticity.acr.all
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction      species 
    ##   <fct>         <chr>              <dbl> <dbl> <chr>         <chr>   
    ## 1 Mahana Low    tp1.2              1.85  2.88  All Responses Acropora
    ## 2 Hilton Medium tp1.2              0.364 0.567 All Responses Acropora
    ## 3 Manava High   tp1.2              0.843 1.31  All Responses Acropora
    ## 4 Mahana Low    tp2.3              2.97  4.64  All Responses Acropora
    ## 5 Hilton Medium tp2.3              2.05  3.19  All Responses Acropora
    ## 6 Manava High   tp2.3              3.65  5.69  All Responses Acropora
    ## 7 Mahana Low    tp3.4              3.70  5.77  All Responses Acropora
    ## 8 Hilton Medium tp3.4              2.42  3.77  All Responses Acropora
    ## 9 Manava High   tp3.4              3.69  5.76  All Responses Acropora

### Porites

``` r
por_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

por_data_afdw<-por_data_afdw[complete.cases(por_data_afdw), ]

scaled_por_afdw<-prcomp(por_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 

por_info<-por_data_afdw[c(2,4)]

por_data<-scaled_por_afdw%>%
  augment(por_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

por.centroids<-por_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate mean distance between all points and the centroid of all
points (spread)

``` r
#calculate mean centroid location
mean.centroid.por <- por.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points
spread.por<-por_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.por[1])^2+(PC2.mean-mean.centroid.por[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #summarize across groups
  summarise(distance=mean(distance)) #summarize across groups

#spread<-acr_data[c(12,13)]
#spread<-as.matrix(dist(spread, mean.centroid, method="euclidean"))
#spread<-mean(spread)
spread.por
```

    ## # A tibble: 1 × 1
    ##   distance
    ##      <dbl>
    ## 1    0.687

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.por<-por.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.por.all<-distance.por%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.por$distance)%>%
  mutate(fraction="All Responses")%>%
  mutate(species="Porites")
plasticity.por.all
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction      species
    ##   <fct>         <chr>              <dbl> <dbl> <chr>         <chr>  
    ## 1 Mahana Low    tp1.2              0.435 0.633 All Responses Porites
    ## 2 Hilton Medium tp1.2              2.11  3.07  All Responses Porites
    ## 3 Manava High   tp1.2              0.131 0.191 All Responses Porites
    ## 4 Mahana Low    tp2.3              0.540 0.786 All Responses Porites
    ## 5 Hilton Medium tp2.3              0.932 1.36  All Responses Porites
    ## 6 Manava High   tp2.3              0.599 0.871 All Responses Porites
    ## 7 Mahana Low    tp3.4              1.09  1.59  All Responses Porites
    ## 8 Hilton Medium tp3.4              2.50  3.63  All Responses Porites
    ## 9 Manava High   tp3.4              2.21  3.21  All Responses Porites

### Pocillopora

``` r
poc_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

poc_data_afdw<-poc_data_afdw[complete.cases(poc_data_afdw), ]
 
scaled_poc_afdw<-prcomp(poc_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 

poc_info<-poc_data_afdw[c(2,4)]

poc_data<-scaled_poc_afdw%>%
  augment(poc_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

poc.centroids<-poc_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate mean distance between all points and the centroid of all
points (spread)

``` r
#calculate mean centroid location
mean.centroid.poc <- poc.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average distance between mean centroid and location of each point using formula for distance between two points
spread.poc<-poc_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.poc[1])^2+(PC2.mean-mean.centroid.poc[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #summarize across groups
  summarise(distance=mean(distance)) #summarize across groups

#spread<-acr_data[c(12,13)]
#spread<-as.matrix(dist(spread, mean.centroid, method="euclidean"))
#spread<-mean(spread)
spread.poc
```

    ## # A tibble: 1 × 1
    ##   distance
    ##      <dbl>
    ## 1    0.622

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.poc<-poc.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.poc.all<-distance.poc%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.poc$distance)%>%
  mutate(fraction="All Responses")%>%
  mutate(species="Pocillopora")
plasticity.poc.all
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction      species    
    ##   <fct>         <chr>              <dbl> <dbl> <chr>         <chr>      
    ## 1 Mahana Low    tp1.2              1.97   3.17 All Responses Pocillopora
    ## 2 Hilton Medium tp1.2              0.936  1.51 All Responses Pocillopora
    ## 3 Manava High   tp1.2              1.05   1.68 All Responses Pocillopora
    ## 4 Mahana Low    tp2.3              2.23   3.59 All Responses Pocillopora
    ## 5 Hilton Medium tp2.3              2.89   4.64 All Responses Pocillopora
    ## 6 Manava High   tp2.3              3.35   5.39 All Responses Pocillopora
    ## 7 Mahana Low    tp3.4              1.77   2.84 All Responses Pocillopora
    ## 8 Hilton Medium tp3.4              3.71   5.97 All Responses Pocillopora
    ## 9 Manava High   tp3.4              4.45   7.16 All Responses Pocillopora

## Holobiont Responses

Generate a list of symbiont and host responses.

``` r
symbiont_responses<-c("cells.mgAFDW", "Sym_AFDW.mg.cm2", "Total_Chl", "Total_Chl_cell", "Am", "AQY", "Ratio_AFDW.mg.cm2") 
holobiont_responses<-c("cre.umol.mgafdw", "Host_AFDW.mg.cm2", "prot_mg.mgafdw", "Rd", "calc.umol.mgAFDW.hr")
```

### Acropora

``` r
acr_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5) %>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Respiration=Rd, Calc=calc.umol.mgAFDW.hr)

acr_data_afdw<-acr_data_afdw[complete.cases(acr_data_afdw), ]

scaled_acr_afdw<-prcomp(acr_data_afdw[c(5:9)], scale=TRUE, center=TRUE) 

acr_info<-acr_data_afdw[c(2,4)]

acr_data<-scaled_acr_afdw%>%
  augment(acr_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

acr.centroids<-acr_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate sd of mean distance between all points and the centroid of
all points (spread)

``` r
#calculate mean centroid location
mean.centroid.acr <- acr.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points
spread.acr<-acr_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.acr[1])^2+(PC2.mean-mean.centroid.acr[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #calculate the sd of distances
  summarise(distance=mean(distance)) #summarize across groups

spread.acr$distance
```

    ## [1] 0.4972175

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.acr<-acr.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.acr.holo<-distance.acr%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.acr$distance)%>%
  mutate(species="Acropora")%>%
  mutate(fraction="Holobiont")
plasticity.acr.holo
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio species  fraction 
    ##   <fct>         <chr>              <dbl> <dbl> <chr>    <chr>    
    ## 1 Mahana Low    tp1.2              0.434 0.872 Acropora Holobiont
    ## 2 Hilton Medium tp1.2              0.388 0.779 Acropora Holobiont
    ## 3 Manava High   tp1.2              0.481 0.968 Acropora Holobiont
    ## 4 Mahana Low    tp2.3              2.69  5.42  Acropora Holobiont
    ## 5 Hilton Medium tp2.3              1.36  2.73  Acropora Holobiont
    ## 6 Manava High   tp2.3              2.57  5.17  Acropora Holobiont
    ## 7 Mahana Low    tp3.4              2.90  5.82  Acropora Holobiont
    ## 8 Hilton Medium tp3.4              1.90  3.82  Acropora Holobiont
    ## 9 Manava High   tp3.4              3.21  6.45  Acropora Holobiont

### Porites

``` r
por_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Porites")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Respiration=Rd, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

por_data_afdw<-por_data_afdw[complete.cases(por_data_afdw), ]

scaled_por_afdw<-prcomp(por_data_afdw[c(5:9)], scale=TRUE, center=TRUE) 

por_info<-por_data_afdw[c(2,4)]

por_data<-scaled_por_afdw%>%
  augment(por_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

por.centroids<-por_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate mean distance between all points and the centroid of all
points (spread)

``` r
#calculate mean centroid location
mean.centroid.por <- por.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points
spread.por<-por_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.por[1])^2+(PC2.mean-mean.centroid.por[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #summarize across groups
  summarise(distance=mean(distance)) #summarize across groups

#spread<-acr_data[c(12,13)]
#spread<-as.matrix(dist(spread, mean.centroid, method="euclidean"))
#spread<-mean(spread)
spread.por
```

    ## # A tibble: 1 × 1
    ##   distance
    ##      <dbl>
    ## 1    0.458

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.por<-por.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.por.holo<-distance.por%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.por$distance)%>%
  mutate(fraction="Holobiont")%>%
  mutate(species="Porites")
plasticity.por.holo
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction  species
    ##   <fct>         <chr>              <dbl> <dbl> <chr>     <chr>  
    ## 1 Mahana Low    tp1.2              0.886  1.93 Holobiont Porites
    ## 2 Hilton Medium tp1.2              2.29   5.00 Holobiont Porites
    ## 3 Manava High   tp1.2              0.710  1.55 Holobiont Porites
    ## 4 Mahana Low    tp2.3              0.649  1.42 Holobiont Porites
    ## 5 Hilton Medium tp2.3              0.538  1.17 Holobiont Porites
    ## 6 Manava High   tp2.3              0.713  1.56 Holobiont Porites
    ## 7 Mahana Low    tp3.4              0.778  1.70 Holobiont Porites
    ## 8 Hilton Medium tp3.4              1.89   4.12 Holobiont Porites
    ## 9 Manava High   tp3.4              2.13   4.65 Holobiont Porites

### Pocillopora

``` r
poc_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Respiration=Rd, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

poc_data_afdw<-poc_data_afdw[complete.cases(poc_data_afdw), ]
 
scaled_poc_afdw<-prcomp(poc_data_afdw[c(5:9)], scale=TRUE, center=TRUE) 

poc_info<-poc_data_afdw[c(2,4)]

poc_data<-scaled_poc_afdw%>%
  augment(poc_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

poc.centroids<-poc_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate mean distance between all points and the centroid of all
points (spread)

``` r
#calculate mean centroid location
mean.centroid.poc <- poc.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average distance between mean centroid and location of each point using formula for distance between two points
spread.poc<-poc_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.poc[1])^2+(PC2.mean-mean.centroid.poc[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #summarize across groups
  summarise(distance=mean(distance)) #summarize across groups

spread.poc
```

    ## # A tibble: 1 × 1
    ##   distance
    ##      <dbl>
    ## 1    0.489

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.poc<-poc.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.poc.holo<-distance.poc%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.poc$distance)%>%
  mutate(fraction="Holobiont")%>%
  mutate(species="Pocillopora")
plasticity.poc.holo
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction  species    
    ##   <fct>         <chr>              <dbl> <dbl> <chr>     <chr>      
    ## 1 Mahana Low    tp1.2              0.616  1.26 Holobiont Pocillopora
    ## 2 Hilton Medium tp1.2              0.848  1.73 Holobiont Pocillopora
    ## 3 Manava High   tp1.2              1.60   3.28 Holobiont Pocillopora
    ## 4 Mahana Low    tp2.3              0.886  1.81 Holobiont Pocillopora
    ## 5 Hilton Medium tp2.3              1.47   3.02 Holobiont Pocillopora
    ## 6 Manava High   tp2.3              2.33   4.77 Holobiont Pocillopora
    ## 7 Mahana Low    tp3.4              0.846  1.73 Holobiont Pocillopora
    ## 8 Hilton Medium tp3.4              2.56   5.25 Holobiont Pocillopora
    ## 9 Manava High   tp3.4              2.87   5.87 Holobiont Pocillopora

## Symbiont Responses

### Acropora

``` r
acr_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell, SH_Ratio=Ratio_AFDW.mg.cm2)

acr_data_afdw<-acr_data_afdw[complete.cases(acr_data_afdw), ]

scaled_acr_afdw<-prcomp(acr_data_afdw[c(5:11)], scale=TRUE, center=TRUE) 

acr_info<-acr_data_afdw[c(2,4)]

acr_data<-scaled_acr_afdw%>%
  augment(acr_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

acr.centroids<-acr_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate sd of mean distance between all points and the centroid of
all points (spread)

``` r
#calculate mean centroid location
mean.centroid.acr <- acr.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points
spread.acr<-acr_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.acr[1])^2+(PC2.mean-mean.centroid.acr[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #calculate the sd of distances
  summarise(distance=mean(distance)) #summarize across groups

spread.acr$distance
```

    ## [1] 0.4983495

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.acr<-acr.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.acr.sym<-distance.acr%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.acr$distance)%>%
  mutate(fraction="Symbiont")%>%
  mutate(species="Acropora")
plasticity.acr.sym
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction species 
    ##   <fct>         <chr>              <dbl> <dbl> <chr>    <chr>   
    ## 1 Mahana Low    tp1.2              1.98   3.97 Symbiont Acropora
    ## 2 Hilton Medium tp1.2              0.529  1.06 Symbiont Acropora
    ## 3 Manava High   tp1.2              0.924  1.85 Symbiont Acropora
    ## 4 Mahana Low    tp2.3              2.18   4.37 Symbiont Acropora
    ## 5 Hilton Medium tp2.3              1.92   3.85 Symbiont Acropora
    ## 6 Manava High   tp2.3              2.90   5.83 Symbiont Acropora
    ## 7 Mahana Low    tp3.4              2.20   4.42 Symbiont Acropora
    ## 8 Hilton Medium tp3.4              1.42   2.85 Symbiont Acropora
    ## 9 Manava High   tp3.4              2.16   4.33 Symbiont Acropora

### Porites

``` r
por_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Porites")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell, SH_Ratio=Ratio_AFDW.mg.cm2)

por_data_afdw<-por_data_afdw[complete.cases(por_data_afdw), ]

scaled_por_afdw<-prcomp(por_data_afdw[c(5:11)], scale=TRUE, center=TRUE) 

por_info<-por_data_afdw[c(2,4)]

por_data<-scaled_por_afdw%>%
  augment(por_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

por.centroids<-por_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate mean distance between all points and the centroid of all
points (spread)

``` r
#calculate mean centroid location
mean.centroid.por <- por.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points
spread.por<-por_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.por[1])^2+(PC2.mean-mean.centroid.por[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #summarize across groups
  summarise(distance=mean(distance)) #summarize across groups

spread.por
```

    ## # A tibble: 1 × 1
    ##   distance
    ##      <dbl>
    ## 1    0.471

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.por<-por.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.por.sym<-distance.por%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.por$distance)%>%
  mutate(fraction="Symbiont")%>%
  mutate(species="Porites")
plasticity.por.sym
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction species
    ##   <fct>         <chr>              <dbl> <dbl> <chr>    <chr>  
    ## 1 Mahana Low    tp1.2              1.50   3.18 Symbiont Porites
    ## 2 Hilton Medium tp1.2              0.912  1.94 Symbiont Porites
    ## 3 Manava High   tp1.2              1.58   3.36 Symbiont Porites
    ## 4 Mahana Low    tp2.3              1.20   2.55 Symbiont Porites
    ## 5 Hilton Medium tp2.3              1.03   2.19 Symbiont Porites
    ## 6 Manava High   tp2.3              0.656  1.39 Symbiont Porites
    ## 7 Mahana Low    tp3.4              1.12   2.38 Symbiont Porites
    ## 8 Hilton Medium tp3.4              1.68   3.56 Symbiont Porites
    ## 9 Manava High   tp3.4              1.09   2.32 Symbiont Porites

### Pocillopora

``` r
poc_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Pocillopora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell, SH_Ratio=Ratio_AFDW.mg.cm2)

poc_data_afdw<-poc_data_afdw[complete.cases(poc_data_afdw), ]
 
scaled_poc_afdw<-prcomp(poc_data_afdw[c(5:11)], scale=TRUE, center=TRUE) 

poc_info<-poc_data_afdw[c(2,4)]

poc_data<-scaled_poc_afdw%>%
  augment(poc_info)%>%
  group_by(timepoint, site_code)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

poc.centroids<-poc_data %>% 
  select(timepoint, site_code, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate mean distance between all points and the centroid of all
points (spread)

``` r
#calculate mean centroid location
mean.centroid.poc <- poc.centroids%>%
  select(PC1.mean, PC2.mean)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean))%>%
  summarise(x.mean = mean(x.mean), 
         y.mean = mean(y.mean))

#calculate average distance between mean centroid and location of each point using formula for distance between two points
spread.poc<-poc_data%>%
  mutate(distance=sqrt((PC1.mean-mean.centroid.poc[1])^2+(PC2.mean-mean.centroid.poc[2])^2))%>%
  summarise(distance=mean(distance$x.mean))%>% #summarize across groups
  summarise(distance=sd(distance))%>% #summarize across groups
  summarise(distance=mean(distance)) #summarize across groups

spread.poc
```

    ## # A tibble: 1 × 1
    ##   distance
    ##      <dbl>
    ## 1    0.278

2 - calculate distance between each time point centroids for each site
(distance)

``` r
distance.poc<-poc.centroids%>%
  arrange(site_code)%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(site_code)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread, generating a plasticity ratio for each
site, higher ratio = more plasticity

``` r
plasticity.poc.sym<-distance.poc%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.poc$distance)%>%
  mutate(fraction="Symbiont")%>%
  mutate(species="Pocillopora")
plasticity.poc.sym
```

    ## # A tibble: 9 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio fraction species    
    ##   <fct>         <chr>              <dbl> <dbl> <chr>    <chr>      
    ## 1 Mahana Low    tp1.2              2.16   7.77 Symbiont Pocillopora
    ## 2 Hilton Medium tp1.2              0.604  2.17 Symbiont Pocillopora
    ## 3 Manava High   tp1.2              1.54   5.52 Symbiont Pocillopora
    ## 4 Mahana Low    tp2.3              1.92   6.91 Symbiont Pocillopora
    ## 5 Hilton Medium tp2.3              2.57   9.24 Symbiont Pocillopora
    ## 6 Manava High   tp2.3              2.64   9.47 Symbiont Pocillopora
    ## 7 Mahana Low    tp3.4              1.31   4.71 Symbiont Pocillopora
    ## 8 Hilton Medium tp3.4              2.59   9.31 Symbiont Pocillopora
    ## 9 Manava High   tp3.4              2.49   8.93 Symbiont Pocillopora

## Generate figure

Bind all data frames together.

``` r
plasticity.pop<-rbind(plasticity.acr.holo, plasticity.por.holo, plasticity.poc.holo, plasticity.acr.sym, plasticity.poc.sym, plasticity.por.sym)
head(plasticity.pop)
```

    ## # A tibble: 6 × 6
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio species  fraction 
    ##   <fct>         <chr>              <dbl> <dbl> <chr>    <chr>    
    ## 1 Mahana Low    tp1.2              0.434 0.872 Acropora Holobiont
    ## 2 Hilton Medium tp1.2              0.388 0.779 Acropora Holobiont
    ## 3 Manava High   tp1.2              0.481 0.968 Acropora Holobiont
    ## 4 Mahana Low    tp2.3              2.69  5.42  Acropora Holobiont
    ## 5 Hilton Medium tp2.3              1.36  2.73  Acropora Holobiont
    ## 6 Manava High   tp2.3              2.57  5.17  Acropora Holobiont

``` r
#plasticity.acr.all, plasticity.poc.all, plasticity.por.all

#reorder site levels 
plasticity.pop$side_code<-as.factor(plasticity.pop$site_code)
plasticity.pop$site_code<-fct_relevel(plasticity.pop$site_code, "Mahana Low", "Hilton Medium", "Manava High")
```

Plot plasticity values with fill for fraction (all responses, holobiont,
or symbiont), site on x axis in increasing nutrient exposure, and
faceted by species.

First, plot panel that shows fraction/species at each site.

``` r
figure1a<-plasticity.pop%>% 
    group_by(species, site_code, fraction)%>%
    summarise(mean=mean(ratio))%>%
    ggplot(aes(x = site_code, y = mean, fill=fraction, group=interaction(species, site_code))) +
    facet_wrap(~species)+
    geom_point(pch = 21, size=5, position = position_jitterdodge(0.5), color="black") + 
    scale_fill_manual(values = c("darkgray", "white"))+
    xlab("Site") + 
    ylab(expression(bold("Plasticity Score")))+
    theme_classic() + 
    #ylim(0.5,3)+
    theme(
      legend.position="right",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); figure1a
```

![](5_plasticity_analysis_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

Second, plot panel that shows fraction/species by time point comparison.

``` r
figure1b<-plasticity.pop%>% 
    group_by(species, comparison, fraction)%>%
    summarise(mean=mean(ratio))%>%
    ggplot(aes(x = comparison, y = mean, fill=fraction, group=interaction(species, comparison))) +
    facet_wrap(~species)+
    geom_point(pch = 21, size=5, position = position_jitterdodge(0.5), color="black") + 
    scale_fill_manual(values = c("darkgray", "white"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    xlab("Time") + 
    ylab(expression(bold("Plasticity Score")))+
    theme_classic() + 
    #ylim(0.5,3)+
    theme(
      legend.position="right",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); figure1b
```

![](5_plasticity_analysis_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

Export figure.

``` r
plasticity_grid<-plot_grid(figure1a, figure1b, ncol=1, nrow=2)

ggsave(filename="Figures/Plasiticity_Figure.pdf", plot=plasticity_grid, dpi=500, width=14, height=8, units="in")
```

# Colony - level analysis

## All responses

Generate centroid locations and spread for each colony as done in script
4_multivariate_analysis.Rmd. This analysis is completed at the level of
the colony instead of the population level.

``` r
data_all<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr, Total_Chl, Total_Chl_cell)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

data_all<-data_all[complete.cases(data_all), ]

scaled_all<-prcomp(data_all[c(5:16)], scale=TRUE, center=TRUE) 

all_info<-data_all[c(1:4)]

all_data<-scaled_all%>%
  augment(all_info)%>%
  group_by(timepoint, site_code, species, colony_id_corr)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

all.centroids<-all_data %>% 
  select(timepoint, site_code, species, colony_id_corr, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code, species, colony_id_corr)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate distance between all points and the centroid of all points
for each colony (spread)

``` r
#calculate mean centroid location
mean.centroid.all <- all.centroids%>%
  group_by(colony_id_corr)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean)) 

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points

all_data$x.mean<-mean.centroid.all$x.mean[match(all_data$colony_id_corr, mean.centroid.all$colony_id_corr)]
all_data$y.mean<-mean.centroid.all$y.mean[match(all_data$colony_id_corr, mean.centroid.all$colony_id_corr)]

#calculate spread as the square root of squared distances between mean PC location of each colony measurement and average PC location for each colony
all_data<-all_data%>%
  mutate(spread=sqrt((PC1.mean-x.mean)^2+(PC2.mean-y.mean)^2))

#calculate average spread for each colony
spread.all<-all_data%>%
  group_by(colony_id_corr)%>%
  summarise(spread.mean=mean(spread), spread.sd=sd(spread), n=length(spread))
```

2 - calculate distance between each time point centroid colony
(distance)

``` r
distance.all<-all.centroids%>%
  arrange(timepoint)%>%
  gather(variable, value, -(timepoint:colony_id_corr)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(colony_id_corr)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2), 
         tp1.3=sqrt((timepoint1_PC1.mean-timepoint3_PC1.mean)^2+(timepoint1_PC2.mean-timepoint3_PC2.mean)^2), 
         tp1.4=sqrt((timepoint1_PC1.mean-timepoint4_PC1.mean)^2+(timepoint1_PC2.mean-timepoint4_PC2.mean)^2), 
         tp2.4=sqrt((timepoint2_PC1.mean-timepoint4_PC1.mean)^2+(timepoint2_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread for each colony, generating a plasticity
score. The higher score = more distance in physiological trajectory
relative to variability (spread).

``` r
#remove colonies that were only observed at one time point and therefore cannot be calculated due to no distance measurements 
plasticity.all<-left_join(distance.all, spread.all)%>%
  #select(site_code, colony_id_corr, species, tp1.2, tp2.3, tp3.4, tp1.3, tp1.4, tp2.4, spread.mean)%>% #keep desired columns
  select(site_code, colony_id_corr, species, tp1.2, tp2.3, tp3.4, spread.mean)%>% #keep desired columns
  relocate(spread.mean, .after=species)%>%
  filter(!spread.mean==0)%>% 
  mutate(distance.mean = mean(c_across(starts_with("tp")), na.rm = TRUE))%>%
  mutate(plasticity=distance.mean/spread.mean)
```

4 - average these scores across colony to generate a mean plasticity
ratio for each species and site with std. error.

Display mean plasticity scores by site:

``` r
all.plasticity.stats<-plasticity.all%>%
  filter(!is.na(plasticity))%>%
  group_by(site_code, species)%>%
  summarise(plasticity.mean=mean(plasticity, na.rm=TRUE), plasticity.sd=sd(plasticity, na.rm=TRUE), n=length(plasticity), plasticity.se=plasticity.sd/sqrt(n))%>%
  mutate(fraction="All Responses")
all.plasticity.stats
```

    ## # A tibble: 9 × 7
    ## # Groups:   site_code [3]
    ##   site_code   species plasticity.mean plasticity.sd     n plasticity.se fraction
    ##   <fct>       <chr>             <dbl>         <dbl> <int>         <dbl> <chr>   
    ## 1 Mahana Low  Acropo…            1.68         0.416     8        0.147  All Res…
    ## 2 Mahana Low  Pocill…            1.65         0.405    11        0.122  All Res…
    ## 3 Mahana Low  Porites            1.80         0.336    13        0.0931 All Res…
    ## 4 Hilton Med… Acropo…            1.37         0.602     7        0.228  All Res…
    ## 5 Hilton Med… Pocill…            1.91         0.461    13        0.128  All Res…
    ## 6 Hilton Med… Porites            1.68         0.298    15        0.0770 All Res…
    ## 7 Manava High Acropo…            2.01         0.570     9        0.190  All Res…
    ## 8 Manava High Pocill…            1.88         0.410    13        0.114  All Res…
    ## 9 Manava High Porites            1.75         0.427    13        0.118  All Res…

## Holobiont responses

Generate centroid locations and spread for each colony as done in script
4_multivariate_analysis.Rmd. This analysis is completed at the level of
the colony instead of the population level.

Generate a list of symbiont and holobiont responses.

``` r
symbiont_responses<-c("cells.mgAFDW", "Sym_AFDW.mg.cm2", "Total_Chl", "Total_Chl_cell", "Am", "AQY", "Ratio_AFDW.mg.cm2") 
holobiont_responses<-c("cre.umol.mgafdw", "Host_AFDW.mg.cm2", "prot_mg.mgafdw", "Rd", "calc.umol.mgAFDW.hr")
```

``` r
data_holo<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

data_holo<-data_holo[complete.cases(data_holo), ]

scaled_holo<-prcomp(data_holo[c(5:9)], scale=TRUE, center=TRUE) 

holo_info<-data_holo[c(1:4)]

holo_data<-scaled_holo%>%
  augment(holo_info)%>%
  group_by(timepoint, site_code, species, colony_id_corr)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

holo.centroids<-holo_data %>% 
  select(timepoint, site_code, species, colony_id_corr, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code, species, colony_id_corr)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate distance between all points and the centroid of all points
for each colony (spread)

``` r
#calculate mean centroid location
mean.centroid.holo <- holo.centroids%>%
  group_by(colony_id_corr)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean)) 

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points

holo_data$x.mean<-mean.centroid.holo$x.mean[match(holo_data$colony_id_corr, mean.centroid.holo$colony_id_corr)]
holo_data$y.mean<-mean.centroid.holo$y.mean[match(holo_data$colony_id_corr, mean.centroid.holo$colony_id_corr)]

#calculate spread as the square root of squared distances between mean PC location of each colony measurement and average PC location for each colony
holo_data<-holo_data%>%
  mutate(spread=sqrt((PC1.mean-x.mean)^2+(PC2.mean-y.mean)^2))

#calculate average spread for each colony
spread.holo<-holo_data%>%
  group_by(colony_id_corr)%>%
  summarise(spread.mean=mean(spread), spread.sd=sd(spread), n=length(spread))
```

2 - calculate distance between each time point centroid colony
(distance)

``` r
distance.holo<-holo.centroids%>%
  arrange(timepoint)%>%
  gather(variable, value, -(timepoint:colony_id_corr)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(colony_id_corr)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2), 
         tp1.3=sqrt((timepoint1_PC1.mean-timepoint3_PC1.mean)^2+(timepoint1_PC2.mean-timepoint3_PC2.mean)^2), 
         tp1.4=sqrt((timepoint1_PC1.mean-timepoint4_PC1.mean)^2+(timepoint1_PC2.mean-timepoint4_PC2.mean)^2), 
         tp2.4=sqrt((timepoint2_PC1.mean-timepoint4_PC1.mean)^2+(timepoint2_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread for each colony, generating a plasticity
score. The higher score = more distance in physiological trajectory
relative to variability (spread).

``` r
#remove colonies that were only observed at one time point and therefore cannot be calculated due to no distance measurements 
plasticity.holo<-left_join(distance.holo, spread.holo)%>%
  #select(site_code, colony_id_corr, species, tp1.2, tp2.3, tp3.4, tp1.3, tp1.4, tp2.4, spread.mean)%>% #keep desired columns
  select(site_code, colony_id_corr, species, tp1.2, tp2.3, tp3.4, spread.mean)%>% #keep desired columns
  relocate(spread.mean, .after=species)%>%
  filter(!spread.mean==0)%>% 
  mutate(distance.mean = mean(c_across(starts_with("tp")), na.rm = TRUE))%>%
  mutate(plasticity=distance.mean/spread.mean)
```

4 - average these scores across colony to generate a mean plasticity
ratio for each species and site with std. error.

Display mean plasticity scores by site:

``` r
holo.plasticity.stats<-plasticity.holo%>%
  filter(!is.na(plasticity))%>%
  group_by(site_code, species)%>%
  summarise(plasticity.mean=mean(plasticity, na.rm=TRUE), plasticity.sd=sd(plasticity, na.rm=TRUE), n=length(plasticity), plasticity.se=plasticity.sd/sqrt(n))%>%
  mutate(fraction="Holobiont")
holo.plasticity.stats
```

    ## # A tibble: 9 × 7
    ## # Groups:   site_code [3]
    ##   site_code   species plasticity.mean plasticity.sd     n plasticity.se fraction
    ##   <fct>       <chr>             <dbl>         <dbl> <int>         <dbl> <chr>   
    ## 1 Mahana Low  Acropo…            1.95         0.107     8        0.0377 Holobio…
    ## 2 Mahana Low  Pocill…            1.69         0.322    13        0.0893 Holobio…
    ## 3 Mahana Low  Porites            1.76         0.328    13        0.0911 Holobio…
    ## 4 Hilton Med… Acropo…            1.60         0.362     7        0.137  Holobio…
    ## 5 Hilton Med… Pocill…            1.69         0.499    15        0.129  Holobio…
    ## 6 Hilton Med… Porites            1.68         0.308    15        0.0796 Holobio…
    ## 7 Manava High Acropo…            2.03         0.507     9        0.169  Holobio…
    ## 8 Manava High Pocill…            1.80         0.502    14        0.134  Holobio…
    ## 9 Manava High Porites            1.70         0.369    13        0.102  Holobio…

## Symbiont responses

Generate centroid locations and spread for each colony as done in script
4_multivariate_analysis.Rmd. This analysis is completed at the level of
the colony instead of the population level.

``` r
data_sym<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW)

data_sym<-data_sym[complete.cases(data_sym), ]

scaled_sym<-prcomp(data_sym[c(5:11)], scale=TRUE, center=TRUE) 

sym_info<-data_sym[c(1:4)]

sym_data<-scaled_sym%>%
  augment(sym_info)%>%
  group_by(timepoint, site_code, species, colony_id_corr)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

sym.centroids<-sym_data %>% 
  select(timepoint, site_code, species, colony_id_corr, PC1.mean, PC2.mean)%>%
  group_by(timepoint, site_code, species, colony_id_corr)%>%
  summarise(PC1.mean = mean(PC1.mean),
         PC2.mean = mean(PC2.mean))
```

1 - calculate distance between all points and the centroid of all points
for each colony (spread)

``` r
#calculate mean centroid location
mean.centroid.sym <- sym.centroids%>%
  group_by(colony_id_corr)%>%
  summarise(x.mean = mean(PC1.mean), 
         y.mean = mean(PC2.mean)) 

#calculate average standard deviation of mean distance between mean centroid and location of each point using formula for distance between two points

sym_data$x.mean<-mean.centroid.sym$x.mean[match(sym_data$colony_id_corr, mean.centroid.sym$colony_id_corr)]
sym_data$y.mean<-mean.centroid.sym$y.mean[match(sym_data$colony_id_corr, mean.centroid.sym$colony_id_corr)]

#calculate spread as the square root of squared distances between mean PC location of each colony measurement and average PC location for each colony
sym_data<-sym_data%>%
  mutate(spread=sqrt((PC1.mean-x.mean)^2+(PC2.mean-y.mean)^2))

#calculate average spread for each colony
spread.sym<-sym_data%>%
  group_by(colony_id_corr)%>%
  summarise(spread.mean=mean(spread), spread.sd=sd(spread), n=length(spread))
```

2 - calculate distance between each time point centroid colony
(distance)

``` r
distance.sym<-sym.centroids%>%
  arrange(timepoint)%>%
  gather(variable, value, -(timepoint:colony_id_corr)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)%>%
  group_by(colony_id_corr)%>%
  mutate(tp1.2=sqrt((timepoint1_PC1.mean-timepoint2_PC1.mean)^2+(timepoint1_PC2.mean-timepoint2_PC2.mean)^2), #calculate distance between tp1 and tp2 centroids, do for each pair of time points
         tp2.3=sqrt((timepoint2_PC1.mean-timepoint3_PC1.mean)^2+(timepoint2_PC2.mean-timepoint3_PC2.mean)^2),
         tp3.4=sqrt((timepoint3_PC1.mean-timepoint4_PC1.mean)^2+(timepoint3_PC2.mean-timepoint4_PC2.mean)^2), 
         tp1.3=sqrt((timepoint1_PC1.mean-timepoint3_PC1.mean)^2+(timepoint1_PC2.mean-timepoint3_PC2.mean)^2), 
         tp1.4=sqrt((timepoint1_PC1.mean-timepoint4_PC1.mean)^2+(timepoint1_PC2.mean-timepoint4_PC2.mean)^2), 
         tp2.4=sqrt((timepoint2_PC1.mean-timepoint4_PC1.mean)^2+(timepoint2_PC2.mean-timepoint4_PC2.mean)^2))
```

3 - divide distance by spread for each colony, generating a plasticity
score. The higher score = more distance in physiological trajectory
relative to variability (spread).

``` r
#remove colonies that were only observed at one time point and therefore cannot be calculated due to no distance measurements 
plasticity.sym<-left_join(distance.sym, spread.sym)%>%
  #select(site_code, colony_id_corr, species, tp1.2, tp2.3, tp3.4, tp1.3, tp1.4, tp2.4, spread.mean)%>% #keep desired columns
  select(site_code, colony_id_corr, species, tp1.2, tp2.3, tp3.4, spread.mean)%>% #keep desired columns
  relocate(spread.mean, .after=species)%>%
  filter(!spread.mean==0)%>% 
  mutate(distance.mean = mean(c_across(starts_with("tp")), na.rm = TRUE))%>%
  mutate(plasticity=distance.mean/spread.mean)
```

4 - average these scores across colony to generate a mean plasticity
ratio for each species and site with std. error.

Display mean plasticity scores by site:

``` r
sym.plasticity.stats<-plasticity.sym%>%
  filter(!is.na(plasticity))%>%
  group_by(site_code, species)%>%
  summarise(plasticity.mean=mean(plasticity, na.rm=TRUE), plasticity.sd=sd(plasticity, na.rm=TRUE), n=length(plasticity), plasticity.se=plasticity.sd/sqrt(n))%>%
  mutate(fraction="Symbiont")
sym.plasticity.stats
```

    ## # A tibble: 9 × 7
    ## # Groups:   site_code [3]
    ##   site_code   species plasticity.mean plasticity.sd     n plasticity.se fraction
    ##   <fct>       <chr>             <dbl>         <dbl> <int>         <dbl> <chr>   
    ## 1 Mahana Low  Acropo…            1.76         0.439    11        0.132  Symbiont
    ## 2 Mahana Low  Pocill…            1.68         0.351    13        0.0974 Symbiont
    ## 3 Mahana Low  Porites            1.62         0.295    14        0.0790 Symbiont
    ## 4 Hilton Med… Acropo…            1.48         0.553     7        0.209  Symbiont
    ## 5 Hilton Med… Pocill…            2.05         0.328    14        0.0876 Symbiont
    ## 6 Hilton Med… Porites            1.76         0.354    15        0.0914 Symbiont
    ## 7 Manava High Acropo…            1.71         0.365     9        0.122  Symbiont
    ## 8 Manava High Pocill…            1.95         0.365    14        0.0976 Symbiont
    ## 9 Manava High Porites            1.68         0.437    15        0.113  Symbiont

## Generate figure

Our data frames of interest are:

`all.plasticity.stats` `holo.plasticity.stats` `symb.plasticity.stats`

Join data frames.

``` r
df <- rbind(all.plasticity.stats, holo.plasticity.stats, sym.plasticity.stats)
```

Generate figure faceted by fraction with site on x axis and species in
color.

``` r
figure2<-df%>% 
    group_by(species, site_code, fraction)%>%
    summarise(mean=mean(plasticity.mean), se=mean(plasticity.se))%>%
    ggplot(aes(x = site_code, y = mean, fill=species, group=interaction(species, site_code))) +
    facet_wrap(~fraction)+
    geom_point(pch = 21, size=5, position = position_dodge(0.15)) + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=interaction(species, site_code)), width=0, color="black", position=position_dodge(0.15))+
    scale_fill_manual(values = c("darkgray", "orange", "purple"))+
    xlab("Site") + 
    ylab(expression(bold("Plasticity Score (Distance / Spread)")))+
    theme_classic() + 
    ylim(0.5,3)+
    theme(
      legend.position="right",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); figure2
```

![](5_plasticity_analysis_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->
