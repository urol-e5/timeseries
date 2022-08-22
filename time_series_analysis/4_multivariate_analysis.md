Multivariate analysis of E5 time series biological data
================
Ariana S Huffmyer, E5 RoL Team
05/17/2022

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
```

# Build trajectory PCA plots

## All Responses

### PCA’s for centroid and trajectories

#### Acropora

Generate PCA using scaled (scaled_acr_afdw) from dataframe
(acr_data_afdw).

``` r
acr_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

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

Examine PERMANOVA results.

``` r
# scale data
vegan <- scale(acr_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = acr_data_afdw, method='eu')
z_acr<-permanova$aov.tab
z_acr
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    334.05 111.350 14.0024 0.27292  0.001 ***
    ## site_code             2     82.97  41.483  5.2166 0.06778  0.001 ***
    ## timepoint:site_code   6     83.33  13.889  1.7465 0.06808  0.011 *  
    ## Residuals            91    723.65   7.952         0.59122           
    ## Total               102   1224.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
acrPCA<-ggplot(acr_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylim(-6,10)+
  xlim(-5,11)+
  ylab("PC2")+
  xlab("PC1")+
  ggtitle(expression(italic("Acropora")))+
  geom_text(x=7, y=-2.5, label=paste("p(Timepoint)=", z_acr$`Pr(>F)`[1]), size=4, color=ifelse(z_acr$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=7, y=-3.25, label=paste("p(Site)=", z_acr$`Pr(>F)`[2]), size=4, color=ifelse(z_acr$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=7, y=-4, label=paste("p(Timepoint x Site)=", z_acr$`Pr(>F)`[3]), size=4, color=ifelse(z_acr$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         title = element_text(size=25, face="bold"), 
         axis.title = element_text(size=18));acrPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
acrPCAcen<-acrPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=acr.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); acrPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Add segments

``` r
#3. add segments
segpoints<-acr.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
acrPCAfull<-acrPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); acrPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Add bi plot with loadings

``` r
#1. make plot with dots
acrArrows<-ggplot2::autoplot(scaled_acr_afdw, data=acr_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=0.5, loadings.label.hjust=-0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  #scale_colour_manual(values=c("blue4", "darkgray", "springgreen4")) +
  #scale_fill_manual(values=c("blue4", "darkgray", "springgreen4")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         title = element_text(size=25, face="bold"),
         axis.title = element_blank());acrArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Asemble final figure with inset biplot.

``` r
#acrFullPCA<-ggdraw(acrPCAfull) + #theme_half_open(12)) +
  #draw_plot(acrArrows, .5, .5, .5, .5); acrFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Acropora.pdf", plot=acrFullPCA, dpi=300, width=12, height=8, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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

#spread<-acr_data[c(12,13)]
#spread<-as.matrix(dist(spread, mean.centroid, method="euclidean"))
#spread<-mean(spread)
spread.acr$distance
```

    ## [1] 0.6828306

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
plasticity.acr<-distance.acr%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.acr$distance)
plasticity.acr
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              0.463 0.678
    ## 2 Mahana Low    tp1.2              1.92  2.82 
    ## 3 Manava High   tp1.2              0.867 1.27 
    ## 4 Hilton Medium tp2.3              1.07  1.56 
    ## 5 Mahana Low    tp2.3              2.04  2.99 
    ## 6 Manava High   tp2.3              3.20  4.69 
    ## 7 Hilton Medium tp3.4              2.28  3.34 
    ## 8 Mahana Low    tp3.4              3.53  5.17 
    ## 9 Manava High   tp3.4              3.81  5.58

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
acr.plasticity.site<-plasticity.acr%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Acropora")
acr.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species 
    ##   <chr>             <dbl> <chr>   
    ## 1 Hilton Medium      1.86 Acropora
    ## 2 Mahana Low         3.66 Acropora
    ## 3 Manava High        3.85 Acropora

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.acr.time<-plasticity.acr%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Acropora")
plasticity.acr.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species 
    ##   <chr>               <dbl> <chr>   
    ## 1 tp1.2                1.59 Acropora
    ## 2 tp2.3                3.08 Acropora
    ## 3 tp3.4                4.69 Acropora

Display mean plasticity score for the species:

``` r
plasticity.acr.species<-plasticity.acr%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Acropora")
plasticity.acr.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species 
    ##          <dbl> <chr>   
    ## 1         3.12 Acropora

#### Porites

Generate PCA using scaled (scaled_por_afdw) from dataframe
(por_data_afdw).

Calculate centroid locations for each site and timepoint and pull PC1
and PC2 locations.

``` r
por_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

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

``` r
# scale data
vegan <- scale(por_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = por_data_afdw, method='eu')
z_por<-permanova$aov.tab
z_por
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    287.88  95.960 11.2516 0.17511  0.001 ***
    ## site_code             2    192.99  96.495 11.3143 0.11739  0.001 ***
    ## timepoint:site_code   6     88.54  14.756  1.7302 0.05385  0.012 *  
    ## Residuals           126   1074.60   8.529         0.65365           
    ## Total               137   1644.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
porPCA<-ggplot(por_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylab("PC2")+
  xlab("PC1")+
  ylim(-6,10)+
  xlim(-5,11)+
  ggtitle(expression(italic("Porites")))+
  geom_text(x=7, y=-2.5, label=paste("p(Timepoint)=", z_por$`Pr(>F)`[1]), size=4, color=ifelse(z_por$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=7, y=-3.25, label=paste("p(Site)=", z_por$`Pr(>F)`[2]), size=4, color=ifelse(z_por$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=7, y=-4, label=paste("p(Timepoint x Site)=", z_por$`Pr(>F)`[3]), size=4, color=ifelse(z_por$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         legend.title = element_blank(), 
         axis.text = element_text(size=18), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         title = element_text(size=25, face="bold"),
         axis.title = element_text(size=18,  face="bold"));porPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
porPCAcen<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position=c(1,0.3)); porPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
#build a plot for the legend
legend_plot<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=TRUE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="bottom")

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  legend_plot + theme(legend.box.margin = margin(1,1,1,1))
)
```

Add segments

``` r
#3. add segments
segpoints<-por.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
porPCAfull<-porPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); porPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Add bi plot with loadings

``` r
#1. make plot with dots
porArrows<-ggplot2::autoplot(scaled_por_afdw, data=por_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=-0.1, loadings.label.hjust=0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-0.5,0.1)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_blank());porArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#porFullPCA<-ggdraw(porPCAfull)+ #+ theme_half_open(12)) +
  #draw_plot(porArrows, .5, .5, .5, .5); porFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Porites.pdf", plot=porFullPCA, dpi=300, width=12, height=8, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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
    ## 1    0.793

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
plasticity.por<-distance.por%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.por$distance)
plasticity.por
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              2.05  2.58 
    ## 2 Mahana Low    tp1.2              0.833 1.05 
    ## 3 Manava High   tp1.2              0.146 0.184
    ## 4 Hilton Medium tp2.3              0.980 1.24 
    ## 5 Mahana Low    tp2.3              0.628 0.792
    ## 6 Manava High   tp2.3              0.819 1.03 
    ## 7 Hilton Medium tp3.4              2.74  3.45 
    ## 8 Mahana Low    tp3.4              1.24  1.56 
    ## 9 Manava High   tp3.4              2.51  3.16

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
por.plasticity.site<-plasticity.por%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Porites")
por.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species
    ##   <chr>             <dbl> <chr>  
    ## 1 Hilton Medium      2.42 Porites
    ## 2 Mahana Low         1.13 Porites
    ## 3 Manava High        1.46 Porites

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.por.time<-plasticity.por%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Porites")
plasticity.por.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species
    ##   <chr>               <dbl> <chr>  
    ## 1 tp1.2                1.27 Porites
    ## 2 tp2.3                1.02 Porites
    ## 3 tp3.4                2.72 Porites

Display mean plasticity score for the species:

``` r
plasticity.por.species<-plasticity.por%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Porites")
plasticity.por.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species
    ##          <dbl> <chr>  
    ## 1         1.67 Porites

#### Pocillopora

Generate PCA using scaled (scaled_poc_afdw) from dataframe
(poc_data_afdw).

``` r
poc_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)

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

``` r
# scale data
vegan <- scale(poc_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = poc_data_afdw, method='eu')
z_poc<-permanova$aov.tab
z_poc
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    411.81 137.269 16.7951 0.25420  0.001 ***
    ## site_code             2     65.34  32.668  3.9969 0.04033  0.001 ***
    ## timepoint:site_code   6    129.38  21.563  2.6383 0.07986  0.001 ***
    ## Residuals           124   1013.48   8.173         0.62560           
    ## Total               135   1620.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
pocPCA<-ggplot(poc_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylim(-6,10)+
  xlim(-5,11)+
  ylab("PC2")+
  xlab("PC1")+
  ggtitle(expression(italic("Pocillopora")))+
  geom_text(x=7, y=-2.5, label=paste("p(Timepoint)=", z_poc$`Pr(>F)`[1]), size=4, color=ifelse(z_poc$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=7, y=-3.25, label=paste("p(Site)=", z_poc$`Pr(>F)`[2]), size=4, color=ifelse(z_poc$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=7, y=-4, label=paste("p(Timepoint x Site)=", z_poc$`Pr(>F)`[3]), size=4, color=ifelse(z_poc$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         legend.title = element_blank(), 
         title = element_text(size=25, face="bold"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18,  face="bold"));pocPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
pocPCAcen<-pocPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=poc.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); pocPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Add segments

``` r
#3. add segments
segpoints<-poc.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
pocPCAfull<-pocPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); pocPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Add bi plot with loadings

``` r
pocArrows<-ggplot2::autoplot(scaled_poc_afdw, data=poc_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=0, loadings.label.hjust=-0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_blank());pocArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#pocFullPCA<-ggdraw(pocPCAfull) + #theme_half_open(12)) +
  #draw_plot(pocArrows, .5, .5, .5, .5); pocFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Pocillopora.pdf", plot=pocFullPCA, dpi=300, width=12, height=8, units="in")
```

Assemble all plots.

``` r
PCA_full_panel<-plot_grid(acrPCAfull, pocPCAfull, porPCAfull, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
PCA_full_panel_legend<-plot_grid(PCA_full_panel, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Full_Panel_PCAs.png", plot=PCA_full_panel_legend, dpi=500, width=20, height=6, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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
    ## 1    0.698

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
plasticity.poc<-distance.poc%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.poc$distance)
plasticity.poc
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2               1.26  1.81
    ## 2 Mahana Low    tp1.2               2.22  3.17
    ## 3 Manava High   tp1.2               1.42  2.03
    ## 4 Hilton Medium tp2.3               2.85  4.08
    ## 5 Mahana Low    tp2.3               1.84  2.63
    ## 6 Manava High   tp2.3               3.38  4.85
    ## 7 Hilton Medium tp3.4               3.88  5.57
    ## 8 Mahana Low    tp3.4               1.73  2.48
    ## 9 Manava High   tp3.4               4.73  6.78

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
poc.plasticity.site<-plasticity.poc%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Pocillopora")
poc.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species    
    ##   <chr>             <dbl> <chr>      
    ## 1 Hilton Medium      3.82 Pocillopora
    ## 2 Mahana Low         2.76 Pocillopora
    ## 3 Manava High        4.55 Pocillopora

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.poc.time<-plasticity.poc%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Pocillopora")
plasticity.poc.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species    
    ##   <chr>               <dbl> <chr>      
    ## 1 tp1.2                2.34 Pocillopora
    ## 2 tp2.3                3.85 Pocillopora
    ## 3 tp3.4                4.94 Pocillopora

Display mean plasticity score for the species:

``` r
plasticity.poc.species<-plasticity.poc%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Pocillopora")
plasticity.poc.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species    
    ##          <dbl> <chr>      
    ## 1         3.71 Pocillopora

### Examine plasticity scores generated from PCA’s

Generate plasticity dataframe from calculations for each species.

``` r
plasticity.site<-bind_rows(acr.plasticity.site, poc.plasticity.site, por.plasticity.site);plasticity.site
```

    ## # A tibble: 9 × 3
    ##   site_code     mean_site species    
    ##   <chr>             <dbl> <chr>      
    ## 1 Hilton Medium      1.86 Acropora   
    ## 2 Mahana Low         3.66 Acropora   
    ## 3 Manava High        3.85 Acropora   
    ## 4 Hilton Medium      3.82 Pocillopora
    ## 5 Mahana Low         2.76 Pocillopora
    ## 6 Manava High        4.55 Pocillopora
    ## 7 Hilton Medium      2.42 Porites    
    ## 8 Mahana Low         1.13 Porites    
    ## 9 Manava High        1.46 Porites

``` r
plasticity.time<-bind_rows(plasticity.acr.time, plasticity.poc.time, plasticity.por.time);plasticity.time
```

    ## # A tibble: 9 × 3
    ##   comparison mean_timepoint species    
    ##   <chr>               <dbl> <chr>      
    ## 1 tp1.2                1.59 Acropora   
    ## 2 tp2.3                3.08 Acropora   
    ## 3 tp3.4                4.69 Acropora   
    ## 4 tp1.2                2.34 Pocillopora
    ## 5 tp2.3                3.85 Pocillopora
    ## 6 tp3.4                4.94 Pocillopora
    ## 7 tp1.2                1.27 Porites    
    ## 8 tp2.3                1.02 Porites    
    ## 9 tp3.4                2.72 Porites

``` r
plasticity.species<-bind_rows(plasticity.acr.species, plasticity.poc.species, plasticity.por.species);plasticity.species
```

    ## # A tibble: 3 × 2
    ##   mean_species species    
    ##          <dbl> <chr>      
    ## 1         3.12 Acropora   
    ## 2         3.71 Pocillopora
    ## 3         1.67 Porites

Plot plasticity scores for site values.

``` r
plasticity_site<-plasticity.site%>%
                group_by(site_code)%>%
                summarise(mean=mean(mean_site))%>%

  ggplot(., aes(x = site_code, y = mean, fill = site_code, group=interaction(site_code))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    #facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_site
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

Plot plasticity scores for timepoint values.

``` r
plasticity_time<-plasticity.time%>%
                group_by(comparison)%>%
                summarise(mean=mean(mean_timepoint))%>%
  
  ggplot(., aes(x = comparison, y = mean, group=interaction(comparison))) +
    geom_point(pch = 21, size=5, fill="gray") + 
    #scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    #facet_wrap(~species) +
    xlab("Timepoint") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_time
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

Plot plasticity scores for species values.

``` r
plasticity_species<-ggplot(plasticity.species, aes(x = species, y = mean_species, fill=species, group=interaction(species))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("darkgray", "orange", "purple"))+
    #scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    #facet_wrap(~species) +
    xlab("Species") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_species
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

Assemble full plasticity figure.

``` r
plasticity_full_panel<-plot_grid(plasticity_species, plasticity_site, plasticity_time, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

ggsave(filename="Figures/Plasticity_Panel.png", plot=plasticity_full_panel, dpi=500, width=10, height=4, units="in")
```

Generate plots for within species.

Plot plasticity scores for site values.

``` r
plasticity_site2<-ggplot(plasticity.site, aes(x = site_code, y = mean_site, fill = site_code, group=interaction(site_code))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_site2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

Plot plasticity scores for timepoint values.

``` r
plasticity_time2<-ggplot(plasticity.time, aes(x = comparison, y = mean_timepoint, group=interaction(comparison))) +
    geom_point(pch = 21, size=5, fill="gray") + 
    #scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    facet_wrap(~species) +
    xlab("Timepoint") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_time2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

Assemble full plasticity figure.

``` r
plasticity_full_panel2<-plot_grid(plasticity_species, plasticity_site2, plasticity_time2, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(0.5,1,1), rel_widths = c(0.5,1,1), label_size = 20, label_y=1, align="vh")

ggsave(filename="Figures/Plasticity_Panel_species_detail.png", plot=plasticity_full_panel2, dpi=500, width=20, height=4, units="in")
```

### Matrix of PCA’s by species and time

Generate biplot showing groupings by site for each species and time
point.

#### Pocillopora

``` r
poc_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint1")

poc_tp1_data<-poc_tp1_data[complete.cases(poc_tp1_data), ]
 
poc_tp1_scaled<-prcomp(poc_tp1_data[c(5:16)], scale=TRUE, center=TRUE) 

pocTP1<-ggplot2::autoplot(poc_tp1_scaled, data=poc_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

``` r
poc_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint2")
 
poc_tp2_data<-poc_tp2_data[complete.cases(poc_tp2_data), ]
 
poc_tp2_scaled<-prcomp(poc_tp2_data[c(5:16)], scale=TRUE, center=TRUE) 

pocTP2<-ggplot2::autoplot(poc_tp2_scaled, data=poc_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

``` r
poc_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint3")
 
poc_tp3_data<-poc_tp3_data[complete.cases(poc_tp3_data), ]
 
poc_tp3_scaled<-prcomp(poc_tp3_data[c(5:16)], scale=TRUE, center=TRUE) 

pocTP3<-ggplot2::autoplot(poc_tp3_scaled, data=poc_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

``` r
poc_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint4")
 
poc_tp4_data<-poc_tp4_data[complete.cases(poc_tp4_data), ]
 
poc_tp4_scaled<-prcomp(poc_tp4_data[c(5:16)], scale=TRUE, center=TRUE) 

pocTP4<-ggplot2::autoplot(poc_tp4_scaled, data=poc_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

#### Acropora

``` r
acr_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint1")

acr_tp1_data<-acr_tp1_data[complete.cases(acr_tp1_data), ]
 
acr_tp1_scaled<-prcomp(acr_tp1_data[c(5:16)], scale=TRUE, center=TRUE) 

acrTP1<-ggplot2::autoplot(acr_tp1_scaled, data=acr_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Acropora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

``` r
acr_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Acropora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint2")
 
acr_tp2_data<-acr_tp2_data[complete.cases(acr_tp2_data), ]
 
acr_tp2_scaled<-prcomp(acr_tp2_data[c(5:16)], scale=TRUE, center=TRUE) 

acrTP2<-ggplot2::autoplot(acr_tp2_scaled, data=acr_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

``` r
acr_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Acropora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint3")
 
acr_tp3_data<-acr_tp3_data[complete.cases(acr_tp3_data), ]
 
acr_tp3_scaled<-prcomp(acr_tp3_data[c(5:16)], scale=TRUE, center=TRUE) 

acrTP3<-ggplot2::autoplot(acr_tp3_scaled, data=acr_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
acr_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Acropora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint4")
 
acr_tp4_data<-acr_tp4_data[complete.cases(acr_tp4_data), ]
 
acr_tp4_scaled<-prcomp(acr_tp4_data[c(5:16)], scale=TRUE, center=TRUE) 

acrTP4<-ggplot2::autoplot(acr_tp4_scaled, data=acr_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

#### Porites

``` r
por_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint1")

por_tp1_data<-por_tp1_data[complete.cases(por_tp1_data), ]
 
por_tp1_scaled<-prcomp(por_tp1_data[c(5:16)], scale=TRUE, center=TRUE) 

porTP1<-ggplot2::autoplot(por_tp1_scaled, data=por_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Porites)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

``` r
por_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Porites")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint2")
 
por_tp2_data<-por_tp2_data[complete.cases(por_tp2_data), ]
 
por_tp2_scaled<-prcomp(por_tp2_data[c(5:16)], scale=TRUE, center=TRUE) 

porTP2<-ggplot2::autoplot(por_tp2_scaled, data=por_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

``` r
por_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Porites")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint3")
 
por_tp3_data<-por_tp3_data[complete.cases(por_tp3_data), ]
 
por_tp3_scaled<-prcomp(por_tp3_data[c(5:16)], scale=TRUE, center=TRUE) 

porTP3<-ggplot2::autoplot(por_tp3_scaled, data=por_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

``` r
por_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Porites")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Chl_a=chla.ug.mgAFDW, Chl_c2=chlc2.ug.mgAFDW, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint4")
 
por_tp4_data<-por_tp4_data[complete.cases(por_tp4_data), ]
 
por_tp4_scaled<-prcomp(por_tp4_data[c(5:16)], scale=TRUE, center=TRUE) 

porTP4<-ggplot2::autoplot(por_tp4_scaled, data=por_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
 # ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

#### Assemble grid plot of all metrics by species and time

``` r
all_grid<-plot_grid(acrTP1, acrTP2, acrTP3, acrTP4, pocTP1, pocTP2, pocTP3, pocTP4, porTP1, porTP2, porTP3, porTP4, ncol=4, nrow=3)
all_grid2<-plot_grid(all_grid, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Matrix_AllResponses.pdf", plot=all_grid2, dpi=500, width=20, height=20, units="in")
```

## Host Responses

Generate a list of symbiont and host responses.

``` r
symbiont_responses<-c("cells.mgAFDW", "Sym_AFDW.mg.cm2", "Total_Chl", "Total_Chl_cell", "Am", "AQY", "Ratio_AFDW.mg.cm2") 
holobiont_responses<-c("cre.umol.mgafdw", "Host_AFDW.mg.cm2", "prot_mg.mgafdw", "Rd", "calc.umol.mgAFDW.hr")
```

### PCA’s for centroid and trajectories

#### Acropora

Generate PCA using scaled (scaled_acr_afdw) from dataframe
(acr_data_afdw).

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

Examine PERMANOVA results.

``` r
# scale data
vegan <- scale(acr_data_afdw[ ,5:9])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = acr_data_afdw, method='eu')
z_acr<-permanova$aov.tab
z_acr
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    184.06  61.353 20.9953 0.35740  0.001 ***
    ## site_code             2     35.97  17.983  6.1539 0.06984  0.001 ***
    ## timepoint:site_code   6     26.13   4.355  1.4902 0.05074  0.102    
    ## Residuals            92    268.85   2.922         0.52203           
    ## Total               103    515.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
acrPCA<-ggplot(acr_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylim(-5,5)+
  xlim(-5,5)+
  ylab("PC2")+
  xlab("PC1")+
  ggtitle(expression(italic("Acropora")))+
  geom_text(x=2.5, y=-2.7, label=paste("p(Timepoint)=", z_acr$`Pr(>F)`[1]), size=4, color=ifelse(z_acr$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=2.5, y=-3.45, label=paste("p(Site)=", z_acr$`Pr(>F)`[2]), size=4, color=ifelse(z_acr$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=2.5, y=-4.2, label=paste("p(Timepoint x Site)=", z_acr$`Pr(>F)`[3]), size=4, color=ifelse(z_acr$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         title = element_text(size=25, face="bold"), 
         axis.title = element_text(size=18));acrPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
acrPCAcen<-acrPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=acr.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); acrPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

Add segments

``` r
#3. add segments
segpoints<-acr.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
acrPCAfull<-acrPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); acrPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

Add bi plot with loadings

``` r
#1. make plot with dots
acrArrows<-ggplot2::autoplot(scaled_acr_afdw, data=acr_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=0.5, loadings.label.hjust=-0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  #scale_colour_manual(values=c("blue4", "darkgray", "springgreen4")) +
  #scale_fill_manual(values=c("blue4", "darkgray", "springgreen4")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         title = element_text(size=25, face="bold"),
         axis.title = element_blank());acrArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

Asemble final figure with inset biplot.

``` r
#acrFullPCA<-ggdraw(acrPCAfull) + #theme_half_open(12)) +
  #draw_plot(acrArrows, .5, .5, .5, .5); acrFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Acropora.pdf", plot=acrFullPCA, dpi=300, width=12, height=8, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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

#spread<-acr_data[c(12,13)]
#spread<-as.matrix(dist(spread, mean.centroid, method="euclidean"))
#spread<-mean(spread)
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
plasticity.acr<-distance.acr%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.acr$distance)
plasticity.acr
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              0.388 0.779
    ## 2 Mahana Low    tp1.2              0.434 0.872
    ## 3 Manava High   tp1.2              0.481 0.968
    ## 4 Hilton Medium tp2.3              1.36  2.73 
    ## 5 Mahana Low    tp2.3              2.69  5.42 
    ## 6 Manava High   tp2.3              2.57  5.17 
    ## 7 Hilton Medium tp3.4              1.90  3.82 
    ## 8 Mahana Low    tp3.4              2.90  5.82 
    ## 9 Manava High   tp3.4              3.21  6.45

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
acr.plasticity.site<-plasticity.acr%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Acropora")
acr.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species 
    ##   <chr>             <dbl> <chr>   
    ## 1 Hilton Medium      2.44 Acropora
    ## 2 Mahana Low         4.04 Acropora
    ## 3 Manava High        4.20 Acropora

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.acr.time<-plasticity.acr%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Acropora")
plasticity.acr.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species 
    ##   <chr>               <dbl> <chr>   
    ## 1 tp1.2               0.873 Acropora
    ## 2 tp2.3               4.44  Acropora
    ## 3 tp3.4               5.37  Acropora

Display mean plasticity score for the species:

``` r
plasticity.acr.species<-plasticity.acr%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Acropora")
plasticity.acr.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species 
    ##          <dbl> <chr>   
    ## 1         3.56 Acropora

#### Porites

Generate PCA using scaled (scaled_por_afdw) from dataframe
(por_data_afdw).

Calculate centroid locations for each site and timepoint and pull PC1
and PC2 locations.

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

``` r
# scale data
vegan <- scale(por_data_afdw[ ,5:8])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = por_data_afdw, method='eu')
z_por<-permanova$aov.tab
z_por
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3     78.00  26.000  8.8182 0.14029  0.001 ***
    ## site_code             2     73.11  36.555 12.3978 0.13149  0.001 ***
    ## timepoint:site_code   6     27.48   4.580  1.5535 0.04943  0.065 .  
    ## Residuals           128    377.41   2.948         0.67879           
    ## Total               139    556.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
porPCA<-ggplot(por_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylab("PC2")+
  xlab("PC1")+
  ylim(-5,5)+
  xlim(-5,5)+
  ggtitle(expression(italic("Porites")))+
  geom_text(x=2, y=-3, label=paste("p(Timepoint)=", z_por$`Pr(>F)`[1]), size=4, color=ifelse(z_por$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=2, y=-3.75, label=paste("p(Site)=", z_por$`Pr(>F)`[2]), size=4, color=ifelse(z_por$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=2, y=-4.5, label=paste("p(Timepoint x Site)=", z_por$`Pr(>F)`[3]), size=4, color=ifelse(z_por$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         legend.title = element_blank(), 
         axis.text = element_text(size=18), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         title = element_text(size=25, face="bold"),
         axis.title = element_text(size=18,  face="bold"));porPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
porPCAcen<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position=c(1,0.3)); porPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-81-1.png)<!-- -->

``` r
#build a plot for the legend
legend_plot<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=TRUE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="bottom")

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  legend_plot + theme(legend.box.margin = margin(1,1,1,1))
)
```

Add segments

``` r
#3. add segments
segpoints<-por.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
porPCAfull<-porPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); porPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

Add bi plot with loadings

``` r
#1. make plot with dots
porArrows<-ggplot2::autoplot(scaled_por_afdw, data=por_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=-0.1, loadings.label.hjust=0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-0.5,0.1)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_blank());porArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-83-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#porFullPCA<-ggdraw(porPCAfull)+ #+ theme_half_open(12)) +
  #draw_plot(porArrows, .5, .5, .5, .5); porFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Porites.pdf", plot=porFullPCA, dpi=300, width=12, height=8, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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
plasticity.por<-distance.por%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.por$distance)
plasticity.por
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              2.29   5.00
    ## 2 Mahana Low    tp1.2              0.886  1.93
    ## 3 Manava High   tp1.2              0.710  1.55
    ## 4 Hilton Medium tp2.3              0.538  1.17
    ## 5 Mahana Low    tp2.3              0.649  1.42
    ## 6 Manava High   tp2.3              0.713  1.56
    ## 7 Hilton Medium tp3.4              1.89   4.12
    ## 8 Mahana Low    tp3.4              0.778  1.70
    ## 9 Manava High   tp3.4              2.13   4.65

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
por.plasticity.site<-plasticity.por%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Porites")
por.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species
    ##   <chr>             <dbl> <chr>  
    ## 1 Hilton Medium      3.43 Porites
    ## 2 Mahana Low         1.68 Porites
    ## 3 Manava High        2.59 Porites

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.por.time<-plasticity.por%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Porites")
plasticity.por.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species
    ##   <chr>               <dbl> <chr>  
    ## 1 tp1.2                2.83 Porites
    ## 2 tp2.3                1.38 Porites
    ## 3 tp3.4                3.49 Porites

Display mean plasticity score for the species:

``` r
plasticity.por.species<-plasticity.por%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Porites")
plasticity.por.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species
    ##          <dbl> <chr>  
    ## 1         2.57 Porites

#### Pocillopora

Generate PCA using scaled (scaled_poc_afdw) from dataframe
(poc_data_afdw).

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

``` r
# scale data
vegan <- scale(poc_data_afdw[ ,5:9])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = poc_data_afdw, method='eu')
z_poc<-permanova$aov.tab
z_poc
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    172.94  57.648 16.6912 0.23854  0.001 ***
    ## site_code             2     28.94  14.470  4.1895 0.03992  0.001 ***
    ## timepoint:site_code   6     60.31  10.052  2.9104 0.08319  0.001 ***
    ## Residuals           134    462.81   3.454         0.63835           
    ## Total               145    725.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
pocPCA<-ggplot(poc_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylim(-5,5)+
  xlim(-5,5)+
  ylab("PC2")+
  xlab("PC1")+
  ggtitle(expression(italic("Pocillopora")))+
  geom_text(x=2, y=-3, label=paste("p(Timepoint)=", z_poc$`Pr(>F)`[1]), size=4, color=ifelse(z_poc$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=2, y=-3.75, label=paste("p(Site)=", z_poc$`Pr(>F)`[2]), size=4, color=ifelse(z_poc$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=2, y=-4.5, label=paste("p(Timepoint x Site)=", z_poc$`Pr(>F)`[3]), size=4, color=ifelse(z_poc$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         legend.title = element_blank(), 
         title = element_text(size=25, face="bold"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18,  face="bold"));pocPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
pocPCAcen<-pocPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=poc.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); pocPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-94-1.png)<!-- -->

Add segments

``` r
#3. add segments
segpoints<-poc.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
pocPCAfull<-pocPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); pocPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-95-1.png)<!-- -->

Add bi plot with loadings

``` r
pocArrows<-ggplot2::autoplot(scaled_poc_afdw, data=poc_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=0, loadings.label.hjust=-0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_blank());pocArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-96-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#pocFullPCA<-ggdraw(pocPCAfull) + #theme_half_open(12)) +
  #draw_plot(pocArrows, .5, .5, .5, .5); pocFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Pocillopora.pdf", plot=pocFullPCA, dpi=300, width=12, height=8, units="in")
```

Assemble all plots.

``` r
PCA_full_panel<-plot_grid(acrPCAfull, pocPCAfull, porPCAfull, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
PCA_full_panel_legend<-plot_grid(PCA_full_panel, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Full_Panel_PCAs_Holobiont.png", plot=PCA_full_panel_legend, dpi=500, width=20, height=6, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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
plasticity.poc<-distance.poc%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.poc$distance)
plasticity.poc
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              0.848  1.73
    ## 2 Mahana Low    tp1.2              0.616  1.26
    ## 3 Manava High   tp1.2              1.60   3.28
    ## 4 Hilton Medium tp2.3              1.47   3.02
    ## 5 Mahana Low    tp2.3              0.886  1.81
    ## 6 Manava High   tp2.3              2.33   4.77
    ## 7 Hilton Medium tp3.4              2.56   5.25
    ## 8 Mahana Low    tp3.4              0.846  1.73
    ## 9 Manava High   tp3.4              2.87   5.87

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
poc.plasticity.site<-plasticity.poc%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Pocillopora")
poc.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species    
    ##   <chr>             <dbl> <chr>      
    ## 1 Hilton Medium      3.33 Pocillopora
    ## 2 Mahana Low         1.60 Pocillopora
    ## 3 Manava High        4.64 Pocillopora

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.poc.time<-plasticity.poc%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Pocillopora")
plasticity.poc.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species    
    ##   <chr>               <dbl> <chr>      
    ## 1 tp1.2                2.09 Pocillopora
    ## 2 tp2.3                3.20 Pocillopora
    ## 3 tp3.4                4.28 Pocillopora

Display mean plasticity score for the species:

``` r
plasticity.poc.species<-plasticity.poc%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Pocillopora")
plasticity.poc.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species    
    ##          <dbl> <chr>      
    ## 1         3.19 Pocillopora

### Examine plasticity scores generated from PCA’s

Generate plasticity dataframe from calculations for each species.

``` r
plasticity.site<-bind_rows(acr.plasticity.site, poc.plasticity.site, por.plasticity.site);plasticity.site
```

    ## # A tibble: 9 × 3
    ##   site_code     mean_site species    
    ##   <chr>             <dbl> <chr>      
    ## 1 Hilton Medium      2.44 Acropora   
    ## 2 Mahana Low         4.04 Acropora   
    ## 3 Manava High        4.20 Acropora   
    ## 4 Hilton Medium      3.33 Pocillopora
    ## 5 Mahana Low         1.60 Pocillopora
    ## 6 Manava High        4.64 Pocillopora
    ## 7 Hilton Medium      3.43 Porites    
    ## 8 Mahana Low         1.68 Porites    
    ## 9 Manava High        2.59 Porites

``` r
plasticity.time<-bind_rows(plasticity.acr.time, plasticity.poc.time, plasticity.por.time);plasticity.time
```

    ## # A tibble: 9 × 3
    ##   comparison mean_timepoint species    
    ##   <chr>               <dbl> <chr>      
    ## 1 tp1.2               0.873 Acropora   
    ## 2 tp2.3               4.44  Acropora   
    ## 3 tp3.4               5.37  Acropora   
    ## 4 tp1.2               2.09  Pocillopora
    ## 5 tp2.3               3.20  Pocillopora
    ## 6 tp3.4               4.28  Pocillopora
    ## 7 tp1.2               2.83  Porites    
    ## 8 tp2.3               1.38  Porites    
    ## 9 tp3.4               3.49  Porites

``` r
plasticity.species<-bind_rows(plasticity.acr.species, plasticity.poc.species, plasticity.por.species);plasticity.species
```

    ## # A tibble: 3 × 2
    ##   mean_species species    
    ##          <dbl> <chr>      
    ## 1         3.56 Acropora   
    ## 2         3.19 Pocillopora
    ## 3         2.57 Porites

Plot plasticity scores for site values.

``` r
plasticity_site<-plasticity.site%>%
                group_by(site_code)%>%
                summarise(mean=mean(mean_site))%>%

  ggplot(., aes(x = site_code, y = mean, fill = site_code, group=interaction(site_code))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    #facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_site
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-106-1.png)<!-- -->

Plot plasticity scores for timepoint values.

``` r
plasticity_time<-plasticity.time%>%
                group_by(comparison)%>%
                summarise(mean=mean(mean_timepoint))%>%
  
  ggplot(., aes(x = comparison, y = mean, group=interaction(comparison))) +
    geom_point(pch = 21, size=5, fill="gray") + 
    #scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    #facet_wrap(~species) +
    xlab("Timepoint") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_time
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-107-1.png)<!-- -->

Plot plasticity scores for species values.

``` r
plasticity_species<-ggplot(plasticity.species, aes(x = species, y = mean_species, fill=species, group=interaction(species))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("darkgray", "orange", "purple"))+
    #scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    #facet_wrap(~species) +
    xlab("Species") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,5)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_species
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-108-1.png)<!-- -->

Assemble full plasticity figure.

``` r
plasticity_full_panel<-plot_grid(plasticity_species, plasticity_site, plasticity_time, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

ggsave(filename="Figures/Plasticity_Panel_Host.png", plot=plasticity_full_panel, dpi=500, width=10, height=4, units="in")
```

Generate plots for within species.

Plot plasticity scores for site values.

``` r
plasticity_site2<-ggplot(plasticity.site, aes(x = site_code, y = mean_site, fill = site_code, group=interaction(site_code))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,6)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_site2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-110-1.png)<!-- -->

Plot plasticity scores for timepoint values.

``` r
plasticity_time2<-ggplot(plasticity.time, aes(x = comparison, y = mean_timepoint, group=interaction(comparison))) +
    geom_point(pch = 21, size=5, fill="gray") + 
    #scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    facet_wrap(~species) +
    xlab("Timepoint") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,6)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_time2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-111-1.png)<!-- -->

Assemble full plasticity figure.

``` r
plasticity_full_panel2<-plot_grid(plasticity_species, plasticity_site2, plasticity_time2, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(0.5,1,1), rel_widths = c(0.5,1,1), label_size = 20, label_y=1, align="vh")

ggsave(filename="Figures/Plasticity_Panel_species_detail_host.png", plot=plasticity_full_panel2, dpi=500, width=20, height=4, units="in")
```

### Matrix of PCA’s by species and time

Generate biplot showing groupings by site for each species and time
point.

#### Pocillopora

``` r
poc_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint1")

poc_tp1_data<-poc_tp1_data[complete.cases(poc_tp1_data), ]
 
poc_tp1_scaled<-prcomp(poc_tp1_data[c(5:9)], scale=TRUE, center=TRUE) 

pocTP1<-ggplot2::autoplot(poc_tp1_scaled, data=poc_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-113-1.png)<!-- -->

``` r
poc_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint2")
 
poc_tp2_data<-poc_tp2_data[complete.cases(poc_tp2_data), ]
 
poc_tp2_scaled<-prcomp(poc_tp2_data[c(5:9)], scale=TRUE, center=TRUE) 

pocTP2<-ggplot2::autoplot(poc_tp2_scaled, data=poc_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-114-1.png)<!-- -->

``` r
poc_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint3")
 
poc_tp3_data<-poc_tp3_data[complete.cases(poc_tp3_data), ]
 
poc_tp3_scaled<-prcomp(poc_tp3_data[c(5:9)], scale=TRUE, center=TRUE) 

pocTP3<-ggplot2::autoplot(poc_tp3_scaled, data=poc_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-115-1.png)<!-- -->

``` r
poc_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint4")
 
poc_tp4_data<-poc_tp4_data[complete.cases(poc_tp4_data), ]
 
poc_tp4_scaled<-prcomp(poc_tp4_data[c(5:9)], scale=TRUE, center=TRUE) 

pocTP4<-ggplot2::autoplot(poc_tp4_scaled, data=poc_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-116-1.png)<!-- -->

#### Acropora

``` r
acr_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint1")

acr_tp1_data<-acr_tp1_data[complete.cases(acr_tp1_data), ]
 
acr_tp1_scaled<-prcomp(acr_tp1_data[c(5:9)], scale=TRUE, center=TRUE) 

acrTP1<-ggplot2::autoplot(acr_tp1_scaled, data=acr_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Acropora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-117-1.png)<!-- -->

``` r
acr_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint2")
 
acr_tp2_data<-acr_tp2_data[complete.cases(acr_tp2_data), ]
 
acr_tp2_scaled<-prcomp(acr_tp2_data[c(5:9)], scale=TRUE, center=TRUE) 

acrTP2<-ggplot2::autoplot(acr_tp2_scaled, data=acr_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-118-1.png)<!-- -->

``` r
acr_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint3")
 
acr_tp3_data<-acr_tp3_data[complete.cases(acr_tp3_data), ]
 
acr_tp3_scaled<-prcomp(acr_tp3_data[c(5:9)], scale=TRUE, center=TRUE) 

acrTP3<-ggplot2::autoplot(acr_tp3_scaled, data=acr_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-119-1.png)<!-- -->

``` r
acr_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint4")
 
acr_tp4_data<-acr_tp4_data[complete.cases(acr_tp4_data), ]
 
acr_tp4_scaled<-prcomp(acr_tp4_data[c(5:9)], scale=TRUE, center=TRUE) 

acrTP4<-ggplot2::autoplot(acr_tp4_scaled, data=acr_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-120-1.png)<!-- -->

#### Porites

``` r
por_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint1")

por_tp1_data<-por_tp1_data[complete.cases(por_tp1_data), ]
 
por_tp1_scaled<-prcomp(por_tp1_data[c(5:9)], scale=TRUE, center=TRUE) 

porTP1<-ggplot2::autoplot(por_tp1_scaled, data=por_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Porites)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-121-1.png)<!-- -->

``` r
por_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Porites")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint2")
 
por_tp2_data<-por_tp2_data[complete.cases(por_tp2_data), ]
 
por_tp2_scaled<-prcomp(por_tp2_data[c(5:9)], scale=TRUE, center=TRUE) 

porTP2<-ggplot2::autoplot(por_tp2_scaled, data=por_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-122-1.png)<!-- -->

``` r
por_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Porites")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint3")
 
por_tp3_data<-por_tp3_data[complete.cases(por_tp3_data), ]
 
por_tp3_scaled<-prcomp(por_tp3_data[c(5:9)], scale=TRUE, center=TRUE) 

porTP3<-ggplot2::autoplot(por_tp3_scaled, data=por_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-123-1.png)<!-- -->

``` r
por_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(holobiont_responses))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2>0)%>%
  #filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")%>%
  filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Host_Biomass=Host_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
  filter(timepoint=="timepoint4")
 
por_tp4_data<-por_tp4_data[complete.cases(por_tp4_data), ]
 
por_tp4_scaled<-prcomp(por_tp4_data[c(5:9)], scale=TRUE, center=TRUE) 

porTP4<-ggplot2::autoplot(por_tp4_scaled, data=por_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
 # ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-124-1.png)<!-- -->

#### Assemble grid plot of all metrics by species and time

``` r
all_grid<-plot_grid(acrTP1, acrTP2, acrTP3, acrTP4, pocTP1, pocTP2, pocTP3, pocTP4, porTP1, porTP2, porTP3, porTP4, ncol=4, nrow=3)
all_grid2<-plot_grid(all_grid, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Matrix_Host_Responses.pdf", plot=all_grid2, dpi=500, width=20, height=20, units="in")
```

## Symbiont Responses

### PCA’s for centroid and trajectories

#### Acropora

Generate PCA using scaled (scaled_acr_afdw) from dataframe
(acr_data_afdw).

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

Examine PERMANOVA results.

``` r
# scale data
vegan <- scale(acr_data_afdw[ ,5:11])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = acr_data_afdw, method='eu')
z_acr<-permanova$aov.tab
z_acr
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    172.28  57.426 11.5172 0.21974  0.001 ***
    ## site_code             2     49.19  24.594  4.9324 0.06274  0.001 ***
    ## timepoint:site_code   6     58.93   9.822  1.9698 0.07517  0.009 ** 
    ## Residuals           101    503.60   4.986         0.64235           
    ## Total               112    784.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
acrPCA<-ggplot(acr_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylim(-5,6)+
  xlim(-5,6)+
  ylab("PC2")+
  xlab("PC1")+
  ggtitle(expression(italic("Acropora")))+
  geom_text(x=3, y=-3, label=paste("p(Timepoint)=", z_acr$`Pr(>F)`[1]), size=4, color=ifelse(z_acr$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=3, y=-3.75, label=paste("p(Site)=", z_acr$`Pr(>F)`[2]), size=4, color=ifelse(z_acr$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=3, y=-4.5, label=paste("p(Timepoint x Site)=", z_acr$`Pr(>F)`[3]), size=4, color=ifelse(z_acr$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         title = element_text(size=25, face="bold"), 
         axis.title = element_text(size=18));acrPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-128-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
acrPCAcen<-acrPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=acr.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); acrPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-129-1.png)<!-- -->

Add segments

``` r
#3. add segments
segpoints<-acr.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
acrPCAfull<-acrPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); acrPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-130-1.png)<!-- -->

Add bi plot with loadings

``` r
#1. make plot with dots
acrArrows<-ggplot2::autoplot(scaled_acr_afdw, data=acr_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=0.5, loadings.label.hjust=-0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  #scale_colour_manual(values=c("blue4", "darkgray", "springgreen4")) +
  #scale_fill_manual(values=c("blue4", "darkgray", "springgreen4")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         title = element_text(size=25, face="bold"),
         axis.title = element_blank());acrArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-131-1.png)<!-- -->

Asemble final figure with inset biplot.

``` r
#acrFullPCA<-ggdraw(acrPCAfull) + #theme_half_open(12)) +
  #draw_plot(acrArrows, .5, .5, .5, .5); acrFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Acropora.pdf", plot=acrFullPCA, dpi=300, width=12, height=8, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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

#spread<-acr_data[c(12,13)]
#spread<-as.matrix(dist(spread, mean.centroid, method="euclidean"))
#spread<-mean(spread)
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
plasticity.acr<-distance.acr%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.acr$distance)
plasticity.acr
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              0.529  1.06
    ## 2 Mahana Low    tp1.2              1.98   3.97
    ## 3 Manava High   tp1.2              0.924  1.85
    ## 4 Hilton Medium tp2.3              1.92   3.85
    ## 5 Mahana Low    tp2.3              2.18   4.37
    ## 6 Manava High   tp2.3              2.90   5.83
    ## 7 Hilton Medium tp3.4              1.42   2.85
    ## 8 Mahana Low    tp3.4              2.20   4.42
    ## 9 Manava High   tp3.4              2.16   4.33

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
acr.plasticity.site<-plasticity.acr%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Acropora")
acr.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species 
    ##   <chr>             <dbl> <chr>   
    ## 1 Hilton Medium      2.59 Acropora
    ## 2 Mahana Low         4.25 Acropora
    ## 3 Manava High        4.00 Acropora

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.acr.time<-plasticity.acr%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Acropora")
plasticity.acr.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species 
    ##   <chr>               <dbl> <chr>   
    ## 1 tp1.2                2.30 Acropora
    ## 2 tp2.3                4.68 Acropora
    ## 3 tp3.4                3.87 Acropora

Display mean plasticity score for the species:

``` r
plasticity.acr.species<-plasticity.acr%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Acropora")
plasticity.acr.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species 
    ##          <dbl> <chr>   
    ## 1         3.61 Acropora

#### Porites

Generate PCA using scaled (scaled_por_afdw) from dataframe
(por_data_afdw).

Calculate centroid locations for each site and timepoint and pull PC1
and PC2 locations.

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

``` r
# scale data
vegan <- scale(por_data_afdw[ ,5:11])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = por_data_afdw, method='eu')
z_por<-permanova$aov.tab
z_por
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    125.81  41.937  7.6005 0.11521  0.001 ***
    ## site_code             2    122.49  61.243 11.0993 0.11217  0.001 ***
    ## timepoint:site_code   6     43.63   7.272  1.3179 0.03996  0.122    
    ## Residuals           145    800.07   5.518         0.73267           
    ## Total               156   1092.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
porPCA<-ggplot(por_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylab("PC2")+
  xlab("PC1")+
  ylim(-5,6)+
  xlim(-5,6)+
  ggtitle(expression(italic("Porites")))+
  geom_text(x=3, y=-3, label=paste("p(Timepoint)=", z_por$`Pr(>F)`[1]), size=4, color=ifelse(z_por$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=3, y=-3.75, label=paste("p(Site)=", z_por$`Pr(>F)`[2]), size=4, color=ifelse(z_por$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=3, y=-4.5, label=paste("p(Timepoint x Site)=", z_por$`Pr(>F)`[3]), size=4, color=ifelse(z_por$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         legend.title = element_blank(), 
         axis.text = element_text(size=18), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         title = element_text(size=25, face="bold"),
         axis.title = element_text(size=18,  face="bold"));porPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-141-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
porPCAcen<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position=c(1,0.3)); porPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-142-1.png)<!-- -->

``` r
#build a plot for the legend
legend_plot<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=TRUE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="bottom")

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  legend_plot + theme(legend.box.margin = margin(1,1,1,1))
)
```

Add segments

``` r
#3. add segments
segpoints<-por.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
porPCAfull<-porPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); porPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-143-1.png)<!-- -->

Add bi plot with loadings

``` r
#1. make plot with dots
porArrows<-ggplot2::autoplot(scaled_por_afdw, data=por_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=-0.1, loadings.label.hjust=0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-0.5,0.1)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_blank());porArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-144-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#porFullPCA<-ggdraw(porPCAfull)+ #+ theme_half_open(12)) +
 # draw_plot(porArrows, .5, .5, .5, .5); porFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Porites.pdf", plot=porFullPCA, dpi=300, width=12, height=8, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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
plasticity.por<-distance.por%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.por$distance)
plasticity.por
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              0.912  1.94
    ## 2 Mahana Low    tp1.2              1.50   3.18
    ## 3 Manava High   tp1.2              1.58   3.36
    ## 4 Hilton Medium tp2.3              1.03   2.19
    ## 5 Mahana Low    tp2.3              1.20   2.55
    ## 6 Manava High   tp2.3              0.656  1.39
    ## 7 Hilton Medium tp3.4              1.68   3.56
    ## 8 Mahana Low    tp3.4              1.12   2.38
    ## 9 Manava High   tp3.4              1.09   2.32

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
por.plasticity.site<-plasticity.por%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Porites")
por.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species
    ##   <chr>             <dbl> <chr>  
    ## 1 Hilton Medium      2.56 Porites
    ## 2 Mahana Low         2.70 Porites
    ## 3 Manava High        2.36 Porites

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.por.time<-plasticity.por%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Porites")
plasticity.por.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species
    ##   <chr>               <dbl> <chr>  
    ## 1 tp1.2                2.83 Porites
    ## 2 tp2.3                2.04 Porites
    ## 3 tp3.4                2.75 Porites

Display mean plasticity score for the species:

``` r
plasticity.por.species<-plasticity.por%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Porites")
plasticity.por.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species
    ##          <dbl> <chr>  
    ## 1         2.54 Porites

#### Pocillopora

Generate PCA using scaled (scaled_poc_afdw) from dataframe
(poc_data_afdw).

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

``` r
# scale data
vegan <- scale(poc_data_afdw[ ,5:11])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis(vegan ~ timepoint*site_code, data = poc_data_afdw, method='eu')
z_poc<-permanova$aov.tab
z_poc
```

    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## timepoint             3    218.61  72.869 14.0694 0.20960  0.001 ***
    ## site_code             2     32.32  16.159  3.1200 0.03099  0.001 ***
    ## timepoint:site_code   6     77.34  12.889  2.4887 0.07415  0.001 ***
    ## Residuals           138    714.74   5.179         0.68527           
    ## Total               149   1043.00                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Assemble plot with background points

``` r
#1. make plot with dots
pocPCA<-ggplot(poc_data, aes(.fittedPC1, .fittedPC2, color=site_code)) + 
  geom_point(size = 4, alpha=0.2, show.legend=FALSE)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  ylim(-5,6)+
  xlim(-5,6)+
  ylab("PC2")+
  xlab("PC1")+
  ggtitle(expression(italic("Pocillopora")))+
  geom_text(x=3, y=-3, label=paste("p(Timepoint)=", z_poc$`Pr(>F)`[1]), size=4, color=ifelse(z_poc$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=3, y=-3.75, label=paste("p(Site)=", z_poc$`Pr(>F)`[2]), size=4, color=ifelse(z_poc$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=3, y=-4.5, label=paste("p(Timepoint x Site)=", z_poc$`Pr(>F)`[3]), size=4, color=ifelse(z_poc$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=18), 
         legend.position="none",
         plot.background = element_blank(),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         legend.title = element_blank(), 
         title = element_text(size=25, face="bold"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18,  face="bold"));pocPCA
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-154-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
pocPCAcen<-pocPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=poc.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); pocPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-155-1.png)<!-- -->

Add segments

``` r
#3. add segments
segpoints<-poc.centroids%>%
  gather(variable, value, -(timepoint:site_code)) %>%
  unite(temp, timepoint, variable) %>%
  spread(temp, value)
  
pocPCAfull<-pocPCAcen + 
  geom_segment(aes(x = timepoint1_PC1.mean, y = timepoint1_PC2.mean, xend = timepoint2_PC1.mean, yend = timepoint2_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint2_PC1.mean, y = timepoint2_PC2.mean, xend = timepoint3_PC1.mean, yend = timepoint3_PC2.mean, colour = site_code), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = timepoint3_PC1.mean, y = timepoint3_PC2.mean, xend = timepoint4_PC1.mean, yend = timepoint4_PC2.mean, colour = site_code), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); pocPCAfull
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-156-1.png)<!-- -->

Add bi plot with loadings

``` r
pocArrows<-ggplot2::autoplot(scaled_poc_afdw, data=poc_data_afdw, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, loadings.label.size=3.5, loadings.label.vjust=0, loadings.label.hjust=-0.1, loadings.label.repel=TRUE, size=4, alpha=0.0) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         plot.background = element_blank(),
         legend.title = element_blank(), 
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_blank());pocArrows
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-157-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#pocFullPCA<-ggdraw(pocPCAfull) + #theme_half_open(12)) +
  #draw_plot(pocArrows, .5, .5, .5, .5); pocFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Pocillopora.pdf", plot=pocFullPCA, dpi=300, width=12, height=8, units="in")
```

Assemble all plots.

``` r
PCA_full_panel<-plot_grid(acrPCAfull, pocPCAfull, porPCAfull, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
PCA_full_panel_legend<-plot_grid(PCA_full_panel, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Full_Panel_PCAs_Symbiont.png", plot=PCA_full_panel_legend, dpi=500, width=20, height=6, units="in")
```

##### Calculate plasticity using centroid travel calculations.

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
plasticity.poc<-distance.poc%>%
  select(site_code, tp1.2, tp2.3, tp3.4)%>% #keep desired columns
  gather(key="comparison", value="distance_time", -site_code)%>% #gather and generate new rows to designate time comparisons
  mutate(ratio=distance_time/spread.poc$distance)
plasticity.poc
```

    ## # A tibble: 9 × 4
    ## # Groups:   site_code [3]
    ##   site_code     comparison distance_time ratio
    ##   <chr>         <chr>              <dbl> <dbl>
    ## 1 Hilton Medium tp1.2              0.604  2.17
    ## 2 Mahana Low    tp1.2              2.16   7.77
    ## 3 Manava High   tp1.2              1.54   5.52
    ## 4 Hilton Medium tp2.3              2.57   9.24
    ## 5 Mahana Low    tp2.3              1.92   6.91
    ## 6 Manava High   tp2.3              2.64   9.47
    ## 7 Hilton Medium tp3.4              2.59   9.31
    ## 8 Mahana Low    tp3.4              1.31   4.71
    ## 9 Manava High   tp3.4              2.49   8.93

4 - average these ratios across sites to generate a mean plasticity
ratio for time, site, and species

Display mean plasticity scores by site:

``` r
poc.plasticity.site<-plasticity.poc%>%
  group_by(site_code)%>%
  summarise(mean_site=mean(ratio))%>%
  mutate(species="Pocillopora")
poc.plasticity.site
```

    ## # A tibble: 3 × 3
    ##   site_code     mean_site species    
    ##   <chr>             <dbl> <chr>      
    ## 1 Hilton Medium      6.91 Pocillopora
    ## 2 Mahana Low         6.47 Pocillopora
    ## 3 Manava High        7.97 Pocillopora

Display mean plasticity scores by timepoint comparison:

``` r
plasticity.poc.time<-plasticity.poc%>%
  group_by(comparison)%>%
  summarise(mean_timepoint=mean(ratio))%>%
  mutate(species="Pocillopora")
plasticity.poc.time
```

    ## # A tibble: 3 × 3
    ##   comparison mean_timepoint species    
    ##   <chr>               <dbl> <chr>      
    ## 1 tp1.2                5.15 Pocillopora
    ## 2 tp2.3                8.54 Pocillopora
    ## 3 tp3.4                7.65 Pocillopora

Display mean plasticity score for the species:

``` r
plasticity.poc.species<-plasticity.poc%>%
  summarise(mean_species=mean(ratio))%>%
  summarise(mean_species=mean(mean_species))%>%
  mutate(species="Pocillopora")
plasticity.poc.species
```

    ## # A tibble: 1 × 2
    ##   mean_species species    
    ##          <dbl> <chr>      
    ## 1         7.12 Pocillopora

### Examine plasticity scores generated from PCA’s

Generate plasticity dataframe from calculations for each species.

``` r
plasticity.site<-bind_rows(acr.plasticity.site, poc.plasticity.site, por.plasticity.site);plasticity.site
```

    ## # A tibble: 9 × 3
    ##   site_code     mean_site species    
    ##   <chr>             <dbl> <chr>      
    ## 1 Hilton Medium      2.59 Acropora   
    ## 2 Mahana Low         4.25 Acropora   
    ## 3 Manava High        4.00 Acropora   
    ## 4 Hilton Medium      6.91 Pocillopora
    ## 5 Mahana Low         6.47 Pocillopora
    ## 6 Manava High        7.97 Pocillopora
    ## 7 Hilton Medium      2.56 Porites    
    ## 8 Mahana Low         2.70 Porites    
    ## 9 Manava High        2.36 Porites

``` r
plasticity.time<-bind_rows(plasticity.acr.time, plasticity.poc.time, plasticity.por.time);plasticity.time
```

    ## # A tibble: 9 × 3
    ##   comparison mean_timepoint species    
    ##   <chr>               <dbl> <chr>      
    ## 1 tp1.2                2.30 Acropora   
    ## 2 tp2.3                4.68 Acropora   
    ## 3 tp3.4                3.87 Acropora   
    ## 4 tp1.2                5.15 Pocillopora
    ## 5 tp2.3                8.54 Pocillopora
    ## 6 tp3.4                7.65 Pocillopora
    ## 7 tp1.2                2.83 Porites    
    ## 8 tp2.3                2.04 Porites    
    ## 9 tp3.4                2.75 Porites

``` r
plasticity.species<-bind_rows(plasticity.acr.species, plasticity.poc.species, plasticity.por.species);plasticity.species
```

    ## # A tibble: 3 × 2
    ##   mean_species species    
    ##          <dbl> <chr>      
    ## 1         3.61 Acropora   
    ## 2         7.12 Pocillopora
    ## 3         2.54 Porites

Plot plasticity scores for site values.

``` r
plasticity_site<-plasticity.site%>%
                group_by(site_code)%>%
                summarise(mean=mean(mean_site))%>%

  ggplot(., aes(x = site_code, y = mean, fill = site_code, group=interaction(site_code))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    #facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,8)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_site
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-167-1.png)<!-- -->

Plot plasticity scores for timepoint values.

``` r
plasticity_time<-plasticity.time%>%
                group_by(comparison)%>%
                summarise(mean=mean(mean_timepoint))%>%
  
  ggplot(., aes(x = comparison, y = mean, group=interaction(comparison))) +
    geom_point(pch = 21, size=5, fill="gray") + 
    #scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    #facet_wrap(~species) +
    xlab("Timepoint") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,8)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_time
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-168-1.png)<!-- -->

Plot plasticity scores for species values.

``` r
plasticity_species<-ggplot(plasticity.species, aes(x = species, y = mean_species, fill=species, group=interaction(species))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("darkgray", "orange", "purple"))+
    #scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    #facet_wrap(~species) +
    xlab("Species") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,8)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_species
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-169-1.png)<!-- -->

Assemble full plasticity figure.

``` r
plasticity_full_panel<-plot_grid(plasticity_species, plasticity_site, plasticity_time, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

ggsave(filename="Figures/Plasticity_Panel_Symbiont.png", plot=plasticity_full_panel, dpi=500, width=10, height=4, units="in")
```

Generate plots for within species.

Plot plasticity scores for site values.

``` r
plasticity_site2<-ggplot(plasticity.site, aes(x = site_code, y = mean_site, fill = site_code, group=interaction(site_code))) +
    geom_point(pch = 21, size=5, position = position_jitterdodge(0)) + 
    scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,8)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_site2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-171-1.png)<!-- -->

Plot plasticity scores for timepoint values.

``` r
plasticity_time2<-ggplot(plasticity.time, aes(x = comparison, y = mean_timepoint, group=interaction(comparison))) +
    geom_point(pch = 21, size=5, fill="gray") + 
    #scale_fill_manual(values = c("#374d7c", "#00cccc", "#ff6633"))+
    scale_x_discrete(labels=c("tp1.2" = "Jan-March", "tp2.3" = "March-Sept", "tp3.4" = "Sept-Nov"))+
    facet_wrap(~species) +
    xlab("Timepoint") + 
    ylab(expression(bold("Plasticity (Distance / Spread)")))+
    theme_classic() + 
    ylim(0,9)+
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=12),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      ); plasticity_time2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-172-1.png)<!-- -->

Assemble full plasticity figure.

``` r
plasticity_full_panel2<-plot_grid(plasticity_species, plasticity_site2, plasticity_time2, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(0.5,1,1), rel_widths = c(0.5,1,1), label_size = 20, label_y=1, align="vh")

ggsave(filename="Figures/Plasticity_Panel_species_detail_symbiont.png", plot=plasticity_full_panel2, dpi=500, width=20, height=4, units="in")
```

### Matrix of PCA’s by species and time

Generate biplot showing groupings by site for each species and time
point.

#### Pocillopora

``` r
poc_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Pocillopora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2)%>%
  filter(timepoint=="timepoint1")

poc_tp1_data<-poc_tp1_data[complete.cases(poc_tp1_data), ]
 
poc_tp1_scaled<-prcomp(poc_tp1_data[c(5:11)], scale=TRUE, center=TRUE) 

pocTP1<-ggplot2::autoplot(poc_tp1_scaled, data=poc_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-174-1.png)<!-- -->

``` r
poc_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint2")
 
poc_tp2_data<-poc_tp2_data[complete.cases(poc_tp2_data), ]
 
poc_tp2_scaled<-prcomp(poc_tp2_data[c(5:11)], scale=TRUE, center=TRUE) 

pocTP2<-ggplot2::autoplot(poc_tp2_scaled, data=poc_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-175-1.png)<!-- -->

``` r
poc_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Pocillopora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint3")
 
poc_tp3_data<-poc_tp3_data[complete.cases(poc_tp3_data), ]
 
poc_tp3_scaled<-prcomp(poc_tp3_data[c(5:11)], scale=TRUE, center=TRUE) 

pocTP3<-ggplot2::autoplot(poc_tp3_scaled, data=poc_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-176-1.png)<!-- -->

``` r
poc_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint4")
 
poc_tp4_data<-poc_tp4_data[complete.cases(poc_tp4_data), ]
 
poc_tp4_scaled<-prcomp(poc_tp4_data[c(5:11)], scale=TRUE, center=TRUE) 

pocTP4<-ggplot2::autoplot(poc_tp4_scaled, data=poc_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));pocTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-177-1.png)<!-- -->

#### Acropora

``` r
acr_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint1")

acr_tp1_data<-acr_tp1_data[complete.cases(acr_tp1_data), ]
 
acr_tp1_scaled<-prcomp(acr_tp1_data[c(5:11)], scale=TRUE, center=TRUE) 

acrTP1<-ggplot2::autoplot(acr_tp1_scaled, data=acr_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Acropora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-178-1.png)<!-- -->

``` r
acr_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint2")
 
acr_tp2_data<-acr_tp2_data[complete.cases(acr_tp2_data), ]
 
acr_tp2_scaled<-prcomp(acr_tp2_data[c(5:11)], scale=TRUE, center=TRUE) 

acrTP2<-ggplot2::autoplot(acr_tp2_scaled, data=acr_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-179-1.png)<!-- -->

``` r
acr_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint3")
 
acr_tp3_data<-acr_tp3_data[complete.cases(acr_tp3_data), ]
 
acr_tp3_scaled<-prcomp(acr_tp3_data[c(5:11)], scale=TRUE, center=TRUE) 

acrTP3<-ggplot2::autoplot(acr_tp3_scaled, data=acr_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-180-1.png)<!-- -->

``` r
acr_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Acropora")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint4")
 
acr_tp4_data<-acr_tp4_data[complete.cases(acr_tp4_data), ]
 
acr_tp4_scaled<-prcomp(acr_tp4_data[c(5:11)], scale=TRUE, center=TRUE) 

acrTP4<-ggplot2::autoplot(acr_tp4_scaled, data=acr_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));acrTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-181-1.png)<!-- -->

#### Porites

``` r
por_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Porites")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint1")

por_tp1_data<-por_tp1_data[complete.cases(por_tp1_data), ]
 
por_tp1_scaled<-prcomp(por_tp1_data[c(5:11)], scale=TRUE, center=TRUE) 

porTP1<-ggplot2::autoplot(por_tp1_scaled, data=por_tp1_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0, loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("January 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  annotate("text", x=-16, y=0, label="italic(Porites)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP1
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-182-1.png)<!-- -->

``` r
por_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Porites")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint2")
 
por_tp2_data<-por_tp2_data[complete.cases(por_tp2_data), ]
 
por_tp2_scaled<-prcomp(por_tp2_data[c(5:11)], scale=TRUE, center=TRUE) 

porTP2<-ggplot2::autoplot(por_tp2_scaled, data=por_tp2_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("March 2020")+
  coord_cartesian(xlim=c(-8, 8), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP2
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-183-1.png)<!-- -->

``` r
por_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Porites")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint3")
 
por_tp3_data<-por_tp3_data[complete.cases(por_tp3_data), ]
 
por_tp3_scaled<-prcomp(por_tp3_data[c(5:11)], scale=TRUE, center=TRUE) 

porTP3<-ggplot2::autoplot(por_tp3_scaled, data=por_tp3_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
  #ggtitle("September 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP3
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-184-1.png)<!-- -->

``` r
por_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, all_of(symbiont_responses))%>%
  #filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(Total_Chl<13)%>%
  filter(species=="Porites")%>%
  #filter(prot_mg.mgafdw<1.5)%>% #remove outliers
  rename(Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Symbiont_Density=cells.mgAFDW, Total_Chl=Total_Chl, Total_Chl_Cell=Total_Chl_cell)%>%
  filter(timepoint=="timepoint4")
 
por_tp4_data<-por_tp4_data[complete.cases(por_tp4_data), ]
 
por_tp4_scaled<-prcomp(por_tp4_data[c(5:11)], scale=TRUE, center=TRUE) 

porTP4<-ggplot2::autoplot(por_tp4_scaled, data=por_tp4_data, loadings=TRUE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=TRUE, frame=FALSE, frame.type="convex", frame.level=0.95, loadings.label.size=5, loadings.label.repel=TRUE, size=4, alpha=1.0,   loadings.arrow = grid::arrow(length = grid::unit(10, "points")), scale=0) + 
  #geom_point(aes(colour=site_code), size=3)+
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
  scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) + 
  theme_classic()+
  #xlim(-.3,.6)+
  #ylim(-.3, .3)+
 # ggtitle("November 2020")+
  coord_cartesian(xlim=c(-8, 9), ylim=c(-5, 5), clip = "off")+
  #annotate("text", x=-8, y=0, label="italic(Pocillopora)", fontface="bold", parse = TRUE, angle=90, size=10)+
   theme(legend.text = element_text(size=18), 
         legend.position="none",
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         plot.background = element_blank(),
         plot.title = element_text(size=22, face="bold", hjust = 0.5), 
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         axis.text = element_text(size=18), 
         axis.title = element_text(size=18));porTP4
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-185-1.png)<!-- -->

#### Assemble grid plot of all metrics by species and time

``` r
all_grid<-plot_grid(acrTP1, acrTP2, acrTP3, acrTP4, pocTP1, pocTP2, pocTP3, pocTP4, porTP1, porTP2, porTP3, porTP4, ncol=4, nrow=3)
all_grid2<-plot_grid(all_grid, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Matrix_Symbiont_Responses.pdf", plot=all_grid2, dpi=500, width=20, height=20, units="in")
```
