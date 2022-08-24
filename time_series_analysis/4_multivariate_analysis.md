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

#reorder site levels 
master$side_code<-as.factor(master$site_code)
master$site_code<-fct_relevel(master$site_code, "Mahana Low", "Hilton Medium", "Manava High")
```

# Build trajectory PCA plots

## All Responses

### PCA’s for centroid and trajectories

#### Acropora

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

Examine PERMANOVA results.

``` r
# scale data
vegan <- scale(acr_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis2(vegan ~ timepoint*site_code, data = acr_data_afdw, method='eu')
z_acr<-permanova
z_acr
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = acr_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   342.87 0.28012 14.4921  0.001 ***
    ## site_code             2    85.26 0.06966  5.4055  0.001 ***
    ## timepoint:site_code   6    78.22 0.06391  1.6531  0.021 *  
    ## Residual             91   717.65 0.58632                   
    ## Total               102  1224.00 1.00000                   
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

#### Porites

Generate PCA using scaled (scaled_por_afdw) from dataframe
(por_data_afdw).

Calculate centroid locations for each site and timepoint and pull PC1
and PC2 locations.

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

``` r
# scale data
vegan <- scale(por_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis2(vegan ~ timepoint*site_code, data = por_data_afdw, method='eu')
z_por<-permanova
z_por
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = por_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   274.09 0.16672 10.5063  0.001 ***
    ## site_code             2   185.52 0.11285 10.6669  0.001 ***
    ## timepoint:site_code   6    88.68 0.05394  1.6996  0.016 *  
    ## Residual            126  1095.71 0.66649                   
    ## Total               137  1644.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
porPCAcen<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position=c(1,0.3)); porPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

#### Pocillopora

Generate PCA using scaled (scaled_poc_afdw) from dataframe
(poc_data_afdw).

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

``` r
# scale data
vegan <- scale(poc_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis2(vegan ~ timepoint*site_code, data = poc_data_afdw, method='eu')
z_poc<-permanova
z_poc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = poc_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   390.79 0.24123 15.4405  0.001 ***
    ## site_code             2    62.30 0.03846  3.6923  0.001 ***
    ## timepoint:site_code   6   120.81 0.07457  2.3866  0.001 ***
    ## Residual            124  1046.11 0.64575                   
    ## Total               135  1620.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
pocPCAcen<-pocPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=poc.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); pocPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Assemble all plots.

``` r
PCA_full_panel<-plot_grid(acrPCAfull, pocPCAfull, porPCAfull, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
PCA_full_panel_legend<-plot_grid(PCA_full_panel, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Full_Panel_PCAs.png", plot=PCA_full_panel_legend, dpi=500, width=20, height=6, units="in")
```

### Matrix of PCA’s by species and time

Generate biplot showing groupings by site for each species and time
point.

#### Pocillopora

``` r
poc_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
poc_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
poc_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
poc_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Pocillopora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

#### Acropora

``` r
acr_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
acr_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Acropora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
acr_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Acropora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
acr_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Acropora")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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
  coord_cartesian(xlim=c(-8, 9), ylim=c(-7, 5), clip = "off")+
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

#### Porites

``` r
por_tp1_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
por_tp2_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Porites")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
por_tp3_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Porites")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
por_tp4_data<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, Total_Chl, Total_Chl_cell, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  #filter(chla.ug.mgAFDW<2.2)%>%
  filter(species=="Porites")%>%
  rename(Host_Biomass=Host_AFDW.mg.cm2, Symbiont_Biomass=Sym_AFDW.mg.cm2, S_H_Biomass_Ratio=Ratio_AFDW.mg.cm2, Host_Protein=prot_mg.mgafdw, Symbiont_Density=cells.mgAFDW, Antiox_Capacity=cre.umol.mgafdw, Calc=calc.umol.mgAFDW.hr)%>%
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

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
permanova<-adonis2(vegan ~ timepoint*site_code, data = acr_data_afdw, method='eu')
z_acr<-permanova
z_acr
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = acr_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   184.06 0.35740 20.9953  0.001 ***
    ## site_code             2    35.97 0.06984  6.1539  0.001 ***
    ## timepoint:site_code   6    26.13 0.05074  1.4902  0.097 .  
    ## Residual             92   268.85 0.52203                   
    ## Total               103   515.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
acrPCAcen<-acrPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=acr.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); acrPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#acrFullPCA<-ggdraw(acrPCAfull) + #theme_half_open(12)) +
  #draw_plot(acrArrows, .5, .5, .5, .5); acrFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Acropora.pdf", plot=acrFullPCA, dpi=300, width=12, height=8, units="in")
```

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
permanova<-adonis2(vegan ~ timepoint*site_code, data = por_data_afdw, method='eu')
z_por<-permanova
z_por
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = por_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3    78.00 0.14029  8.8182  0.001 ***
    ## site_code             2    73.11 0.13149 12.3978  0.001 ***
    ## timepoint:site_code   6    27.48 0.04943  1.5535  0.062 .  
    ## Residual            128   377.41 0.67879                   
    ## Total               139   556.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
porPCAcen<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position=c(1,0.3)); porPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

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
permanova<-adonis2(vegan ~ timepoint*site_code, data = poc_data_afdw, method='eu')
z_poc<-permanova
z_poc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = poc_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   172.94 0.23854 16.6912  0.001 ***
    ## site_code             2    28.94 0.03992  4.1895  0.001 ***
    ## timepoint:site_code   6    60.31 0.08319  2.9104  0.001 ***
    ## Residual            134   462.81 0.63835                   
    ## Total               145   725.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
pocPCAcen<-pocPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=poc.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); pocPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

Assemble all plots.

``` r
PCA_full_panel<-plot_grid(acrPCAfull, pocPCAfull, porPCAfull, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
PCA_full_panel_legend<-plot_grid(PCA_full_panel, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Full_Panel_PCAs_Holobiont.png", plot=PCA_full_panel_legend, dpi=500, width=20, height=6, units="in")
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

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
permanova<-adonis2(vegan ~ timepoint*site_code, data = acr_data_afdw, method='eu')
z_acr<-permanova
z_acr
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = acr_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   172.28 0.21974 11.5172  0.001 ***
    ## site_code             2    49.19 0.06274  4.9324  0.001 ***
    ## timepoint:site_code   6    58.93 0.07517  1.9698  0.009 ** 
    ## Residual            101   503.60 0.64235                   
    ## Total               112   784.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
acrPCAcen<-acrPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=acr.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); acrPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

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
permanova<-adonis2(vegan ~ timepoint*site_code, data = por_data_afdw, method='eu')
z_por<-permanova
z_por
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = por_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   125.81 0.11521  7.6005  0.001 ***
    ## site_code             2   122.49 0.11217 11.0993  0.001 ***
    ## timepoint:site_code   6    43.63 0.03996  1.3179  0.126    
    ## Residual            145   800.07 0.73267                   
    ## Total               156  1092.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
porPCAcen<-porPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=por.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position=c(1,0.3)); porPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->

Assemble final figure with inset biplot.

``` r
#porFullPCA<-ggdraw(porPCAfull)+ #+ theme_half_open(12)) +
 # draw_plot(porArrows, .5, .5, .5, .5); porFullPCA #x, y, w, h

#ggsave(filename="Figures/FullPCA_Porites.pdf", plot=porFullPCA, dpi=300, width=12, height=8, units="in")
```

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
permanova<-adonis2(vegan ~ timepoint*site_code, data = poc_data_afdw, method='eu')
z_poc<-permanova
z_poc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ timepoint * site_code, data = poc_data_afdw, method = "eu")
    ##                      Df SumOfSqs      R2       F Pr(>F)    
    ## timepoint             3   218.61 0.20960 14.0694  0.001 ***
    ## site_code             2    32.32 0.03099  3.1200  0.003 ** 
    ## timepoint:site_code   6    77.34 0.07415  2.4887  0.001 ***
    ## Residual            138   714.74 0.68527                   
    ## Total               149  1043.00 1.00000                   
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-84-1.png)<!-- -->

Add centroids

``` r
#2. add centroids 
pocPCAcen<-pocPCA + geom_point(aes(x=PC1.mean, y=PC2.mean, colour=site_code), data=poc.centroids, size=3, show.legend=FALSE) + 
  scale_colour_manual(values=c("#374d7c", "#00cccc", "#ff6633"))+
  theme(legend.position="none"); pocPCAcen
```

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-85-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-86-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-87-1.png)<!-- -->

Assemble all plots.

``` r
PCA_full_panel<-plot_grid(acrPCAfull, pocPCAfull, porPCAfull, labels = c("", "", ""), ncol=3, nrow=1, rel_heights= c(1,1,1), rel_widths = c(1,1,1), label_size = 20, label_y=1, align="vh")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
PCA_full_panel_legend<-plot_grid(PCA_full_panel, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Full_Panel_PCAs_Symbiont.png", plot=PCA_full_panel_legend, dpi=500, width=20, height=6, units="in")
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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-89-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-90-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-91-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-92-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-94-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-95-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-96-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-97-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-98-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-99-1.png)<!-- -->

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

![](4_multivariate_analysis_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->

#### Assemble grid plot of all metrics by species and time

``` r
all_grid<-plot_grid(acrTP1, acrTP2, acrTP3, acrTP4, pocTP1, pocTP2, pocTP3, pocTP4, porTP1, porTP2, porTP3, porTP4, ncol=4, nrow=3)
all_grid2<-plot_grid(all_grid, legend, rel_heights = c(4, .2), ncol=1, nrow=2)

ggsave(filename="Figures/Matrix_Symbiont_Responses.pdf", plot=all_grid2, dpi=500, width=20, height=20, units="in")
```
