Examining normalizers and preliminary PCA plots
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

# Load dataframe

Load in master dataframe generated from 1_assemble_data.Rmd.

``` r
master<-read.csv("Output/master_timeseries.csv")

#reorder site levels 
master$side_code<-as.factor(master$site_code)
master$site_code<-fct_relevel(master$site_code, "Mahana Low", "Hilton Medium", "Manava High")
```

# Examine multivariate profiles for each normalizer

## All Responses (host and symbiont)

### PCA’s of each normalizer

#### PCAs: Surface Area normalizers

##### Combined PCA’s

Plot PCA using prcomp to center and scale.

``` r
pca_data_cm2<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, prot_mg.cm2, cells.cm2, cre.umol.cm2, Am, AQY, Rd, calc.umol.cm2.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)

pca_data_cm2<-pca_data_cm2[complete.cases(pca_data_cm2), ]

scaled_data_cm2<-prcomp(pca_data_cm2[c(5:16)], scale=TRUE, center=TRUE) 
```

Groupings by species:

``` r
pcaSpecies_cm2<-ggplot2::autoplot(scaled_data_cm2, data=pca_data_cm2, frame.colour="species", loadings=FALSE,  colour="species", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("orange", "darkgray", "purple")) +
  scale_fill_manual(values=c("orange", "darkgray", "purple")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaSpecies_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

PERMANOVA of the effect of species

``` r
# scale data
vegan_cm2 <- scale(pca_data_cm2[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ species, data = pca_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ species, data = pca_data_cm2, method = "eu")
    ##           Df SumOfSqs      R2     F Pr(>F)    
    ## species    2   2200.9 0.48522 177.2  0.001 ***
    ## Residual 376   2335.1 0.51478                 
    ## Total    378   4536.0 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by timepoint:

``` r
pcaTimepoint_cm2<-ggplot2::autoplot(scaled_data_cm2, data=pca_data_cm2, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaTimepoint_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ timepoint, data = pca_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ timepoint, data = pca_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3      299 0.06591 8.8208  0.001 ***
    ## Residual  375     4237 0.93409                  
    ## Total     378     4536 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
pcaSite_cm2<-ggplot2::autoplot(scaled_data_cm2, data=pca_data_cm2, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaSite_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ site_code, data = pca_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ site_code, data = pca_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## site_code   2       89 0.01962 3.7621  0.005 **
    ## Residual  376     4447 0.98038                 
    ## Total     378     4536 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generate a plot of all PCA’s for surface area normalization.

``` r
pcaplots_cm2<-plot_grid(pcaSpecies_cm2, pcaSite_cm2, pcaTimepoint_cm2, labels = c("Species: Surface Area", "Site: Surface Area", "Timepoint: Surface Area"), ncol=3, nrow=1, rel_widths = c(1, 1, 1), label_size = 20, label_x = 0.01)

#ggsave(filename="Figures/All_PCAs_cm2.pdf", plot=pcaplots_cm2, dpi=300, width=24, height=5, units="in")
```

##### Species-specific PCA’s

Next generate PCA for site and timepoint for each species separately.

Separate by species.

``` r
#porites
por_data_cm2<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, prot_mg.cm2, cells.cm2, cre.umol.cm2, Am, AQY, Rd, calc.umol.cm2.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Porites")

por_data_cm2<-por_data_cm2[complete.cases(por_data_cm2), ]

scaled_por_cm2<-prcomp(por_data_cm2[c(5:16)], scale=TRUE, center=TRUE) 

#acropora
acr_data_cm2<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, prot_mg.cm2, cells.cm2, cre.umol.cm2, Am, AQY, Rd, calc.umol.cm2.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")

acr_data_cm2<-acr_data_cm2[complete.cases(acr_data_cm2), ]

scaled_acr_cm2<-prcomp(acr_data_cm2[c(5:16)], scale=TRUE, center=TRUE) 

#pocillopora
poc_data_cm2<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.cm2, chlc2.ug.cm2, prot_mg.cm2, cells.cm2, cre.umol.cm2, Am, AQY, Rd, calc.umol.cm2.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")

poc_data_cm2<-poc_data_cm2[complete.cases(poc_data_cm2), ]

scaled_poc_cm2<-prcomp(poc_data_cm2[c(5:16)], scale=TRUE, center=TRUE) 
```

###### Porites PCA’s

Groupings by timepoint:

``` r
porTimepoint_cm2<-ggplot2::autoplot(scaled_por_cm2, data=por_data_cm2, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.7,0.5)+
  #ylim(-0.7, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));porTimepoint_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_cm2 <- scale(por_data_cm2[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ timepoint, data = por_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ timepoint, data = por_data_cm2, method = "eu")
    ##            Df SumOfSqs     R2      F Pr(>F)    
    ## timepoint   3   240.61 0.1453 7.6498  0.001 ***
    ## Residual  135  1415.39 0.8547                  
    ## Total     138  1656.00 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
porSite_cm2<-ggplot2::autoplot(scaled_por_cm2, data=por_data_cm2, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.7,0.5)+
  #ylim(-0.7, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));porSite_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ site_code, data = por_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ site_code, data = por_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2    F Pr(>F)    
    ## site_code   2   225.14 0.13595 10.7  0.001 ***
    ## Residual  136  1430.86 0.86405                
    ## Total     138  1656.00 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###### Acropora PCA’s

Groupings by timepoint:

``` r
acrTimepoint_cm2<-ggplot2::autoplot(scaled_acr_cm2, data=acr_data_cm2, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));acrTimepoint_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_cm2 <- scale(acr_data_cm2[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ timepoint, data = acr_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ timepoint, data = acr_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3    327.6 0.26505 12.021  0.001 ***
    ## Residual  100    908.4 0.73495                  
    ## Total     103   1236.0 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
acrSite_cm2<-ggplot2::autoplot(scaled_acr_cm2, data=acr_data_cm2, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));acrSite_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ site_code, data = acr_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ site_code, data = acr_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## site_code   2   120.53 0.09752 5.4568  0.001 ***
    ## Residual  101  1115.47 0.90248                  
    ## Total     103  1236.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###### Pocillopora PCA’s

Groupings by timepoint:

``` r
pocTimepoint_cm2<-ggplot2::autoplot(scaled_poc_cm2, data=poc_data_cm2, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.6,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pocTimepoint_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_cm2 <- scale(poc_data_cm2[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ timepoint, data = poc_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ timepoint, data = poc_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3   312.86 0.19312 10.531  0.001 ***
    ## Residual  132  1307.14 0.80688                  
    ## Total     135  1620.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
pocSite_cm2<-ggplot2::autoplot(scaled_poc_cm2, data=poc_data_cm2, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pocSite_cm2
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_cm2 ~ site_code, data = poc_data_cm2, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_cm2 ~ site_code, data = poc_data_cm2, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## site_code   2    89.23 0.05508 3.8762  0.001 ***
    ## Residual  133  1530.77 0.94492                  
    ## Total     135  1620.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generate a plot of all PCA’s for all Species

``` r
species_plots_cm2<-plot_grid(porSite_cm2, porTimepoint_cm2, acrSite_cm2, acrTimepoint_cm2, pocSite_cm2, pocTimepoint_cm2, labels = c("Porites: Surface Area", "Porites: Surface Area", "Acropora: Surface Area", "Acropora: Surface Area", "Pocillopora: Surface Area", "Pocillopora: Surface Area"), ncol=2, nrow=3, rel_widths = c(1, 1, 1, 1, 1, 1), label_size = 20, label_x = 0.01, label_y=1.0)

#ggsave(filename="Figures/Species_PCAs_cm2.pdf", plot=species_plots_cm2, dpi=300, width=18, height=15, units="in")
```

#### PCAs: Protein normalizers

##### Combined PCA’s

Plot PCA using prcomp to center and scale.

``` r
pca_data_prot<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgprot, chlc2.ug.mgprot, prot_mg.cm2, cells.mgprot, cre.umol.mgprot, Am, AQY, Rd, calc.umol.mgprot.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)

pca_data_prot<-pca_data_prot[complete.cases(pca_data_prot), ]

scaled_data_prot<-prcomp(pca_data_prot[c(5:16)], scale=TRUE, center=TRUE) 
```

Groupings by species:

``` r
pcaSpecies_prot<-ggplot2::autoplot(scaled_data_prot, data=pca_data_prot, frame.colour="species", loadings=FALSE,  colour="species", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("orange", "darkgray", "purple")) +
  scale_fill_manual(values=c("orange", "darkgray", "purple")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaSpecies_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

PERMANOVA of the effect of species

``` r
# scale data
vegan_prot <- scale(pca_data_prot[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ species, data = pca_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ species, data = pca_data_prot, method = "eu")
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## species    2   1958.4 0.43175 142.84  0.001 ***
    ## Residual 376   2577.6 0.56825                  
    ## Total    378   4536.0 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by timepoint:

``` r
pcaTimepoint_prot<-ggplot2::autoplot(scaled_data_prot, data=pca_data_prot, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaTimepoint_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ timepoint, data = pca_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ timepoint, data = pca_data_prot, method = "eu")
    ##            Df SumOfSqs      R2     F Pr(>F)    
    ## timepoint   3    396.4 0.08739 11.97  0.001 ***
    ## Residual  375   4139.6 0.91261                 
    ## Total     378   4536.0 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
pcaSite_prot<-ggplot2::autoplot(scaled_data_prot, data=pca_data_prot, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaSite_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ site_code, data = pca_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ site_code, data = pca_data_prot, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## site_code   2     77.9 0.01717 3.2842  0.006 **
    ## Residual  376   4458.1 0.98283                 
    ## Total     378   4536.0 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generate a plot of all PCA’s for protein normalization.

``` r
pcaplots_prot<-plot_grid(pcaSpecies_prot, pcaSite_prot, pcaTimepoint_prot, labels = c("Species: Protein", "Site: Protein", "Timepoint: Protein"), ncol=3, nrow=1, rel_widths = c(1, 1, 1), label_size = 20, label_x = 0.02)

#ggsave(filename="Figures/All_PCAs_prot.pdf", plot=pcaplots_prot, dpi=300, width=24, height=5, units="in")
```

##### Species-specific PCA’s

Next generate PCA for site and timepoint for each species separately.

Separate by species.

``` r
#porites
por_data_prot<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgprot, chlc2.ug.mgprot, prot_mg.cm2, cells.mgprot, cre.umol.mgprot, Am, AQY, Rd, calc.umol.mgprot.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Porites")

por_data_prot<-por_data_prot[complete.cases(por_data_prot), ]

scaled_por_prot<-prcomp(por_data_prot[c(5:16)], scale=TRUE, center=TRUE) 

#acropora
acr_data_prot<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgprot, chlc2.ug.mgprot, prot_mg.cm2, cells.mgprot, cre.umol.mgprot, Am, AQY, Rd, calc.umol.mgprot.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Acropora")

acr_data_prot<-acr_data_prot[complete.cases(acr_data_prot), ]

scaled_acr_prot<-prcomp(acr_data_prot[c(5:16)], scale=TRUE, center=TRUE) 

#pocillopora
poc_data_prot<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgprot, chlc2.ug.mgprot, prot_mg.cm2, cells.mgprot, cre.umol.mgprot, Am, AQY, Rd, calc.umol.mgprot.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")

poc_data_prot<-poc_data_prot[complete.cases(poc_data_prot), ]

scaled_poc_prot<-prcomp(poc_data_prot[c(5:16)], scale=TRUE, center=TRUE) 
```

###### Porites PCA’s

Groupings by timepoint:

``` r
porTimepoint_prot<-ggplot2::autoplot(scaled_por_prot, data=por_data_prot, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.7,0.5)+
  #ylim(-0.7, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));porTimepoint_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_prot <- scale(por_data_prot[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ timepoint, data = por_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ timepoint, data = por_data_prot, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3   237.89 0.14365 7.5488  0.001 ***
    ## Residual  135  1418.11 0.85635                  
    ## Total     138  1656.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
porSite_prot<-ggplot2::autoplot(scaled_por_prot, data=por_data_prot, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.7,0.5)+
  #ylim(-0.7, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));porSite_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ site_code, data = por_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ site_code, data = por_data_prot, method = "eu")
    ##            Df SumOfSqs     R2      F Pr(>F)    
    ## site_code   2   142.91 0.0863 6.4224  0.001 ***
    ## Residual  136  1513.09 0.9137                  
    ## Total     138  1656.00 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###### Acropora PCA’s

Groupings by timepoint:

``` r
acrTimepoint_prot<-ggplot2::autoplot(scaled_acr_prot, data=acr_data_prot, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));acrTimepoint_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_prot <- scale(acr_data_prot[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ timepoint, data = acr_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ timepoint, data = acr_data_prot, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3   325.79 0.26358 11.931  0.001 ***
    ## Residual  100   910.21 0.73642                  
    ## Total     103  1236.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
acrSite_prot<-ggplot2::autoplot(scaled_acr_prot, data=acr_data_prot, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));acrSite_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ site_code, data = acr_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ site_code, data = acr_data_prot, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## site_code   2   101.01 0.08172 4.4943  0.001 ***
    ## Residual  101  1134.99 0.91828                  
    ## Total     103  1236.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###### Pocillopora PCA’s

Groupings by timepoint:

``` r
pocTimepoint_prot<-ggplot2::autoplot(scaled_poc_prot, data=poc_data_prot, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.6,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pocTimepoint_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_prot <- scale(poc_data_prot[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ timepoint, data = poc_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ timepoint, data = poc_data_prot, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3   311.07 0.19202 10.457  0.001 ***
    ## Residual  132  1308.93 0.80798                  
    ## Total     135  1620.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
pocSite_prot<-ggplot2::autoplot(scaled_poc_prot, data=poc_data_prot, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pocSite_prot
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_prot ~ site_code, data = poc_data_prot, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_prot ~ site_code, data = poc_data_prot, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## site_code   2    84.54 0.05219 3.6615  0.001 ***
    ## Residual  133  1535.46 0.94781                  
    ## Total     135  1620.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generate a plot of all PCA’s for all Species

``` r
species_plots_prot<-plot_grid(porSite_prot, porTimepoint_prot, acrSite_prot, acrTimepoint_prot, pocSite_prot, pocTimepoint_prot, labels = c("Porites: Protein", "Porites: Protein", "Acropora: Protein", "Acropora: Protein", "Pocillopora: Protein", "Pocillopora: Protein"), ncol=2, nrow=3, rel_widths = c(1, 1, 1, 1, 1, 1), label_size = 20, label_x = 0.01, label_y=1.0)

#ggsave(filename="Figures/Species_PCAs_prot.pdf", plot=species_plots_prot, dpi=300, width=18, height=15, units="in")
```

#### PCAs: Biomass normalizers

##### Combined PCA’s

Plot PCA using prcomp to center and scale.

``` r
pca_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(prot_mg.mgafdw<1.5) #remove outliers


pca_data_afdw<-pca_data_afdw[complete.cases(pca_data_afdw), ]

scaled_data_afdw<-prcomp(pca_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 
```

``` r
# scale data
vegan <- scale(pca_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
permanova<-adonis2(vegan ~ species, data = pca_data_afdw, method='eu')
z_species<-permanova
z_species
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan ~ species, data = pca_data_afdw, method = "eu")
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## species    2   1762.1 0.39054 119.83  0.001 ***
    ## Residual 374   2749.9 0.60946                  
    ## Total    376   4512.0 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by species:

``` r
pcaSpecies_afdw<-ggplot2::autoplot(scaled_data_afdw, data=pca_data_afdw, frame.colour="species", loadings=FALSE,  colour="species", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm', alpha=0.5) + 
  scale_colour_manual(values=c("darkgray", "orange", "purple")) +
  scale_fill_manual(values=c("darkgray", "orange", "purple")) + 
  theme_classic()+
  geom_text(x=0.15, y=-0.1, label=paste("p(Species)=", z_species$`Pr(>F)`[1]), size=4, color=ifelse(z_species$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
  #geom_text(x=0.3, y=-0.15, label=paste("p(Species)=", z_species$`Pr(>F)`[1]), size=6, color="black") + 
  xlab("PC1")+
  ylab("PC2")+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_blank(), 
        axis.text = element_text(size=18, color="black"), 
        axis.title = element_text(size=18, color="black"));pcaSpecies_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

``` r
ggsave(filename="Figures/Species_Overall_PCA.png", plot=pcaSpecies_afdw, dpi=300, width=7, height=5, units="in")
```

PERMANOVA of the effect of species

``` r
# scale data
vegan_afdw <- scale(pca_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ species, data = pca_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ species, data = pca_data_afdw, method = "eu")
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## species    2   1762.1 0.39054 119.83  0.001 ***
    ## Residual 374   2749.9 0.60946                  
    ## Total    376   4512.0 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by timepoint:

``` r
pcaTimepoint_afdw<-ggplot2::autoplot(scaled_data_afdw, data=pca_data_afdw, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaTimepoint_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ timepoint, data = pca_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ timepoint, data = pca_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3    517.9 0.11478 16.122  0.001 ***
    ## Residual  373   3994.1 0.88522                  
    ## Total     376   4512.0 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
pcaSite_afdw<-ggplot2::autoplot(scaled_data_afdw, data=pca_data_afdw, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.3,0.5)+
  #ylim(-0.3, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pcaSite_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ site_code, data = pca_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ site_code, data = pca_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## site_code   2       90 0.01995 3.8067  0.001 ***
    ## Residual  374     4422 0.98005                  
    ## Total     376     4512 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generate a plot of all PCA’s for afdw normalization.

``` r
pcaplots_afdw<-plot_grid(pcaSpecies_afdw, pcaSite_afdw, pcaTimepoint_afdw, labels = c("Species: Biomass", "Site: Biomass", "Timepoint: Biomass"), ncol=3, nrow=1, rel_widths = c(1, 1, 1), label_size = 20, label_x = 0.02)

#ggsave(filename="Figures/All_PCAs_afdw.pdf", plot=pcaplots_afdw, dpi=300, width=24, height=5, units="in")
```

##### Species-specific PCA’s

Next generate PCA for site and timepoint for each species separately.

Separate by species.

``` r
#porites
por_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Porites")

por_data_afdw<-por_data_afdw[complete.cases(por_data_afdw), ]

scaled_por_afdw<-prcomp(por_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 

#acropora
acr_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(species=="Acropora")

acr_data_afdw<-acr_data_afdw[complete.cases(acr_data_afdw), ]

scaled_acr_afdw<-prcomp(acr_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 

#pocillopora
poc_data_afdw<-master%>%
  select(colony_id_corr, timepoint, species, site_code, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chla.ug.mgAFDW, chlc2.ug.mgAFDW, prot_mg.mgafdw, cells.mgAFDW, cre.umol.mgafdw, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(species=="Pocillopora")

poc_data_afdw<-poc_data_afdw[complete.cases(poc_data_afdw), ]

scaled_poc_afdw<-prcomp(poc_data_afdw[c(5:16)], scale=TRUE, center=TRUE) 
```

###### Porites PCA’s

Groupings by timepoint:

``` r
porTimepoint_afdw<-ggplot2::autoplot(scaled_por_afdw, data=por_data_afdw, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.7,0.5)+
  #ylim(-0.7, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));porTimepoint_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_afdw <- scale(por_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ timepoint, data = por_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ timepoint, data = por_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3   253.23 0.15292 8.1236  0.001 ***
    ## Residual  135  1402.77 0.84708                  
    ## Total     138  1656.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
porSite_afdw<-ggplot2::autoplot(scaled_por_afdw, data=por_data_afdw, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
 # xlim(-0.7,0.5)+
 # ylim(-0.7, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));porSite_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ site_code, data = por_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ site_code, data = por_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## site_code   2   186.68 0.11273 8.6394  0.001 ***
    ## Residual  136  1469.32 0.88727                  
    ## Total     138  1656.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###### Acropora PCA’s

Groupings by timepoint:

``` r
acrTimepoint_afdw<-ggplot2::autoplot(scaled_acr_afdw, data=acr_data_afdw, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));acrTimepoint_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_afdw <- scale(acr_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ timepoint, data = acr_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ timepoint, data = acr_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## timepoint   3   315.59 0.25533 11.429  0.001 ***
    ## Residual  100   920.41 0.74467                  
    ## Total     103  1236.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
acrSite_afdw<-ggplot2::autoplot(scaled_acr_afdw, data=acr_data_afdw, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));acrSite_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ site_code, data = acr_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ site_code, data = acr_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2    F Pr(>F)    
    ## site_code   2    79.26 0.06412 3.46  0.001 ***
    ## Residual  101  1156.74 0.93588                
    ## Total     103  1236.00 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

###### Pocillopora PCA’s

Groupings by timepoint:

``` r
pocTimepoint_afdw<-ggplot2::autoplot(scaled_poc_afdw, data=poc_data_afdw, frame.colour="timepoint", loadings=FALSE,  colour="timepoint", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) +
  scale_fill_manual(values=c("darkturquoise", "darkgray", "mediumvioletred", "black"), labels=c("Jan2020", "March2020", "Sept2020", "Nov2020")) + 
  theme_classic()+
 # xlim(-0.6,0.5)+
  #ylim(-0.5, 0.8)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pocTimepoint_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

PERMANOVA of the effect of timepoint

``` r
# scale data
vegan_afdw <- scale(poc_data_afdw[ ,5:16])

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ timepoint, data = poc_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ timepoint, data = poc_data_afdw, method = "eu")
    ##            Df SumOfSqs     R2      F Pr(>F)    
    ## timepoint   3   411.81 0.2542 14.997  0.001 ***
    ## Residual  132  1208.19 0.7458                  
    ## Total     135  1620.00 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Groupings by site:

``` r
pocSite_afdw<-ggplot2::autoplot(scaled_poc_afdw, data=poc_data_afdw, frame.colour="site_code", loadings=FALSE,  colour="site_code", loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, loadings.label.size=5, loadings.label.vjust=-1, size=5, frame.type = 'norm') + 
  scale_colour_manual(values=c("blue4", "darkgray", "springgreen3")) +
  scale_fill_manual(values=c("blue4", "darkgray", "springgreen3")) + 
  theme_classic()+
  #xlim(-0.5,0.5)+
  #ylim(-0.5, 0.5)+
   theme(legend.text = element_text(size=18), 
         legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18), 
        axis.title = element_text(size=18,  face="bold"));pocSite_afdw
```

![](3_normalizer_pca_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

PERMANOVA of the effect of site

``` r
# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(vegan_afdw ~ site_code, data = poc_data_afdw, method='eu')
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = vegan_afdw ~ site_code, data = poc_data_afdw, method = "eu")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## site_code   2    77.35 0.04774 3.3342  0.002 **
    ## Residual  133  1542.65 0.95226                 
    ## Total     135  1620.00 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generate a plot of all PCA’s for all Species

``` r
species_plots_afdw<-plot_grid(porSite_afdw, porTimepoint_afdw, acrSite_afdw, acrTimepoint_afdw, pocSite_afdw, pocTimepoint_afdw, labels = c("Porites: Biomass", "Porites: Biomass", "Acropora: Biomass", "Acropora: Biomass", "Pocillopora: Biomass", "Pocillopora: Biomass"), ncol=2, nrow=3, rel_widths = c(1, 1, 1, 1, 1, 1), label_size = 20, label_x = 0.01, label_y=1.0)

#ggsave(filename="Figures/Species_PCAs_afdw.pdf", plot=species_plots_afdw, dpi=300, width=18, height=15, units="in")
```

#### Assemble all PCA plots

Generate plot of all combined PCA’s

``` r
All_PCAs<-plot_grid(pcaSpecies_prot, pcaSite_prot, pcaTimepoint_prot, pcaSpecies_afdw, pcaSite_afdw, pcaTimepoint_afdw, pcaSpecies_cm2, pcaSite_cm2, pcaTimepoint_cm2, labels = c("Protein", "Protein", "Protein", "Biomass", "Biomass", "Biomass","Surface Area", "Surface Area", "Surface Area"), ncol=3, nrow=3, rel_widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1), label_size = 20, label_x = 0.06, label_y=1.0)

ggsave(filename="Figures/All_PCAs.pdf", plot=All_PCAs, dpi=300, width=24, height=15, units="in")
```

Generate plot of all species specific PCA’s

``` r
All_PCAs_species<-plot_grid(porSite_prot, porTimepoint_prot, porSite_afdw, porTimepoint_afdw, porSite_cm2, porTimepoint_cm2, 
                    acrSite_prot, acrTimepoint_prot, acrSite_afdw, acrTimepoint_afdw, acrSite_cm2, acrTimepoint_cm2,
                    pocSite_prot, pocTimepoint_prot, pocSite_afdw, pocTimepoint_afdw, pocSite_cm2, pocTimepoint_cm2,
                    labels = c("Porites: Protein", "Porites: Protein", "Porites: Biomass", "Porites: Biomass", "Porites: SA", "Porites: SA",
                               "Acropora: Protein", "Acropora: Protein", "Acropora: Biomass", "Acropora: Biomass", "Acropora: SA", "Acropora: SA",
                               "Pocillopora: Protein", "Pocillopora: Protein", "Pocillopora: Biomass", "Pocillopora: Biomass", "Pocillopora: SA", "Pocillopora: SA"), ncol=6, nrow=3, rel_widths = c(1), label_size = 20, label_x = 0.06, label_y=1.0)

ggsave(filename="Figures/All_Species_PCAs.pdf", plot=All_PCAs_species, dpi=300, width=48, height=15, units="in")
```
