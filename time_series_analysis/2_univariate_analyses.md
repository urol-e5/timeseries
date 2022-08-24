Univariate analysis of E5 time series biological data
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

# Univariate Responses

Plot univariate response plots and then generate panel of all plots at
the end of this section. Individual plots will not be displayed.

## Antioxidant capacity

### Plot: CRE umol mgprot

View by site and timepoint for each species.

``` r
tacplot_prot<-master %>%
  filter(!is.na(cre.umol.mgprot)) %>% 
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, cre.umol.mgprot)%>%
  
  ggplot(., aes(x = site_code, y = cre.umol.mgprot, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    ylab(expression(bold(paste("Antioxidant Capacity (", mu, "mol CRE mg protein"^-1,")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: CRE umol mgafdw

View by site and timepoint for each species.

``` r
tacplot_afdw<-master %>%
  filter(!is.na(cre.umol.mgafdw)) %>%
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, cre.umol.mgafdw)%>%
  
  ggplot(., aes(x = site_code, y = cre.umol.mgafdw, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    ylab(expression(bold(paste("Antioxidant Capacity (", mu, "mol CRE mg afdw"^-1,")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: CRE umol cm2

View by site and timepoint for each species.

``` r
tacplot_cm2<-master %>%
  filter(!is.na(cre.umol.cm2)) %>%
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, cre.umol.cm2)%>%
  
  ggplot(., aes(x = site_code, y = cre.umol.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
     scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    ylab(expression(bold(paste("Antioxidant Capacity (", mu, "mol CRE mg afdw"^-1,")"))))+
    theme_classic() + 
    theme(
      legend.position="right",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots for each normalization together.

``` r
tac_figure<-plot_grid(tacplot_prot, tacplot_afdw, tacplot_cm2, labels = c("", "", ""), ncol=3, nrow=1, rel_widths = c(.8, .8, 1), label_size = 20, label_x = 0.1);tac_figure
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave(filename="Figures/TAC_Figure.pdf", plot=tac_figure, dpi=300, width=24, height=5, units="in")
```

### Analysis

Build a mixed model for univariate analysis and examine data
distribution of biomass normalized antioxidant capacity.

`antiox_model<-lmer(cre.umol.mgafdw~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
antiox_model<-lmer(cre.umol.mgafdw~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(antiox_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    ## 107 218 
    ##  99 207

``` r
#hist(residuals(antiox_model))
```

Residuals are not normally distributed. Attempt with log transformation.

`antiox_model<-lmer(log(1+cre.umol.mgafdw)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
antiox_model<-lmer(log(1+cre.umol.mgafdw)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(antiox_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## 218 107 
    ## 207  99

``` r
#hist(residuals(antiox_model))
```

Generate a Type III Anova of model.

``` r
anova(antiox_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
    ## timepoint              0.27598 0.09199     3 281.77  57.0562 < 2.2e-16 ***
    ## site                   0.04320 0.02160     2 149.57  13.3983 4.441e-06 ***
    ## species                0.81874 0.40937     2 145.61 253.9005 < 2.2e-16 ***
    ## timepoint:site         0.02060 0.00343     6 280.88   2.1290   0.05020 .  
    ## timepoint:species      0.13378 0.02230     6 278.09  13.8293 9.338e-14 ***
    ## site:species           0.05596 0.01399     4 144.77   8.6768 2.632e-06 ***
    ## timepoint:site:species 0.03285 0.00274    12 277.56   1.6981   0.06681 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(antiox_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + cre.umol.mgafdw) ~ timepoint * site * species + (1 |  
    ##     colony_id_corr)
    ##    Data: master
    ## 
    ## REML criterion at convergence: -1104.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1453 -0.3846 -0.0599  0.3447  4.2274 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001030 0.03210 
    ##  Residual                   0.001612 0.04015 
    ## Number of obs: 403, groups:  colony_id_corr, 138
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                        4.102e-02  1.469e-02
    ## timepointtimepoint2                               -5.741e-02  1.886e-02
    ## timepointtimepoint3                                5.144e-02  3.340e-02
    ## timepointtimepoint4                                5.743e-02  2.111e-02
    ## siteMahana                                         1.339e-02  2.059e-02
    ## siteManava                                         4.891e-04  1.973e-02
    ## speciesPocillopora                                -1.748e-02  2.030e-02
    ## speciesPorites                                     1.443e-01  2.256e-02
    ## timepointtimepoint2:siteMahana                    -7.980e-03  3.116e-02
    ## timepointtimepoint3:siteMahana                     1.909e-02  3.844e-02
    ## timepointtimepoint4:siteMahana                    -3.066e-02  2.798e-02
    ## timepointtimepoint2:siteManava                     3.914e-02  2.806e-02
    ## timepointtimepoint3:siteManava                     8.460e-02  3.783e-02
    ## timepointtimepoint4:siteManava                    -2.012e-02  2.678e-02
    ## timepointtimepoint2:speciesPocillopora             4.533e-02  2.452e-02
    ## timepointtimepoint3:speciesPocillopora             3.973e-03  3.705e-02
    ## timepointtimepoint4:speciesPocillopora            -1.787e-02  2.609e-02
    ## timepointtimepoint2:speciesPorites                 9.737e-02  2.711e-02
    ## timepointtimepoint3:speciesPorites                 3.428e-02  3.805e-02
    ## timepointtimepoint4:speciesPorites                -5.862e-02  2.824e-02
    ## siteMahana:speciesPocillopora                     -1.241e-02  2.908e-02
    ## siteManava:speciesPocillopora                      2.094e-03  2.821e-02
    ## siteMahana:speciesPorites                         -7.241e-02  3.072e-02
    ## siteManava:speciesPorites                          3.102e-02  3.014e-02
    ## timepointtimepoint2:siteMahana:speciesPocillopora  1.996e-02  3.862e-02
    ## timepointtimepoint3:siteMahana:speciesPocillopora -3.869e-02  4.490e-02
    ## timepointtimepoint4:siteMahana:speciesPocillopora  5.111e-03  3.625e-02
    ## timepointtimepoint2:siteManava:speciesPocillopora -4.646e-02  3.612e-02
    ## timepointtimepoint3:siteManava:speciesPocillopora -9.930e-02  4.454e-02
    ## timepointtimepoint4:siteManava:speciesPocillopora  1.070e-02  3.510e-02
    ## timepointtimepoint2:siteMahana:speciesPorites     -7.363e-03  4.042e-02
    ## timepointtimepoint3:siteMahana:speciesPorites     -1.995e-02  4.603e-02
    ## timepointtimepoint4:siteMahana:speciesPorites      6.062e-02  3.786e-02
    ## timepointtimepoint2:siteManava:speciesPorites      2.570e-03  3.815e-02
    ## timepointtimepoint3:siteManava:speciesPorites     -4.254e-02  4.523e-02
    ## timepointtimepoint4:siteManava:speciesPorites      3.705e-02  3.687e-02
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                        3.042e+02   2.792 0.005565
    ## timepointtimepoint2                                3.089e+02  -3.043 0.002541
    ## timepointtimepoint3                                3.011e+02   1.540 0.124630
    ## timepointtimepoint4                                2.886e+02   2.721 0.006904
    ## siteMahana                                         3.178e+02   0.650 0.516180
    ## siteManava                                         3.010e+02   0.025 0.980239
    ## speciesPocillopora                                 3.069e+02  -0.861 0.389916
    ## speciesPorites                                     3.409e+02   6.394 5.31e-10
    ## timepointtimepoint2:siteMahana                     3.023e+02  -0.256 0.798065
    ## timepointtimepoint3:siteMahana                     3.107e+02   0.497 0.619834
    ## timepointtimepoint4:siteMahana                     3.080e+02  -1.096 0.273966
    ## timepointtimepoint2:siteManava                     3.016e+02   1.395 0.164043
    ## timepointtimepoint3:siteManava                     2.977e+02   2.236 0.026077
    ## timepointtimepoint4:siteManava                     2.840e+02  -0.751 0.453172
    ## timepointtimepoint2:speciesPocillopora             2.888e+02   1.849 0.065480
    ## timepointtimepoint3:speciesPocillopora             2.938e+02   0.107 0.914675
    ## timepointtimepoint4:speciesPocillopora             2.776e+02  -0.685 0.494035
    ## timepointtimepoint2:speciesPorites                 2.921e+02   3.591 0.000386
    ## timepointtimepoint3:speciesPorites                 2.946e+02   0.901 0.368355
    ## timepointtimepoint4:speciesPorites                 2.822e+02  -2.076 0.038812
    ## siteMahana:speciesPocillopora                      3.197e+02  -0.427 0.669795
    ## siteManava:speciesPocillopora                      3.081e+02   0.074 0.940887
    ## siteMahana:speciesPorites                          3.362e+02  -2.357 0.018982
    ## siteManava:speciesPorites                          3.308e+02   1.029 0.304118
    ## timepointtimepoint2:siteMahana:speciesPocillopora  2.890e+02   0.517 0.605605
    ## timepointtimepoint3:siteMahana:speciesPocillopora  2.982e+02  -0.862 0.389576
    ## timepointtimepoint4:siteMahana:speciesPocillopora  2.899e+02   0.141 0.887962
    ## timepointtimepoint2:siteManava:speciesPocillopora  2.870e+02  -1.286 0.199370
    ## timepointtimepoint3:siteManava:speciesPocillopora  2.885e+02  -2.229 0.026558
    ## timepointtimepoint4:siteManava:speciesPocillopora  2.742e+02   0.305 0.760640
    ## timepointtimepoint2:siteMahana:speciesPorites      2.904e+02  -0.182 0.855580
    ## timepointtimepoint3:siteMahana:speciesPorites      2.989e+02  -0.433 0.664988
    ## timepointtimepoint4:siteMahana:speciesPorites      2.925e+02   1.601 0.110433
    ## timepointtimepoint2:siteManava:speciesPorites      2.894e+02   0.067 0.946344
    ## timepointtimepoint3:siteManava:speciesPorites      2.892e+02  -0.941 0.347708
    ## timepointtimepoint4:siteManava:speciesPorites      2.767e+02   1.005 0.315864
    ##                                                      
    ## (Intercept)                                       ** 
    ## timepointtimepoint2                               ** 
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                               ** 
    ## siteMahana                                           
    ## siteManava                                           
    ## speciesPocillopora                                   
    ## speciesPorites                                    ***
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                    *  
    ## timepointtimepoint4:siteManava                       
    ## timepointtimepoint2:speciesPocillopora            .  
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                ***
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                *  
    ## siteMahana:speciesPocillopora                        
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                         *  
    ## siteManava:speciesPorites                            
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora    
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora    
    ## timepointtimepoint3:siteManava:speciesPocillopora *  
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites        
    ## timepointtimepoint3:siteMahana:speciesPorites        
    ## timepointtimepoint4:siteMahana:speciesPorites        
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites        
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, run a model for site and time effects
within each species.

*Acropora*:

`antiox_model_acr<-lmer(log(1+cre.umol.mgafdw)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))`

``` r
antiox_model_acr<-lmer(log(1+cre.umol.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))
qqPlot(residuals(antiox_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ## 218 215 
    ##  63  60

Generate a Type III Anova of model.

``` r
anova(antiox_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      0.140750 0.046917     3 80.208 37.7194 2.625e-15 ***
    ## site           0.007394 0.003697     2 59.382  2.9724   0.05885 .  
    ## timepoint:site 0.022442 0.003740     6 79.584  3.0071   0.01062 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(antiox_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + cre.umol.mgafdw) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Acropora")
    ## 
    ## REML criterion at convergence: -324.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2683 -0.3726 -0.0835  0.2970  5.3691 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0002775 0.01666 
    ##  Residual                   0.0012438 0.03527 
    ## Number of obs: 107, groups:  colony_id_corr, 48
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     0.0408642  0.0112351 93.1264758   3.637
    ## timepointtimepoint2            -0.0571109  0.0160229 79.6070126  -3.564
    ## timepointtimepoint3             0.0527673  0.0283164 83.5042341   1.863
    ## timepointtimepoint4             0.0565369  0.0181059 76.7962616   3.123
    ## siteMahana                      0.0143213  0.0158545 93.8310837   0.903
    ## siteManava                      0.0003999  0.0150733 93.0795925   0.027
    ## timepointtimepoint2:siteMahana -0.0081308  0.0265196 79.9725342  -0.307
    ## timepointtimepoint3:siteMahana  0.0168448  0.0324498 84.3289205   0.519
    ## timepointtimepoint4:siteMahana -0.0298523  0.0237479 80.3775155  -1.257
    ## timepointtimepoint2:siteManava  0.0442807  0.0238787 80.0601914   1.854
    ## timepointtimepoint3:siteManava  0.0840321  0.0321672 81.9741951   2.612
    ## timepointtimepoint4:siteManava -0.0187899  0.0230431 75.4599620  -0.815
    ##                                Pr(>|t|)    
    ## (Intercept)                    0.000452 ***
    ## timepointtimepoint2            0.000620 ***
    ## timepointtimepoint3            0.065909 .  
    ## timepointtimepoint4            0.002528 ** 
    ## siteMahana                     0.368681    
    ## siteManava                     0.978890    
    ## timepointtimepoint2:siteMahana 0.759950    
    ## timepointtimepoint3:siteMahana 0.605049    
    ## timepointtimepoint4:siteMahana 0.212376    
    ## timepointtimepoint2:siteManava 0.067364 .  
    ## timepointtimepoint3:siteManava 0.010693 *  
    ## timepointtimepoint4:siteManava 0.417396    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.601                                                    
    ## tmpnttmpnt3      -0.327  0.251                                             
    ## tmpnttmpnt4      -0.510  0.365  0.203                                      
    ## siteMahana       -0.709  0.426  0.232  0.361                               
    ## siteManava       -0.745  0.448  0.244  0.380  0.528                        
    ## tmpnttmpnt2:stMh  0.363 -0.604 -0.152 -0.221 -0.507 -0.271                 
    ## tmpnttmpnt3:stMh  0.285 -0.219 -0.873 -0.177 -0.422 -0.213  0.269          
    ## tmpnttmpnt4:stMh  0.389 -0.278 -0.155 -0.762 -0.574 -0.290  0.351          
    ## tmpnttmpnt2:stMn  0.403 -0.671 -0.168 -0.245 -0.286 -0.537  0.405          
    ## tmpnttmpnt3:stMn  0.288 -0.221 -0.880 -0.179 -0.204 -0.388  0.133          
    ## tmpnttmpnt4:stMn  0.401 -0.287 -0.160 -0.786 -0.284 -0.540  0.173          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.291                                            
    ## tmpnttmpnt2:stMn  0.147            0.187                           
    ## tmpnttmpnt3:stMn  0.768            0.136            0.268          
    ## tmpnttmpnt4:stMn  0.139            0.599            0.347          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.257

*Pocillopora*:

`antiox_model_poc<-lmer(log(1+cre.umol.mgafdw)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))`

``` r
antiox_model_poc<-lmer(log(1+cre.umol.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))
qqPlot(residuals(antiox_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## 353 363 
    ## 127 137

``` r
#hist(residuals(antiox_model_poc))
```

Generate a Type III Anova of model.

``` r
anova(antiox_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq   Mean Sq NumDF   DenDF F value  Pr(>F)    
    ## timepoint      0.069336 0.0231119     3 108.525 37.6135 < 2e-16 ***
    ## site           0.001272 0.0006360     2  39.426  1.0351 0.36466    
    ## timepoint:site 0.006856 0.0011426     6 108.377  1.8596 0.09433 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(antiox_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + cre.umol.mgafdw) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -606.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2773 -0.5355 -0.1382  0.3644  4.3232 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 2.529e-05 0.005029
    ##  Residual                   6.145e-04 0.024788
    ## Number of obs: 153, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      0.023876   0.007013 140.719092   3.405
    ## timepointtimepoint2             -0.012021   0.009567 106.429498  -1.257
    ## timepointtimepoint3              0.054044   0.009750 108.704272   5.543
    ## timepointtimepoint4              0.039229   0.009405 104.318435   4.171
    ## siteMahana                       0.001183   0.010358 140.803717   0.114
    ## siteManava                       0.001920   0.010122 140.743328   0.190
    ## timepointtimepoint2:siteMahana   0.011368   0.013857 109.040954   0.820
    ## timepointtimepoint3:siteMahana  -0.018809   0.014101 109.169296  -1.334
    ## timepointtimepoint4:siteMahana  -0.025293   0.014014 108.706268  -1.805
    ## timepointtimepoint2:siteManava  -0.006899   0.013811 109.068527  -0.500
    ## timepointtimepoint3:siteManava  -0.014208   0.014248 110.592686  -0.997
    ## timepointtimepoint4:siteManava  -0.009408   0.013839 107.147862  -0.680
    ##                                Pr(>|t|)    
    ## (Intercept)                    0.000864 ***
    ## timepointtimepoint2            0.211684    
    ## timepointtimepoint3            2.10e-07 ***
    ## timepointtimepoint4            6.29e-05 ***
    ## siteMahana                     0.909241    
    ## siteManava                     0.849814    
    ## timepointtimepoint2:siteMahana 0.413792    
    ## timepointtimepoint3:siteMahana 0.185023    
    ## timepointtimepoint4:siteMahana 0.073871 .  
    ## timepointtimepoint2:siteManava 0.618422    
    ## timepointtimepoint3:siteManava 0.320841    
    ## timepointtimepoint4:siteManava 0.498058    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.708                                                    
    ## tmpnttmpnt3      -0.695  0.510                                             
    ## tmpnttmpnt4      -0.720  0.528  0.518                                      
    ## siteMahana       -0.677  0.479  0.471  0.488                               
    ## siteManava       -0.693  0.491  0.482  0.499  0.469                        
    ## tmpnttmpnt2:stMh  0.489 -0.690 -0.352 -0.365 -0.724 -0.339                 
    ## tmpnttmpnt3:stMh  0.481 -0.352 -0.691 -0.358 -0.711 -0.333  0.532          
    ## tmpnttmpnt4:stMh  0.483 -0.354 -0.348 -0.671 -0.716 -0.335  0.535          
    ## tmpnttmpnt2:stMn  0.491 -0.693 -0.353 -0.366 -0.332 -0.709  0.478          
    ## tmpnttmpnt3:stMn  0.476 -0.349 -0.684 -0.355 -0.322 -0.687  0.241          
    ## tmpnttmpnt4:stMn  0.489 -0.359 -0.352 -0.680 -0.331 -0.707  0.248          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.527                                            
    ## tmpnttmpnt2:stMn  0.244            0.245                           
    ## tmpnttmpnt3:stMn  0.473            0.238            0.503          
    ## tmpnttmpnt4:stMn  0.244            0.456            0.518          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.503

*Porites*:

`antiox_model_por<-lmer(log(1+cre.umol.mgafdw)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))`

``` r
antiox_model_por<-lmer(log(1+cre.umol.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))
qqPlot(residuals(antiox_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## 107 273 
    ##  24  68

``` r
#hist(residuals(antiox_model_por))
```

Generate a Type III Anova of model.

``` r
anova(antiox_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      0.204764 0.068255     3 95.015  23.436 1.936e-11 ***
    ## site           0.070881 0.035441     2 42.574  12.169 6.608e-05 ***
    ## timepoint:site 0.023188 0.003865     6 94.960   1.327    0.2528    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(antiox_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + cre.umol.mgafdw) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Porites")
    ## 
    ## REML criterion at convergence: -306.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.51585 -0.62220 -0.01774  0.57557  2.85135 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.002789 0.05281 
    ##  Residual                   0.002912 0.05397 
    ## Number of obs: 143, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     1.848e-01  2.444e-02  1.212e+02   7.560
    ## timepointtimepoint2             4.048e-02  2.633e-02  9.696e+01   1.537
    ## timepointtimepoint3             8.621e-02  2.461e-02  9.595e+01   3.503
    ## timepointtimepoint4            -8.111e-04  2.536e-02  9.644e+01  -0.032
    ## siteMahana                     -5.906e-02  3.274e-02  1.144e+02  -1.804
    ## siteManava                      2.997e-02  3.272e-02  1.144e+02   0.916
    ## timepointtimepoint2:siteMahana -1.625e-02  3.478e-02  9.620e+01  -0.467
    ## timepointtimepoint3:siteMahana  7.517e-04  3.417e-02  9.567e+01   0.022
    ## timepointtimepoint4:siteMahana  3.101e-02  3.447e-02  9.637e+01   0.900
    ## timepointtimepoint2:siteManava  4.357e-02  3.494e-02  9.685e+01   1.247
    ## timepointtimepoint3:siteManava  4.308e-02  3.345e-02  9.507e+01   1.288
    ## timepointtimepoint4:siteManava  1.732e-02  3.420e-02  9.490e+01   0.506
    ##                                Pr(>|t|)    
    ## (Intercept)                    8.56e-12 ***
    ## timepointtimepoint2            0.127504    
    ## timepointtimepoint3            0.000701 ***
    ## timepointtimepoint4            0.974546    
    ## siteMahana                     0.073849 .  
    ## siteManava                     0.361675    
    ## timepointtimepoint2:siteMahana 0.641389    
    ## timepointtimepoint3:siteMahana 0.982497    
    ## timepointtimepoint4:siteMahana 0.370473    
    ## timepointtimepoint2:siteManava 0.215473    
    ## timepointtimepoint3:siteManava 0.200967    
    ## timepointtimepoint4:siteManava 0.613744    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.637                                                    
    ## tmpnttmpnt3      -0.684  0.632                                             
    ## tmpnttmpnt4      -0.663  0.607  0.658                                      
    ## siteMahana       -0.747  0.475  0.511  0.495                               
    ## siteManava       -0.747  0.476  0.511  0.495  0.558                        
    ## tmpnttmpnt2:stMh  0.482 -0.757 -0.479 -0.459 -0.611 -0.360                 
    ## tmpnttmpnt3:stMh  0.493 -0.455 -0.720 -0.474 -0.618 -0.368  0.574          
    ## tmpnttmpnt4:stMh  0.488 -0.446 -0.484 -0.736 -0.618 -0.364  0.571          
    ## tmpnttmpnt2:stMn  0.480 -0.754 -0.476 -0.457 -0.358 -0.613  0.571          
    ## tmpnttmpnt3:stMn  0.503 -0.465 -0.736 -0.484 -0.376 -0.634  0.352          
    ## tmpnttmpnt4:stMn  0.491 -0.450 -0.488 -0.741 -0.367 -0.614  0.341          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.600                                            
    ## tmpnttmpnt2:stMn  0.343            0.336                           
    ## tmpnttmpnt3:stMn  0.530            0.356            0.594          
    ## tmpnttmpnt4:stMn  0.352            0.545            0.571          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.607

## Symbiont Densities

### Plot: Surface Area

View by site and timepoint for each species.

``` r
symbplot_cm2<-master %>%
  filter(!is.na(cells.cm2)) %>%
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, cells.cm2)%>%
  
  ggplot(., aes(x = site_code, y = cells.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Symbiont Cells cm"^-2))))+
    theme_classic() + 
    theme(
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Protein

View by site and timepoint for each species.

``` r
symbplot_prot<-master %>%
  filter(!is.na(cells.mgprot)) %>%
   filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, cells.mgprot)%>%
  
  ggplot(., aes(x = site_code, y = cells.mgprot, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Symbiont Cells mg protein"^-1))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Biomass

View by site and timepoint for each species.

``` r
symbplot_afdw<-master %>%
  filter(!is.na(cells.mgAFDW)) %>%
    filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, cells.mgAFDW)%>%
  
  ggplot(., aes(x = site_code, y = cells.mgAFDW, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Symbiont Cells mg afdw"^-1))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots for each normalization together.

``` r
symbcells_figure<-plot_grid(symbplot_prot, symbplot_afdw, symbplot_cm2, labels = c("", "", ""), ncol=3, nrow=1, rel_widths = c(.8, .8, 1), label_size = 20, label_x = 0.1)

ggsave(filename="Figures/CellDensity_Figure.pdf", plot=symbcells_figure, dpi=300, width=24, height=5, units="in")
```

### Analysis

Build a mixed model for univariate analysis and examine data
distribution.

`density_model<-lmer(cells.mgAFDW~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
density_model<-lmer(cells.mgAFDW~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(cells.mgAFDW>0))
qqPlot(residuals(density_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

    ## 421 215 
    ## 407 206

Residuals are not normally distributed. Attempt with log transformation.

`density_model<-lmer(log(cells.mgAFDW)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
density_model<-lmer(log(cells.mgAFDW)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(cells.mgAFDW>0))
qqPlot(residuals(density_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

    ## 167 160 
    ## 161 154

Residuals are improved, but we need to revisit this model to improve
fit.

Generate a Type III Anova of model.

``` r
anova(density_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
    ## timepoint              55.504 18.5013     3 304.91 154.1478 < 2.2e-16 ***
    ## site                    2.683  1.3415     2 128.13  11.1772 3.355e-05 ***
    ## species                54.347 27.1736     2 122.19 226.4034 < 2.2e-16 ***
    ## timepoint:site          4.684  0.7806     6 302.83   6.5036 1.806e-06 ***
    ## timepoint:species       7.658  1.2763     6 297.83  10.6334 1.049e-10 ***
    ## site:species            2.275  0.5688     4 120.26   4.7391  0.001383 ** 
    ## timepoint:site:species  2.740  0.2284    12 296.20   1.9027  0.033649 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of time and site within
each species.

*Acropora*

`density_model_acr<-lmer(log(cells.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(subset=c(cells.mgAFDW>0 & species=="Acropora"))`

``` r
density_model_acr<-lmer(log(cells.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(cells.mgAFDW>0 & species=="Acropora"))
qqPlot(residuals(density_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

    ## 215 421 
    ##  61 112

Generate a Type III Anova of model.

``` r
anova(density_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      27.2401  9.0800     3   102 74.2792 < 2.2e-16 ***
    ## site            0.1990  0.0995     2   102  0.8142  0.445872    
    ## timepoint:site  2.8186  0.4698     6   102  3.8429  0.001699 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

*Pocillopora*

`density_model_poc<-lmer(log(cells.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))`

``` r
density_model_poc<-lmer(log(cells.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(cells.mgAFDW>0 & species=="Pocillopora"))
qqPlot(residuals(density_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

    ## 167  76 
    ##  70  34

Generate a Type III Anova of model.

``` r
anova(density_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      22.3771  7.4590     3 104.546 48.8462 < 2.2e-16 ***
    ## site            0.4200  0.2100     2  35.554  1.3752 0.2659206    
    ## timepoint:site  3.8412  0.6402     6 104.449  4.1924 0.0008103 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

*Porites*

`density_model_por<-lmer(log(cells.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))`

``` r
density_model_por<-lmer(log(cells.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(cells.mgAFDW>0 & species=="Porites"))
qqPlot(residuals(density_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

    ## 178 296 
    ##  36  92

Generate a Type III Anova of model.

``` r
anova(density_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      12.4289  4.1430     3 105.366 48.5681 < 2.2e-16 ***
    ## site            3.5955  1.7978     2  36.864 21.0752  7.89e-07 ***
    ## timepoint:site  1.6595  0.2766     6 105.276  3.2423   0.00579 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Biomass

### Plot Host Biomass

View by site and timepoint for each species in the host.

``` r
hAFDWplot<-master %>%
  filter(!is.na(Host_AFDW.mg.cm2)) %>%
  filter(!is.na(species))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  select(colony_id_corr, species, site_code, timepoint, Host_AFDW.mg.cm2)%>%
  
  ggplot(., aes(x = site_code, y = Host_AFDW.mg.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Host Biomass (mg AFDW cm"^-2, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot Symbiont biomass

View by site and timepoint for each species in the symbiont

``` r
sAFDWplot<-master %>%
  filter(!is.na(Sym_AFDW.mg.cm2)) %>%
  filter(!is.na(species))%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  select(colony_id_corr, species, site_code, timepoint, Sym_AFDW.mg.cm2)%>%
  
  ggplot(., aes(x = site_code, y = Sym_AFDW.mg.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Symbiont Biomass (mg AFDW cm"^-2, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot S:H Ratio

Calculate biomass as a ratio of host to symbiont and plot.

``` r
#master <- master %>% 
  #mutate(Ratio_AFDW.mg.cm2=Sym_AFDW.mg.cm2/(Sym_AFDW.mg.cm2+Host_AFDW.mg.cm2))

ratioAFDWplot<-master %>%
  filter(!is.na(Host_AFDW.mg.cm2)) %>%
  filter(!is.na(species))%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  select(colony_id_corr, species, site_code, timepoint, Ratio_AFDW.mg.cm2)%>%
  
  ggplot(., aes(x = site_code, y = Ratio_AFDW.mg.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Symbiont : Holobiont Biomass"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots together.

``` r
afdw_figure<-plot_grid(hAFDWplot, sAFDWplot, ratioAFDWplot, labels = c("", "", ""), ncol=3, nrow=1, rel_widths = c(.8, .8, 1), label_size = 20, label_x = 0.1)

ggsave(filename="Figures/Biomass_Figure.pdf", plot=afdw_figure, dpi=300, width=24, height=5, units="in")
```

### Analysis

Build a mixed model for univariate analysis and examine data
distribution.

#### Host Biomass:

`hAFDW_model<-lmer(Host_AFDW.mg.cm2~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0))`

``` r
hAFDW_model<-lmer(Host_AFDW.mg.cm2~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0))
qqPlot(residuals(hAFDW_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

    ##  88 447 
    ##  86 436

Residuals are not normally distributed. Attempt with log transformation.

`hAFDW_model<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0))`

``` r
hAFDW_model<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0))
qqPlot(residuals(hAFDW_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

    ##  91 439 
    ##  89 429

Residuals are improved.

Generate a Type III Anova of model.

``` r
anova(hAFDW_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF   F value    Pr(>F)    
    ## timepoint               1.408   0.469     3 327.69   11.6985 2.664e-07 ***
    ## site                    0.966   0.483     2 139.22   12.0418 1.501e-05 ***
    ## species                88.762  44.381     2 130.36 1106.0647 < 2.2e-16 ***
    ## timepoint:site          0.843   0.140     6 325.05    3.5014  0.002272 ** 
    ## timepoint:species       0.456   0.076     6 320.06    1.8926  0.081582 .  
    ## site:species            0.349   0.087     4 127.92    2.1769  0.075210 .  
    ## timepoint:site:species  0.623   0.052    12 317.98    1.2942  0.220301    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(hAFDW_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Host_AFDW.mg.cm2) ~ timepoint * site * species + (1 |  
    ##     colony_id_corr)
    ##    Data: master
    ##  Subset: c(Host_AFDW.mg.cm2 > 0)
    ## 
    ## REML criterion at convergence: -49
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6070 -0.5238  0.0300  0.5556  3.4919 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001514 0.03891 
    ##  Residual                   0.040125 0.20031 
    ## Number of obs: 436, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         0.624586   0.056591
    ## timepointtimepoint2                                -0.070285   0.087435
    ## timepointtimepoint3                                -0.037026   0.154070
    ## timepointtimepoint4                                 0.062195   0.099508
    ## siteMahana                                          0.107289   0.080024
    ## siteManava                                          0.163037   0.077316
    ## speciesPocillopora                                  0.318284   0.077320
    ## speciesPorites                                      1.387143   0.078588
    ## timepointtimepoint2:siteMahana                     -0.029306   0.132464
    ## timepointtimepoint3:siteMahana                     -0.099857   0.176019
    ## timepointtimepoint4:siteMahana                     -0.036484   0.128274
    ## timepointtimepoint2:siteManava                      0.058610   0.124263
    ## timepointtimepoint3:siteManava                     -0.172119   0.175941
    ## timepointtimepoint4:siteManava                     -0.034701   0.127577
    ## timepointtimepoint2:speciesPocillopora             -0.024822   0.114858
    ## timepointtimepoint3:speciesPocillopora             -0.017837   0.171794
    ## timepointtimepoint4:speciesPocillopora              0.169023   0.123499
    ## timepointtimepoint2:speciesPorites                 -0.261977   0.117844
    ## timepointtimepoint3:speciesPorites                 -0.144946   0.171129
    ## timepointtimepoint4:speciesPorites                 -0.042793   0.126004
    ## siteMahana:speciesPocillopora                      -0.217321   0.111273
    ## siteManava:speciesPocillopora                      -0.070119   0.107377
    ## siteMahana:speciesPorites                          -0.380122   0.112157
    ## siteManava:speciesPorites                          -0.255111   0.108293
    ## timepointtimepoint2:siteMahana:speciesPocillopora   0.208143   0.169913
    ## timepointtimepoint3:siteMahana:speciesPocillopora   0.187541   0.207238
    ## timepointtimepoint4:siteMahana:speciesPocillopora  -0.150025   0.168110
    ## timepointtimepoint2:siteManava:speciesPocillopora  -0.009893   0.162900
    ## timepointtimepoint3:siteManava:speciesPocillopora  -0.022624   0.207574
    ## timepointtimepoint4:siteManava:speciesPocillopora  -0.006294   0.165535
    ## timepointtimepoint2:siteMahana:speciesPorites       0.364555   0.171945
    ## timepointtimepoint3:siteMahana:speciesPorites       0.321444   0.208095
    ## timepointtimepoint4:siteMahana:speciesPorites       0.003822   0.169960
    ## timepointtimepoint2:siteManava:speciesPorites       0.237505   0.165019
    ## timepointtimepoint3:siteManava:speciesPorites       0.343045   0.205618
    ## timepointtimepoint4:siteManava:speciesPorites       0.072952   0.168205
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       398.985817  11.037  < 2e-16
    ## timepointtimepoint2                               345.366271  -0.804 0.422033
    ## timepointtimepoint3                               384.141718  -0.240 0.810210
    ## timepointtimepoint4                               343.191132   0.625 0.532368
    ## siteMahana                                        399.246415   1.341 0.180776
    ## siteManava                                        399.023072   2.109 0.035594
    ## speciesPocillopora                                398.874207   4.116 4.68e-05
    ## speciesPorites                                    399.040537  17.651  < 2e-16
    ## timepointtimepoint2:siteMahana                    344.184571  -0.221 0.825040
    ## timepointtimepoint3:siteMahana                    381.790933  -0.567 0.570836
    ## timepointtimepoint4:siteMahana                    353.596044  -0.284 0.776253
    ## timepointtimepoint2:siteManava                    345.374551   0.472 0.637468
    ## timepointtimepoint3:siteManava                    375.355187  -0.978 0.328566
    ## timepointtimepoint4:siteManava                    334.739219  -0.272 0.785789
    ## timepointtimepoint2:speciesPocillopora            322.341066  -0.216 0.829037
    ## timepointtimepoint3:speciesPocillopora            370.846232  -0.104 0.917360
    ## timepointtimepoint4:speciesPocillopora            322.971847   1.369 0.172068
    ## timepointtimepoint2:speciesPorites                328.670644  -2.223 0.026890
    ## timepointtimepoint3:speciesPorites                370.539033  -0.847 0.397543
    ## timepointtimepoint4:speciesPorites                327.661184  -0.340 0.734361
    ## siteMahana:speciesPocillopora                     399.139594  -1.953 0.051512
    ## siteManava:speciesPocillopora                     398.890438  -0.653 0.514119
    ## siteMahana:speciesPorites                         399.208167  -3.389 0.000771
    ## siteManava:speciesPorites                         398.978762  -2.356 0.018968
    ## timepointtimepoint2:siteMahana:speciesPocillopora 324.407032   1.225 0.221464
    ## timepointtimepoint3:siteMahana:speciesPocillopora 362.142850   0.905 0.366092
    ## timepointtimepoint4:siteMahana:speciesPocillopora 329.747769  -0.892 0.372817
    ## timepointtimepoint2:siteManava:speciesPocillopora 322.483853  -0.061 0.951612
    ## timepointtimepoint3:siteManava:speciesPocillopora 357.468834  -0.109 0.913271
    ## timepointtimepoint4:siteManava:speciesPocillopora 316.811619  -0.038 0.969694
    ## timepointtimepoint2:siteMahana:speciesPorites     327.331188   2.120 0.034742
    ## timepointtimepoint3:siteMahana:speciesPorites     362.889632   1.545 0.123290
    ## timepointtimepoint4:siteMahana:speciesPorites     332.253964   0.022 0.982074
    ## timepointtimepoint2:siteManava:speciesPorites     325.717219   1.439 0.151037
    ## timepointtimepoint3:siteManava:speciesPorites     356.326313   1.668 0.096122
    ## timepointtimepoint4:siteManava:speciesPorites     320.665335   0.434 0.664790
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                                  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## siteMahana                                           
    ## siteManava                                        *  
    ## speciesPocillopora                                ***
    ## speciesPorites                                    ***
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                *  
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## siteMahana:speciesPocillopora                     .  
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                         ***
    ## siteManava:speciesPorites                         *  
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora    
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora    
    ## timepointtimepoint3:siteManava:speciesPocillopora    
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites     *  
    ## timepointtimepoint3:siteMahana:speciesPorites        
    ## timepointtimepoint4:siteMahana:speciesPorites        
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites     .  
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of site and time within
each species.

*Acropora*

`hAFDW_model_acr<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0 & species=="Acropora"))`

``` r
hAFDW_model_acr<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0 & species=="Acropora"))
qqPlot(residuals(hAFDW_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

    ## 315 320 
    ##  83  88

Generate a Type III Anova of model.

``` r
anova(hAFDW_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
    ## timepoint      0.283771 0.094590     3 86.423  3.4813 0.01930 *
    ## site           0.194289 0.097145     2 43.480  3.5753 0.03652 *
    ## timepoint:site 0.086992 0.014499     6 84.743  0.5336 0.78129  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(hAFDW_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Host_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Host_AFDW.mg.cm2 > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: -48.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.41422 -0.60879 -0.04488  0.45923  2.86388 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001047 0.03236 
    ##  Residual                   0.027171 0.16484 
    ## Number of obs: 114, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.62458    0.04659 101.84169  13.407   <2e-16
    ## timepointtimepoint2             -0.07028    0.07196  82.98367  -0.977   0.3316
    ## timepointtimepoint3             -0.03710    0.12682  95.76578  -0.293   0.7705
    ## timepointtimepoint4              0.06216    0.08190  82.28297   0.759   0.4500
    ## siteMahana                       0.10730    0.06588 101.90069   1.629   0.1064
    ## siteManava                       0.16302    0.06365 101.85036   2.561   0.0119
    ## timepointtimepoint2:siteMahana  -0.02941    0.10902  82.60257  -0.270   0.7880
    ## timepointtimepoint3:siteMahana  -0.09984    0.14488  94.95078  -0.689   0.4924
    ## timepointtimepoint4:siteMahana  -0.03649    0.10557  85.58698  -0.346   0.7305
    ## timepointtimepoint2:siteManava   0.05868    0.10227  82.97298   0.574   0.5677
    ## timepointtimepoint3:siteManava  -0.17205    0.14481  92.72816  -1.188   0.2378
    ## timepointtimepoint4:siteManava  -0.03467    0.10499  79.66163  -0.330   0.7421
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                        
    ## siteManava                     *  
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.629                                                    
    ## tmpnttmpnt3      -0.354  0.235                                             
    ## tmpnttmpnt4      -0.548  0.357  0.202                                      
    ## siteMahana       -0.707  0.445  0.250  0.387                               
    ## siteManava       -0.732  0.460  0.259  0.401  0.518                        
    ## tmpnttmpnt2:stMh  0.415 -0.660 -0.155 -0.236 -0.585 -0.304                 
    ## tmpnttmpnt3:stMh  0.310 -0.205 -0.875 -0.177 -0.442 -0.227  0.271          
    ## tmpnttmpnt4:stMh  0.425 -0.277 -0.157 -0.776 -0.607 -0.311  0.366          
    ## tmpnttmpnt2:stMn  0.442 -0.704 -0.165 -0.251 -0.313 -0.603  0.464          
    ## tmpnttmpnt3:stMn  0.310 -0.205 -0.876 -0.177 -0.219 -0.424  0.136          
    ## tmpnttmpnt4:stMn  0.427 -0.279 -0.157 -0.780 -0.302 -0.585  0.184          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.279                                            
    ## tmpnttmpnt2:stMn  0.145            0.195                           
    ## tmpnttmpnt3:stMn  0.767            0.137            0.270          
    ## tmpnttmpnt4:stMn  0.138            0.605            0.366          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.258

*Pocillopora*

`hAFDW_model_poc<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0 & species=="Pocillopora"))`

``` r
hAFDW_model_poc<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0 & species=="Pocillopora"))
qqPlot(residuals(hAFDW_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

    ##  56 351 
    ##  16 126

Generate a Type III Anova of model.

``` r
anova(hAFDW_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      1.24961 0.41654     3 112.878 20.9783 7.136e-11 ***
    ## site           0.28615 0.14307     2  40.171  7.2057  0.002117 ** 
    ## timepoint:site 0.83511 0.13918     6 112.801  7.0098 2.347e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(hAFDW_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Host_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Host_AFDW.mg.cm2 > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -109.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1567 -0.4548  0.0945  0.5712  2.5886 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.003899 0.06244 
    ##  Residual                   0.019856 0.14091 
    ## Number of obs: 163, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.94287    0.03980 140.35548  23.693  < 2e-16
    ## timepointtimepoint2             -0.09300    0.05248 110.38394  -1.772   0.0791
    ## timepointtimepoint3             -0.05701    0.05363 111.90648  -1.063   0.2901
    ## timepointtimepoint4              0.23122    0.05145 108.95300   4.494 1.75e-05
    ## siteMahana                      -0.10355    0.05831 142.21208  -1.776   0.0779
    ## siteManava                       0.09292    0.05628 140.35548   1.651   0.1010
    ## timepointtimepoint2:siteMahana   0.17025    0.07506 111.69419   2.268   0.0252
    ## timepointtimepoint3:siteMahana   0.08328    0.07722 112.32423   1.079   0.2831
    ## timepointtimepoint4:siteMahana  -0.18970    0.07668 112.00509  -2.474   0.0149
    ## timepointtimepoint2:siteManava   0.04931    0.07421 110.39036   0.664   0.5078
    ## timepointtimepoint3:siteManava  -0.19122    0.07790 114.10635  -2.455   0.0156
    ## timepointtimepoint4:siteManava  -0.04096    0.07435 110.95362  -0.551   0.5828
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            .  
    ## timepointtimepoint3               
    ## timepointtimepoint4            ***
    ## siteMahana                     .  
    ## siteManava                        
    ## timepointtimepoint2:siteMahana *  
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana *  
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava *  
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.634                                                    
    ## tmpnttmpnt3      -0.620  0.470                                             
    ## tmpnttmpnt4      -0.646  0.490  0.480                                      
    ## siteMahana       -0.682  0.433  0.423  0.441                               
    ## siteManava       -0.707  0.448  0.439  0.457  0.483                        
    ## tmpnttmpnt2:stMh  0.443 -0.699 -0.328 -0.343 -0.658 -0.313                 
    ## tmpnttmpnt3:stMh  0.431 -0.326 -0.695 -0.333 -0.637 -0.305  0.494          
    ## tmpnttmpnt4:stMh  0.434 -0.329 -0.322 -0.671 -0.641 -0.307  0.498          
    ## tmpnttmpnt2:stMn  0.448 -0.707 -0.332 -0.347 -0.306 -0.634  0.494          
    ## tmpnttmpnt3:stMn  0.427 -0.323 -0.688 -0.330 -0.291 -0.604  0.226          
    ## tmpnttmpnt4:stMn  0.447 -0.339 -0.332 -0.692 -0.305 -0.633  0.237          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.487                                            
    ## tmpnttmpnt2:stMn  0.231            0.233                           
    ## tmpnttmpnt3:stMn  0.478            0.222            0.457          
    ## tmpnttmpnt4:stMn  0.231            0.464            0.479          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.462

*Porites*

`hAFDW_model_por<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0 & species=="Porites"))`

``` r
hAFDW_model_por<-lmer(log(1+Host_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Host_AFDW.mg.cm2>0 & species=="Porites"))
qqPlot(residuals(hAFDW_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

    ##  91 439 
    ##  11 152

Generate a Type III Anova of model.

``` r
anova(hAFDW_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
    ## timepoint      0.48562 0.16187     3   147  2.3337 0.076412 . 
    ## site           0.88675 0.44337     2   147  6.3919 0.002179 **
    ## timepoint:site 0.68469 0.11411     6   147  1.6451 0.138687   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(hAFDW_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Host_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Host_AFDW.mg.cm2 > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: 55.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.83140 -0.58068 -0.00518  0.62543  2.76146 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.06936  0.2634  
    ## Number of obs: 159, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      2.01298    0.07039 147.00000  28.598  < 2e-16
    ## timepointtimepoint2             -0.33393    0.10361 147.00000  -3.223  0.00156
    ## timepointtimepoint3             -0.18322    0.09787 147.00000  -1.872  0.06319
    ## timepointtimepoint4              0.01857    0.10144 147.00000   0.183  0.85498
    ## siteMahana                      -0.27282    0.10144 147.00000  -2.689  0.00799
    ## siteManava                      -0.09332    0.09787 147.00000  -0.954  0.34189
    ## timepointtimepoint2:siteMahana   0.33565    0.14386 147.00000   2.333  0.02099
    ## timepointtimepoint3:siteMahana   0.22350    0.14567 147.00000   1.534  0.12711
    ## timepointtimepoint4:siteMahana  -0.03111    0.14631 147.00000  -0.213  0.83190
    ## timepointtimepoint2:siteManava   0.29874    0.14253 147.00000   2.096  0.03779
    ## timepointtimepoint3:siteManava   0.17337    0.13978 147.00000   1.240  0.21685
    ## timepointtimepoint4:siteManava   0.03861    0.14386 147.00000   0.268  0.78879
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            ** 
    ## timepointtimepoint3            .  
    ## timepointtimepoint4               
    ## siteMahana                     ** 
    ## siteManava                        
    ## timepointtimepoint2:siteMahana *  
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava *  
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.679                                                    
    ## tmpnttmpnt3      -0.719  0.489                                             
    ## tmpnttmpnt4      -0.694  0.471  0.499                                      
    ## siteMahana       -0.694  0.471  0.499  0.481                               
    ## siteManava       -0.719  0.489  0.517  0.499  0.499                        
    ## tmpnttmpnt2:stMh  0.489 -0.720 -0.352 -0.340 -0.705 -0.352                 
    ## tmpnttmpnt3:stMh  0.483 -0.328 -0.672 -0.335 -0.696 -0.348  0.491          
    ## tmpnttmpnt4:stMh  0.481 -0.327 -0.346 -0.693 -0.693 -0.346  0.489          
    ## tmpnttmpnt2:stMn  0.494 -0.727 -0.355 -0.343 -0.343 -0.687  0.524          
    ## tmpnttmpnt3:stMn  0.504 -0.342 -0.700 -0.349 -0.349 -0.700  0.246          
    ## tmpnttmpnt4:stMn  0.489 -0.332 -0.352 -0.705 -0.340 -0.680  0.239          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.483                                            
    ## tmpnttmpnt2:stMn  0.239            0.238                           
    ## tmpnttmpnt3:stMn  0.470            0.242            0.481          
    ## tmpnttmpnt4:stMn  0.236            0.489            0.467          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.476          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

#### Symbiont Biomass:

`sAFDW_model<-lmer(Sym_AFDW.mg.cm2~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0))`

``` r
sAFDW_model<-lmer(Sym_AFDW.mg.cm2~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0))
qqPlot(residuals(sAFDW_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

    ## 195 186 
    ## 191 182

Residuals are not normally distributed. Attempt with log transformation.

`sAFDW_model<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0))`

``` r
sAFDW_model<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0))
qqPlot(residuals(sAFDW_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

    ## 195 186 
    ## 191 182

Residuals are improved.

Generate a Type III Anova of model.

``` r
anova(sAFDW_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
    ## timepoint               1.341  0.4470     3 328.28   7.2618 9.917e-05 ***
    ## site                    0.768  0.3840     2 156.10   6.2393  0.002473 ** 
    ## species                51.169 25.5847     2 149.80 415.6608 < 2.2e-16 ***
    ## timepoint:site          0.867  0.1444     6 326.46   2.3464  0.031156 *  
    ## timepoint:species       0.327  0.0545     6 322.13   0.8850  0.506029    
    ## site:species            2.754  0.6885     4 147.56  11.1853 5.937e-08 ***
    ## timepoint:site:species  0.345  0.0288    12 320.70   0.4674  0.932824    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(sAFDW_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Sym_AFDW.mg.cm2) ~ timepoint * site * species + (1 |  
    ##     colony_id_corr)
    ##    Data: master
    ##  Subset: c(Sym_AFDW.mg.cm2 > 0)
    ## 
    ## REML criterion at convergence: 176.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7950 -0.3592  0.0330  0.4235  8.9453 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.01411  0.1188  
    ##  Residual                   0.06155  0.2481  
    ## Number of obs: 437, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         0.508196   0.076128
    ## timepointtimepoint2                                -0.103546   0.111225
    ## timepointtimepoint3                                 0.071014   0.198553
    ## timepointtimepoint4                                 0.011549   0.126165
    ## siteMahana                                          0.277962   0.107422
    ## siteManava                                          0.092552   0.104001
    ## speciesPocillopora                                 -0.005207   0.105699
    ## speciesPorites                                      1.158311   0.104112
    ## timepointtimepoint2:siteMahana                     -0.097026   0.168207
    ## timepointtimepoint3:siteMahana                     -0.128169   0.226859
    ## timepointtimepoint4:siteMahana                     -0.279131   0.163588
    ## timepointtimepoint2:siteManava                      0.057460   0.154611
    ## timepointtimepoint3:siteManava                     -0.056093   0.225743
    ## timepointtimepoint4:siteManava                      0.007596   0.161154
    ## timepointtimepoint2:speciesPocillopora             -0.007416   0.148682
    ## timepointtimepoint3:speciesPocillopora              0.201499   0.220684
    ## timepointtimepoint4:speciesPocillopora              0.110361   0.156389
    ## timepointtimepoint2:speciesPorites                 -0.044533   0.144606
    ## timepointtimepoint3:speciesPorites                 -0.045470   0.218243
    ## timepointtimepoint4:speciesPorites                 -0.046317   0.156389
    ## siteMahana:speciesPocillopora                      -0.116174   0.150664
    ## siteManava:speciesPocillopora                      -0.071275   0.145729
    ## siteMahana:speciesPorites                          -0.586657   0.148239
    ## siteManava:speciesPorites                           0.042992   0.144581
    ## timepointtimepoint2:siteMahana:speciesPocillopora   0.051476   0.216738
    ## timepointtimepoint3:siteMahana:speciesPocillopora  -0.126352   0.266832
    ## timepointtimepoint4:siteMahana:speciesPocillopora   0.061713   0.215736
    ## timepointtimepoint2:siteManava:speciesPocillopora   0.069495   0.205378
    ## timepointtimepoint3:siteManava:speciesPocillopora  -0.074476   0.263864
    ## timepointtimepoint4:siteManava:speciesPocillopora   0.003541   0.208442
    ## timepointtimepoint2:siteMahana:speciesPorites       0.154809   0.213045
    ## timepointtimepoint3:siteMahana:speciesPorites       0.177209   0.263274
    ## timepointtimepoint4:siteMahana:speciesPorites       0.056302   0.211995
    ## timepointtimepoint2:siteManava:speciesPorites      -0.068288   0.201622
    ## timepointtimepoint3:siteManava:speciesPorites       0.053363   0.260968
    ## timepointtimepoint4:siteManava:speciesPorites       0.018564   0.210774
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       380.337742   6.676 8.74e-11
    ## timepointtimepoint2                               355.911394  -0.931   0.3525
    ## timepointtimepoint3                               367.996814   0.358   0.7208
    ## timepointtimepoint4                               344.254753   0.092   0.9271
    ## siteMahana                                        385.345814   2.588   0.0100
    ## siteManava                                        380.214203   0.890   0.3741
    ## speciesPocillopora                                380.436435  -0.049   0.9607
    ## speciesPorites                                    377.766138  11.126  < 2e-16
    ## timepointtimepoint2:siteMahana                    349.835493  -0.577   0.5644
    ## timepointtimepoint3:siteMahana                    370.815411  -0.565   0.5724
    ## timepointtimepoint4:siteMahana                    358.640766  -1.706   0.0888
    ## timepointtimepoint2:siteManava                    347.050939   0.372   0.7104
    ## timepointtimepoint3:siteManava                    361.622437  -0.248   0.8039
    ## timepointtimepoint4:siteManava                    336.234620   0.047   0.9624
    ## timepointtimepoint2:speciesPocillopora            336.324455  -0.050   0.9602
    ## timepointtimepoint3:speciesPocillopora            357.458231   0.913   0.3618
    ## timepointtimepoint4:speciesPocillopora            328.073448   0.706   0.4809
    ## timepointtimepoint2:speciesPorites                332.722285  -0.308   0.7583
    ## timepointtimepoint3:speciesPorites                356.706242  -0.208   0.8351
    ## timepointtimepoint4:speciesPorites                328.064010  -0.296   0.7673
    ## siteMahana:speciesPocillopora                     383.551317  -0.771   0.4411
    ## siteManava:speciesPocillopora                     379.029617  -0.489   0.6251
    ## siteMahana:speciesPorites                         381.281188  -3.958 9.03e-05
    ## siteManava:speciesPorites                         377.601280   0.297   0.7664
    ## timepointtimepoint2:siteMahana:speciesPocillopora 333.373268   0.238   0.8124
    ## timepointtimepoint3:siteMahana:speciesPocillopora 354.326743  -0.474   0.6361
    ## timepointtimepoint4:siteMahana:speciesPocillopora 338.420559   0.286   0.7750
    ## timepointtimepoint2:siteManava:speciesPocillopora 328.687706   0.338   0.7353
    ## timepointtimepoint3:siteManava:speciesPocillopora 347.797321  -0.282   0.7779
    ## timepointtimepoint4:siteManava:speciesPocillopora 321.718605   0.017   0.9865
    ## timepointtimepoint2:siteMahana:speciesPorites     331.010297   0.727   0.4680
    ## timepointtimepoint3:siteMahana:speciesPorites     354.380243   0.673   0.5013
    ## timepointtimepoint4:siteMahana:speciesPorites     336.120078   0.266   0.7907
    ## timepointtimepoint2:siteManava:speciesPorites     326.162719  -0.339   0.7351
    ## timepointtimepoint3:siteManava:speciesPorites     346.892735   0.204   0.8381
    ## timepointtimepoint4:siteManava:speciesPorites     322.951853   0.088   0.9299
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                                  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## siteMahana                                        *  
    ## siteManava                                           
    ## speciesPocillopora                                   
    ## speciesPorites                                    ***
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                    .  
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## siteMahana:speciesPocillopora                        
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                         ***
    ## siteManava:speciesPorites                            
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora    
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora    
    ## timepointtimepoint3:siteManava:speciesPocillopora    
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites        
    ## timepointtimepoint3:siteMahana:speciesPorites        
    ## timepointtimepoint4:siteMahana:speciesPorites        
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites        
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of site and time within
each species.

*Acropora*

`sAFDW_model_acr<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0 & species=="Acropora"))`

``` r
sAFDW_model_acr<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0 & species=="Acropora"))
qqPlot(residuals(sAFDW_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

    ## 215   7 
    ##  62   7

Generate a Type III Anova of model.

``` r
anova(sAFDW_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq  Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## timepoint      0.26993 0.089976     3 90.217  4.2761 0.007181 **
    ## site           0.23316 0.116578     2 51.752  5.5403 0.006603 **
    ## timepoint:site 0.33092 0.055153     6 88.910  2.6211 0.021874 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(sAFDW_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Sym_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Sym_AFDW.mg.cm2 > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: -74.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.77336 -0.52304 -0.02092  0.64541  2.54673 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001139 0.03375 
    ##  Residual                   0.021042 0.14506 
    ## Number of obs: 115, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     5.069e-01  4.130e-02  1.028e+02  12.274
    ## timepointtimepoint2            -9.849e-02  6.349e-02  8.795e+01  -1.551
    ## timepointtimepoint3             9.516e-02  1.121e-01  9.753e+01   0.849
    ## timepointtimepoint4             1.447e-02  7.225e-02  8.708e+01   0.200
    ## siteMahana                      2.775e-01  5.840e-02  1.029e+02   4.752
    ## siteManava                      9.585e-02  5.643e-02  1.028e+02   1.699
    ## timepointtimepoint2:siteMahana -1.011e-01  9.618e-02  8.747e+01  -1.051
    ## timepointtimepoint3:siteMahana -1.505e-01  1.281e-01  9.703e+01  -1.175
    ## timepointtimepoint4:siteMahana -2.805e-01  9.318e-02  8.994e+01  -3.011
    ## timepointtimepoint2:siteManava  5.476e-02  8.850e-02  8.675e+01   0.619
    ## timepointtimepoint3:siteManava -8.042e-02  1.279e-01  9.517e+01  -0.629
    ## timepointtimepoint4:siteManava -5.139e-04  9.259e-02  8.484e+01  -0.006
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## timepointtimepoint2             0.12443    
    ## timepointtimepoint3             0.39803    
    ## timepointtimepoint4             0.84176    
    ## siteMahana                     6.55e-06 ***
    ## siteManava                      0.09240 .  
    ## timepointtimepoint2:siteMahana  0.29617    
    ## timepointtimepoint3:siteMahana  0.24269    
    ## timepointtimepoint4:siteMahana  0.00338 ** 
    ## timepointtimepoint2:siteManava  0.53771    
    ## timepointtimepoint3:siteManava  0.53115    
    ## timepointtimepoint4:siteManava  0.99558    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.625                                                    
    ## tmpnttmpnt3      -0.350  0.235                                             
    ## tmpnttmpnt4      -0.543  0.356  0.201                                      
    ## siteMahana       -0.707  0.442  0.247  0.384                               
    ## siteManava       -0.732  0.457  0.256  0.397  0.518                        
    ## tmpnttmpnt2:stMh  0.412 -0.660 -0.155 -0.235 -0.580 -0.302                 
    ## tmpnttmpnt3:stMh  0.306 -0.206 -0.875 -0.176 -0.439 -0.224  0.271          
    ## tmpnttmpnt4:stMh  0.421 -0.276 -0.156 -0.775 -0.603 -0.308  0.365          
    ## tmpnttmpnt2:stMn  0.448 -0.717 -0.169 -0.256 -0.317 -0.610  0.474          
    ## tmpnttmpnt3:stMn  0.306 -0.206 -0.876 -0.176 -0.217 -0.420  0.136          
    ## tmpnttmpnt4:stMn  0.423 -0.278 -0.157 -0.780 -0.299 -0.580  0.184          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.279                                            
    ## tmpnttmpnt2:stMn  0.148            0.198                           
    ## tmpnttmpnt3:stMn  0.767            0.136            0.276          
    ## tmpnttmpnt4:stMn  0.137            0.605            0.372          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.257

*Pocillopora*

`sAFDW_model_poc<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0 & species=="Pocillopora"))`

``` r
sAFDW_model_poc<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0 & species=="Pocillopora"))
qqPlot(residuals(sAFDW_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

    ## 254  76 
    ##  94  35

Generate a Type III Anova of model.

``` r
anova(sAFDW_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      1.00837 0.33612     3 107.113 16.0125 1.138e-08 ***
    ## site           0.01185 0.00593     2  40.756  0.2823  0.755527    
    ## timepoint:site 0.44062 0.07344     6 107.000  3.4984  0.003374 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(sAFDW_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Sym_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Sym_AFDW.mg.cm2 > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -78.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8207 -0.4873  0.0510  0.6122  1.9208 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.009551 0.09773 
    ##  Residual                   0.020991 0.14488 
    ## Number of obs: 157, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.50644    0.04643 119.64360  10.908  < 2e-16
    ## timepointtimepoint2             -0.11496    0.05790 107.76710  -1.985  0.04964
    ## timepointtimepoint3              0.26647    0.05646 106.57128   4.720 7.21e-06
    ## timepointtimepoint4              0.11846    0.05402 103.81610   2.193  0.03055
    ## siteMahana                       0.15884    0.06683 121.02020   2.377  0.01903
    ## siteManava                       0.01782    0.06474 117.43874   0.275  0.78356
    ## timepointtimepoint2:siteMahana  -0.04207    0.08016 107.35238  -0.525  0.60081
    ## timepointtimepoint3:siteMahana  -0.24939    0.08230 106.11865  -3.030  0.00307
    ## timepointtimepoint4:siteMahana  -0.21153    0.08252 107.66161  -2.563  0.01174
    ## timepointtimepoint2:siteManava   0.13324    0.07919 105.92681   1.683  0.09541
    ## timepointtimepoint3:siteManava  -0.12210    0.08009 106.74943  -1.525  0.13031
    ## timepointtimepoint4:siteManava   0.01871    0.07735 104.87354   0.242  0.80934
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            *  
    ## timepointtimepoint3            ***
    ## timepointtimepoint4            *  
    ## siteMahana                     *  
    ## siteManava                        
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana ** 
    ## timepointtimepoint4:siteMahana *  
    ## timepointtimepoint2:siteManava .  
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.567                                                    
    ## tmpnttmpnt3      -0.581  0.474                                             
    ## tmpnttmpnt4      -0.606  0.488  0.499                                      
    ## siteMahana       -0.695  0.394  0.404  0.421                               
    ## siteManava       -0.717  0.407  0.417  0.434  0.498                        
    ## tmpnttmpnt2:stMh  0.410 -0.722 -0.342 -0.352 -0.597 -0.294                 
    ## tmpnttmpnt3:stMh  0.398 -0.325 -0.686 -0.342 -0.570 -0.286  0.479          
    ## tmpnttmpnt4:stMh  0.396 -0.319 -0.327 -0.655 -0.576 -0.284  0.481          
    ## tmpnttmpnt2:stMn  0.415 -0.731 -0.346 -0.356 -0.288 -0.570  0.528          
    ## tmpnttmpnt3:stMn  0.410 -0.334 -0.705 -0.352 -0.284 -0.564  0.241          
    ## tmpnttmpnt4:stMn  0.423 -0.340 -0.349 -0.698 -0.294 -0.583  0.246          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.465                                            
    ## tmpnttmpnt2:stMn  0.238            0.233                           
    ## tmpnttmpnt3:stMn  0.484            0.230            0.464          
    ## tmpnttmpnt4:stMn  0.239            0.457            0.477          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.480

*Porites*

`sAFDW_model_poc<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0 & species=="Porites"))`

``` r
sAFDW_model_por<-lmer(log(1+Sym_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Sym_AFDW.mg.cm2>0 & species=="Porites"))
qqPlot(residuals(sAFDW_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

    ## 195 186 
    ##  55  46

Generate a Type III Anova of model.

``` r
anova(sAFDW_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      0.8232 0.27441     3 115.439  2.1932   0.09261 .  
    ## site           3.3242 1.66210     2  41.784 13.2846 3.423e-05 ***
    ## timepoint:site 0.4873 0.08122     6 115.386  0.6492   0.69067    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(sAFDW_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Sym_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Sym_AFDW.mg.cm2 > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: 174.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9510 -0.5072  0.0388  0.4482  6.2336 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.03089  0.1757  
    ##  Residual                   0.12512  0.3537  
    ## Number of obs: 165, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      1.666507   0.101981 138.208767  16.341
    ## timepointtimepoint2             -0.147636   0.131766 113.683882  -1.120
    ## timepointtimepoint3              0.025544   0.129159 112.379657   0.198
    ## timepointtimepoint4             -0.034725   0.131766 113.683882  -0.264
    ## siteMahana                      -0.308492   0.146677 139.032293  -2.103
    ## siteManava                       0.135544   0.144223 138.208767   0.940
    ## timepointtimepoint2:siteMahana   0.057137   0.186435 114.298866   0.306
    ## timepointtimepoint3:siteMahana   0.048683   0.190503 114.494553   0.256
    ## timepointtimepoint4:siteMahana  -0.223228   0.192280 115.068518  -1.161
    ## timepointtimepoint2:siteManava  -0.011271   0.184511 113.044640  -0.061
    ## timepointtimepoint3:siteManava  -0.002172   0.186707 114.276470  -0.012
    ## timepointtimepoint4:siteManava   0.027249   0.193747 116.514170   0.141
    ##                                Pr(>|t|)    
    ## (Intercept)                      <2e-16 ***
    ## timepointtimepoint2              0.2649    
    ## timepointtimepoint3              0.8436    
    ## timepointtimepoint4              0.7926    
    ## siteMahana                       0.0372 *  
    ## siteManava                       0.3489    
    ## timepointtimepoint2:siteMahana   0.7598    
    ## timepointtimepoint3:siteMahana   0.7988    
    ## timepointtimepoint4:siteMahana   0.2481    
    ## timepointtimepoint2:siteManava   0.9514    
    ## timepointtimepoint3:siteManava   0.9907    
    ## timepointtimepoint4:siteManava   0.8884    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.621                                                    
    ## tmpnttmpnt3      -0.633  0.490                                             
    ## tmpnttmpnt4      -0.621  0.480  0.490                                      
    ## siteMahana       -0.695  0.432  0.440  0.432                               
    ## siteManava       -0.707  0.439  0.448  0.439  0.492                        
    ## tmpnttmpnt2:stMh  0.439 -0.707 -0.346 -0.339 -0.636 -0.310                 
    ## tmpnttmpnt3:stMh  0.429 -0.332 -0.678 -0.332 -0.618 -0.304  0.486          
    ## tmpnttmpnt4:stMh  0.425 -0.329 -0.336 -0.685 -0.613 -0.301  0.482          
    ## tmpnttmpnt2:stMn  0.443 -0.714 -0.350 -0.343 -0.308 -0.627  0.505          
    ## tmpnttmpnt3:stMn  0.438 -0.339 -0.692 -0.339 -0.305 -0.620  0.240          
    ## tmpnttmpnt4:stMn  0.422 -0.326 -0.333 -0.680 -0.294 -0.597  0.231          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.478                                            
    ## tmpnttmpnt2:stMn  0.237            0.235                           
    ## tmpnttmpnt3:stMn  0.469            0.232            0.484          
    ## tmpnttmpnt4:stMn  0.226            0.466            0.466          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.467

#### S:H Biomass:

`shAFDW_model<-lmer(Ratio_AFDW.mg.cm2~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0))`

``` r
shAFDW_model<-lmer(Ratio_AFDW.mg.cm2~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0))
qqPlot(residuals(shAFDW_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

    ## 195 254 
    ## 186 242

Residuals are not normally distributed. Attempt with log transformation.

`shAFDW_model<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0))`

``` r
shAFDW_model<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0))
qqPlot(residuals(shAFDW_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

    ## 195 254 
    ## 186 242

Residuals are improved.

Generate a Type III Anova of model.

``` r
anova(shAFDW_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                          Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint              0.259969 0.086656     3 305.25 33.0855 < 2.2e-16 ***
    ## site                   0.004618 0.002309     2 143.34  0.8817 0.4163259    
    ## species                0.209225 0.104613     2 138.79 39.9412 1.990e-14 ***
    ## timepoint:site         0.059435 0.009906     6 303.81  3.7821 0.0011969 ** 
    ## timepoint:species      0.067047 0.011175     6 299.84  4.2664 0.0003836 ***
    ## site:species           0.094543 0.023636     4 137.04  9.0241 1.685e-06 ***
    ## timepoint:site:species 0.066830 0.005569    12 298.66  2.1263 0.0153546 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(shAFDW_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Ratio_AFDW.mg.cm2) ~ timepoint * site * species + (1 |  
    ##     colony_id_corr)
    ##    Data: master
    ##  Subset: c(Ratio_AFDW.mg.cm2 > 0)
    ## 
    ## REML criterion at convergence: -1032.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.5320 -0.4517  0.0540  0.4963  4.6452 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0009742 0.03121 
    ##  Residual                   0.0026192 0.05118 
    ## Number of obs: 427, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         0.358846   0.016552
    ## timepointtimepoint2                                -0.024333   0.023289
    ## timepointtimepoint3                                 0.043498   0.041647
    ## timepointtimepoint4                                -0.009035   0.026325
    ## siteMahana                                          0.063178   0.023304
    ## siteManava                                         -0.019870   0.022615
    ## speciesPocillopora                                 -0.112053   0.022982
    ## speciesPorites                                     -0.012532   0.022982
    ## timepointtimepoint2:siteMahana                     -0.005771   0.035156
    ## timepointtimepoint3:siteMahana                     -0.013396   0.047640
    ## timepointtimepoint4:siteMahana                     -0.089953   0.034268
    ## timepointtimepoint2:siteManava                      0.001002   0.032973
    ## timepointtimepoint3:siteManava                      0.039851   0.047276
    ## timepointtimepoint4:siteManava                      0.014714   0.033569
    ## timepointtimepoint2:speciesPocillopora              0.007525   0.030975
    ## timepointtimepoint3:speciesPocillopora              0.069156   0.046166
    ## timepointtimepoint4:speciesPocillopora              0.003759   0.032510
    ## timepointtimepoint2:speciesPorites                  0.029697   0.030967
    ## timepointtimepoint3:speciesPorites                 -0.015101   0.045808
    ## timepointtimepoint4:speciesPorites                 -0.017788   0.033004
    ## siteMahana:speciesPocillopora                       0.039187   0.032714
    ## siteManava:speciesPocillopora                       0.012335   0.031705
    ## siteMahana:speciesPorites                          -0.086571   0.032706
    ## siteManava:speciesPorites                           0.054429   0.031705
    ## timepointtimepoint2:siteMahana:speciesPocillopora  -0.059564   0.045118
    ## timepointtimepoint3:siteMahana:speciesPocillopora  -0.106900   0.055795
    ## timepointtimepoint4:siteMahana:speciesPocillopora   0.039602   0.044961
    ## timepointtimepoint2:siteManava:speciesPocillopora   0.041515   0.043419
    ## timepointtimepoint3:siteManava:speciesPocillopora  -0.016936   0.055289
    ## timepointtimepoint4:siteManava:speciesPocillopora  -0.003838   0.043274
    ## timepointtimepoint2:siteMahana:speciesPorites      -0.028018   0.045106
    ## timepointtimepoint3:siteMahana:speciesPorites      -0.007972   0.055567
    ## timepointtimepoint4:siteMahana:speciesPorites       0.060417   0.044730
    ## timepointtimepoint2:siteManava:speciesPorites      -0.018772   0.043217
    ## timepointtimepoint3:siteManava:speciesPorites      -0.060281   0.054596
    ## timepointtimepoint4:siteManava:speciesPorites      -0.005783   0.044138
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       348.360186  21.680  < 2e-16
    ## timepointtimepoint2                               338.785062  -1.045  0.29685
    ## timepointtimepoint3                               339.987334   1.044  0.29702
    ## timepointtimepoint4                               321.093214  -0.343  0.73167
    ## siteMahana                                        358.159946   2.711  0.00703
    ## siteManava                                        347.325097  -0.879  0.38021
    ## speciesPocillopora                                347.676045  -4.876 1.65e-06
    ## speciesPorites                                    347.678725  -0.545  0.58591
    ## timepointtimepoint2:siteMahana                    329.777051  -0.164  0.86971
    ## timepointtimepoint3:siteMahana                    346.022186  -0.281  0.77874
    ## timepointtimepoint4:siteMahana                    339.359514  -2.625  0.00906
    ## timepointtimepoint2:siteManava                    328.141556   0.030  0.97579
    ## timepointtimepoint3:siteManava                    334.167207   0.843  0.39986
    ## timepointtimepoint4:siteManava                    313.174958   0.438  0.66145
    ## timepointtimepoint2:speciesPocillopora            315.902893   0.243  0.80821
    ## timepointtimepoint3:speciesPocillopora            329.601314   1.498  0.13510
    ## timepointtimepoint4:speciesPocillopora            304.763611   0.116  0.90802
    ## timepointtimepoint2:speciesPorites                315.353823   0.959  0.33830
    ## timepointtimepoint3:speciesPorites                329.023578  -0.330  0.74187
    ## timepointtimepoint4:speciesPorites                306.190948  -0.539  0.59031
    ## siteMahana:speciesPocillopora                     354.030213   1.198  0.23177
    ## siteManava:speciesPocillopora                     344.741003   0.389  0.69748
    ## siteMahana:speciesPorites                         354.363214  -2.647  0.00848
    ## siteManava:speciesPorites                         344.742436   1.717  0.08693
    ## timepointtimepoint2:siteMahana:speciesPocillopora 311.780173  -1.320  0.18774
    ## timepointtimepoint3:siteMahana:speciesPocillopora 329.119090  -1.916  0.05624
    ## timepointtimepoint4:siteMahana:speciesPocillopora 316.925697   0.881  0.37909
    ## timepointtimepoint2:siteManava:speciesPocillopora 308.366940   0.956  0.33975
    ## timepointtimepoint3:siteManava:speciesPocillopora 320.965737  -0.306  0.75956
    ## timepointtimepoint4:siteManava:speciesPocillopora 298.496322  -0.089  0.92939
    ## timepointtimepoint2:siteMahana:speciesPorites     311.251752  -0.621  0.53495
    ## timepointtimepoint3:siteMahana:speciesPorites     329.853330  -0.143  0.88601
    ## timepointtimepoint4:siteMahana:speciesPorites     316.308354   1.351  0.17776
    ## timepointtimepoint2:siteManava:speciesPorites     307.747694  -0.434  0.66433
    ## timepointtimepoint3:siteManava:speciesPorites     320.126170  -1.104  0.27037
    ## timepointtimepoint4:siteManava:speciesPorites     300.361198  -0.131  0.89584
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                                  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## siteMahana                                        ** 
    ## siteManava                                           
    ## speciesPocillopora                                ***
    ## speciesPorites                                       
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                    ** 
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## siteMahana:speciesPocillopora                        
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                         ** 
    ## siteManava:speciesPorites                         .  
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora .  
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora    
    ## timepointtimepoint3:siteManava:speciesPocillopora    
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites        
    ## timepointtimepoint3:siteMahana:speciesPorites        
    ## timepointtimepoint4:siteMahana:speciesPorites        
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites        
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of time and site within
species.

*Acropora*

`shAFDW_model_acr<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0 & species=="Acropora"))`

``` r
shAFDW_model_acr<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0 & species=="Acropora"))
qqPlot(residuals(shAFDW_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

    ## [1] 15  7

Generate a Type III Anova of model.

``` r
anova(shAFDW_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq   Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      0.076426 0.0254753     3 85.752  8.4170 5.764e-05 ***
    ## site           0.024952 0.0124762     2 46.796  4.1221   0.02244 *  
    ## timepoint:site 0.042436 0.0070727     6 84.321  2.3368   0.03892 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(shAFDW_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Ratio_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Ratio_AFDW.mg.cm2 > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: -269.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4822 -0.4154  0.0847  0.5533  2.2850 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0002339 0.01529 
    ##  Residual                   0.0030267 0.05502 
    ## Number of obs: 114, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     3.577e-01  1.583e-02  1.015e+02  22.594
    ## timepointtimepoint2            -2.018e-02  2.417e-02  8.359e+01  -0.835
    ## timepointtimepoint3             5.272e-02  4.278e-02  9.387e+01   1.232
    ## timepointtimepoint4            -1.002e-02  2.749e-02  8.212e+01  -0.364
    ## siteMahana                      6.467e-02  2.238e-02  1.017e+02   2.889
    ## siteManava                     -1.858e-02  2.163e-02  1.015e+02  -0.859
    ## timepointtimepoint2:siteMahana -1.675e-02  3.660e-02  8.279e+01  -0.458
    ## timepointtimepoint3:siteMahana -2.734e-02  4.886e-02  9.349e+01  -0.560
    ## timepointtimepoint4:siteMahana -9.256e-02  3.548e-02  8.569e+01  -2.608
    ## timepointtimepoint2:siteManava  4.238e-04  3.434e-02  8.300e+01   0.012
    ## timepointtimepoint3:siteManava  2.965e-02  4.879e-02  9.106e+01   0.608
    ## timepointtimepoint4:siteManava  1.321e-02  3.522e-02  7.960e+01   0.375
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## timepointtimepoint2             0.40617    
    ## timepointtimepoint3             0.22089    
    ## timepointtimepoint4             0.71653    
    ## siteMahana                      0.00472 ** 
    ## siteManava                      0.39249    
    ## timepointtimepoint2:siteMahana  0.64838    
    ## timepointtimepoint3:siteMahana  0.57710    
    ## timepointtimepoint4:siteMahana  0.01073 *  
    ## timepointtimepoint2:siteManava  0.99018    
    ## timepointtimepoint3:siteManava  0.54490    
    ## timepointtimepoint4:siteManava  0.70863    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.619                                                    
    ## tmpnttmpnt3      -0.344  0.236                                             
    ## tmpnttmpnt4      -0.535  0.355  0.199                                      
    ## siteMahana       -0.707  0.438  0.243  0.378                               
    ## siteManava       -0.732  0.453  0.252  0.392  0.518                        
    ## tmpnttmpnt2:stMh  0.408 -0.660 -0.156 -0.234 -0.573 -0.299                 
    ## tmpnttmpnt3:stMh  0.301 -0.206 -0.876 -0.174 -0.434 -0.220  0.272          
    ## tmpnttmpnt4:stMh  0.414 -0.275 -0.154 -0.775 -0.597 -0.303  0.364          
    ## tmpnttmpnt2:stMn  0.435 -0.704 -0.166 -0.250 -0.308 -0.592  0.465          
    ## tmpnttmpnt3:stMn  0.302 -0.207 -0.877 -0.175 -0.213 -0.413  0.137          
    ## tmpnttmpnt4:stMn  0.418 -0.277 -0.156 -0.781 -0.295 -0.572  0.183          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.279                                            
    ## tmpnttmpnt2:stMn  0.145            0.194                           
    ## tmpnttmpnt3:stMn  0.768            0.135            0.272          
    ## tmpnttmpnt4:stMn  0.136            0.605            0.365          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.256

*Pocillopora*

`shAFDW_model_poc<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0 & species=="Pocillopora"))`

``` r
shAFDW_model_poc<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0 & species=="Pocillopora"))
qqPlot(residuals(shAFDW_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->

    ## 254 240 
    ##  93  80

Generate a Type III Anova of model.

``` r
anova(shAFDW_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      0.245129 0.081710     3 106.919 30.4598 2.598e-14 ***
    ## site           0.024031 0.012016     2  41.175  4.4792 0.0173742 *  
    ## timepoint:site 0.082312 0.013719     6 106.787  5.1141 0.0001185 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(shAFDW_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Ratio_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Ratio_AFDW.mg.cm2 > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -379.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.5269 -0.4834  0.0612  0.5382  1.9987 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.000856 0.02926 
    ##  Residual                   0.002683 0.05179 
    ## Number of obs: 155, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      0.246566   0.015836 127.782531  15.570
    ## timepointtimepoint2             -0.016471   0.020646 107.291768  -0.798
    ## timepointtimepoint3              0.113008   0.020144 105.897125   5.610
    ## timepointtimepoint4             -0.005048   0.019301 102.826115  -0.262
    ## siteMahana                       0.102773   0.022807 128.593719   4.506
    ## siteManava                      -0.007308   0.022061 126.129599  -0.331
    ## timepointtimepoint2:siteMahana  -0.065854   0.028591 106.590796  -2.303
    ## timepointtimepoint3:siteMahana  -0.120532   0.029371 105.378199  -4.104
    ## timepointtimepoint4:siteMahana  -0.050724   0.029426 107.074599  -1.724
    ## timepointtimepoint2:siteManava   0.042019   0.028566 105.759521   1.471
    ## timepointtimepoint3:siteManava   0.022392   0.028984 106.787117   0.773
    ## timepointtimepoint4:siteManava   0.010300   0.027625 103.945604   0.373
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## timepointtimepoint2              0.4268    
    ## timepointtimepoint3            1.63e-07 ***
    ## timepointtimepoint4              0.7942    
    ## siteMahana                     1.47e-05 ***
    ## siteManava                       0.7410    
    ## timepointtimepoint2:siteMahana   0.0232 *  
    ## timepointtimepoint3:siteMahana 8.03e-05 ***
    ## timepointtimepoint4:siteMahana   0.0876 .  
    ## timepointtimepoint2:siteManava   0.1443    
    ## timepointtimepoint3:siteManava   0.4415    
    ## timepointtimepoint4:siteManava   0.7100    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.594                                                    
    ## tmpnttmpnt3      -0.609  0.473                                             
    ## tmpnttmpnt4      -0.634  0.488  0.499                                      
    ## siteMahana       -0.694  0.413  0.423  0.440                               
    ## siteManava       -0.718  0.427  0.437  0.455  0.498                        
    ## tmpnttmpnt2:stMh  0.429 -0.722 -0.342 -0.352 -0.624 -0.308                 
    ## tmpnttmpnt3:stMh  0.417 -0.325 -0.686 -0.342 -0.598 -0.300  0.480          
    ## tmpnttmpnt4:stMh  0.416 -0.320 -0.327 -0.656 -0.602 -0.298  0.481          
    ## tmpnttmpnt2:stMn  0.430 -0.723 -0.342 -0.352 -0.298 -0.592  0.522          
    ## tmpnttmpnt3:stMn  0.423 -0.329 -0.695 -0.347 -0.294 -0.583  0.237          
    ## tmpnttmpnt4:stMn  0.443 -0.341 -0.349 -0.699 -0.307 -0.611  0.246          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.466                                            
    ## tmpnttmpnt2:stMn  0.235            0.231                           
    ## tmpnttmpnt3:stMn  0.477            0.228            0.451          
    ## tmpnttmpnt4:stMn  0.239            0.458            0.472          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.473

*Porites*

`shAFDW_model_por<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0 & species=="Porites"))`

``` r
shAFDW_model_por<-lmer(log(1+Ratio_AFDW.mg.cm2)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Ratio_AFDW.mg.cm2>0 & species=="Porites"))
qqPlot(residuals(shAFDW_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

    ## 195  88 
    ##  52   8

Generate a Type III Anova of model.

``` r
anova(shAFDW_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq   Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      0.045464 0.0151548     3 107.275  6.4703 0.0004562 ***
    ## site           0.039853 0.0199264     2  41.337  8.5076 0.0008046 ***
    ## timepoint:site 0.009323 0.0015538     6 107.210  0.6634 0.6793134    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(shAFDW_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Ratio_AFDW.mg.cm2) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Ratio_AFDW.mg.cm2 > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: -388.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0599 -0.4008  0.0754  0.4509  4.6024 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001591 0.03989 
    ##  Residual                   0.002342 0.04840 
    ## Number of obs: 158, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      0.345475   0.016609 107.320284  20.800
    ## timepointtimepoint2              0.006628   0.019374 107.516625   0.342
    ## timepointtimepoint3              0.029236   0.018054 104.782839   1.619
    ## timepointtimepoint4             -0.026018   0.018876 106.535200  -1.378
    ## siteMahana                      -0.023179   0.023865 109.476964  -0.971
    ## siteManava                       0.035398   0.023197 104.820443   1.526
    ## timepointtimepoint2:siteMahana  -0.034429   0.026815 107.375358  -1.284
    ## timepointtimepoint3:siteMahana  -0.019915   0.027135 107.129263  -0.734
    ## timepointtimepoint4:siteMahana  -0.028956   0.027279 107.305532  -1.061
    ## timepointtimepoint2:siteManava  -0.018718   0.026482 106.248695  -0.707
    ## timepointtimepoint3:siteManava  -0.019918   0.025869 105.718393  -0.770
    ## timepointtimepoint4:siteManava   0.009614   0.027198 107.443782   0.353
    ##                                Pr(>|t|)    
    ## (Intercept)                      <2e-16 ***
    ## timepointtimepoint2               0.733    
    ## timepointtimepoint3               0.108    
    ## timepointtimepoint4               0.171    
    ## siteMahana                        0.334    
    ## siteManava                        0.130    
    ## timepointtimepoint2:siteMahana    0.202    
    ## timepointtimepoint3:siteMahana    0.465    
    ## timepointtimepoint4:siteMahana    0.291    
    ## timepointtimepoint2:siteManava    0.481    
    ## timepointtimepoint3:siteManava    0.443    
    ## timepointtimepoint4:siteManava    0.724    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.530                                                    
    ## tmpnttmpnt3      -0.566  0.488                                             
    ## tmpnttmpnt4      -0.543  0.463  0.500                                      
    ## siteMahana       -0.696  0.369  0.394  0.378                               
    ## siteManava       -0.716  0.380  0.405  0.389  0.498                        
    ## tmpnttmpnt2:stMh  0.383 -0.722 -0.352 -0.335 -0.560 -0.274                 
    ## tmpnttmpnt3:stMh  0.377 -0.324 -0.665 -0.332 -0.547 -0.270  0.488          
    ## tmpnttmpnt4:stMh  0.376 -0.320 -0.346 -0.692 -0.544 -0.269  0.483          
    ## tmpnttmpnt2:stMn  0.388 -0.732 -0.357 -0.339 -0.270 -0.532  0.529          
    ## tmpnttmpnt3:stMn  0.395 -0.340 -0.698 -0.349 -0.275 -0.543  0.246          
    ## tmpnttmpnt4:stMn  0.377 -0.321 -0.347 -0.694 -0.262 -0.517  0.232          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.491                                            
    ## tmpnttmpnt2:stMn  0.237            0.234                           
    ## tmpnttmpnt3:stMn  0.464            0.241            0.476          
    ## tmpnttmpnt4:stMn  0.231            0.480            0.450          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.474

## Chlorophyll

### Plot: Surface area

View by site and timepoint for each species in Chl a.

``` r
chlaplot_cm2<-master %>%
  filter(!is.na(chla.ug.cm2)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chla.ug.cm2)%>%
  
  ggplot(., aes(x = site_code, y = chla.ug.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll a (", mu, "g cm"^-2, ")"))))+
    theme_classic() + 
    theme(
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

View by site and timepoint for each species in Chl c2.

``` r
chlc2plot_cm2<-master %>%
  filter(!is.na(chlc2.ug.cm2)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chlc2.ug.cm2)%>%
  
  ggplot(., aes(x = site_code, y = chlc2.ug.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll c2 (", mu, "g cm"^-2, ")"))))+
    theme_classic() + 
    theme(
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Protein

View by site and timepoint for each species in Chl a.

``` r
chlaplot_prot<-master %>%
  filter(!is.na(chla.ug.mgprot)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chla.ug.mgprot)%>%
  
  ggplot(., aes(x = site_code, y = chla.ug.mgprot, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll a (", mu, "g mg protein"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

View by site and timepoint for each species in Chl c2.

``` r
chlc2plot_prot<-master %>%
  filter(!is.na(chlc2.ug.mgprot)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chlc2.ug.mgprot)%>%
  
  ggplot(., aes(x = site_code, y = chlc2.ug.mgprot, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll c2 (", mu, "g mg protein"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Biomass

View by site and timepoint for each species in Chl a.

``` r
chlaplot_afdw<-master %>%
  filter(!is.na(chla.ug.mgAFDW)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chla.ug.mgAFDW)%>%
  
  ggplot(., aes(x = site_code, y = chla.ug.mgAFDW, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll a (", mu, "g mg afdw"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

View by site and timepoint for each species in Chl c2.

``` r
chlc2plot_afdw<-master %>%
  filter(!is.na(chlc2.ug.mgAFDW)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chlc2.ug.mgAFDW)%>%
  
  ggplot(., aes(x = site_code, y = chlc2.ug.mgAFDW, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll c2 (", mu, "g mg afdw"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Per cell

View by site and timepoint for each species in Chl a.

``` r
chlaplot_cell<-master %>%
  filter(!is.na(chla.ug.cell)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chla.ug.cell)%>%
  
  ggplot(., aes(x = site_code, y = chla.ug.cell, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll a (", mu, "g cell"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

View by site and timepoint for each species in Chl c2.

``` r
chlc2plot_cell<-master %>%
  filter(!is.na(chlc2.ug.cell)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, chlc2.ug.cell)%>%
  
  ggplot(., aes(x = site_code, y = chlc2.ug.cell, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Chlorophyll c2 (", mu, "g cell"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots for each normalization together.

``` r
chl_figure<-plot_grid(chlaplot_prot, chlaplot_afdw, chlaplot_cell, chlaplot_cm2,chlc2plot_prot, chlc2plot_afdw, chlc2plot_cell, chlc2plot_cm2,  labels = c("", "", ""), ncol=4, nrow=2, rel_widths = c(.8, .8, .8, 1, .8, .8, .8, 1), label_size = 20, label_x = 0.1)

ggsave(filename="Figures/Chlorophyll_Figure.pdf", plot=chl_figure, dpi=300, width=32, height=10, units="in")
```

### Plot: Total chl

View by site and timepoint for each species for total chl per unit AFDW.

``` r
chlplot_total<-master %>%
  filter(!is.na(Total_Chl)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, Total_Chl)%>%
  
  ggplot(., aes(x = site_code, y = Total_Chl, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Total Chlorophyll (", mu, "g mg afdw"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

View by site and timepoint for each species for total chl per cell.

``` r
chlplot_total_cell<-master %>%
  filter(!is.na(Total_Chl_cell)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, Total_Chl_cell)%>%
  
  ggplot(., aes(x = site_code, y = Total_Chl_cell, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Total Chlorophyll (", mu, "g cell"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Analysis

Build a mixed model for univariate analysis and examine data
distribution.

#### Chlorophyll a:

`chla_model<-lmer(chla.ug.mgAFDW~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
chla_model<-lmer(chla.ug.mgAFDW~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(chla_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

    ## 109 218 
    ## 104 210

Residuals are not normally distributed. Attempt with log transformation.

`chla_model<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0))`

``` r
chla_model<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0))
qqPlot(residuals(chla_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

    ## 109 218 
    ## 104 209

Generate a Type III Anova of model.

``` r
anova(chla_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint              5.8930 1.96433     3 306.02 46.1930 < 2.2e-16 ***
    ## site                   1.0917 0.54584     2 112.15 12.8360 9.555e-06 ***
    ## species                3.5914 1.79568     2 104.54 42.2270 3.614e-14 ***
    ## timepoint:site         2.4282 0.40470     6 302.31  9.5170 1.383e-09 ***
    ## timepoint:species      0.6016 0.10027     6 296.02  2.3579  0.030625 *  
    ## site:species           0.7428 0.18570     4 102.22  4.3668  0.002665 ** 
    ## timepoint:site:species 0.5409 0.04507    12 293.06  1.0599  0.394111    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chla_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## log(1 + chla.ug.mgAFDW) ~ timepoint * site * species + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chla.ug.mgAFDW > 0)
    ## 
    ## REML criterion at convergence: -20.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2989 -0.5185 -0.0788  0.4204  6.5727 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.002152 0.04639 
    ##  Residual                   0.042524 0.20621 
    ## Number of obs: 433, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         0.929535   0.061002
    ## timepointtimepoint2                                 0.265650   0.091873
    ## timepointtimepoint3                                 0.154037   0.160613
    ## timepointtimepoint4                                -0.178917   0.103913
    ## siteMahana                                          0.436467   0.084585
    ## siteManava                                          0.232480   0.081845
    ## speciesPocillopora                                 -0.229521   0.084590
    ## speciesPorites                                     -0.196218   0.083133
    ## timepointtimepoint2:siteMahana                     -0.364504   0.137759
    ## timepointtimepoint3:siteMahana                     -0.411703   0.183098
    ## timepointtimepoint4:siteMahana                     -0.361437   0.133357
    ## timepointtimepoint2:siteManava                     -0.354425   0.129375
    ## timepointtimepoint3:siteManava                     -0.014424   0.182957
    ## timepointtimepoint4:siteManava                     -0.192499   0.132552
    ## timepointtimepoint2:speciesPocillopora             -0.062919   0.121576
    ## timepointtimepoint3:speciesPocillopora              0.132155   0.179892
    ## timepointtimepoint4:speciesPocillopora             -0.012117   0.130090
    ## timepointtimepoint2:speciesPorites                 -0.163927   0.122746
    ## timepointtimepoint3:speciesPorites                 -0.028427   0.177983
    ## timepointtimepoint4:speciesPorites                  0.040943   0.130910
    ## siteMahana:speciesPocillopora                      -0.083108   0.118421
    ## siteManava:speciesPocillopora                       0.124546   0.114504
    ## siteMahana:speciesPorites                          -0.286039   0.117385
    ## siteManava:speciesPorites                           0.137240   0.113432
    ## timepointtimepoint2:siteMahana:speciesPocillopora  -0.148879   0.177323
    ## timepointtimepoint3:siteMahana:speciesPocillopora  -0.140082   0.215996
    ## timepointtimepoint4:siteMahana:speciesPocillopora   0.232090   0.175389
    ## timepointtimepoint2:siteManava:speciesPocillopora  -0.047161   0.170174
    ## timepointtimepoint3:siteManava:speciesPocillopora  -0.239233   0.216298
    ## timepointtimepoint4:siteManava:speciesPocillopora  -0.002482   0.172690
    ## timepointtimepoint2:siteMahana:speciesPorites       0.061685   0.178128
    ## timepointtimepoint3:siteMahana:speciesPorites       0.130028   0.215863
    ## timepointtimepoint4:siteMahana:speciesPorites       0.254711   0.175999
    ## timepointtimepoint2:siteManava:speciesPorites       0.076750   0.171012
    ## timepointtimepoint3:siteManava:speciesPorites      -0.224832   0.213262
    ## timepointtimepoint4:siteManava:speciesPorites       0.007546   0.174127
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       394.987665  15.238  < 2e-16
    ## timepointtimepoint2                               339.310032   2.891  0.00408
    ## timepointtimepoint3                               389.849950   0.959  0.33813
    ## timepointtimepoint4                               316.894093  -1.722  0.08608
    ## siteMahana                                        395.377499   5.160 3.92e-07
    ## siteManava                                        394.896458   2.840  0.00474
    ## speciesPocillopora                                395.185696  -2.713  0.00695
    ## speciesPorites                                    394.927111  -2.360  0.01875
    ## timepointtimepoint2:siteMahana                    330.067077  -2.646  0.00854
    ## timepointtimepoint3:siteMahana                    385.010961  -2.249  0.02511
    ## timepointtimepoint4:siteMahana                    332.521699  -2.710  0.00707
    ## timepointtimepoint2:siteManava                    332.085166  -2.740  0.00649
    ## timepointtimepoint3:siteManava                    378.460286  -0.079  0.93720
    ## timepointtimepoint4:siteManava                    308.722860  -1.452  0.14745
    ## timepointtimepoint2:speciesPocillopora            310.996294  -0.518  0.60515
    ## timepointtimepoint3:speciesPocillopora            373.576990   0.735  0.46302
    ## timepointtimepoint4:speciesPocillopora            298.192882  -0.093  0.92585
    ## timepointtimepoint2:speciesPorites                313.226577  -1.335  0.18268
    ## timepointtimepoint3:speciesPorites                374.358530  -0.160  0.87319
    ## timepointtimepoint4:speciesPorites                300.141608   0.313  0.75468
    ## siteMahana:speciesPocillopora                     395.332147  -0.702  0.48322
    ## siteManava:speciesPocillopora                     394.850569   1.088  0.27739
    ## siteMahana:speciesPorites                         395.211427  -2.437  0.01526
    ## siteManava:speciesPorites                         394.695248   1.210  0.22705
    ## timepointtimepoint2:siteMahana:speciesPocillopora 306.083259  -0.840  0.40179
    ## timepointtimepoint3:siteMahana:speciesPocillopora 359.273470  -0.649  0.51705
    ## timepointtimepoint4:siteMahana:speciesPocillopora 306.288526   1.323  0.18673
    ## timepointtimepoint2:siteManava:speciesPocillopora 304.248380  -0.277  0.78186
    ## timepointtimepoint3:siteManava:speciesPocillopora 353.614657  -1.106  0.26946
    ## timepointtimepoint4:siteManava:speciesPocillopora 290.930187  -0.014  0.98854
    ## timepointtimepoint2:siteMahana:speciesPorites     307.193522   0.346  0.72936
    ## timepointtimepoint3:siteMahana:speciesPorites     360.647838   0.602  0.54731
    ## timepointtimepoint4:siteMahana:speciesPorites     307.473756   1.447  0.14885
    ## timepointtimepoint2:siteManava:speciesPorites     305.469187   0.449  0.65390
    ## timepointtimepoint3:siteManava:speciesPorites     353.257478  -1.054  0.29249
    ## timepointtimepoint4:siteManava:speciesPorites     293.265566   0.043  0.96546
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                               ** 
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                               .  
    ## siteMahana                                        ***
    ## siteManava                                        ** 
    ## speciesPocillopora                                ** 
    ## speciesPorites                                    *  
    ## timepointtimepoint2:siteMahana                    ** 
    ## timepointtimepoint3:siteMahana                    *  
    ## timepointtimepoint4:siteMahana                    ** 
    ## timepointtimepoint2:siteManava                    ** 
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## siteMahana:speciesPocillopora                        
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                         *  
    ## siteManava:speciesPorites                            
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora    
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora    
    ## timepointtimepoint3:siteManava:speciesPocillopora    
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites        
    ## timepointtimepoint3:siteMahana:speciesPorites        
    ## timepointtimepoint4:siteMahana:speciesPorites        
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites        
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of site and time within
each species.

*Acropora*

`chla_model_acr<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0 & species=="Acropora"))`

``` r
chla_model_acr<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0 & species=="Acropora"))
qqPlot(residuals(chla_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

    ## 218 421 
    ##  63 111

Generate a Type III Anova of model.

``` r
anova(chla_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)    
    ## timepoint      2.75985 0.91995     3   101 13.5451 1.7e-07 ***
    ## site           0.27663 0.13832     2   101  2.0365 0.13580    
    ## timepoint:site 1.00481 0.16747     6   101  2.4658 0.02881 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chla_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chla.ug.mgAFDW) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chla.ug.mgAFDW > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: 40.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7863 -0.5977 -0.1531  0.5029  3.7289 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.06792  0.2606  
    ## Number of obs: 113, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.92972    0.07523 101.00000  12.358  < 2e-16
    ## timepointtimepoint2              0.26461    0.11492 101.00000   2.303   0.0234
    ## timepointtimepoint3              0.15769    0.19904 101.00000   0.792   0.4301
    ## timepointtimepoint4             -0.18381    0.13030 101.00000  -1.411   0.1614
    ## siteMahana                       0.43580    0.10433 101.00000   4.177 6.28e-05
    ## siteManava                       0.23167    0.10093 101.00000   2.295   0.0238
    ## timepointtimepoint2:siteMahana  -0.37114    0.17248 101.00000  -2.152   0.0338
    ## timepointtimepoint3:siteMahana  -0.41809    0.22723 101.00000  -1.840   0.0687
    ## timepointtimepoint4:siteMahana  -0.35883    0.16692 101.00000  -2.150   0.0340
    ## timepointtimepoint2:siteManava  -0.35160    0.16194 101.00000  -2.171   0.0323
    ## timepointtimepoint3:siteManava  -0.02123    0.22736 101.00000  -0.093   0.9258
    ## timepointtimepoint4:siteManava  -0.18959    0.16638 101.00000  -1.140   0.2572
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            *  
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                     ***
    ## siteManava                     *  
    ## timepointtimepoint2:siteMahana *  
    ## timepointtimepoint3:siteMahana .  
    ## timepointtimepoint4:siteMahana *  
    ## timepointtimepoint2:siteManava *  
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.655                                                    
    ## tmpnttmpnt3      -0.378  0.247                                             
    ## tmpnttmpnt4      -0.577  0.378  0.218                                      
    ## siteMahana       -0.721  0.472  0.273  0.416                               
    ## siteManava       -0.745  0.488  0.282  0.430  0.537                        
    ## tmpnttmpnt2:stMh  0.436 -0.666 -0.165 -0.252 -0.605 -0.325                 
    ## tmpnttmpnt3:stMh  0.331 -0.217 -0.876 -0.191 -0.459 -0.247  0.278          
    ## tmpnttmpnt4:stMh  0.451 -0.295 -0.170 -0.781 -0.625 -0.336  0.378          
    ## tmpnttmpnt2:stMn  0.465 -0.710 -0.176 -0.268 -0.335 -0.623  0.473          
    ## tmpnttmpnt3:stMn  0.331 -0.217 -0.875 -0.191 -0.239 -0.444  0.144          
    ## tmpnttmpnt4:stMn  0.452 -0.296 -0.171 -0.783 -0.326 -0.607  0.197          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.287                                            
    ## tmpnttmpnt2:stMn  0.154            0.209                           
    ## tmpnttmpnt3:stMn  0.767            0.149            0.277          
    ## tmpnttmpnt4:stMn  0.150            0.611            0.378          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.269          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

*Pocillopora*

`chla_model_poc<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0 & species=="Pocillopora"))`

``` r
chla_model_poc<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0 & species=="Pocillopora"))
qqPlot(residuals(chla_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-76-1.png)<!-- -->

    ## 173  50 
    ##  76  10

Generate a Type III Anova of model.

``` r
anova(chla_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      2.48237 0.82746     3 110.02 22.9023 1.357e-11 ***
    ## site           0.48318 0.24159     2  36.53  6.6867  0.003353 ** 
    ## timepoint:site 1.69643 0.28274     6 109.94  7.8256 5.022e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chla_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chla.ug.mgAFDW) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chla.ug.mgAFDW > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -34.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1491 -0.5406 -0.0246  0.4595  3.5457 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001701 0.04125 
    ##  Residual                   0.036130 0.19008 
    ## Number of obs: 161, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.70001    0.05393 148.46607  12.980  < 2e-16
    ## timepointtimepoint2              0.20275    0.07338 111.13824   2.763 0.006707
    ## timepointtimepoint3              0.28609    0.07467 109.00485   3.831 0.000214
    ## timepointtimepoint4             -0.19103    0.07214 108.83737  -2.648 0.009294
    ## siteMahana                       0.35326    0.07627 148.43581   4.632 7.86e-06
    ## siteManava                       0.35703    0.07369 148.29811   4.845 3.16e-06
    ## timepointtimepoint2:siteMahana  -0.51329    0.10290 110.09832  -4.988 2.29e-06
    ## timepointtimepoint3:siteMahana  -0.55162    0.10561 109.10920  -5.223 8.49e-07
    ## timepointtimepoint4:siteMahana  -0.12928    0.10499 110.32764  -1.231 0.220804
    ## timepointtimepoint2:siteManava  -0.40157    0.10189 108.94743  -3.941 0.000144
    ## timepointtimepoint3:siteManava  -0.25369    0.10633 111.35611  -2.386 0.018728
    ## timepointtimepoint4:siteManava  -0.19511    0.10202 108.93189  -1.913 0.058438
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            ** 
    ## timepointtimepoint3            ***
    ## timepointtimepoint4            ** 
    ## siteMahana                     ***
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ***
    ## timepointtimepoint3:siteMahana ***
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava ***
    ## timepointtimepoint3:siteManava *  
    ## timepointtimepoint4:siteManava .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.707                                                    
    ## tmpnttmpnt3      -0.692  0.509                                             
    ## tmpnttmpnt4      -0.718  0.528  0.518                                      
    ## siteMahana       -0.707  0.500  0.490  0.508                               
    ## siteManava       -0.732  0.517  0.507  0.526  0.517                        
    ## tmpnttmpnt2:stMh  0.504 -0.713 -0.363 -0.377 -0.712 -0.369                 
    ## tmpnttmpnt3:stMh  0.490 -0.360 -0.707 -0.366 -0.692 -0.358  0.513          
    ## tmpnttmpnt4:stMh  0.494 -0.363 -0.356 -0.687 -0.697 -0.361  0.517          
    ## tmpnttmpnt2:stMn  0.509 -0.720 -0.366 -0.380 -0.360 -0.693  0.514          
    ## tmpnttmpnt3:stMn  0.486 -0.357 -0.702 -0.363 -0.344 -0.663  0.255          
    ## tmpnttmpnt4:stMn  0.508 -0.373 -0.366 -0.707 -0.359 -0.692  0.266          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.504                                            
    ## tmpnttmpnt2:stMn  0.259            0.261                           
    ## tmpnttmpnt3:stMn  0.497            0.250            0.479          
    ## tmpnttmpnt4:stMn  0.259            0.486            0.501          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.481

*Porites*

`chla_model_por<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0 & species=="Porites"))`

``` r
chla_model_por<-lmer(log(1+chla.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chla.ug.mgAFDW>0 & species=="Porites"))
qqPlot(residuals(chla_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

    ## 109  88 
    ##  28   8

Generate a Type III Anova of model.

``` r
anova(chla_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      1.15466 0.38489     3 95.100  12.375 6.661e-07 ***
    ## site           0.99446 0.49723     2 27.129  15.988 2.586e-05 ***
    ## timepoint:site 0.49507 0.08251     6 94.980   2.653    0.0201 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chla_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chla.ug.mgAFDW) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chla.ug.mgAFDW > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: -44
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0991 -0.5207 -0.0360  0.3962  7.0965 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.004814 0.06938 
    ##  Residual                   0.031101 0.17636 
    ## Number of obs: 159, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.73519    0.05058 138.44177  14.536  < 2e-16
    ## timepointtimepoint2              0.09761    0.06995  97.18832   1.396  0.16604
    ## timepointtimepoint3              0.12374    0.06566  90.76824   1.885  0.06267
    ## timepointtimepoint4             -0.13925    0.06833  94.85972  -2.038  0.04435
    ## siteMahana                       0.14886    0.07287 139.04630   2.043  0.04295
    ## siteManava                       0.36785    0.07037 137.29564   5.227 6.27e-07
    ## timepointtimepoint2:siteMahana  -0.29901    0.09694  95.52424  -3.085  0.00267
    ## timepointtimepoint3:siteMahana  -0.27628    0.09813  95.25889  -2.815  0.00592
    ## timepointtimepoint4:siteMahana  -0.10400    0.09860  95.71258  -1.055  0.29420
    ## timepointtimepoint2:siteManava  -0.27296    0.09593  94.14886  -2.845  0.00544
    ## timepointtimepoint3:siteManava  -0.22775    0.09388  92.21337  -2.426  0.01722
    ## timepointtimepoint4:siteManava  -0.17474    0.09692  95.25997  -1.803  0.07455
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3            .  
    ## timepointtimepoint4            *  
    ## siteMahana                     *  
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ** 
    ## timepointtimepoint3:siteMahana ** 
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava ** 
    ## timepointtimepoint3:siteManava *  
    ## timepointtimepoint4:siteManava .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.634                                                    
    ## tmpnttmpnt3      -0.674  0.488                                             
    ## tmpnttmpnt4      -0.648  0.468  0.499                                      
    ## siteMahana       -0.694  0.440  0.468  0.450                               
    ## siteManava       -0.719  0.455  0.484  0.466  0.499                        
    ## tmpnttmpnt2:stMh  0.457 -0.722 -0.352 -0.337 -0.661 -0.329                 
    ## tmpnttmpnt3:stMh  0.451 -0.327 -0.669 -0.334 -0.651 -0.324  0.490          
    ## tmpnttmpnt4:stMh  0.449 -0.324 -0.346 -0.693 -0.648 -0.323  0.486          
    ## tmpnttmpnt2:stMn  0.462 -0.729 -0.356 -0.341 -0.321 -0.639  0.526          
    ## tmpnttmpnt3:stMn  0.471 -0.341 -0.699 -0.349 -0.327 -0.652  0.246          
    ## tmpnttmpnt4:stMn  0.457 -0.330 -0.352 -0.705 -0.317 -0.632  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.486                                            
    ## tmpnttmpnt2:stMn  0.238            0.236                           
    ## tmpnttmpnt3:stMn  0.468            0.242            0.479          
    ## tmpnttmpnt4:stMn  0.236            0.489            0.463          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.478

#### Chlorophyll c2:

Build a mixed model for univariate analysis and examine data
distribution.

`chlc2_model<-lmer(chlc2.ug.mgAFDW~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0))`

``` r
chlc2_model<-lmer(chlc2.ug.mgAFDW~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0))
qqPlot(residuals(chlc2_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->

    ## 109 218 
    ## 104 209

Residuals are not normally distributed. Attempt with log transformation.

`chlc2_model<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0))`

``` r
chlc2_model<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0))
qqPlot(residuals(chlc2_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-81-1.png)<!-- -->

    ## 109 218 
    ## 104 209

Return to check residual distribution in this model.

Generate a Type III Anova of model.

``` r
anova(chlc2_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint              11.2548  3.7516     3   397 43.9910 < 2.2e-16 ***
    ## species                 1.7706  0.8853     2   397 10.3810 4.032e-05 ***
    ## site                    2.9639  1.4820     2   397 17.3774 5.826e-08 ***
    ## timepoint:species       0.3157  0.0526     6   397  0.6169    0.7168    
    ## timepoint:site          6.5801  1.0967     6   397 12.8597 2.726e-13 ***
    ## species:site            0.2064  0.0516     4   397  0.6051    0.6591    
    ## timepoint:species:site  1.2330  0.1027    12   397  1.2048    0.2772    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chlc2_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chlc2.ug.mgAFDW) ~ timepoint * species * site + (1 |  
    ##     colony_id_corr)
    ##    Data: master
    ##  Subset: c(chlc2.ug.mgAFDW > 0)
    ## 
    ## REML criterion at convergence: 237.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1759 -0.3557 -0.0915  0.2969  7.3659 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.000   
    ##  Residual                   0.08528  0.292   
    ## Number of obs: 433, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                        3.755e-01  8.430e-02
    ## timepointtimepoint2                                1.052e-01  1.288e-01
    ## timepointtimepoint3                                1.573e-01  2.230e-01
    ## timepointtimepoint4                                3.103e-02  1.460e-01
    ## speciesPocillopora                                -6.921e-02  1.169e-01
    ## speciesPorites                                    -3.545e-02  1.149e-01
    ## siteMahana                                         8.355e-01  1.169e-01
    ## siteManava                                         5.614e-01  1.131e-01
    ## timepointtimepoint2:speciesPocillopora            -6.839e-02  1.710e-01
    ## timepointtimepoint3:speciesPocillopora             6.570e-02  2.507e-01
    ## timepointtimepoint4:speciesPocillopora            -1.061e-01  1.832e-01
    ## timepointtimepoint2:speciesPorites                -1.263e-01  1.726e-01
    ## timepointtimepoint3:speciesPorites                -1.810e-02  2.480e-01
    ## timepointtimepoint4:speciesPorites                -7.614e-02  1.843e-01
    ## timepointtimepoint2:siteMahana                    -8.547e-01  1.933e-01
    ## timepointtimepoint3:siteMahana                    -6.040e-01  2.546e-01
    ## timepointtimepoint4:siteMahana                    -7.837e-01  1.870e-01
    ## timepointtimepoint2:siteManava                    -6.425e-01  1.815e-01
    ## timepointtimepoint3:siteManava                     5.588e-02  2.548e-01
    ## timepointtimepoint4:siteManava                    -6.174e-01  1.864e-01
    ## speciesPocillopora:siteMahana                     -2.625e-01  1.637e-01
    ## speciesPorites:siteMahana                         -4.203e-01  1.622e-01
    ## speciesPocillopora:siteManava                     -1.944e-01  1.582e-01
    ## speciesPorites:siteManava                          9.936e-04  1.567e-01
    ## timepointtimepoint2:speciesPocillopora:siteMahana  2.078e-01  2.495e-01
    ## timepointtimepoint3:speciesPocillopora:siteMahana -1.732e-02  3.018e-01
    ## timepointtimepoint4:speciesPocillopora:siteMahana  3.786e-01  2.468e-01
    ## timepointtimepoint2:speciesPorites:siteMahana      3.805e-01  2.506e-01
    ## timepointtimepoint3:speciesPorites:siteMahana      2.601e-01  3.015e-01
    ## timepointtimepoint4:speciesPorites:siteMahana      4.806e-01  2.476e-01
    ## timepointtimepoint2:speciesPocillopora:siteManava  2.527e-01  2.395e-01
    ## timepointtimepoint3:speciesPocillopora:siteManava -1.372e-01  3.024e-01
    ## timepointtimepoint4:speciesPocillopora:siteManava  3.006e-01  2.434e-01
    ## timepointtimepoint2:speciesPorites:siteManava      1.087e-01  2.406e-01
    ## timepointtimepoint3:speciesPorites:siteManava     -4.308e-01  2.982e-01
    ## timepointtimepoint4:speciesPorites:siteManava      8.466e-02  2.454e-01
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                        3.970e+02   4.454 1.10e-05
    ## timepointtimepoint2                                3.970e+02   0.817 0.414447
    ## timepointtimepoint3                                3.970e+02   0.705 0.480981
    ## timepointtimepoint4                                3.970e+02   0.213 0.831801
    ## speciesPocillopora                                 3.970e+02  -0.592 0.554148
    ## speciesPorites                                     3.970e+02  -0.309 0.757809
    ## siteMahana                                         3.970e+02   7.147 4.28e-12
    ## siteManava                                         3.970e+02   4.963 1.03e-06
    ## timepointtimepoint2:speciesPocillopora             3.970e+02  -0.400 0.689364
    ## timepointtimepoint3:speciesPocillopora             3.970e+02   0.262 0.793418
    ## timepointtimepoint4:speciesPocillopora             3.970e+02  -0.579 0.562813
    ## timepointtimepoint2:speciesPorites                 3.970e+02  -0.732 0.464666
    ## timepointtimepoint3:speciesPorites                 3.970e+02  -0.073 0.941868
    ## timepointtimepoint4:speciesPorites                 3.970e+02  -0.413 0.679772
    ## timepointtimepoint2:siteMahana                     3.970e+02  -4.422 1.26e-05
    ## timepointtimepoint3:siteMahana                     3.970e+02  -2.372 0.018171
    ## timepointtimepoint4:siteMahana                     3.970e+02  -4.190 3.45e-05
    ## timepointtimepoint2:siteManava                     3.970e+02  -3.541 0.000446
    ## timepointtimepoint3:siteManava                     3.970e+02   0.219 0.826503
    ## timepointtimepoint4:siteManava                     3.970e+02  -3.312 0.001013
    ## speciesPocillopora:siteMahana                      3.970e+02  -1.604 0.109505
    ## speciesPorites:siteMahana                          3.970e+02  -2.591 0.009931
    ## speciesPocillopora:siteManava                      3.970e+02  -1.228 0.219991
    ## speciesPorites:siteManava                          3.970e+02   0.006 0.994945
    ## timepointtimepoint2:speciesPocillopora:siteMahana  3.970e+02   0.833 0.405373
    ## timepointtimepoint3:speciesPocillopora:siteMahana  3.970e+02  -0.057 0.954266
    ## timepointtimepoint4:speciesPocillopora:siteMahana  3.970e+02   1.534 0.125823
    ## timepointtimepoint2:speciesPorites:siteMahana      3.970e+02   1.518 0.129696
    ## timepointtimepoint3:speciesPorites:siteMahana      3.970e+02   0.863 0.388873
    ## timepointtimepoint4:speciesPorites:siteMahana      3.970e+02   1.941 0.052978
    ## timepointtimepoint2:speciesPocillopora:siteManava  3.970e+02   1.055 0.292012
    ## timepointtimepoint3:speciesPocillopora:siteManava  3.970e+02  -0.454 0.650417
    ## timepointtimepoint4:speciesPocillopora:siteManava  3.970e+02   1.235 0.217555
    ## timepointtimepoint2:speciesPorites:siteManava      3.970e+02   0.452 0.651807
    ## timepointtimepoint3:speciesPorites:siteManava      3.970e+02  -1.445 0.149348
    ## timepointtimepoint4:speciesPorites:siteManava      3.970e+02   0.345 0.730256
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                                  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## speciesPocillopora                                   
    ## speciesPorites                                       
    ## siteMahana                                        ***
    ## siteManava                                        ***
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## timepointtimepoint2:siteMahana                    ***
    ## timepointtimepoint3:siteMahana                    *  
    ## timepointtimepoint4:siteMahana                    ***
    ## timepointtimepoint2:siteManava                    ***
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                    ** 
    ## speciesPocillopora:siteMahana                        
    ## speciesPorites:siteMahana                         ** 
    ## speciesPocillopora:siteManava                        
    ## speciesPorites:siteManava                            
    ## timepointtimepoint2:speciesPocillopora:siteMahana    
    ## timepointtimepoint3:speciesPocillopora:siteMahana    
    ## timepointtimepoint4:speciesPocillopora:siteMahana    
    ## timepointtimepoint2:speciesPorites:siteMahana        
    ## timepointtimepoint3:speciesPorites:siteMahana        
    ## timepointtimepoint4:speciesPorites:siteMahana     .  
    ## timepointtimepoint2:speciesPocillopora:siteManava    
    ## timepointtimepoint3:speciesPocillopora:siteManava    
    ## timepointtimepoint4:speciesPocillopora:siteManava    
    ## timepointtimepoint2:speciesPorites:siteManava        
    ## timepointtimepoint3:speciesPorites:siteManava        
    ## timepointtimepoint4:speciesPorites:siteManava        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

Since species is significant, test for effects of site and time within
species.

*Acropora*

`chlc2_model_acr<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0 & species=="Acropora))`

``` r
chlc2_model_acr<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0 & species=="Acropora"))
qqPlot(residuals(chlc2_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-83-1.png)<!-- -->

    ## 218 421 
    ##  63 111

Generate a Type III Anova of model.

``` r
anova(chlc2_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      4.2589 1.41965     3 84.374 11.3654 2.466e-06 ***
    ## site           1.0598 0.52990     2 36.810  4.2423 0.0219825 *  
    ## timepoint:site 3.6619 0.61032     6 81.855  4.8861 0.0002587 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chlc2_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chlc2.ug.mgAFDW) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chlc2.ug.mgAFDW > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: 103.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6018 -0.3471 -0.1160  0.3766  3.4504 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001137 0.03372 
    ##  Residual                   0.124910 0.35343 
    ## Number of obs: 113, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.37535    0.10249 100.99107   3.662 0.000400
    ## timepointtimepoint2              0.10576    0.15615  82.33380   0.677 0.500111
    ## timepointtimepoint3              0.15728    0.27097  99.31887   0.580 0.562944
    ## timepointtimepoint4              0.03165    0.17698  76.94102   0.179 0.858534
    ## siteMahana                       0.83615    0.14212 100.99376   5.883 5.26e-08
    ## siteManava                       0.56126    0.13750 100.99064   4.082 8.95e-05
    ## timepointtimepoint2:siteMahana  -0.85303    0.23433  80.36848  -3.640 0.000480
    ## timepointtimepoint3:siteMahana  -0.60332    0.30925  97.68292  -1.951 0.053931
    ## timepointtimepoint4:siteMahana  -0.78487    0.22679  80.99134  -3.461 0.000862
    ## timepointtimepoint2:siteManava  -0.64307    0.22002  81.12971  -2.923 0.004493
    ## timepointtimepoint3:siteManava   0.05624    0.30935  95.91189   0.182 0.856115
    ## timepointtimepoint4:siteManava  -0.61694    0.22593  74.74814  -2.731 0.007879
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                     ***
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ***
    ## timepointtimepoint3:siteMahana .  
    ## timepointtimepoint4:siteMahana ***
    ## timepointtimepoint2:siteManava ** 
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.652                                                    
    ## tmpnttmpnt3      -0.377  0.249                                             
    ## tmpnttmpnt4      -0.574  0.377  0.218                                      
    ## siteMahana       -0.721  0.470  0.272  0.414                               
    ## siteManava       -0.745  0.486  0.281  0.428  0.537                        
    ## tmpnttmpnt2:stMh  0.435 -0.666 -0.166 -0.252 -0.602 -0.324                 
    ## tmpnttmpnt3:stMh  0.330 -0.218 -0.876 -0.191 -0.457 -0.246  0.278          
    ## tmpnttmpnt4:stMh  0.448 -0.295 -0.170 -0.780 -0.622 -0.334  0.377          
    ## tmpnttmpnt2:stMn  0.463 -0.710 -0.176 -0.268 -0.334 -0.621  0.473          
    ## tmpnttmpnt3:stMn  0.330 -0.218 -0.876 -0.191 -0.238 -0.442  0.145          
    ## tmpnttmpnt4:stMn  0.450 -0.296 -0.171 -0.783 -0.324 -0.603  0.197          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.287                                            
    ## tmpnttmpnt2:stMn  0.155            0.209                           
    ## tmpnttmpnt3:stMn  0.768            0.149            0.278          
    ## tmpnttmpnt4:stMn  0.150            0.611            0.378          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.269

*Pocillopora*

`chlc2_model_poc<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0 & species=="Pocillopora))`

``` r
chlc2_model_poc<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0 & species=="Pocillopora"))
qqPlot(residuals(chlc2_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-85-1.png)<!-- -->

    ## 50 43 
    ## 10  3

Generate a Type III Anova of model.

``` r
anova(chlc2_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      3.7698 1.25660     3   149 22.1419 6.505e-12 ***
    ## site           0.9549 0.47744     2   149  8.4127 0.0003455 ***
    ## timepoint:site 2.3385 0.38976     6   149  6.8677 1.878e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chlc2_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chlc2.ug.mgAFDW) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chlc2.ug.mgAFDW > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: 26.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1624 -0.4052 -0.1152  0.3581  5.0793 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.05675  0.2382  
    ## Number of obs: 161, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.30627    0.06607 149.00000   4.635 7.72e-06
    ## timepointtimepoint2              0.03681    0.09176 149.00000   0.401  0.68889
    ## timepointtimepoint3              0.22304    0.09344 149.00000   2.387  0.01824
    ## timepointtimepoint4             -0.07507    0.09027 149.00000  -0.832  0.40694
    ## siteMahana                       0.57301    0.09344 149.00000   6.132 7.41e-09
    ## siteManava                       0.36699    0.09027 149.00000   4.065 7.74e-05
    ## timepointtimepoint2:siteMahana  -0.64691    0.12872 149.00000  -5.026 1.42e-06
    ## timepointtimepoint3:siteMahana  -0.62128    0.13214 149.00000  -4.702 5.83e-06
    ## timepointtimepoint4:siteMahana  -0.40510    0.13132 149.00000  -3.085  0.00243
    ## timepointtimepoint2:siteManava  -0.38983    0.12750 149.00000  -3.057  0.00265
    ## timepointtimepoint3:siteManava  -0.08128    0.13294 149.00000  -0.611  0.54186
    ## timepointtimepoint4:siteManava  -0.31679    0.12766 149.00000  -2.481  0.01420
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3            *  
    ## timepointtimepoint4               
    ## siteMahana                     ***
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ***
    ## timepointtimepoint3:siteMahana ***
    ## timepointtimepoint4:siteMahana ** 
    ## timepointtimepoint2:siteManava ** 
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.720                                                    
    ## tmpnttmpnt3      -0.707  0.509                                             
    ## tmpnttmpnt4      -0.732  0.527  0.518                                      
    ## siteMahana       -0.707  0.509  0.500  0.518                               
    ## siteManava       -0.732  0.527  0.518  0.536  0.518                        
    ## tmpnttmpnt2:stMh  0.513 -0.713 -0.363 -0.376 -0.726 -0.376                 
    ## tmpnttmpnt3:stMh  0.500 -0.360 -0.707 -0.366 -0.707 -0.366  0.513          
    ## tmpnttmpnt4:stMh  0.503 -0.362 -0.356 -0.687 -0.712 -0.368  0.517          
    ## tmpnttmpnt2:stMn  0.518 -0.720 -0.366 -0.379 -0.366 -0.708  0.513          
    ## tmpnttmpnt3:stMn  0.497 -0.358 -0.703 -0.364 -0.351 -0.679  0.255          
    ## tmpnttmpnt4:stMn  0.518 -0.373 -0.366 -0.707 -0.366 -0.707  0.266          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.503                                            
    ## tmpnttmpnt2:stMn  0.259            0.261                           
    ## tmpnttmpnt3:stMn  0.497            0.250            0.481          
    ## tmpnttmpnt4:stMn  0.259            0.486            0.501          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.480          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

*Porites*

`chlc2_model_por<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0 & species=="Porites))`

``` r
chlc2_model_por<-lmer(log(1+chlc2.ug.mgAFDW)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(chlc2.ug.mgAFDW>0 & species=="Porites"))
qqPlot(residuals(chlc2_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-87-1.png)<!-- -->

    ## 109  87 
    ##  28   7

Generate a Type III Anova of model.

``` r
anova(chlc2_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      3.6146 1.20486     3   147 13.9791 4.564e-08 ***
    ## site           1.1335 0.56675     2   147  6.5756  0.001840 ** 
    ## timepoint:site 1.5475 0.25791     6   147  2.9924  0.008689 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chlc2_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + chlc2.ug.mgAFDW) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(chlc2.ug.mgAFDW > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: 87.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5646 -0.3353 -0.0583  0.2260  7.3270 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.08619  0.2936  
    ## Number of obs: 159, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.34003    0.07846 147.00000   4.334 2.70e-05
    ## timepointtimepoint2             -0.02110    0.11549 147.00000  -0.183 0.855274
    ## timepointtimepoint3              0.13923    0.10910 147.00000   1.276 0.203896
    ## timepointtimepoint4             -0.04510    0.11308 147.00000  -0.399 0.690563
    ## siteMahana                       0.41525    0.11308 147.00000   3.672 0.000336
    ## siteManava                       0.56237    0.10910 147.00000   5.155 8.06e-07
    ## timepointtimepoint2:siteMahana  -0.47422    0.16036 147.00000  -2.957 0.003618
    ## timepointtimepoint3:siteMahana  -0.34385    0.16238 147.00000  -2.118 0.035896
    ## timepointtimepoint4:siteMahana  -0.30309    0.16309 147.00000  -1.858 0.065109
    ## timepointtimepoint2:siteManava  -0.53385    0.15888 147.00000  -3.360 0.000992
    ## timepointtimepoint3:siteManava  -0.37493    0.15582 147.00000  -2.406 0.017359
    ## timepointtimepoint4:siteManava  -0.53275    0.16036 147.00000  -3.322 0.001127
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                     ***
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ** 
    ## timepointtimepoint3:siteMahana *  
    ## timepointtimepoint4:siteMahana .  
    ## timepointtimepoint2:siteManava ***
    ## timepointtimepoint3:siteManava *  
    ## timepointtimepoint4:siteManava ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.679                                                    
    ## tmpnttmpnt3      -0.719  0.489                                             
    ## tmpnttmpnt4      -0.694  0.471  0.499                                      
    ## siteMahana       -0.694  0.471  0.499  0.481                               
    ## siteManava       -0.719  0.489  0.517  0.499  0.499                        
    ## tmpnttmpnt2:stMh  0.489 -0.720 -0.352 -0.340 -0.705 -0.352                 
    ## tmpnttmpnt3:stMh  0.483 -0.328 -0.672 -0.335 -0.696 -0.348  0.491          
    ## tmpnttmpnt4:stMh  0.481 -0.327 -0.346 -0.693 -0.693 -0.346  0.489          
    ## tmpnttmpnt2:stMn  0.494 -0.727 -0.355 -0.343 -0.343 -0.687  0.524          
    ## tmpnttmpnt3:stMn  0.504 -0.342 -0.700 -0.349 -0.349 -0.700  0.246          
    ## tmpnttmpnt4:stMn  0.489 -0.332 -0.352 -0.705 -0.340 -0.680  0.239          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.483                                            
    ## tmpnttmpnt2:stMn  0.239            0.238                           
    ## tmpnttmpnt3:stMn  0.470            0.242            0.481          
    ## tmpnttmpnt4:stMn  0.236            0.489            0.467          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.476          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

#### Total Chlorophyll (a + c2):

Build a mixed model for univariate analysis and examine data
distribution.

`chltotal_model<-lmer(Total_Chl~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0))`

``` r
chltotal_model<-lmer(Total_Chl~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0))
qqPlot(residuals(chltotal_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-89-1.png)<!-- -->

    ## 109 218 
    ## 104 209

Residuals are not normally distributed. Attempt with log transformation.

`chltotal_model<-lmer(log(1+Total_Chl)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0))`

``` r
chltotal_model<-lmer(log(1+Total_Chl)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0))
qqPlot(residuals(chltotal_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-90-1.png)<!-- -->

    ## 109 218 
    ## 104 209

Generate a Type III Anova of model.

``` r
anova(chltotal_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint              11.7178  3.9059     3   397 45.2337 < 2.2e-16 ***
    ## site                    2.8718  1.4359     2   397 16.6287 1.161e-07 ***
    ## species                 5.2448  2.6224     2   397 30.3696 5.336e-13 ***
    ## timepoint:site          6.1919  1.0320     6   397 11.9512 2.401e-12 ***
    ## timepoint:species       0.6818  0.1136     6   397  1.3159   0.24876    
    ## site:species            0.8461  0.2115     4   397  2.4497   0.04572 *  
    ## timepoint:site:species  1.0626  0.0886    12   397  1.0255   0.42418    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotal_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl) ~ timepoint * site * species + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl > 0)
    ## 
    ## REML criterion at convergence: 242
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4249 -0.5018 -0.0624  0.4266  6.7987 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.08635  0.2939  
    ## Number of obs: 433, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         1.093535   0.084828
    ## timepointtimepoint2                                 0.271835   0.129577
    ## timepointtimepoint3                                 0.209067   0.224434
    ## timepointtimepoint4                                -0.141790   0.146927
    ## siteMahana                                          0.743750   0.117635
    ## siteManava                                          0.463407   0.113809
    ## speciesPocillopora                                 -0.235001   0.117635
    ## speciesPorites                                     -0.182802   0.115601
    ## timepointtimepoint2:siteMahana                     -0.685543   0.194484
    ## timepointtimepoint3:siteMahana                     -0.620104   0.256219
    ## timepointtimepoint4:siteMahana                     -0.645647   0.188217
    ## timepointtimepoint2:siteManava                     -0.601395   0.182594
    ## timepointtimepoint3:siteManava                     -0.014878   0.256362
    ## timepointtimepoint4:siteManava                     -0.465546   0.187601
    ## timepointtimepoint2:speciesPocillopora             -0.079960   0.172048
    ## timepointtimepoint3:speciesPocillopora              0.146107   0.252300
    ## timepointtimepoint4:speciesPocillopora             -0.065284   0.184354
    ## timepointtimepoint2:speciesPorites                 -0.197031   0.173649
    ## timepointtimepoint3:speciesPorites                 -0.032726   0.249590
    ## timepointtimepoint4:speciesPorites                 -0.001579   0.185466
    ## siteMahana:speciesPocillopora                      -0.163384   0.164690
    ## siteManava:speciesPocillopora                       0.019089   0.159221
    ## siteMahana:speciesPorites                          -0.404052   0.163243
    ## siteManava:speciesPorites                           0.114620   0.157724
    ## timepointtimepoint2:siteMahana:speciesPocillopora  -0.065680   0.251064
    ## timepointtimepoint3:siteMahana:speciesPocillopora  -0.143004   0.303673
    ## timepointtimepoint4:siteMahana:speciesPocillopora   0.356717   0.248319
    ## timepointtimepoint2:siteManava:speciesPocillopora   0.068844   0.240988
    ## timepointtimepoint3:siteManava:speciesPocillopora  -0.244399   0.304324
    ## timepointtimepoint4:siteManava:speciesPocillopora   0.153813   0.244932
    ## timepointtimepoint2:siteMahana:speciesPorites       0.180445   0.252164
    ## timepointtimepoint3:siteMahana:speciesPorites       0.219437   0.303422
    ## timepointtimepoint4:siteMahana:speciesPorites       0.414365   0.249146
    ## timepointtimepoint2:siteManava:speciesPorites       0.114866   0.242134
    ## timepointtimepoint3:siteManava:speciesPorites      -0.358796   0.300075
    ## timepointtimepoint4:siteManava:speciesPorites       0.056444   0.246894
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       397.000000  12.891  < 2e-16
    ## timepointtimepoint2                               397.000000   2.098 0.036549
    ## timepointtimepoint3                               397.000000   0.932 0.352145
    ## timepointtimepoint4                               397.000000  -0.965 0.335112
    ## siteMahana                                        397.000000   6.323 6.93e-10
    ## siteManava                                        397.000000   4.072 5.63e-05
    ## speciesPocillopora                                397.000000  -1.998 0.046430
    ## speciesPorites                                    397.000000  -1.581 0.114603
    ## timepointtimepoint2:siteMahana                    397.000000  -3.525 0.000473
    ## timepointtimepoint3:siteMahana                    397.000000  -2.420 0.015959
    ## timepointtimepoint4:siteMahana                    397.000000  -3.430 0.000666
    ## timepointtimepoint2:siteManava                    397.000000  -3.294 0.001078
    ## timepointtimepoint3:siteManava                    397.000000  -0.058 0.953751
    ## timepointtimepoint4:siteManava                    397.000000  -2.482 0.013493
    ## timepointtimepoint2:speciesPocillopora            397.000000  -0.465 0.642364
    ## timepointtimepoint3:speciesPocillopora            397.000000   0.579 0.562849
    ## timepointtimepoint4:speciesPocillopora            397.000000  -0.354 0.723434
    ## timepointtimepoint2:speciesPorites                397.000000  -1.135 0.257205
    ## timepointtimepoint3:speciesPorites                397.000000  -0.131 0.895746
    ## timepointtimepoint4:speciesPorites                397.000000  -0.009 0.993213
    ## siteMahana:speciesPocillopora                     397.000000  -0.992 0.321767
    ## siteManava:speciesPocillopora                     397.000000   0.120 0.904629
    ## siteMahana:speciesPorites                         397.000000  -2.475 0.013734
    ## siteManava:speciesPorites                         397.000000   0.727 0.467833
    ## timepointtimepoint2:siteMahana:speciesPocillopora 397.000000  -0.262 0.793762
    ## timepointtimepoint3:siteMahana:speciesPocillopora 397.000000  -0.471 0.637961
    ## timepointtimepoint4:siteMahana:speciesPocillopora 397.000000   1.437 0.151640
    ## timepointtimepoint2:siteManava:speciesPocillopora 397.000000   0.286 0.775276
    ## timepointtimepoint3:siteManava:speciesPocillopora 397.000000  -0.803 0.422404
    ## timepointtimepoint4:siteManava:speciesPocillopora 397.000000   0.628 0.530376
    ## timepointtimepoint2:siteMahana:speciesPorites     397.000000   0.716 0.474668
    ## timepointtimepoint3:siteMahana:speciesPorites     397.000000   0.723 0.469978
    ## timepointtimepoint4:siteMahana:speciesPorites     397.000000   1.663 0.097074
    ## timepointtimepoint2:siteManava:speciesPorites     397.000000   0.474 0.635481
    ## timepointtimepoint3:siteManava:speciesPorites     397.000000  -1.196 0.232532
    ## timepointtimepoint4:siteManava:speciesPorites     397.000000   0.229 0.819284
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                               *  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## siteMahana                                        ***
    ## siteManava                                        ***
    ## speciesPocillopora                                *  
    ## speciesPorites                                       
    ## timepointtimepoint2:siteMahana                    ***
    ## timepointtimepoint3:siteMahana                    *  
    ## timepointtimepoint4:siteMahana                    ***
    ## timepointtimepoint2:siteManava                    ** 
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                    *  
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## siteMahana:speciesPocillopora                        
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                         *  
    ## siteManava:speciesPorites                            
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora    
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora    
    ## timepointtimepoint3:siteManava:speciesPocillopora    
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites        
    ## timepointtimepoint3:siteMahana:speciesPorites        
    ## timepointtimepoint4:siteMahana:speciesPorites     .  
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites        
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

Since species is significant, test for effects of time and site within
species.

*Acropora*

`chltotal_model_acr<-lmer(log(1+Total_Chl)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0 & species=="Acropora"))`

``` r
chltotal_model_acr<-lmer(log(1+Total_Chl)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0 & species=="Acropora"))
qqPlot(residuals(chltotal_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-92-1.png)<!-- -->

    ## 218 421 
    ##  63 111

Generate a Type III Anova of model.

``` r
anova(chltotal_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)    
    ## timepoint      4.5957 1.53190     3   101 12.8885 3.39e-07 ***
    ## site           0.8249 0.41244     2   101  3.4700 0.034875 *  
    ## timepoint:site 2.7880 0.46466     6   101  3.9094 0.001491 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotal_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: 97.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6404 -0.5666 -0.1683  0.4349  3.5370 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.0000   0.0000  
    ##  Residual                   0.1189   0.3448  
    ## Number of obs: 113, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      1.09353    0.09952 101.00000  10.988  < 2e-16
    ## timepointtimepoint2              0.27184    0.15202 101.00000   1.788 0.076757
    ## timepointtimepoint3              0.20907    0.26331 101.00000   0.794 0.429064
    ## timepointtimepoint4             -0.14179    0.17238 101.00000  -0.823 0.412701
    ## siteMahana                       0.74375    0.13801 101.00000   5.389 4.65e-07
    ## siteManava                       0.46341    0.13352 101.00000   3.471 0.000765
    ## timepointtimepoint2:siteMahana  -0.68554    0.22817 101.00000  -3.004 0.003355
    ## timepointtimepoint3:siteMahana  -0.62010    0.30060 101.00000  -2.063 0.041691
    ## timepointtimepoint4:siteMahana  -0.64565    0.22082 101.00000  -2.924 0.004268
    ## timepointtimepoint2:siteManava  -0.60140    0.21422 101.00000  -2.807 0.005996
    ## timepointtimepoint3:siteManava  -0.01488    0.30077 101.00000  -0.049 0.960647
    ## timepointtimepoint4:siteManava  -0.46555    0.22010 101.00000  -2.115 0.036876
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            .  
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                     ***
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ** 
    ## timepointtimepoint3:siteMahana *  
    ## timepointtimepoint4:siteMahana ** 
    ## timepointtimepoint2:siteManava ** 
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.655                                                    
    ## tmpnttmpnt3      -0.378  0.247                                             
    ## tmpnttmpnt4      -0.577  0.378  0.218                                      
    ## siteMahana       -0.721  0.472  0.273  0.416                               
    ## siteManava       -0.745  0.488  0.282  0.430  0.537                        
    ## tmpnttmpnt2:stMh  0.436 -0.666 -0.165 -0.252 -0.605 -0.325                 
    ## tmpnttmpnt3:stMh  0.331 -0.217 -0.876 -0.191 -0.459 -0.247  0.278          
    ## tmpnttmpnt4:stMh  0.451 -0.295 -0.170 -0.781 -0.625 -0.336  0.378          
    ## tmpnttmpnt2:stMn  0.465 -0.710 -0.176 -0.268 -0.335 -0.623  0.473          
    ## tmpnttmpnt3:stMn  0.331 -0.217 -0.875 -0.191 -0.239 -0.444  0.144          
    ## tmpnttmpnt4:stMn  0.452 -0.296 -0.171 -0.783 -0.326 -0.607  0.197          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.287                                            
    ## tmpnttmpnt2:stMn  0.154            0.209                           
    ## tmpnttmpnt3:stMn  0.767            0.149            0.277          
    ## tmpnttmpnt4:stMn  0.150            0.611            0.378          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.269          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

*Pocillopora*

`chltotal_model_poc<-lmer(log(1+Total_Chl)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0 & species=="Pocillopora"))`

``` r
chltotal_model_poc<-lmer(log(1+Total_Chl)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0 & species=="Pocillopora"))
qqPlot(residuals(chltotal_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-94-1.png)<!-- -->

    ##  50 173 
    ##  10  76

Generate a Type III Anova of model.

``` r
anova(chltotal_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      4.8368 1.61226     3   149 22.9655 2.815e-12 ***
    ## site           1.1714 0.58568     2   149  8.3427 0.0003679 ***
    ## timepoint:site 3.2347 0.53912     6   149  7.6794 3.361e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotal_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: 58.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7072 -0.5111 -0.0675  0.4762  4.1418 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.0000   0.000   
    ##  Residual                   0.0702   0.265   
    ## Number of obs: 161, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.85853    0.07349 149.00000  11.683  < 2e-16
    ## timepointtimepoint2              0.19188    0.10205 149.00000   1.880 0.062038
    ## timepointtimepoint3              0.35517    0.10393 149.00000   3.418 0.000815
    ## timepointtimepoint4             -0.20707    0.10040 149.00000  -2.062 0.040900
    ## siteMahana                       0.58037    0.10393 149.00000   5.584 1.08e-07
    ## siteManava                       0.48250    0.10040 149.00000   4.806 3.73e-06
    ## timepointtimepoint2:siteMahana  -0.75122    0.14316 149.00000  -5.247 5.21e-07
    ## timepointtimepoint3:siteMahana  -0.76311    0.14697 149.00000  -5.192 6.71e-07
    ## timepointtimepoint4:siteMahana  -0.28893    0.14605 149.00000  -1.978 0.049742
    ## timepointtimepoint2:siteManava  -0.53255    0.14181 149.00000  -3.755 0.000248
    ## timepointtimepoint3:siteManava  -0.25928    0.14786 149.00000  -1.754 0.081570
    ## timepointtimepoint4:siteManava  -0.31173    0.14199 149.00000  -2.195 0.029677
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            .  
    ## timepointtimepoint3            ***
    ## timepointtimepoint4            *  
    ## siteMahana                     ***
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ***
    ## timepointtimepoint3:siteMahana ***
    ## timepointtimepoint4:siteMahana *  
    ## timepointtimepoint2:siteManava ***
    ## timepointtimepoint3:siteManava .  
    ## timepointtimepoint4:siteManava *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.720                                                    
    ## tmpnttmpnt3      -0.707  0.509                                             
    ## tmpnttmpnt4      -0.732  0.527  0.518                                      
    ## siteMahana       -0.707  0.509  0.500  0.518                               
    ## siteManava       -0.732  0.527  0.518  0.536  0.518                        
    ## tmpnttmpnt2:stMh  0.513 -0.713 -0.363 -0.376 -0.726 -0.376                 
    ## tmpnttmpnt3:stMh  0.500 -0.360 -0.707 -0.366 -0.707 -0.366  0.513          
    ## tmpnttmpnt4:stMh  0.503 -0.362 -0.356 -0.687 -0.712 -0.368  0.517          
    ## tmpnttmpnt2:stMn  0.518 -0.720 -0.366 -0.379 -0.366 -0.708  0.513          
    ## tmpnttmpnt3:stMn  0.497 -0.358 -0.703 -0.364 -0.351 -0.679  0.255          
    ## tmpnttmpnt4:stMn  0.518 -0.373 -0.366 -0.707 -0.366 -0.707  0.266          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.503                                            
    ## tmpnttmpnt2:stMn  0.259            0.261                           
    ## tmpnttmpnt3:stMn  0.497            0.250            0.481          
    ## tmpnttmpnt4:stMn  0.259            0.486            0.501          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.480          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

*Porites*

`chltotal_model_poc<-lmer(log(1+Total_Chl)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0 & species=="Porites"))`

``` r
chltotal_model_por<-lmer(log(1+Total_Chl)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl>0 & species=="Porites"))
qqPlot(residuals(chltotal_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-96-1.png)<!-- -->

    ## 109  88 
    ##  28   8

Generate a Type III Anova of model.

``` r
anova(chltotal_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      3.0968 1.03226     3   147 12.8422 1.687e-07 ***
    ## site           2.0147 1.00734     2   147 12.5322 9.431e-06 ***
    ## timepoint:site 1.4164 0.23607     6   147  2.9369  0.009788 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotal_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: 77.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5498 -0.4355 -0.0216  0.3014  7.0466 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.08038  0.2835  
    ## Number of obs: 159, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.91073    0.07577 147.00000  12.019  < 2e-16
    ## timepointtimepoint2              0.07480    0.11153 147.00000   0.671  0.50347
    ## timepointtimepoint3              0.17634    0.10536 147.00000   1.674  0.09631
    ## timepointtimepoint4             -0.14337    0.10920 147.00000  -1.313  0.19126
    ## siteMahana                       0.33970    0.10920 147.00000   3.111  0.00224
    ## siteManava                       0.57803    0.10536 147.00000   5.486 1.75e-07
    ## timepointtimepoint2:siteMahana  -0.50510    0.15486 147.00000  -3.262  0.00138
    ## timepointtimepoint3:siteMahana  -0.40067    0.15681 147.00000  -2.555  0.01163
    ## timepointtimepoint4:siteMahana  -0.23128    0.15750 147.00000  -1.468  0.14411
    ## timepointtimepoint2:siteManava  -0.48653    0.15343 147.00000  -3.171  0.00185
    ## timepointtimepoint3:siteManava  -0.37367    0.15047 147.00000  -2.483  0.01414
    ## timepointtimepoint4:siteManava  -0.40910    0.15486 147.00000  -2.642  0.00914
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3            .  
    ## timepointtimepoint4               
    ## siteMahana                     ** 
    ## siteManava                     ***
    ## timepointtimepoint2:siteMahana ** 
    ## timepointtimepoint3:siteMahana *  
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava ** 
    ## timepointtimepoint3:siteManava *  
    ## timepointtimepoint4:siteManava ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.679                                                    
    ## tmpnttmpnt3      -0.719  0.489                                             
    ## tmpnttmpnt4      -0.694  0.471  0.499                                      
    ## siteMahana       -0.694  0.471  0.499  0.481                               
    ## siteManava       -0.719  0.489  0.517  0.499  0.499                        
    ## tmpnttmpnt2:stMh  0.489 -0.720 -0.352 -0.340 -0.705 -0.352                 
    ## tmpnttmpnt3:stMh  0.483 -0.328 -0.672 -0.335 -0.696 -0.348  0.491          
    ## tmpnttmpnt4:stMh  0.481 -0.327 -0.346 -0.693 -0.693 -0.346  0.489          
    ## tmpnttmpnt2:stMn  0.494 -0.727 -0.355 -0.343 -0.343 -0.687  0.524          
    ## tmpnttmpnt3:stMn  0.504 -0.342 -0.700 -0.349 -0.349 -0.700  0.246          
    ## tmpnttmpnt4:stMn  0.489 -0.332 -0.352 -0.705 -0.340 -0.680  0.239          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.483                                            
    ## tmpnttmpnt2:stMn  0.239            0.238                           
    ## tmpnttmpnt3:stMn  0.470            0.242            0.481          
    ## tmpnttmpnt4:stMn  0.236            0.489            0.467          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.476          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

#### Total Chlorophyll per cell (a + c2):

Build a mixed model for univariate analysis and examine data
distribution.

`chltotalcell_model<-lmer(Total_Chl_cell~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0))`

``` r
chltotalcell_model<-lmer(Total_Chl_cell~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0))
qqPlot(residuals(chltotalcell_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-98-1.png)<!-- -->

    ## 49 76 
    ## 48 72

Residuals are not normally distributed. Attempt with log transformation.

`chltotalcell_model<-lmer(log(1+Total_Chl_cell)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0))`

``` r
chltotalcell_model<-lmer(log(1+Total_Chl_cell)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0))
qqPlot(residuals(chltotalcell_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-99-1.png)<!-- -->

    ## 49 76 
    ## 48 72

Return to check distribution of residuals.

Generate a Type III Anova of model.

``` r
anova(chltotalcell_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                            Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint              1.6493e-09 5.4976e-10     3 35.561 13.3411 5.444e-06 ***
    ## site                   7.2250e-10 3.6124e-10     2 36.512  8.7663 0.0007775 ***
    ## species                5.8367e-09 2.9183e-09     2 36.152 70.8199 3.128e-13 ***
    ## timepoint:site         6.6950e-10 1.1159e-10     6 35.502  2.7079 0.0285950 *  
    ## timepoint:species      5.0250e-10 8.3750e-11     6 35.411  2.0324 0.0869746 .  
    ## site:species           1.2550e-10 3.1370e-11     4 36.045  0.7613 0.5573710    
    ## timepoint:site:species 3.1810e-10 2.6510e-11    12 35.360  0.6433 0.7910120    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotalcell_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## log(1 + Total_Chl_cell) ~ timepoint * site * species + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl_cell > 0)
    ## 
    ## REML criterion at convergence: -8399.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0963 -0.4586 -0.0885  0.2482  7.7372 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev. 
    ##  colony_id_corr (Intercept) 2.990e-12 1.729e-06
    ##  Residual                   4.121e-11 6.419e-06
    ## Number of obs: 440, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                        3.786e-06  1.918e-06
    ## timepointtimepoint2                                1.236e-07  2.872e-06
    ## timepointtimepoint3                                7.973e-06  5.037e-06
    ## timepointtimepoint4                                4.020e-06  3.245e-06
    ## siteMahana                                         4.555e-06  2.660e-06
    ## siteManava                                         3.494e-06  2.574e-06
    ## speciesPocillopora                                 5.311e-06  2.712e-06
    ## speciesPorites                                     1.197e-06  2.574e-06
    ## timepointtimepoint2:siteMahana                    -4.029e-06  4.304e-06
    ## timepointtimepoint3:siteMahana                    -2.765e-06  5.739e-06
    ## timepointtimepoint4:siteMahana                    -2.082e-06  4.167e-06
    ## timepointtimepoint2:siteManava                    -3.574e-06  3.966e-06
    ## timepointtimepoint3:siteManava                    -3.209e-06  5.732e-06
    ## timepointtimepoint4:siteManava                    -4.338e-06  4.137e-06
    ## timepointtimepoint2:speciesPocillopora             5.365e-06  3.832e-06
    ## timepointtimepoint3:speciesPocillopora            -7.796e-07  5.659e-06
    ## timepointtimepoint4:speciesPocillopora             3.146e-06  4.093e-06
    ## timepointtimepoint2:speciesPorites                 4.370e-07  3.735e-06
    ## timepointtimepoint3:speciesPorites                -4.832e-06  5.556e-06
    ## timepointtimepoint4:speciesPorites                 1.587e-06  4.029e-06
    ## siteMahana:speciesPocillopora                      4.329e-06  3.798e-06
    ## siteManava:speciesPocillopora                      3.620e-06  3.640e-06
    ## siteMahana:speciesPorites                          3.734e-06  3.601e-06
    ## siteManava:speciesPorites                          1.857e-06  3.538e-06
    ## timepointtimepoint2:siteMahana:speciesPocillopora -5.933e-06  5.584e-06
    ## timepointtimepoint3:siteMahana:speciesPocillopora -6.570e-07  6.800e-06
    ## timepointtimepoint4:siteMahana:speciesPocillopora  2.662e-06  5.552e-06
    ## timepointtimepoint2:siteManava:speciesPocillopora -1.669e-06  5.259e-06
    ## timepointtimepoint3:siteManava:speciesPocillopora -1.085e-06  6.787e-06
    ## timepointtimepoint4:siteManava:speciesPocillopora -1.005e-06  5.411e-06
    ## timepointtimepoint2:siteMahana:speciesPorites     -3.346e-06  5.452e-06
    ## timepointtimepoint3:siteMahana:speciesPorites     -9.254e-07  6.683e-06
    ## timepointtimepoint4:siteMahana:speciesPorites     -6.177e-06  5.412e-06
    ## timepointtimepoint2:siteManava:speciesPorites     -1.909e-06  5.189e-06
    ## timepointtimepoint3:siteManava:speciesPorites     -1.855e-06  6.655e-06
    ## timepointtimepoint4:siteManava:speciesPorites     -3.246e-06  5.389e-06
    ##                                                           df t value Pr(>|t|)  
    ## (Intercept)                                        3.775e+01   1.974   0.0558 .
    ## timepointtimepoint2                                3.617e+01   0.043   0.9659  
    ## timepointtimepoint3                                3.709e+01   1.583   0.1219  
    ## timepointtimepoint4                                3.578e+01   1.239   0.2235  
    ## siteMahana                                         3.772e+01   1.713   0.0950 .
    ## siteManava                                         3.775e+01   1.358   0.1827  
    ## speciesPocillopora                                 3.771e+01   1.958   0.0576 .
    ## speciesPorites                                     3.777e+01   0.465   0.6446  
    ## timepointtimepoint2:siteMahana                     3.600e+01  -0.936   0.3554  
    ## timepointtimepoint3:siteMahana                     3.697e+01  -0.482   0.6328  
    ## timepointtimepoint4:siteMahana                     3.604e+01  -0.500   0.6204  
    ## timepointtimepoint2:siteManava                     3.596e+01  -0.901   0.3735  
    ## timepointtimepoint3:siteManava                     3.682e+01  -0.560   0.5790  
    ## timepointtimepoint4:siteManava                     3.563e+01  -1.048   0.3015  
    ## timepointtimepoint2:speciesPocillopora             3.573e+01   1.400   0.1701  
    ## timepointtimepoint3:speciesPocillopora             3.673e+01  -0.138   0.8912  
    ## timepointtimepoint4:speciesPocillopora             3.549e+01   0.769   0.4472  
    ## timepointtimepoint2:speciesPorites                 3.562e+01   0.117   0.9075  
    ## timepointtimepoint3:speciesPorites                 3.673e+01  -0.870   0.3900  
    ## timepointtimepoint4:speciesPorites                 3.543e+01   0.394   0.6961  
    ## siteMahana:speciesPocillopora                      3.771e+01   1.140   0.2616  
    ## siteManava:speciesPocillopora                      3.774e+01   0.995   0.3262  
    ## siteMahana:speciesPorites                          3.776e+01   1.037   0.3064  
    ## siteManava:speciesPorites                          3.777e+01   0.525   0.6027  
    ## timepointtimepoint2:siteMahana:speciesPocillopora  3.564e+01  -1.062   0.2952  
    ## timepointtimepoint3:siteMahana:speciesPocillopora  3.648e+01  -0.097   0.9236  
    ## timepointtimepoint4:siteMahana:speciesPocillopora  3.562e+01   0.479   0.6346  
    ## timepointtimepoint2:siteManava:speciesPocillopora  3.553e+01  -0.317   0.7528  
    ## timepointtimepoint3:siteManava:speciesPocillopora  3.637e+01  -0.160   0.8738  
    ## timepointtimepoint4:siteManava:speciesPocillopora  3.535e+01  -0.186   0.8538  
    ## timepointtimepoint2:siteMahana:speciesPorites      3.554e+01  -0.614   0.5433  
    ## timepointtimepoint3:siteMahana:speciesPorites      3.649e+01  -0.138   0.8906  
    ## timepointtimepoint4:siteMahana:speciesPorites      3.560e+01  -1.141   0.2614  
    ## timepointtimepoint2:siteManava:speciesPorites      3.546e+01  -0.368   0.7152  
    ## timepointtimepoint3:siteManava:speciesPorites      3.636e+01  -0.279   0.7820  
    ## timepointtimepoint4:siteManava:speciesPorites      3.533e+01  -0.602   0.5508  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of time and set within
species.

*Acropora*

`chltotalcell_model_acr<-lmer(log(1+Total_Chl_cell)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0 & species=="Acropora"))`

``` r
chltotalcell_model_acr<-lmer(log(1+Total_Chl_cell)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0 & species=="Acropora"))
qqPlot(residuals(chltotalcell_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-101-1.png)<!-- -->

    ##  16 218 
    ##  16  64

Generate a Type III Anova of model.

``` r
anova(chltotalcell_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                    Sum Sq    Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## timepoint      6.0291e-10 2.0097e-10     3 7.8637 14.9640 0.001289 **
    ## site           6.7710e-11 3.3854e-11     2 6.6787  2.5207 0.152900   
    ## timepoint:site 7.7480e-11 1.2913e-11     6 7.8664  0.9615 0.505753   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotalcell_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl_cell) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl_cell > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: -2228.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7031 -0.4517 -0.0883  0.3702  3.5747 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev. 
    ##  colony_id_corr (Intercept) 1.417e-12 1.190e-06
    ##  Residual                   1.343e-11 3.665e-06
    ## Number of obs: 114, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     3.794e-06  1.111e-06  7.639e+00   3.414
    ## timepointtimepoint2             1.217e-07  1.649e-06  7.865e+00   0.074
    ## timepointtimepoint3             7.970e-06  2.905e-06  7.795e+00   2.744
    ## timepointtimepoint4             4.050e-06  1.860e-06  7.881e+00   2.177
    ## siteMahana                      4.565e-06  1.540e-06  7.648e+00   2.963
    ## siteManava                      3.482e-06  1.491e-06  7.637e+00   2.335
    ## timepointtimepoint2:siteMahana -3.990e-06  2.470e-06  7.874e+00  -1.616
    ## timepointtimepoint3:siteMahana -2.783e-06  3.308e-06  7.807e+00  -0.841
    ## timepointtimepoint4:siteMahana -2.140e-06  2.392e-06  7.872e+00  -0.895
    ## timepointtimepoint2:siteManava -3.532e-06  2.276e-06  7.875e+00  -1.552
    ## timepointtimepoint3:siteManava -3.176e-06  3.301e-06  7.824e+00  -0.962
    ## timepointtimepoint4:siteManava -4.353e-06  2.370e-06  7.884e+00  -1.836
    ##                                Pr(>|t|)   
    ## (Intercept)                     0.00982 **
    ## timepointtimepoint2             0.94301   
    ## timepointtimepoint3             0.02593 * 
    ## timepointtimepoint4             0.06163 . 
    ## siteMahana                      0.01898 * 
    ## siteManava                      0.04925 * 
    ## timepointtimepoint2:siteMahana  0.14545   
    ## timepointtimepoint3:siteMahana  0.42522   
    ## timepointtimepoint4:siteMahana  0.39749   
    ## timepointtimepoint2:siteManava  0.15986   
    ## timepointtimepoint3:siteManava  0.36474   
    ## timepointtimepoint4:siteManava  0.10419   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.631                                                    
    ## tmpnttmpnt3      -0.364  0.260                                             
    ## tmpnttmpnt4      -0.542  0.373  0.217                                      
    ## siteMahana       -0.721  0.455  0.262  0.391                               
    ## siteManava       -0.745  0.470  0.271  0.404  0.538                        
    ## tmpnttmpnt2:stMh  0.421 -0.668 -0.174 -0.249 -0.576 -0.314                 
    ## tmpnttmpnt3:stMh  0.319 -0.228 -0.878 -0.190 -0.443 -0.238  0.286          
    ## tmpnttmpnt4:stMh  0.421 -0.290 -0.169 -0.778 -0.598 -0.314  0.373          
    ## tmpnttmpnt2:stMn  0.457 -0.725 -0.189 -0.270 -0.330 -0.607  0.484          
    ## tmpnttmpnt3:stMn  0.320 -0.229 -0.880 -0.191 -0.231 -0.422  0.153          
    ## tmpnttmpnt4:stMn  0.425 -0.293 -0.170 -0.785 -0.307 -0.572  0.195          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.291                                            
    ## tmpnttmpnt2:stMn  0.166            0.210                           
    ## tmpnttmpnt3:stMn  0.773            0.148            0.293          
    ## tmpnttmpnt4:stMn  0.149            0.610            0.380          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.267

*Pocillopora*

`chltotalcell_model_poc<-lmer(log(1+Total_Chl_cell)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0 & species=="Pocillopora"))`

``` r
chltotalcell_model_poc<-lmer(log(1+Total_Chl_cell)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0 & species=="Pocillopora"))
qqPlot(residuals(chltotalcell_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-103-1.png)<!-- -->

    ## 49 76 
    ##  9 33

Generate a Type III Anova of model.

``` r
anova(chltotalcell_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                    Sum Sq    Mean Sq NumDF  DenDF F value  Pr(>F)  
    ## timepoint      9.3484e-10 3.1161e-10     3 20.299  4.0484 0.02091 *
    ## site           6.1445e-10 3.0722e-10     2 33.008  3.9914 0.02802 *
    ## timepoint:site 6.5307e-10 1.0884e-10     6 20.278  1.4141 0.25734  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotalcell_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl_cell) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl_cell > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -2961.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5302 -0.5738 -0.1957  0.3135  5.4996 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev. 
    ##  colony_id_corr (Intercept) 8.554e-12 2.925e-06
    ##  Residual                   7.697e-11 8.773e-06
    ## Number of obs: 159, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     9.077e-06  2.664e-06  4.928e+01   3.408
    ## timepointtimepoint2             5.516e-06  3.474e-06  2.090e+01   1.588
    ## timepointtimepoint3             7.188e-06  3.530e-06  2.041e+01   2.036
    ## timepointtimepoint4             7.186e-06  3.414e-06  2.036e+01   2.105
    ## siteMahana                      9.008e-06  3.768e-06  4.965e+01   2.391
    ## siteManava                      7.135e-06  3.577e-06  5.113e+01   1.995
    ## timepointtimepoint2:siteMahana -1.009e-05  4.872e-06  2.072e+01  -2.072
    ## timepointtimepoint3:siteMahana -3.535e-06  4.993e-06  2.050e+01  -0.708
    ## timepointtimepoint4:siteMahana  4.954e-07  5.022e-06  2.051e+01   0.099
    ## timepointtimepoint2:siteManava -5.271e-06  4.726e-06  2.002e+01  -1.115
    ## timepointtimepoint3:siteManava -4.414e-06  4.975e-06  2.067e+01  -0.887
    ## timepointtimepoint4:siteManava -5.405e-06  4.773e-06  2.017e+01  -1.132
    ##                                Pr(>|t|)   
    ## (Intercept)                     0.00131 **
    ## timepointtimepoint2             0.12735   
    ## timepointtimepoint3             0.05490 . 
    ## timepointtimepoint4             0.04792 * 
    ## siteMahana                      0.02064 * 
    ## siteManava                      0.05143 . 
    ## timepointtimepoint2:siteMahana  0.05096 . 
    ## timepointtimepoint3:siteMahana  0.48685   
    ## timepointtimepoint4:siteMahana  0.92238   
    ## timepointtimepoint2:siteManava  0.27791   
    ## timepointtimepoint3:siteManava  0.38518   
    ## timepointtimepoint4:siteManava  0.27071   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.706                                                    
    ## tmpnttmpnt3      -0.691  0.530                                             
    ## tmpnttmpnt4      -0.717  0.551  0.539                                      
    ## siteMahana       -0.707  0.499  0.489  0.507                               
    ## siteManava       -0.745  0.526  0.515  0.534  0.526                        
    ## tmpnttmpnt2:stMh  0.503 -0.713 -0.378 -0.393 -0.712 -0.375                 
    ## tmpnttmpnt3:stMh  0.489 -0.375 -0.707 -0.381 -0.691 -0.364  0.535          
    ## tmpnttmpnt4:stMh  0.488 -0.374 -0.367 -0.680 -0.686 -0.363  0.531          
    ## tmpnttmpnt2:stMn  0.519 -0.735 -0.390 -0.405 -0.367 -0.690  0.524          
    ## tmpnttmpnt3:stMn  0.490 -0.376 -0.710 -0.383 -0.347 -0.653  0.268          
    ## tmpnttmpnt4:stMn  0.513 -0.394 -0.386 -0.715 -0.363 -0.683  0.281          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.519                                            
    ## tmpnttmpnt2:stMn  0.276            0.275                           
    ## tmpnttmpnt3:stMn  0.502            0.260            0.495          
    ## tmpnttmpnt4:stMn  0.273            0.486            0.517          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.493

*Porites*

`chltotalcell_model_por<-lmer(log(1+Total_Chl_cell)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0 & species=="Porites"))`

``` r
chltotalcell_model_por<-lmer(log(1+Total_Chl_cell)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(Total_Chl_cell>0 & species=="Porites"))
qqPlot(residuals(chltotalcell_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-105-1.png)<!-- -->

    ## 109  87 
    ##  31   9

Generate a Type III Anova of model.

``` r
anova(chltotalcell_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                    Sum Sq    Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## timepoint      4.8104e-10 1.6035e-10     3 13.298  6.5781 0.005835 **
    ## site           3.5270e-10 1.7635e-10     2 13.298  7.2345 0.007482 **
    ## timepoint:site 3.8719e-10 6.4531e-11     6 13.298  2.6473 0.065342 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(chltotalcell_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + Total_Chl_cell) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(Total_Chl_cell > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: -3316.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0934 -0.3865 -0.0978  0.2288  5.7417 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev. 
    ##  colony_id_corr (Intercept) 0.000e+00 0.000e+00
    ##  Residual                   2.438e-11 4.937e-06
    ## Number of obs: 167, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     4.983e-06  1.275e-06  1.330e+01   3.909
    ## timepointtimepoint2             5.715e-07  1.835e-06  1.330e+01   0.311
    ## timepointtimepoint3             3.141e-06  1.803e-06  1.330e+01   1.742
    ## timepointtimepoint4             5.614e-06  1.835e-06  1.330e+01   3.060
    ## siteMahana                      8.289e-06  1.803e-06  1.330e+01   4.598
    ## siteManava                      5.351e-06  1.803e-06  1.330e+01   2.968
    ## timepointtimepoint2:siteMahana -7.386e-06  2.572e-06  1.330e+01  -2.871
    ## timepointtimepoint3:siteMahana -3.768e-06  2.628e-06  1.330e+01  -1.434
    ## timepointtimepoint4:siteMahana -8.344e-06  2.650e-06  1.330e+01  -3.149
    ## timepointtimepoint2:siteManava -5.494e-06  2.572e-06  1.330e+01  -2.136
    ## timepointtimepoint3:siteManava -5.187e-06  2.598e-06  1.330e+01  -1.997
    ## timepointtimepoint4:siteManava -7.733e-06  2.650e-06  1.330e+01  -2.918
    ##                                Pr(>|t|)    
    ## (Intercept)                    0.001723 ** 
    ## timepointtimepoint2            0.760253    
    ## timepointtimepoint3            0.104549    
    ## timepointtimepoint4            0.008919 ** 
    ## siteMahana                     0.000472 ***
    ## siteManava                     0.010659 *  
    ## timepointtimepoint2:siteMahana 0.012859 *  
    ## timepointtimepoint3:siteMahana 0.174750    
    ## timepointtimepoint4:siteMahana 0.007506 ** 
    ## timepointtimepoint2:siteManava 0.051825 .  
    ## timepointtimepoint3:siteManava 0.066759 .  
    ## timepointtimepoint4:siteManava 0.011746 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.695                                                    
    ## tmpnttmpnt3      -0.707  0.491                                             
    ## tmpnttmpnt4      -0.695  0.483  0.491                                      
    ## siteMahana       -0.707  0.491  0.500  0.491                               
    ## siteManava       -0.707  0.491  0.500  0.491  0.500                        
    ## tmpnttmpnt2:stMh  0.496 -0.713 -0.350 -0.344 -0.701 -0.350                 
    ## tmpnttmpnt3:stMh  0.485 -0.337 -0.686 -0.337 -0.686 -0.343  0.481          
    ## tmpnttmpnt4:stMh  0.481 -0.334 -0.340 -0.692 -0.680 -0.340  0.477          
    ## tmpnttmpnt2:stMn  0.496 -0.713 -0.350 -0.344 -0.350 -0.701  0.509          
    ## tmpnttmpnt3:stMn  0.491 -0.341 -0.694 -0.341 -0.347 -0.694  0.243          
    ## tmpnttmpnt4:stMn  0.481 -0.334 -0.340 -0.692 -0.340 -0.680  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.467                                            
    ## tmpnttmpnt2:stMn  0.240            0.238                           
    ## tmpnttmpnt3:stMn  0.476            0.236            0.486          
    ## tmpnttmpnt4:stMn  0.233            0.479            0.477          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.472          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

## Protein

### Plot: Surface Area

View by site and timepoint for each species.

``` r
protplot_cm2<-master %>%
  filter(!is.na(prot_mg.cm2)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, prot_mg.cm2)%>%
  
  ggplot(., aes(x = site_code, y = prot_mg.cm2, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Protein (mg cm"^-2, ")"))))+
    theme_classic() + 
    theme(
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: AFDW

View by site and timepoint for each species.

``` r
protplot_afdw<-master %>%
  filter(!is.na(prot_mg.mgafdw)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, prot_mg.mgafdw)%>%
  
  ggplot(., aes(x = site_code, y = prot_mg.mgafdw, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Protein (mg mg afdw"^-1, ")"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots for each normalization together.

``` r
prot_figure<-plot_grid(protplot_afdw, protplot_cm2, labels = c("", "", ""), ncol=2, nrow=1, rel_widths = c(.8, 1), label_size = 20, label_x = 0.1)

ggsave(filename="Figures/Protein_Figure.pdf", plot=prot_figure, dpi=300, width=16, height=5, units="in")
```

### Analysis

Build a mixed model for univariate analysis and examine data
distribution.

`prot_model<-lmer(prot_mg.mgafdw~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_ug.mgAFDW>0))`

``` r
prot_model<-lmer(prot_mg.mgafdw~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_mg.mgafdw>0))
qqPlot(residuals(prot_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-110-1.png)<!-- -->

    ## 218 109 
    ## 206 101

Residuals are not normally distributed. Attempt with log transformation.

`prot_model<-lmer(log(1+prot_mg.mgafdw)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_ug.mgafdw>0))`

``` r
prot_model<-lmer(log(1+prot_mg.mgafdw)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_mg.mgafdw>0))
qqPlot(residuals(prot_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-111-1.png)<!-- -->

    ## 218 109 
    ## 206 101

Generate a Type III Anova of model.

``` r
anova(prot_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint              0.86170 0.287232     3 299.62 55.2077 < 2.2e-16 ***
    ## species                0.45716 0.228580     2 123.43 43.9345 3.899e-15 ***
    ## site                   0.04766 0.023829     2 128.23  4.5800 0.0119887 *  
    ## timepoint:species      0.14482 0.024137     6 293.12  4.6393 0.0001598 ***
    ## timepoint:site         0.16951 0.028252     6 297.78  5.4302 2.400e-05 ***
    ## species:site           0.20114 0.050284     4 121.80  9.6649 8.023e-07 ***
    ## timepoint:species:site 0.13994 0.011662    12 291.77  2.2414 0.0101727 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(prot_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## log(1 + prot_mg.mgafdw) ~ timepoint * species * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(prot_mg.mgafdw > 0)
    ## 
    ## REML criterion at convergence: -780.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8872 -0.4131 -0.1314  0.3402  7.4991 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001713 0.04139 
    ##  Residual                   0.005203 0.07213 
    ## Number of obs: 431, groups:  colony_id_corr, 139
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         0.378361   0.023910
    ## timepointtimepoint2                                -0.080540   0.033118
    ## timepointtimepoint3                                -0.012647   0.058662
    ## timepointtimepoint4                                -0.148638   0.037327
    ## speciesPocillopora                                 -0.136466   0.032137
    ## speciesPorites                                      0.046924   0.033676
    ## siteMahana                                          0.041729   0.033043
    ## siteManava                                         -0.055531   0.032083
    ## timepointtimepoint2:speciesPocillopora              0.111798   0.042654
    ## timepointtimepoint3:speciesPocillopora              0.068907   0.064786
    ## timepointtimepoint4:speciesPocillopora              0.083209   0.045684
    ## timepointtimepoint2:speciesPorites                 -0.009535   0.044737
    ## timepointtimepoint3:speciesPorites                 -0.008438   0.065087
    ## timepointtimepoint4:speciesPorites                 -0.023755   0.047494
    ## timepointtimepoint2:siteMahana                      0.038259   0.049657
    ## timepointtimepoint3:siteMahana                      0.028303   0.067019
    ## timepointtimepoint4:siteMahana                      0.011008   0.048356
    ## timepointtimepoint2:siteManava                      0.041198   0.046615
    ## timepointtimepoint3:siteManava                      0.244914   0.066561
    ## timepointtimepoint4:siteManava                      0.039256   0.047457
    ## speciesPocillopora:siteMahana                       0.014322   0.046012
    ## speciesPorites:siteMahana                          -0.129036   0.046686
    ## speciesPocillopora:siteManava                       0.078268   0.044501
    ## speciesPorites:siteManava                           0.123705   0.045308
    ## timepointtimepoint2:speciesPocillopora:siteMahana  -0.095081   0.063129
    ## timepointtimepoint3:speciesPocillopora:siteMahana  -0.090526   0.078100
    ## timepointtimepoint4:speciesPocillopora:siteMahana  -0.002542   0.062651
    ## timepointtimepoint2:speciesPorites:siteMahana       0.017083   0.064253
    ## timepointtimepoint3:speciesPorites:siteMahana      -0.052634   0.078657
    ## timepointtimepoint4:speciesPorites:siteMahana       0.116717   0.063669
    ## timepointtimepoint2:speciesPocillopora:siteManava  -0.086113   0.060396
    ## timepointtimepoint3:speciesPocillopora:siteManava  -0.259784   0.077857
    ## timepointtimepoint4:speciesPocillopora:siteManava  -0.043021   0.061108
    ## timepointtimepoint2:speciesPorites:siteManava      -0.016314   0.061645
    ## timepointtimepoint3:speciesPorites:siteManava      -0.239297   0.077351
    ## timepointtimepoint4:speciesPorites:siteManava      -0.008879   0.062538
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       352.451277  15.824  < 2e-16
    ## timepointtimepoint2                               330.275916  -2.432 0.015549
    ## timepointtimepoint3                               338.971204  -0.216 0.829441
    ## timepointtimepoint4                               312.071355  -3.982 8.51e-05
    ## speciesPocillopora                                346.967863  -4.246 2.79e-05
    ## speciesPorites                                    361.401451   1.393 0.164365
    ## siteMahana                                        361.534415   1.263 0.207445
    ## siteManava                                        351.242087  -1.731 0.084359
    ## timepointtimepoint2:speciesPocillopora            304.499061   2.621 0.009206
    ## timepointtimepoint3:speciesPocillopora            326.738750   1.064 0.288291
    ## timepointtimepoint4:speciesPocillopora            295.220101   1.821 0.069558
    ## timepointtimepoint2:speciesPorites                310.696950  -0.213 0.831352
    ## timepointtimepoint3:speciesPorites                326.954420  -0.130 0.896934
    ## timepointtimepoint4:speciesPorites                300.295749  -0.500 0.617321
    ## timepointtimepoint2:siteMahana                    324.305418   0.770 0.441585
    ## timepointtimepoint3:siteMahana                    345.507135   0.422 0.673061
    ## timepointtimepoint4:siteMahana                    334.168978   0.228 0.820066
    ## timepointtimepoint2:siteManava                    322.387148   0.884 0.377463
    ## timepointtimepoint3:siteManava                    332.559945   3.680 0.000272
    ## timepointtimepoint4:siteManava                    305.536229   0.827 0.408767
    ## speciesPocillopora:siteMahana                     358.621654   0.311 0.755787
    ## speciesPorites:siteMahana                         362.664115  -2.764 0.006002
    ## speciesPocillopora:siteManava                     348.498501   1.759 0.079492
    ## speciesPorites:siteManava                         354.265286   2.730 0.006643
    ## timepointtimepoint2:speciesPocillopora:siteMahana 303.913319  -1.506 0.133074
    ## timepointtimepoint3:speciesPocillopora:siteMahana 327.003079  -1.159 0.247257
    ## timepointtimepoint4:speciesPocillopora:siteMahana 309.321075  -0.041 0.967663
    ## timepointtimepoint2:speciesPorites:siteMahana     306.556782   0.266 0.790519
    ## timepointtimepoint3:speciesPorites:siteMahana     327.586150  -0.669 0.503867
    ## timepointtimepoint4:speciesPorites:siteMahana     311.192539   1.833 0.067729
    ## timepointtimepoint2:speciesPocillopora:siteManava 299.908883  -1.426 0.154964
    ## timepointtimepoint3:speciesPocillopora:siteManava 317.846017  -3.337 0.000948
    ## timepointtimepoint4:speciesPocillopora:siteManava 290.990908  -0.704 0.481991
    ## timepointtimepoint2:speciesPorites:siteManava     302.965177  -0.265 0.791468
    ## timepointtimepoint3:speciesPorites:siteManava     317.119355  -3.094 0.002153
    ## timepointtimepoint4:speciesPorites:siteManava     293.892068  -0.142 0.887201
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                               *  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                               ***
    ## speciesPocillopora                                ***
    ## speciesPorites                                       
    ## siteMahana                                           
    ## siteManava                                        .  
    ## timepointtimepoint2:speciesPocillopora            ** 
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora            .  
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                    ***
    ## timepointtimepoint4:siteManava                       
    ## speciesPocillopora:siteMahana                        
    ## speciesPorites:siteMahana                         ** 
    ## speciesPocillopora:siteManava                     .  
    ## speciesPorites:siteManava                         ** 
    ## timepointtimepoint2:speciesPocillopora:siteMahana    
    ## timepointtimepoint3:speciesPocillopora:siteMahana    
    ## timepointtimepoint4:speciesPocillopora:siteMahana    
    ## timepointtimepoint2:speciesPorites:siteMahana        
    ## timepointtimepoint3:speciesPorites:siteMahana        
    ## timepointtimepoint4:speciesPorites:siteMahana     .  
    ## timepointtimepoint2:speciesPocillopora:siteManava    
    ## timepointtimepoint3:speciesPocillopora:siteManava ***
    ## timepointtimepoint4:speciesPocillopora:siteManava    
    ## timepointtimepoint2:speciesPorites:siteManava        
    ## timepointtimepoint3:speciesPorites:siteManava     ** 
    ## timepointtimepoint4:speciesPorites:siteManava        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of site and time within
each species.

*Acropora*

`prot_model_acr<-lmer(log(1+prot_mg.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_ug.mgafdw>0 & species=="Acropora"))`

``` r
prot_model_acr<-lmer(log(1+prot_mg.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_mg.mgafdw>0 & species=="Acropora"))
qqPlot(residuals(prot_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-113-1.png)<!-- -->

    ## 218 215 
    ##  63  60

Generate a Type III Anova of model.

``` r
anova(prot_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      0.47259 0.157529     3 88.328 14.9515 5.956e-08 ***
    ## site           0.03382 0.016911     2 52.627  1.6050   0.21056    
    ## timepoint:site 0.17730 0.029549     6 87.181  2.8046   0.01527 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(prot_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + prot_mg.mgafdw) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(prot_mg.mgafdw > 0 & species == "Acropora")
    ## 
    ## REML criterion at convergence: -141
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3014 -0.4414 -0.0442  0.2986  5.8555 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0007401 0.0272  
    ##  Residual                   0.0105360 0.1026  
    ## Number of obs: 113, groups:  colony_id_corr, 49
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     3.777e-01  3.065e-02  1.007e+02  12.326
    ## timepointtimepoint2            -7.829e-02  4.576e-02  8.541e+01  -1.711
    ## timepointtimepoint3             8.195e-05  8.007e-02  9.463e+01   0.001
    ## timepointtimepoint4            -1.538e-01  5.186e-02  8.418e+01  -2.966
    ## siteMahana                      4.007e-02  4.249e-02  1.008e+02   0.943
    ## siteManava                     -5.551e-02  4.111e-02  1.007e+02  -1.350
    ## timepointtimepoint2:siteMahana  2.722e-02  6.870e-02  8.550e+01   0.396
    ## timepointtimepoint3:siteMahana  1.380e-02  9.135e-02  9.439e+01   0.151
    ## timepointtimepoint4:siteMahana  1.739e-02  6.660e-02  8.758e+01   0.261
    ## timepointtimepoint2:siteManava  4.971e-02  6.451e-02  8.563e+01   0.771
    ## timepointtimepoint3:siteManava  2.343e-01  9.123e-02  9.250e+01   2.568
    ## timepointtimepoint4:siteManava  4.474e-02  6.613e-02  8.262e+01   0.677
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## timepointtimepoint2             0.09069 .  
    ## timepointtimepoint3             0.99919    
    ## timepointtimepoint4             0.00393 ** 
    ## siteMahana                      0.34792    
    ## siteManava                      0.18001    
    ## timepointtimepoint2:siteMahana  0.69291    
    ## timepointtimepoint3:siteMahana  0.88020    
    ## timepointtimepoint4:siteMahana  0.79464    
    ## timepointtimepoint2:siteManava  0.44306    
    ## timepointtimepoint3:siteManava  0.01182 *  
    ## timepointtimepoint4:siteManava  0.50061    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.636                                                    
    ## tmpnttmpnt3      -0.358  0.248                                             
    ## tmpnttmpnt4      -0.553  0.373  0.212                                      
    ## siteMahana       -0.721  0.458  0.258  0.399                               
    ## siteManava       -0.745  0.474  0.267  0.412  0.538                        
    ## tmpnttmpnt2:stMh  0.423 -0.666 -0.165 -0.248 -0.583 -0.316                 
    ## tmpnttmpnt3:stMh  0.314 -0.218 -0.877 -0.185 -0.442 -0.234  0.279          
    ## tmpnttmpnt4:stMh  0.430 -0.290 -0.165 -0.779 -0.606 -0.321  0.373          
    ## tmpnttmpnt2:stMn  0.451 -0.709 -0.176 -0.264 -0.325 -0.603  0.472          
    ## tmpnttmpnt3:stMn  0.314 -0.218 -0.878 -0.186 -0.227 -0.423  0.145          
    ## tmpnttmpnt4:stMn  0.433 -0.292 -0.166 -0.784 -0.313 -0.583  0.195          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.286                                            
    ## tmpnttmpnt2:stMn  0.154            0.206                           
    ## tmpnttmpnt3:stMn  0.769            0.145            0.279          
    ## tmpnttmpnt4:stMn  0.145            0.611            0.374          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.264

*Pocillopora*

`prot_model_poc<-lmer(log(1+prot_mg.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_ug.mgafdw>0 & species=="Pocillopora"))`

``` r
prot_model_poc<-lmer(log(1+prot_mg.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_mg.mgafdw>0 & species=="Pocillopora"))
qqPlot(residuals(prot_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-115-1.png)<!-- -->

    ## 173 157 
    ##  75  60

Generate a Type III Anova of model.

``` r
anova(prot_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      0.184139 0.061380     3   149 33.3678 < 2.2e-16 ***
    ## site           0.025378 0.012689     2   149  6.8980  0.001364 ** 
    ## timepoint:site 0.031783 0.005297     6   149  2.8797  0.011025 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(prot_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + prot_mg.mgafdw) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(prot_mg.mgafdw > 0 & species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -484.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9048 -0.5754 -0.2107  0.3614  5.7366 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.000000 0.00000 
    ##  Residual                   0.001839 0.04289 
    ## Number of obs: 161, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      0.241894   0.011074 149.000000  21.844
    ## timepointtimepoint2              0.032055   0.015938 149.000000   2.011
    ## timepointtimepoint3              0.056200   0.016252 149.000000   3.458
    ## timepointtimepoint4             -0.065429   0.015661 149.000000  -4.178
    ## siteMahana                       0.056639   0.016611 149.000000   3.410
    ## siteManava                       0.023095   0.015938 149.000000   1.449
    ## timepointtimepoint2:siteMahana  -0.058208   0.023021 149.000000  -2.529
    ## timepointtimepoint3:siteMahana  -0.061669   0.023642 149.000000  -2.609
    ## timepointtimepoint4:siteMahana   0.008603   0.023491 149.000000   0.366
    ## timepointtimepoint2:siteManava  -0.045293   0.022733 149.000000  -1.992
    ## timepointtimepoint3:siteManava  -0.020183   0.023722 149.000000  -0.851
    ## timepointtimepoint4:siteManava  -0.008664   0.022763 149.000000  -0.381
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## timepointtimepoint2            0.046105 *  
    ## timepointtimepoint3            0.000710 ***
    ## timepointtimepoint4            4.99e-05 ***
    ## siteMahana                     0.000837 ***
    ## siteManava                     0.149422    
    ## timepointtimepoint2:siteMahana 0.012496 *  
    ## timepointtimepoint3:siteMahana 0.010020 *  
    ## timepointtimepoint4:siteMahana 0.714733    
    ## timepointtimepoint2:siteManava 0.048158 *  
    ## timepointtimepoint3:siteManava 0.396253    
    ## timepointtimepoint4:siteManava 0.704029    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.695                                                    
    ## tmpnttmpnt3      -0.681  0.473                                             
    ## tmpnttmpnt4      -0.707  0.491  0.482                                      
    ## siteMahana       -0.667  0.463  0.454  0.471                               
    ## siteManava       -0.695  0.483  0.473  0.491  0.463                        
    ## tmpnttmpnt2:stMh  0.481 -0.692 -0.328 -0.340 -0.722 -0.334                 
    ## tmpnttmpnt3:stMh  0.468 -0.325 -0.687 -0.331 -0.703 -0.325  0.507          
    ## tmpnttmpnt4:stMh  0.471 -0.328 -0.321 -0.667 -0.707 -0.328  0.510          
    ## tmpnttmpnt2:stMn  0.487 -0.701 -0.332 -0.344 -0.325 -0.701  0.485          
    ## tmpnttmpnt3:stMn  0.467 -0.324 -0.685 -0.330 -0.311 -0.672  0.225          
    ## tmpnttmpnt4:stMn  0.486 -0.338 -0.331 -0.688 -0.324 -0.700  0.234          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.497                                            
    ## tmpnttmpnt2:stMn  0.228            0.230                           
    ## tmpnttmpnt3:stMn  0.471            0.220            0.471          
    ## tmpnttmpnt4:stMn  0.228            0.459            0.491          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.470          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

*Porites*

`prot_model_por<-lmer(log(1+prot_mg.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_ug.mgafdw>0 & species=="Porites"))`

``` r
prot_model_por<-lmer(log(1+prot_mg.mgafdw)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(prot_mg.mgafdw>0 & species=="Porites"))
qqPlot(residuals(prot_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-117-1.png)<!-- -->

    ## 109  93 
    ##  26  13

Generate a Type III Anova of model.

``` r
anova(prot_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      0.289692 0.096564     3 102.000  19.370 5.123e-10 ***
    ## site           0.113534 0.056767     2  37.564  11.387 0.0001362 ***
    ## timepoint:site 0.090508 0.015085     6 101.930   3.026 0.0091749 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(prot_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + prot_mg.mgafdw) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(prot_mg.mgafdw > 0 & species == "Porites")
    ## 
    ## REML criterion at convergence: -268.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1191 -0.4810 -0.1088  0.3148  4.5303 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.004402 0.06635 
    ##  Residual                   0.004985 0.07061 
    ## Number of obs: 157, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.42459    0.02701 104.81081  15.719  < 2e-16
    ## timepointtimepoint2             -0.09089    0.02974 104.04497  -3.056  0.00285
    ## timepointtimepoint3             -0.02039    0.02772 101.07533  -0.736  0.46364
    ## timepointtimepoint4             -0.17095    0.02897 102.94050  -5.900 4.68e-08
    ## siteMahana                      -0.08683    0.03774 101.43486  -2.301  0.02346
    ## siteManava                       0.06886    0.03682  96.83640   1.870  0.06444
    ## timepointtimepoint2:siteMahana   0.05637    0.04022 103.09731   1.401  0.16412
    ## timepointtimepoint3:siteMahana  -0.02014    0.04056 102.29325  -0.497  0.62054
    ## timepointtimepoint4:siteMahana   0.12856    0.04084 102.72229   3.148  0.00215
    ## timepointtimepoint2:siteManava   0.02648    0.03973 102.04637   0.667  0.50655
    ## timepointtimepoint3:siteManava   0.01036    0.03874 101.14914   0.267  0.78972
    ## timepointtimepoint4:siteManava   0.03401    0.04014 102.45858   0.847  0.39884
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2            ** 
    ## timepointtimepoint3               
    ## timepointtimepoint4            ***
    ## siteMahana                     *  
    ## siteManava                     .  
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana ** 
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.550                                                    
    ## tmpnttmpnt3      -0.582  0.536                                             
    ## tmpnttmpnt4      -0.562  0.513  0.548                                      
    ## siteMahana       -0.716  0.394  0.417  0.402                               
    ## siteManava       -0.734  0.404  0.427  0.412  0.525                        
    ## tmpnttmpnt2:stMh  0.407 -0.739 -0.397 -0.379 -0.556 -0.299                 
    ## tmpnttmpnt3:stMh  0.398 -0.367 -0.683 -0.374 -0.541 -0.292  0.512          
    ## tmpnttmpnt4:stMh  0.399 -0.364 -0.388 -0.709 -0.540 -0.292  0.508          
    ## tmpnttmpnt2:stMn  0.412 -0.748 -0.401 -0.384 -0.295 -0.529  0.553          
    ## tmpnttmpnt3:stMn  0.417 -0.384 -0.716 -0.392 -0.298 -0.539  0.284          
    ## tmpnttmpnt4:stMn  0.406 -0.370 -0.395 -0.722 -0.290 -0.522  0.274          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.515                                            
    ## tmpnttmpnt2:stMn  0.274            0.273                           
    ## tmpnttmpnt3:stMn  0.489            0.278            0.502          
    ## tmpnttmpnt4:stMn  0.270            0.512            0.484          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.507

## PI Curve Parameters

### Plot: Am

View by site and timepoint for each species.

``` r
am_plot<-master %>%
  filter(!is.na(AQY)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, AQY, Am, Rd)%>%
  
  ggplot(., aes(x = site_code, y = Am, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Am"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: AQY

View by site and timepoint for each species.

``` r
aqy_plot<-master %>%
  filter(!is.na(AQY)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, AQY, Am, Rd)%>%
  
  ggplot(., aes(x = site_code, y = AQY, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("AQY"))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Rd

View by site and timepoint for each species.

``` r
rd_plot<-master %>%
  filter(!is.na(AQY)) %>%
  filter(!is.na(species)) %>%
  select(colony_id_corr, species, site_code, timepoint, AQY, Am, Rd)%>%
  
  ggplot(., aes(x = site_code, y = Rd, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    ylab(expression(bold(paste("Respiration"))))+
    theme_classic() + 
    theme(
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots together.

``` r
pi_figure<-plot_grid(am_plot, aqy_plot, rd_plot, labels = c("", "", ""), ncol=3, nrow=1, rel_widths = c(.8, .8, 1), label_size = 20, label_x = 0.1)

ggsave(filename="Figures/PI_Figure.pdf", plot=pi_figure, dpi=300, width=22, height=5, units="in")
```

### Analysis of Am

Build a mixed model for univariate analysis and examine data
distribution.

`am_model<-lmer(Am~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
am_model<-lmer(Am~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(am_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-123-1.png)<!-- -->

    ## 415 334 
    ## 414 334

Residuals are not normally distributed. Attempt with log transformation.

`am_model<-lmer(log(Am)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
am_model<-lmer(log(Am)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(am_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-124-1.png)<!-- -->

    ## [1] 334  56

Generate a Type III Anova of model.

``` r
anova(am_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
    ## timepoint               4.534   1.511     3 337.77  19.7291 8.287e-12 ***
    ## species                65.636  32.818     2 131.91 428.4050 < 2.2e-16 ***
    ## site                    0.860   0.430     2 142.02   5.6102  0.004519 ** 
    ## timepoint:species       0.911   0.152     6 329.28   1.9811  0.067845 .  
    ## timepoint:site          0.457   0.076     6 334.85   0.9950  0.428488    
    ## species:site            0.993   0.248     4 129.25   3.2393  0.014348 *  
    ## timepoint:species:site  1.095   0.091    12 327.06   1.1912  0.287927    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(am_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Am) ~ timepoint * species * site + (1 | colony_id_corr)
    ##    Data: master
    ## 
    ## REML criterion at convergence: 207.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6121 -0.5139  0.0468  0.5446  3.9066 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001602 0.04002 
    ##  Residual                   0.076605 0.27678 
    ## Number of obs: 447, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                        -0.231806   0.077561
    ## timepointtimepoint2                                 0.280986   0.120465
    ## timepointtimepoint3                                -0.184618   0.211758
    ## timepointtimepoint4                                -0.209823   0.137110
    ## speciesPocillopora                                  0.010977   0.105969
    ## speciesPorites                                      0.820814   0.105969
    ## siteMahana                                          0.434039   0.107708
    ## siteManava                                          0.250807   0.105967
    ## timepointtimepoint2:speciesPocillopora             -0.112904   0.158423
    ## timepointtimepoint3:speciesPocillopora              0.179487   0.236340
    ## timepointtimepoint4:speciesPocillopora              0.145246   0.170333
    ## timepointtimepoint2:speciesPorites                 -0.168742   0.158423
    ## timepointtimepoint3:speciesPorites                  0.136772   0.234639
    ## timepointtimepoint4:speciesPorites                  0.164879   0.171421
    ## timepointtimepoint2:siteMahana                     -0.455256   0.181379
    ## timepointtimepoint3:siteMahana                     -0.309954   0.241103
    ## timepointtimepoint4:siteMahana                     -0.068725   0.175405
    ## timepointtimepoint2:siteManava                     -0.226974   0.167994
    ## timepointtimepoint3:siteManava                     -0.241244   0.241975
    ## timepointtimepoint4:siteManava                     -0.055069   0.175858
    ## speciesPocillopora:siteMahana                      -0.378159   0.151096
    ## speciesPorites:siteMahana                          -0.355786   0.148420
    ## speciesPocillopora:siteManava                       0.009015   0.147162
    ## speciesPorites:siteManava                          -0.136348   0.147162
    ## timepointtimepoint2:speciesPocillopora:siteMahana   0.311361   0.233452
    ## timepointtimepoint3:speciesPocillopora:siteMahana   0.218797   0.284516
    ## timepointtimepoint4:speciesPocillopora:siteMahana  -0.004093   0.232121
    ## timepointtimepoint2:speciesPorites:siteMahana       0.514797   0.231729
    ## timepointtimepoint3:speciesPorites:siteMahana       0.144970   0.282597
    ## timepointtimepoint4:speciesPorites:siteMahana       0.036258   0.229933
    ## timepointtimepoint2:speciesPocillopora:siteManava  -0.044987   0.221409
    ## timepointtimepoint3:speciesPocillopora:siteManava  -0.132388   0.284751
    ## timepointtimepoint4:speciesPocillopora:siteManava  -0.108352   0.228376
    ## timepointtimepoint2:speciesPorites:siteManava       0.351701   0.221409
    ## timepointtimepoint3:speciesPorites:siteManava       0.338999   0.282456
    ## timepointtimepoint4:speciesPorites:siteManava       0.132124   0.230279
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       410.656433  -2.989  0.00297
    ## timepointtimepoint2                               356.018143   2.333  0.02023
    ## timepointtimepoint3                               397.437591  -0.872  0.38382
    ## timepointtimepoint4                               355.095245  -1.530  0.12683
    ## speciesPocillopora                                410.620713   0.104  0.91755
    ## speciesPorites                                    410.620713   7.746 7.50e-14
    ## siteMahana                                        410.729813   4.030 6.65e-05
    ## siteManava                                        410.670648   2.367  0.01840
    ## timepointtimepoint2:speciesPocillopora            333.335131  -0.713  0.47654
    ## timepointtimepoint3:speciesPocillopora            384.332598   0.759  0.44805
    ## timepointtimepoint4:speciesPocillopora            334.677429   0.853  0.39443
    ## timepointtimepoint2:speciesPorites                333.334593  -1.065  0.28759
    ## timepointtimepoint3:speciesPorites                383.720842   0.583  0.56030
    ## timepointtimepoint4:speciesPorites                336.108604   0.962  0.33682
    ## timepointtimepoint2:siteMahana                    358.051031  -2.510  0.01251
    ## timepointtimepoint3:siteMahana                    395.048055  -1.286  0.19935
    ## timepointtimepoint4:siteMahana                    361.173768  -0.392  0.69543
    ## timepointtimepoint2:siteManava                    352.809657  -1.351  0.17753
    ## timepointtimepoint3:siteManava                    388.665251  -0.997  0.31939
    ## timepointtimepoint4:siteManava                    346.588860  -0.313  0.75436
    ## speciesPocillopora:siteMahana                     410.700800  -2.503  0.01271
    ## speciesPorites:siteMahana                         410.661849  -2.397  0.01697
    ## speciesPocillopora:siteManava                     410.627155   0.061  0.95118
    ## speciesPorites:siteManava                         410.627155  -0.927  0.35472
    ## timepointtimepoint2:speciesPocillopora:siteMahana 337.303516   1.334  0.18319
    ## timepointtimepoint3:speciesPocillopora:siteMahana 375.345228   0.769  0.44237
    ## timepointtimepoint4:speciesPocillopora:siteMahana 340.067350  -0.018  0.98594
    ## timepointtimepoint2:speciesPorites:siteMahana     335.571158   2.222  0.02698
    ## timepointtimepoint3:speciesPorites:siteMahana     375.651045   0.513  0.60826
    ## timepointtimepoint4:speciesPorites:siteMahana     339.511367   0.158  0.87480
    ## timepointtimepoint2:speciesPocillopora:siteManava 330.096649  -0.203  0.83912
    ## timepointtimepoint3:speciesPocillopora:siteManava 370.253241  -0.465  0.64226
    ## timepointtimepoint4:speciesPocillopora:siteManava 328.489887  -0.474  0.63550
    ## timepointtimepoint2:speciesPorites:siteManava     330.096371   1.588  0.11314
    ## timepointtimepoint3:speciesPorites:siteManava     369.204183   1.200  0.23084
    ## timepointtimepoint4:speciesPorites:siteManava     330.472965   0.574  0.56652
    ##                                                      
    ## (Intercept)                                       ** 
    ## timepointtimepoint2                               *  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## speciesPocillopora                                   
    ## speciesPorites                                    ***
    ## siteMahana                                        ***
    ## siteManava                                        *  
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## timepointtimepoint2:siteMahana                    *  
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## speciesPocillopora:siteMahana                     *  
    ## speciesPorites:siteMahana                         *  
    ## speciesPocillopora:siteManava                        
    ## speciesPorites:siteManava                            
    ## timepointtimepoint2:speciesPocillopora:siteMahana    
    ## timepointtimepoint3:speciesPocillopora:siteMahana    
    ## timepointtimepoint4:speciesPocillopora:siteMahana    
    ## timepointtimepoint2:speciesPorites:siteMahana     *  
    ## timepointtimepoint3:speciesPorites:siteMahana        
    ## timepointtimepoint4:speciesPorites:siteMahana        
    ## timepointtimepoint2:speciesPocillopora:siteManava    
    ## timepointtimepoint3:speciesPocillopora:siteManava    
    ## timepointtimepoint4:speciesPocillopora:siteManava    
    ## timepointtimepoint2:speciesPorites:siteManava        
    ## timepointtimepoint3:speciesPorites:siteManava        
    ## timepointtimepoint4:speciesPorites:siteManava        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of time and site within
species.

*Acropora*

`am_model<-lmer(log(Am)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))`

``` r
am_model_acr<-lmer(log(Am)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))
qqPlot(residuals(am_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-126-1.png)<!-- -->

    ## 334 333 
    ## 103 102

Generate a Type III Anova of model.

``` r
anova(am_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      2.37550 0.79183     3 84.277 12.7822 5.807e-07 ***
    ## site           0.46616 0.23308     2 49.962  3.7625   0.03005 *  
    ## timepoint:site 0.57786 0.09631     6 83.112  1.5547   0.17078    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(am_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Am) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Acropora")
    ## 
    ## REML criterion at convergence: 46.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6691 -0.4407  0.0903  0.4592  3.8759 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.01041  0.1020  
    ##  Residual                   0.06195  0.2489  
    ## Number of obs: 116, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                     -0.23516    0.07451 101.82701  -3.156  0.00210
    ## timepointtimepoint2              0.29672    0.11075  84.60449   2.679  0.00886
    ## timepointtimepoint3             -0.20340    0.19727  91.01112  -1.031  0.30523
    ## timepointtimepoint4             -0.21149    0.12579  81.57705  -1.681  0.09654
    ## siteMahana                       0.43312    0.10336 102.46731   4.190 5.92e-05
    ## siteManava                       0.25473    0.10180 101.84257   2.502  0.01393
    ## timepointtimepoint2:siteMahana  -0.48271    0.16669  83.53035  -2.896  0.00483
    ## timepointtimepoint3:siteMahana  -0.28910    0.22454  91.58103  -1.287  0.20116
    ## timepointtimepoint4:siteMahana  -0.06657    0.16140  84.48689  -0.412  0.68103
    ## timepointtimepoint2:siteManava  -0.23445    0.15411  82.10016  -1.521  0.13202
    ## timepointtimepoint3:siteManava  -0.21284    0.22450  88.51681  -0.948  0.34569
    ## timepointtimepoint4:siteManava  -0.05352    0.16084  78.98018  -0.333  0.74020
    ##                                   
    ## (Intercept)                    ** 
    ## timepointtimepoint2            ** 
    ## timepointtimepoint3               
    ## timepointtimepoint4            .  
    ## siteMahana                     ***
    ## siteManava                     *  
    ## timepointtimepoint2:siteMahana ** 
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.597                                                    
    ## tmpnttmpnt3      -0.325  0.239                                             
    ## tmpnttmpnt4      -0.509  0.352  0.195                                      
    ## siteMahana       -0.721  0.430  0.234  0.367                               
    ## siteManava       -0.732  0.437  0.238  0.372  0.528                        
    ## tmpnttmpnt2:stMh  0.397 -0.664 -0.158 -0.234 -0.544 -0.290                 
    ## tmpnttmpnt3:stMh  0.285 -0.210 -0.879 -0.171 -0.410 -0.209  0.269          
    ## tmpnttmpnt4:stMh  0.396 -0.274 -0.152 -0.779 -0.567 -0.290  0.352          
    ## tmpnttmpnt2:stMn  0.429 -0.719 -0.171 -0.253 -0.309 -0.581  0.477          
    ## tmpnttmpnt3:stMn  0.285 -0.210 -0.879 -0.171 -0.206 -0.392  0.139          
    ## tmpnttmpnt4:stMn  0.398 -0.275 -0.153 -0.782 -0.287 -0.546  0.183          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.273                                            
    ## tmpnttmpnt2:stMn  0.151            0.197                           
    ## tmpnttmpnt3:stMn  0.772            0.134            0.281          
    ## tmpnttmpnt4:stMn  0.134            0.610            0.367          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.252

*Pocillopora*

`am_model_poc<-lmer(log(Am)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))`

``` r
am_model_poc<-lmer(log(Am)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))
qqPlot(residuals(am_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-128-1.png)<!-- -->

    ##  56 352 
    ##  16 128

Generate a Type III Anova of model.

``` r
anova(am_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      1.12774 0.37591     3 114.357  6.3913 0.0004831 ***
    ## site           0.11516 0.05758     2  40.163  0.9790 0.3844659    
    ## timepoint:site 0.59153 0.09859     6 114.279  1.6762 0.1330818    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(am_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Am) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Pocillopora")
    ## 
    ## REML criterion at convergence: 46.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.9885 -0.4961  0.0079  0.5418  2.3759 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.006783 0.08236 
    ##  Residual                   0.058816 0.24252 
    ## Number of obs: 164, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                    -2.208e-01  6.613e-02  1.475e+02  -3.339
    ## timepointtimepoint2             1.675e-01  9.026e-02  1.120e+02   1.856
    ## timepointtimepoint3            -3.341e-04  9.218e-02  1.138e+02  -0.004
    ## timepointtimepoint4            -6.458e-02  8.856e-02  1.104e+02  -0.729
    ## siteMahana                      6.460e-02  9.699e-02  1.484e+02   0.666
    ## siteManava                      2.598e-01  9.352e-02  1.475e+02   2.778
    ## timepointtimepoint2:siteMahana -1.520e-01  1.290e-01  1.133e+02  -1.178
    ## timepointtimepoint3:siteMahana -1.020e-01  1.327e-01  1.141e+02  -0.768
    ## timepointtimepoint4:siteMahana -7.912e-02  1.336e-01  1.151e+02  -0.592
    ## timepointtimepoint2:siteManava -2.714e-01  1.264e-01  1.112e+02  -2.146
    ## timepointtimepoint3:siteManava -3.779e-01  1.319e-01  1.150e+02  -2.864
    ## timepointtimepoint4:siteManava -1.637e-01  1.278e-01  1.125e+02  -1.280
    ##                                Pr(>|t|)   
    ## (Intercept)                     0.00106 **
    ## timepointtimepoint2             0.06608 . 
    ## timepointtimepoint3             0.99711   
    ## timepointtimepoint4             0.46740   
    ## siteMahana                      0.50643   
    ## siteManava                      0.00618 **
    ## timepointtimepoint2:siteMahana  0.24111   
    ## timepointtimepoint3:siteMahana  0.44390   
    ## timepointtimepoint4:siteMahana  0.55496   
    ## timepointtimepoint2:siteManava  0.03402 * 
    ## timepointtimepoint3:siteManava  0.00497 **
    ## timepointtimepoint4:siteManava  0.20308   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.657                                                    
    ## tmpnttmpnt3      -0.643  0.471                                             
    ## tmpnttmpnt4      -0.670  0.491  0.480                                      
    ## siteMahana       -0.682  0.448  0.439  0.457                               
    ## siteManava       -0.707  0.465  0.455  0.473  0.482                        
    ## tmpnttmpnt2:stMh  0.460 -0.700 -0.329 -0.343 -0.679 -0.325                 
    ## tmpnttmpnt3:stMh  0.447 -0.327 -0.695 -0.334 -0.659 -0.316  0.495          
    ## tmpnttmpnt4:stMh  0.444 -0.325 -0.318 -0.663 -0.654 -0.314  0.492          
    ## tmpnttmpnt2:stMn  0.469 -0.714 -0.336 -0.350 -0.320 -0.663  0.499          
    ## tmpnttmpnt3:stMn  0.449 -0.329 -0.699 -0.336 -0.306 -0.636  0.230          
    ## tmpnttmpnt4:stMn  0.464 -0.340 -0.333 -0.693 -0.316 -0.656  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.480                                            
    ## tmpnttmpnt2:stMn  0.233            0.232                           
    ## tmpnttmpnt3:stMn  0.485            0.222            0.470          
    ## tmpnttmpnt4:stMn  0.231            0.459            0.485          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.468

*Porites*

`am_model_por<-lmer(log(Am)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))`

``` r
am_model_por<-lmer(log(Am)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))
qqPlot(residuals(am_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-130-1.png)<!-- -->

    ## 107 198 
    ##  29  59

Generate a Type III Anova of model.

``` r
anova(am_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
    ## timepoint      1.46332 0.48777     3   155  5.1688 0.001971 **
    ## site           1.09079 0.54540     2   155  5.7794 0.003795 **
    ## timepoint:site 0.32376 0.05396     6   155  0.5718 0.752343   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(am_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Am) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Porites")
    ## 
    ## REML criterion at convergence: 105.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.79518 -0.59976  0.06642  0.71022  2.26242 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.00000  0.0000  
    ##  Residual                   0.09437  0.3072  
    ## Number of obs: 167, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                      0.58901    0.07932 155.00000   7.426 6.99e-12
    ## timepointtimepoint2              0.11137    0.11416 155.00000   0.976    0.331
    ## timepointtimepoint3             -0.04785    0.11217 155.00000  -0.427    0.670
    ## timepointtimepoint4             -0.04534    0.11416 155.00000  -0.397    0.692
    ## siteMahana                       0.07825    0.11217 155.00000   0.698    0.486
    ## siteManava                       0.11446    0.11217 155.00000   1.020    0.309
    ## timepointtimepoint2:siteMahana   0.06041    0.16004 155.00000   0.377    0.706
    ## timepointtimepoint3:siteMahana  -0.16330    0.16352 155.00000  -0.999    0.319
    ## timepointtimepoint4:siteMahana  -0.03039    0.16488 155.00000  -0.184    0.854
    ## timepointtimepoint2:siteManava   0.12560    0.16004 155.00000   0.785    0.434
    ## timepointtimepoint3:siteManava   0.09743    0.16166 155.00000   0.603    0.548
    ## timepointtimepoint4:siteManava   0.07803    0.16488 155.00000   0.473    0.637
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                        
    ## siteManava                        
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.695                                                    
    ## tmpnttmpnt3      -0.707  0.491                                             
    ## tmpnttmpnt4      -0.695  0.483  0.491                                      
    ## siteMahana       -0.707  0.491  0.500  0.491                               
    ## siteManava       -0.707  0.491  0.500  0.491  0.500                        
    ## tmpnttmpnt2:stMh  0.496 -0.713 -0.350 -0.344 -0.701 -0.350                 
    ## tmpnttmpnt3:stMh  0.485 -0.337 -0.686 -0.337 -0.686 -0.343  0.481          
    ## tmpnttmpnt4:stMh  0.481 -0.334 -0.340 -0.692 -0.680 -0.340  0.477          
    ## tmpnttmpnt2:stMn  0.496 -0.713 -0.350 -0.344 -0.350 -0.701  0.509          
    ## tmpnttmpnt3:stMn  0.491 -0.341 -0.694 -0.341 -0.347 -0.694  0.243          
    ## tmpnttmpnt4:stMn  0.481 -0.334 -0.340 -0.692 -0.340 -0.680  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.467                                            
    ## tmpnttmpnt2:stMn  0.240            0.238                           
    ## tmpnttmpnt3:stMn  0.476            0.236            0.486          
    ## tmpnttmpnt4:stMn  0.233            0.479            0.477          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.472          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

### Analysis of AQY

Build a mixed model for univariate analysis and examine data
distribution.

`aqy_model<-lmer(AQY~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
aqy_model<-lmer(AQY~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(aqy_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-132-1.png)<!-- -->

    ## [1] 281  79

Residuals are not normally distributed. Attempt with log transformation.

`aqy_model<-lmer(log(AQY)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
aqy_model<-lmer(log(AQY)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(aqy_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-133-1.png)<!-- -->

    ## [1] 334 180

Generate a Type II Anova of model.

``` r
anova(aqy_model, type="II")
```

    ## Type II Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
    ## timepoint               5.851   1.950     3 312.35  13.3936 3.028e-08 ***
    ## species                89.853  44.926     2 117.32 308.5228 < 2.2e-16 ***
    ## site                    0.321   0.160     2 110.30   1.1019   0.33587    
    ## timepoint:species       2.449   0.408     6 318.72   2.8035   0.01133 *  
    ## timepoint:site          0.974   0.162     6 312.76   1.1149   0.35334    
    ## species:site            1.911   0.478     4 120.16   3.2808   0.01366 *  
    ## timepoint:species:site  1.760   0.147    12 321.22   1.0069   0.44207    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aqy_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(AQY) ~ timepoint * species * site + (1 | colony_id_corr)
    ##    Data: master
    ## 
    ## REML criterion at convergence: 468.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3884 -0.5780  0.0372  0.5565  3.6905 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.001825 0.04272 
    ##  Residual                   0.145618 0.38160 
    ## Number of obs: 447, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                        -6.110598   0.106497
    ## timepointtimepoint2                                -0.065435   0.165844
    ## timepointtimepoint3                                 0.361498   0.291136
    ## timepointtimepoint4                                -0.207863   0.188763
    ## speciesPocillopora                                  0.090012   0.145503
    ## speciesPorites                                      1.181590   0.145503
    ## siteMahana                                          0.291852   0.147894
    ## siteManava                                          0.127450   0.145502
    ## timepointtimepoint2:speciesPocillopora             -0.147683   0.218224
    ## timepointtimepoint3:speciesPocillopora             -0.248182   0.325096
    ## timepointtimepoint4:speciesPocillopora              0.177000   0.234621
    ## timepointtimepoint2:speciesPorites                 -0.132178   0.218224
    ## timepointtimepoint3:speciesPorites                 -0.191789   0.322763
    ## timepointtimepoint4:speciesPorites                 -0.122055   0.236112
    ## timepointtimepoint2:siteMahana                     -0.295721   0.249690
    ## timepointtimepoint3:siteMahana                     -0.655468   0.331516
    ## timepointtimepoint4:siteMahana                      0.023280   0.241446
    ## timepointtimepoint2:siteManava                     -0.086146   0.231295
    ## timepointtimepoint3:siteManava                     -0.656662   0.332796
    ## timepointtimepoint4:siteManava                      0.047292   0.242161
    ## speciesPocillopora:siteMahana                      -0.319229   0.207469
    ## speciesPorites:siteMahana                          -0.500351   0.203793
    ## speciesPocillopora:siteManava                       0.143830   0.202064
    ## speciesPorites:siteManava                          -0.186851   0.202064
    ## timepointtimepoint2:speciesPocillopora:siteMahana   0.233567   0.321543
    ## timepointtimepoint3:speciesPocillopora:siteMahana   0.862177   0.391477
    ## timepointtimepoint4:speciesPocillopora:siteMahana  -0.081060   0.319688
    ## timepointtimepoint2:speciesPorites:siteMahana       0.259285   0.319184
    ## timepointtimepoint3:speciesPorites:siteMahana       0.716593   0.388833
    ## timepointtimepoint4:speciesPorites:siteMahana       0.261314   0.316679
    ## timepointtimepoint2:speciesPocillopora:siteManava  -0.130515   0.305009
    ## timepointtimepoint3:speciesPocillopora:siteManava   0.596536   0.391858
    ## timepointtimepoint4:speciesPocillopora:siteManava  -0.089564   0.314619
    ## timepointtimepoint2:speciesPorites:siteManava       0.225222   0.305009
    ## timepointtimepoint3:speciesPorites:siteManava       0.424829   0.388711
    ## timepointtimepoint4:speciesPorites:siteManava      -0.008168   0.317225
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       410.861900 -57.378  < 2e-16
    ## timepointtimepoint2                               351.125958  -0.395   0.6934
    ## timepointtimepoint3                               397.282955   1.242   0.2151
    ## timepointtimepoint4                               350.833732  -1.101   0.2716
    ## speciesPocillopora                                410.847752   0.619   0.5365
    ## speciesPorites                                    410.847752   8.121 5.46e-15
    ## siteMahana                                        410.891391   1.973   0.0491
    ## siteManava                                        410.867970   0.876   0.3816
    ## timepointtimepoint2:speciesPocillopora            327.291601  -0.677   0.4990
    ## timepointtimepoint3:speciesPocillopora            383.142842  -0.763   0.4457
    ## timepointtimepoint4:speciesPocillopora            329.091272   0.754   0.4511
    ## timepointtimepoint2:speciesPorites                327.291250  -0.606   0.5451
    ## timepointtimepoint3:speciesPorites                382.462205  -0.594   0.5527
    ## timepointtimepoint4:speciesPorites                330.650063  -0.517   0.6055
    ## timepointtimepoint2:siteMahana                    353.730879  -1.184   0.2371
    ## timepointtimepoint3:siteMahana                    394.438689  -1.977   0.0487
    ## timepointtimepoint4:siteMahana                    357.048084   0.096   0.9232
    ## timepointtimepoint2:siteManava                    348.104114  -0.372   0.7098
    ## timepointtimepoint3:siteManava                    387.743472  -1.973   0.0492
    ## timepointtimepoint4:siteManava                    341.836555   0.195   0.8453
    ## speciesPocillopora:siteMahana                     410.879986  -1.539   0.1247
    ## speciesPorites:siteMahana                         410.864279  -2.455   0.0145
    ## speciesPocillopora:siteManava                     410.850545   0.712   0.4770
    ## speciesPorites:siteManava                         410.850544  -0.925   0.3557
    ## timepointtimepoint2:speciesPocillopora:siteMahana 331.737167   0.726   0.4681
    ## timepointtimepoint3:speciesPocillopora:siteMahana 373.116414   2.202   0.0283
    ## timepointtimepoint4:speciesPocillopora:siteMahana 334.722272  -0.254   0.8000
    ## timepointtimepoint2:speciesPorites:siteMahana     329.884674   0.812   0.4172
    ## timepointtimepoint3:speciesPorites:siteMahana     373.419270   1.843   0.0661
    ## timepointtimepoint4:speciesPorites:siteMahana     334.072789   0.825   0.4099
    ## timepointtimepoint2:speciesPocillopora:siteManava 324.098007  -0.428   0.6690
    ## timepointtimepoint3:speciesPocillopora:siteManava 367.802644   1.522   0.1288
    ## timepointtimepoint4:speciesPocillopora:siteManava 322.616479  -0.285   0.7761
    ## timepointtimepoint2:speciesPorites:siteManava     324.097828   0.738   0.4608
    ## timepointtimepoint3:speciesPorites:siteManava     366.641017   1.093   0.2751
    ## timepointtimepoint4:speciesPorites:siteManava     324.765488  -0.026   0.9795
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                                  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                                  
    ## speciesPocillopora                                   
    ## speciesPorites                                    ***
    ## siteMahana                                        *  
    ## siteManava                                           
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                   
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                    *  
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                    *  
    ## timepointtimepoint4:siteManava                       
    ## speciesPocillopora:siteMahana                        
    ## speciesPorites:siteMahana                         *  
    ## speciesPocillopora:siteManava                        
    ## speciesPorites:siteManava                            
    ## timepointtimepoint2:speciesPocillopora:siteMahana    
    ## timepointtimepoint3:speciesPocillopora:siteMahana *  
    ## timepointtimepoint4:speciesPocillopora:siteMahana    
    ## timepointtimepoint2:speciesPorites:siteMahana        
    ## timepointtimepoint3:speciesPorites:siteMahana     .  
    ## timepointtimepoint4:speciesPorites:siteMahana        
    ## timepointtimepoint2:speciesPocillopora:siteManava    
    ## timepointtimepoint3:speciesPocillopora:siteManava    
    ## timepointtimepoint4:speciesPocillopora:siteManava    
    ## timepointtimepoint2:speciesPorites:siteManava        
    ## timepointtimepoint3:speciesPorites:siteManava        
    ## timepointtimepoint4:speciesPorites:siteManava        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effect of site and time within
species.

*Acropora*

`aqy_model_acr<-lmer(log(AQY)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))`

``` r
aqy_model_acr<-lmer(log(AQY)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))
qqPlot(residuals(aqy_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-135-1.png)<!-- -->

    ## 334 328 
    ## 103  97

Generate a Type III Anova of model.

``` r
anova(aqy_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
    ## timepoint      0.74545 0.248482     3 71.974  2.6507 0.05518 .
    ## site           0.08729 0.043645     2 48.014  0.4656 0.63057  
    ## timepoint:site 0.64921 0.108202     6 71.308  1.1543 0.34053  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aqy_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(AQY) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Acropora")
    ## 
    ## REML criterion at convergence: 109.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4627 -0.5023  0.0327  0.4672  3.0306 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.04820  0.2195  
    ##  Residual                   0.09374  0.3062  
    ## Number of obs: 116, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                Estimate Std. Error       df t value Pr(>|t|)
    ## (Intercept)                    -6.12462    0.10376 91.15240 -59.028   <2e-16
    ## timepointtimepoint2            -0.04635    0.14102 77.08870  -0.329   0.7433
    ## timepointtimepoint3             0.24944    0.25204 74.72474   0.990   0.3255
    ## timepointtimepoint4            -0.24029    0.15884 70.34454  -1.513   0.1348
    ## siteMahana                      0.31657    0.14331 94.10673   2.209   0.0296
    ## siteManava                      0.13994    0.14181 90.53675   0.987   0.3264
    ## timepointtimepoint2:siteMahana -0.33041    0.21148 73.89322  -1.562   0.1225
    ## timepointtimepoint3:siteMahana -0.57096    0.28768 77.20345  -1.985   0.0507
    ## timepointtimepoint4:siteMahana  0.04021    0.20507 74.63434   0.196   0.8451
    ## timepointtimepoint2:siteManava -0.07856    0.19506 72.58916  -0.403   0.6883
    ## timepointtimepoint3:siteManava -0.50493    0.28582 73.16175  -1.767   0.0815
    ## timepointtimepoint4:siteManava  0.09287    0.20221 67.94549   0.459   0.6475
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                     *  
    ## siteManava                        
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana .  
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava .  
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.536                                                    
    ## tmpnttmpnt3      -0.279  0.247                                             
    ## tmpnttmpnt4      -0.439  0.345  0.189                                      
    ## siteMahana       -0.724  0.388  0.202  0.318                               
    ## siteManava       -0.732  0.392  0.204  0.321  0.530                        
    ## tmpnttmpnt2:stMh  0.357 -0.667 -0.165 -0.230 -0.485 -0.261                 
    ## tmpnttmpnt3:stMh  0.244 -0.217 -0.876 -0.165 -0.367 -0.179  0.280          
    ## tmpnttmpnt4:stMh  0.340 -0.267 -0.146 -0.775 -0.507 -0.249  0.349          
    ## tmpnttmpnt2:stMn  0.387 -0.723 -0.179 -0.249 -0.280 -0.516  0.482          
    ## tmpnttmpnt3:stMn  0.246 -0.218 -0.882 -0.167 -0.178 -0.338  0.145          
    ## tmpnttmpnt4:stMn  0.345 -0.271 -0.148 -0.786 -0.250 -0.475  0.181          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.281                                            
    ## tmpnttmpnt2:stMn  0.157            0.193                           
    ## tmpnttmpnt3:stMn  0.773            0.129            0.292          
    ## tmpnttmpnt4:stMn  0.130            0.608            0.360          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.247

*Pocillopora*

`aqy_model_poc<-lmer(log(AQY)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))`

``` r
aqy_model_poc<-lmer(log(AQY)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))
qqPlot(residuals(aqy_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-137-1.png)<!-- -->

    ##  56 144 
    ##  16  49

Generate a Type III Anova of model.

``` r
anova(aqy_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      4.6742 1.55807     3   152 10.9477 1.498e-06 ***
    ## site           1.3674 0.68371     2   152  4.8041  0.009482 ** 
    ## timepoint:site 0.5733 0.09555     6   152  0.6714  0.672915    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aqy_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(AQY) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Pocillopora")
    ## 
    ## REML criterion at convergence: 166.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3994 -0.6031  0.0844  0.6758  2.6285 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.0000   0.0000  
    ##  Residual                   0.1423   0.3773  
    ## Number of obs: 164, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                     -6.02059    0.09741 152.00000 -61.809   <2e-16
    ## timepointtimepoint2             -0.21272    0.14019 152.00000  -1.517   0.1313
    ## timepointtimepoint3              0.11199    0.14295 152.00000   0.783   0.4346
    ## timepointtimepoint4             -0.03086    0.13775 152.00000  -0.224   0.8230
    ## siteMahana                      -0.02737    0.14295 152.00000  -0.191   0.8484
    ## siteManava                       0.27128    0.13775 152.00000   1.969   0.0507
    ## timepointtimepoint2:siteMahana  -0.06256    0.20022 152.00000  -0.312   0.7551
    ## timepointtimepoint3:siteMahana   0.20785    0.20574 152.00000   1.010   0.3140
    ## timepointtimepoint4:siteMahana  -0.05826    0.20703 152.00000  -0.281   0.7788
    ## timepointtimepoint2:siteManava  -0.21706    0.19654 152.00000  -1.104   0.2712
    ## timepointtimepoint3:siteManava  -0.05944    0.20441 152.00000  -0.291   0.7716
    ## timepointtimepoint4:siteManava  -0.04253    0.19852 152.00000  -0.214   0.8306
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4               
    ## siteMahana                        
    ## siteManava                     .  
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.695                                                    
    ## tmpnttmpnt3      -0.681  0.473                                             
    ## tmpnttmpnt4      -0.707  0.491  0.482                                      
    ## siteMahana       -0.681  0.473  0.464  0.482                               
    ## siteManava       -0.707  0.491  0.482  0.500  0.482                        
    ## tmpnttmpnt2:stMh  0.486 -0.700 -0.331 -0.344 -0.714 -0.344                 
    ## tmpnttmpnt3:stMh  0.473 -0.329 -0.695 -0.335 -0.695 -0.335  0.496          
    ## tmpnttmpnt4:stMh  0.470 -0.327 -0.321 -0.665 -0.690 -0.333  0.493          
    ## tmpnttmpnt2:stMn  0.496 -0.713 -0.338 -0.350 -0.338 -0.701  0.499          
    ## tmpnttmpnt3:stMn  0.477 -0.331 -0.699 -0.337 -0.325 -0.674  0.232          
    ## tmpnttmpnt4:stMn  0.491 -0.341 -0.334 -0.694 -0.334 -0.694  0.239          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.480                                            
    ## tmpnttmpnt2:stMn  0.235            0.233                           
    ## tmpnttmpnt3:stMn  0.486            0.224            0.472          
    ## tmpnttmpnt4:stMn  0.232            0.462            0.486          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.468          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

*Porites*

`aqy_model_por<-lmer(log(AQY)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))`

``` r
aqy_model_por<-lmer(log(AQY)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))
qqPlot(residuals(aqy_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-139-1.png)<!-- -->

    ## 180 281 
    ##  41  83

Generate a Type III Anova of model.

``` r
anova(aqy_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      2.80259 0.93420     3   155  5.9448 0.0007289 ***
    ## site           0.47448 0.23724     2   155  1.5097 0.2242046    
    ## timepoint:site 1.29268 0.21545     6   155  1.3710 0.2296000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aqy_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(AQY) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Porites")
    ## 
    ## REML criterion at convergence: 184.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3018 -0.5786  0.0305  0.6003  2.6132 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.0000   0.0000  
    ##  Residual                   0.1571   0.3964  
    ## Number of obs: 167, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                     -4.92901    0.10235 155.00000 -48.157   <2e-16
    ## timepointtimepoint2             -0.19824    0.14731 155.00000  -1.346   0.1803
    ## timepointtimepoint3              0.16971    0.14475 155.00000   1.172   0.2428
    ## timepointtimepoint4             -0.33009    0.14731 155.00000  -2.241   0.0265
    ## siteMahana                      -0.20850    0.14475 155.00000  -1.440   0.1518
    ## siteManava                      -0.05940    0.14475 155.00000  -0.410   0.6821
    ## timepointtimepoint2:siteMahana  -0.03580    0.20653 155.00000  -0.173   0.8626
    ## timepointtimepoint3:siteMahana   0.06267    0.21101 155.00000   0.297   0.7669
    ## timepointtimepoint4:siteMahana   0.28630    0.21277 155.00000   1.346   0.1804
    ## timepointtimepoint2:siteManava   0.13971    0.20653 155.00000   0.676   0.4998
    ## timepointtimepoint3:siteManava  -0.23182    0.20861 155.00000  -1.111   0.2682
    ## timepointtimepoint4:siteManava   0.03985    0.21277 155.00000   0.187   0.8517
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4            *  
    ## siteMahana                        
    ## siteManava                        
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.695                                                    
    ## tmpnttmpnt3      -0.707  0.491                                             
    ## tmpnttmpnt4      -0.695  0.483  0.491                                      
    ## siteMahana       -0.707  0.491  0.500  0.491                               
    ## siteManava       -0.707  0.491  0.500  0.491  0.500                        
    ## tmpnttmpnt2:stMh  0.496 -0.713 -0.350 -0.344 -0.701 -0.350                 
    ## tmpnttmpnt3:stMh  0.485 -0.337 -0.686 -0.337 -0.686 -0.343  0.481          
    ## tmpnttmpnt4:stMh  0.481 -0.334 -0.340 -0.692 -0.680 -0.340  0.477          
    ## tmpnttmpnt2:stMn  0.496 -0.713 -0.350 -0.344 -0.350 -0.701  0.509          
    ## tmpnttmpnt3:stMn  0.491 -0.341 -0.694 -0.341 -0.347 -0.694  0.243          
    ## tmpnttmpnt4:stMn  0.481 -0.334 -0.340 -0.692 -0.340 -0.680  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.467                                            
    ## tmpnttmpnt2:stMn  0.240            0.238                           
    ## tmpnttmpnt3:stMn  0.476            0.236            0.486          
    ## tmpnttmpnt4:stMn  0.233            0.479            0.477          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.472          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

### Analysis of Rd

Build a mixed model for univariate analysis and examine data
distribution.

`rd_model<-lmer(Rd~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
rd_model<-lmer(Rd~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(rd_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-141-1.png)<!-- -->

    ## [1] 334  79

Residuals are not normally distributed. Attempt with log transformation.

`rd_model<-lmer(log(Rd)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
rd_model<-lmer(log(Rd)~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(rd_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-142-1.png)<!-- -->

    ## 383 328 
    ## 382 328

Generate a Type II Anova of model.

``` r
anova(rd_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
    ## timepoint               7.949  2.6496     3 333.47  26.6368 1.810e-15 ***
    ## species                43.369 21.6845     2 127.42 217.9990 < 2.2e-16 ***
    ## site                    1.061  0.5304     2 137.06   5.3318   0.00589 ** 
    ## timepoint:species       4.125  0.6875     6 324.73   6.9116 6.350e-07 ***
    ## timepoint:site          1.289  0.2148     6 330.50   2.1591   0.04660 *  
    ## species:site            2.682  0.6706     4 124.84   6.7419 6.032e-05 ***
    ## timepoint:species:site  1.162  0.0968    12 322.49   0.9734   0.47414    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(rd_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Rd) ~ timepoint * species * site + (1 | colony_id_corr)
    ##    Data: master
    ## 
    ## REML criterion at convergence: 320.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.7715 -0.4827  0.0318  0.5666  4.0040 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.003434 0.0586  
    ##  Residual                   0.099471 0.3154  
    ## Number of obs: 447, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                    Estimate Std. Error
    ## (Intercept)                                        -0.95821    0.08896
    ## timepointtimepoint2                                -0.19669    0.13759
    ## timepointtimepoint3                                -0.11694    0.24235
    ## timepointtimepoint4                                -0.76258    0.15659
    ## speciesPocillopora                                 -0.05159    0.12155
    ## speciesPorites                                      0.80784    0.12155
    ## siteMahana                                          0.40868    0.12354
    ## siteManava                                          0.12758    0.12155
    ## timepointtimepoint2:speciesPocillopora              0.23439    0.18078
    ## timepointtimepoint3:speciesPocillopora              0.01734    0.27027
    ## timepointtimepoint4:speciesPocillopora              0.71056    0.19438
    ## timepointtimepoint2:speciesPorites                  0.05958    0.18078
    ## timepointtimepoint3:speciesPorites                 -0.08134    0.26832
    ## timepointtimepoint4:speciesPorites                  0.45684    0.19563
    ## timepointtimepoint2:siteMahana                     -0.02519    0.20718
    ## timepointtimepoint3:siteMahana                     -0.38553    0.27589
    ## timepointtimepoint4:siteMahana                      0.23441    0.20038
    ## timepointtimepoint2:siteManava                      0.12100    0.19185
    ## timepointtimepoint3:siteManava                     -0.12831    0.27678
    ## timepointtimepoint4:siteManava                     -0.01851    0.20078
    ## speciesPocillopora:siteMahana                      -0.36943    0.17331
    ## speciesPorites:siteMahana                          -0.45159    0.17024
    ## speciesPocillopora:siteManava                      -0.04568    0.16880
    ## speciesPorites:siteManava                          -0.07296    0.16880
    ## timepointtimepoint2:speciesPocillopora:siteMahana  -0.01914    0.26644
    ## timepointtimepoint3:speciesPocillopora:siteMahana   0.44341    0.32522
    ## timepointtimepoint4:speciesPocillopora:siteMahana  -0.22954    0.26495
    ## timepointtimepoint2:speciesPorites:siteMahana       0.11150    0.26446
    ## timepointtimepoint3:speciesPorites:siteMahana       0.26811    0.32303
    ## timepointtimepoint4:speciesPorites:siteMahana      -0.35428    0.26245
    ## timepointtimepoint2:speciesPocillopora:siteManava  -0.17886    0.25263
    ## timepointtimepoint3:speciesPocillopora:siteManava   0.58561    0.32541
    ## timepointtimepoint4:speciesPocillopora:siteManava   0.01673    0.26056
    ## timepointtimepoint2:speciesPorites:siteManava      -0.02234    0.25263
    ## timepointtimepoint3:speciesPorites:siteManava       0.24754    0.32278
    ## timepointtimepoint4:speciesPorites:siteManava      -0.01584    0.26275
    ##                                                          df t value Pr(>|t|)
    ## (Intercept)                                       410.04200 -10.771  < 2e-16
    ## timepointtimepoint2                               353.79841  -1.430 0.153741
    ## timepointtimepoint3                               394.92829  -0.483 0.629686
    ## timepointtimepoint4                               351.76954  -4.870 1.69e-06
    ## speciesPocillopora                                409.94022  -0.424 0.671468
    ## speciesPorites                                    409.94022   6.646 9.62e-11
    ## siteMahana                                        410.24643   3.308 0.001022
    ## siteManava                                        410.07771   1.050 0.294496
    ## timepointtimepoint2:speciesPocillopora            329.82946   1.297 0.195709
    ## timepointtimepoint3:speciesPocillopora            381.02499   0.064 0.948876
    ## timepointtimepoint4:speciesPocillopora            330.62595   3.656 0.000298
    ## timepointtimepoint2:speciesPorites                329.82859   0.330 0.741936
    ## timepointtimepoint3:speciesPorites                380.41082  -0.303 0.761942
    ## timepointtimepoint4:speciesPorites                332.04226   2.335 0.020128
    ## timepointtimepoint2:siteMahana                    355.25577  -0.122 0.903298
    ## timepointtimepoint3:siteMahana                    392.80417  -1.397 0.163073
    ## timepointtimepoint4:siteMahana                    358.54736   1.170 0.242849
    ## timepointtimepoint2:siteManava                    349.83804   0.631 0.528640
    ## timepointtimepoint3:siteManava                    385.71694  -0.464 0.643219
    ## timepointtimepoint4:siteManava                    342.81823  -0.092 0.926608
    ## speciesPocillopora:siteMahana                     410.16282  -2.132 0.033632
    ## speciesPorites:siteMahana                         410.05484  -2.653 0.008297
    ## speciesPocillopora:siteManava                     409.95595  -0.271 0.786839
    ## speciesPorites:siteManava                         409.95595  -0.432 0.665825
    ## timepointtimepoint2:speciesPocillopora:siteMahana 333.55880  -0.072 0.942775
    ## timepointtimepoint3:speciesPocillopora:siteMahana 372.09001   1.363 0.173576
    ## timepointtimepoint4:speciesPocillopora:siteMahana 336.33763  -0.866 0.386921
    ## timepointtimepoint2:speciesPorites:siteMahana     331.79125   0.422 0.673577
    ## timepointtimepoint3:speciesPorites:siteMahana     372.44983   0.830 0.407083
    ## timepointtimepoint4:speciesPorites:siteMahana     335.84991  -1.350 0.177956
    ## timepointtimepoint2:speciesPocillopora:siteManava 326.12421  -0.708 0.479449
    ## timepointtimepoint3:speciesPocillopora:siteManava 366.49750   1.800 0.072747
    ## timepointtimepoint4:speciesPocillopora:siteManava 324.13092   0.064 0.948855
    ## timepointtimepoint2:speciesPorites:siteManava     326.12376  -0.088 0.929597
    ## timepointtimepoint3:speciesPorites:siteManava     365.45049   0.767 0.443626
    ## timepointtimepoint4:speciesPorites:siteManava     326.09768  -0.060 0.951966
    ##                                                      
    ## (Intercept)                                       ***
    ## timepointtimepoint2                                  
    ## timepointtimepoint3                                  
    ## timepointtimepoint4                               ***
    ## speciesPocillopora                                   
    ## speciesPorites                                    ***
    ## siteMahana                                        ** 
    ## siteManava                                           
    ## timepointtimepoint2:speciesPocillopora               
    ## timepointtimepoint3:speciesPocillopora               
    ## timepointtimepoint4:speciesPocillopora            ***
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                   
    ## timepointtimepoint4:speciesPorites                *  
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                       
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## speciesPocillopora:siteMahana                     *  
    ## speciesPorites:siteMahana                         ** 
    ## speciesPocillopora:siteManava                        
    ## speciesPorites:siteManava                            
    ## timepointtimepoint2:speciesPocillopora:siteMahana    
    ## timepointtimepoint3:speciesPocillopora:siteMahana    
    ## timepointtimepoint4:speciesPocillopora:siteMahana    
    ## timepointtimepoint2:speciesPorites:siteMahana        
    ## timepointtimepoint3:speciesPorites:siteMahana        
    ## timepointtimepoint4:speciesPorites:siteMahana        
    ## timepointtimepoint2:speciesPocillopora:siteManava    
    ## timepointtimepoint3:speciesPocillopora:siteManava .  
    ## timepointtimepoint4:speciesPocillopora:siteManava    
    ## timepointtimepoint2:speciesPorites:siteManava        
    ## timepointtimepoint3:speciesPorites:siteManava        
    ## timepointtimepoint4:speciesPorites:siteManava        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, test for effects of site and time within
species.

*Acropora*

`rd_model_acr<-lmer(log(Rd)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))`

``` r
rd_model_acr<-lmer(log(Rd)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))
qqPlot(residuals(rd_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-144-1.png)<!-- -->

    ## 328 334 
    ##  97 103

Generate a Type II Anova of model.

``` r
anova(rd_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      7.8278 2.60927     3 89.862 21.3594 1.566e-10 ***
    ## site           1.8657 0.93287     2 47.361  7.6365  0.001335 ** 
    ## timepoint:site 0.9407 0.15678     6 88.291  1.2834  0.273151    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(rd_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Rd) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Acropora")
    ## 
    ## REML criterion at convergence: 106
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.9653 -0.5064 -0.0597  0.5714  3.6095 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.004358 0.06602 
    ##  Residual                   0.122161 0.34951 
    ## Number of obs: 116, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                     -0.95822    0.09864 103.87514  -9.714 2.97e-16
    ## timepointtimepoint2             -0.19658    0.15251  86.84937  -1.289  0.20082
    ## timepointtimepoint3             -0.11720    0.26866  98.63828  -0.436  0.66361
    ## timepointtimepoint4             -0.76262    0.17357  86.25906  -4.394 3.15e-05
    ## siteMahana                       0.40872    0.13698 103.91454   2.984  0.00355
    ## siteManava                       0.12764    0.13477 103.88215   0.947  0.34581
    ## timepointtimepoint2:siteMahana  -0.02521    0.22964  87.24116  -0.110  0.91283
    ## timepointtimepoint3:siteMahana  -0.38528    0.30584  98.01161  -1.260  0.21075
    ## timepointtimepoint4:siteMahana   0.23428    0.22211  88.16404   1.055  0.29441
    ## timepointtimepoint2:siteManava   0.12092    0.21265  85.73231   0.569  0.57107
    ## timepointtimepoint3:siteManava  -0.12798    0.30682  95.90550  -0.417  0.67753
    ## timepointtimepoint4:siteManava  -0.01849    0.22254  83.77726  -0.083  0.93397
    ##                                   
    ## (Intercept)                    ***
    ## timepointtimepoint2               
    ## timepointtimepoint3               
    ## timepointtimepoint4            ***
    ## siteMahana                     ** 
    ## siteManava                        
    ## timepointtimepoint2:siteMahana    
    ## timepointtimepoint3:siteMahana    
    ## timepointtimepoint4:siteMahana    
    ## timepointtimepoint2:siteManava    
    ## timepointtimepoint3:siteManava    
    ## timepointtimepoint4:siteManava    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.629                                                    
    ## tmpnttmpnt3      -0.355  0.235                                             
    ## tmpnttmpnt4      -0.549  0.357  0.202                                      
    ## siteMahana       -0.720  0.453  0.255  0.395                               
    ## siteManava       -0.732  0.461  0.260  0.402  0.527                        
    ## tmpnttmpnt2:stMh  0.418 -0.664 -0.156 -0.237 -0.579 -0.306                 
    ## tmpnttmpnt3:stMh  0.312 -0.206 -0.878 -0.178 -0.436 -0.228  0.264          
    ## tmpnttmpnt4:stMh  0.429 -0.279 -0.158 -0.781 -0.600 -0.314  0.358          
    ## tmpnttmpnt2:stMn  0.451 -0.717 -0.168 -0.256 -0.325 -0.616  0.476          
    ## tmpnttmpnt3:stMn  0.311 -0.205 -0.876 -0.177 -0.224 -0.425  0.136          
    ## tmpnttmpnt4:stMn  0.428 -0.279 -0.158 -0.780 -0.308 -0.586  0.185          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.272                                            
    ## tmpnttmpnt2:stMn  0.148            0.200                           
    ## tmpnttmpnt3:stMn  0.769            0.138            0.275          
    ## tmpnttmpnt4:stMn  0.138            0.609            0.373          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.258

*Pocillopora*

`rd_model_poc<-lmer(log(Rd)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))`

``` r
rd_model_poc<-lmer(log(Rd)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))
qqPlot(residuals(rd_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-146-1.png)<!-- -->

    ##  56 242 
    ##  16  87

Generate a Type II Anova of model.

``` r
anova(rd_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
    ## timepoint      0.30074 0.10025     3 111.992  1.4328 0.237021   
    ## site           0.56140 0.28070     2  38.691  4.0119 0.026102 * 
    ## timepoint:site 1.26658 0.21110     6 111.919  3.0171 0.009016 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(rd_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Rd) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Pocillopora")
    ## 
    ## REML criterion at convergence: 81.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.95559 -0.50833  0.01336  0.54517  2.14407 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.01440  0.1200  
    ##  Residual                   0.06997  0.2645  
    ## Number of obs: 164, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     -1.009802   0.074995 139.816327 -13.465
    ## timepointtimepoint2              0.038708   0.098512 109.725493   0.393
    ## timepointtimepoint3             -0.094255   0.100688 111.270283  -0.936
    ## timepointtimepoint4             -0.052014   0.096586 108.274545  -0.539
    ## siteMahana                       0.049588   0.109874 141.900036   0.451
    ## siteManava                       0.081905   0.106060 139.816327   0.772
    ## timepointtimepoint2:siteMahana  -0.055678   0.140915 111.081529  -0.395
    ## timepointtimepoint3:siteMahana   0.042084   0.144973 111.707440   0.290
    ## timepointtimepoint4:siteMahana  -0.010120   0.146059 112.561902  -0.069
    ## timepointtimepoint2:siteManava  -0.058869   0.137962 109.013776  -0.427
    ## timepointtimepoint3:siteManava   0.455485   0.144204 112.559810   3.159
    ## timepointtimepoint4:siteManava   0.004465   0.139571 110.312813   0.032
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## timepointtimepoint2             0.69513    
    ## timepointtimepoint3             0.35124    
    ## timepointtimepoint4             0.59132    
    ## siteMahana                      0.65245    
    ## siteManava                      0.44127    
    ## timepointtimepoint2:siteMahana  0.69352    
    ## timepointtimepoint3:siteMahana  0.77214    
    ## timepointtimepoint4:siteMahana  0.94488    
    ## timepointtimepoint2:siteManava  0.67044    
    ## timepointtimepoint3:siteManava  0.00204 ** 
    ## timepointtimepoint4:siteManava  0.97454    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.631                                                    
    ## tmpnttmpnt3      -0.618  0.470                                             
    ## tmpnttmpnt4      -0.644  0.490  0.480                                      
    ## siteMahana       -0.683  0.431  0.422  0.440                               
    ## siteManava       -0.707  0.446  0.437  0.455  0.483                        
    ## tmpnttmpnt2:stMh  0.441 -0.699 -0.328 -0.343 -0.656 -0.312                 
    ## tmpnttmpnt3:stMh  0.429 -0.326 -0.695 -0.333 -0.634 -0.303  0.494          
    ## tmpnttmpnt4:stMh  0.426 -0.324 -0.317 -0.661 -0.630 -0.301  0.491          
    ## tmpnttmpnt2:stMn  0.451 -0.714 -0.335 -0.350 -0.308 -0.638  0.499          
    ## tmpnttmpnt3:stMn  0.431 -0.328 -0.698 -0.335 -0.294 -0.610  0.229          
    ## tmpnttmpnt4:stMn  0.446 -0.339 -0.332 -0.692 -0.304 -0.630  0.237          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.480                                            
    ## tmpnttmpnt2:stMn  0.233            0.231                           
    ## tmpnttmpnt3:stMn  0.485            0.221            0.469          
    ## tmpnttmpnt4:stMn  0.231            0.458            0.484          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.469

*Porites*

`rd_model_por<-lmer(log(Rd)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))`

``` r
rd_model_por<-lmer(log(Rd)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))
qqPlot(residuals(rd_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-148-1.png)<!-- -->

    ## 383 282 
    ## 118  84

Generate a Type II Anova of model.

``` r
anova(rd_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      2.95880 0.98627     3   155  9.3459 1.024e-05 ***
    ## site           0.87634 0.43817     2   155  4.1521   0.01752 *  
    ## timepoint:site 0.36391 0.06065     6   155  0.5747   0.75003    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(rd_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(Rd) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Porites")
    ## 
    ## REML criterion at convergence: 122.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.7767 -0.5007  0.0491  0.5755  2.1571 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance Std.Dev.
    ##  colony_id_corr (Intercept) 0.0000   0.0000  
    ##  Residual                   0.1055   0.3249  
    ## Number of obs: 167, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                     -0.15037    0.08388 155.00000  -1.793   0.0750
    ## timepointtimepoint2             -0.13830    0.12072 155.00000  -1.146   0.2537
    ## timepointtimepoint3             -0.19828    0.11862 155.00000  -1.672   0.0966
    ## timepointtimepoint4             -0.30554    0.12072 155.00000  -2.531   0.0124
    ## siteMahana                      -0.04292    0.11862 155.00000  -0.362   0.7180
    ## siteManava                       0.05463    0.11862 155.00000   0.461   0.6458
    ## timepointtimepoint2:siteMahana   0.08750    0.16924 155.00000   0.517   0.6059
    ## timepointtimepoint3:siteMahana  -0.11563    0.17292 155.00000  -0.669   0.5047
    ## timepointtimepoint4:siteMahana  -0.11827    0.17436 155.00000  -0.678   0.4986
    ## timepointtimepoint2:siteManava   0.09986    0.16924 155.00000   0.590   0.5560
    ## timepointtimepoint3:siteManava   0.11824    0.17095 155.00000   0.692   0.4902
    ## timepointtimepoint4:siteManava  -0.03542    0.17436 155.00000  -0.203   0.8393
    ##                                 
    ## (Intercept)                    .
    ## timepointtimepoint2             
    ## timepointtimepoint3            .
    ## timepointtimepoint4            *
    ## siteMahana                      
    ## siteManava                      
    ## timepointtimepoint2:siteMahana  
    ## timepointtimepoint3:siteMahana  
    ## timepointtimepoint4:siteMahana  
    ## timepointtimepoint2:siteManava  
    ## timepointtimepoint3:siteManava  
    ## timepointtimepoint4:siteManava  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.695                                                    
    ## tmpnttmpnt3      -0.707  0.491                                             
    ## tmpnttmpnt4      -0.695  0.483  0.491                                      
    ## siteMahana       -0.707  0.491  0.500  0.491                               
    ## siteManava       -0.707  0.491  0.500  0.491  0.500                        
    ## tmpnttmpnt2:stMh  0.496 -0.713 -0.350 -0.344 -0.701 -0.350                 
    ## tmpnttmpnt3:stMh  0.485 -0.337 -0.686 -0.337 -0.686 -0.343  0.481          
    ## tmpnttmpnt4:stMh  0.481 -0.334 -0.340 -0.692 -0.680 -0.340  0.477          
    ## tmpnttmpnt2:stMn  0.496 -0.713 -0.350 -0.344 -0.350 -0.701  0.509          
    ## tmpnttmpnt3:stMn  0.491 -0.341 -0.694 -0.341 -0.347 -0.694  0.243          
    ## tmpnttmpnt4:stMn  0.481 -0.334 -0.340 -0.692 -0.340 -0.680  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.467                                            
    ## tmpnttmpnt2:stMn  0.240            0.238                           
    ## tmpnttmpnt3:stMn  0.476            0.236            0.486          
    ## tmpnttmpnt4:stMn  0.233            0.479            0.477          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.472          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

## Calcification

### Plot: Calcification mgprot

View by site and timepoint for each species.

``` r
calcplot_prot<-master %>%
  filter(!is.na(calc.umol.mgprot.hr)) %>% 
  filter(calc.umol.mgprot.hr<3)%>%
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, calc.umol.mgprot.hr)%>%
  
  ggplot(., aes(x = site_code, y = calc.umol.mgprot.hr, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    ylab(expression(bold(paste(mu,"mol Ca", CO[3], " mg protein"^-1, " hr"^-1))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: Calcification mgafdw

View by site and timepoint for each species.

``` r
calcplot_afdw<-master %>%
  filter(!is.na(calc.umol.mgAFDW.hr)) %>% 
  filter(calc.umol.mgAFDW.hr>-.06)%>%
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, calc.umol.mgAFDW.hr)%>%
  
  ggplot(., aes(x = site_code, y = calc.umol.mgAFDW.hr, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    ylab(expression(bold(paste(mu,"mol Ca", CO[3], " mg AFDW"^-1, " hr"^-1))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

### Plot: CRE umol cm2

View by site and timepoint for each species.

``` r
calcplot_cm2<-master %>%
  filter(!is.na(calc.umol.cm2.hr)) %>% 
  filter(!is.na(species))%>%
  select(colony_id_corr, species, site_code, timepoint, calc.umol.cm2.hr)%>%
  
  ggplot(., aes(x = site_code, y = calc.umol.cm2.hr, fill = timepoint, group=interaction(site_code,timepoint))) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    scale_fill_manual(values = brewer.pal(4,"YlGnBu"), limits=c("timepoint1", "timepoint2", "timepoint3", "timepoint4"), labels=c("January 2020", "March 2020", "September 2020", "November 2020"))+
    labs(fill = "Month") +
    facet_wrap(~species) +
    xlab("Site") + 
    scale_x_discrete(labels=c("Mahana Low" = "Mahana \nLow", "Hilton Medium" = "Hilton \nMedium",
                              "Manava High" = "Manava \nHigh"))+
    ylab(expression(bold(paste(mu,"mol Ca", CO[3], " cm"^-2, " hr"^-1))))+
    theme_classic() + 
    theme(
      legend.position="none",
      legend.title=element_text(face="bold", size=14),
      legend.text=element_text(size=14),
      axis.title=element_text(face="bold", size=14),
      axis.text=element_text(size=10, color="black"), 
      strip.text.x=element_text(face="italic", size=14)
      )
```

Join all plots for each normalization together.

``` r
calc_figure<-plot_grid(calcplot_prot, calcplot_afdw, calcplot_cm2, labels = c("", "", ""), ncol=3, nrow=1, rel_widths = c(.8, .8, 1), label_size = 20, label_x = 0.1);calc_figure
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-153-1.png)<!-- -->

``` r
ggsave(filename="Figures/Calcification_Figure.pdf", plot=calc_figure, dpi=300, width=24, height=5, units="in")
```

### Analysis

Build a mixed model for univariate analysis and examine data
distribution of biomass normalized antioxidant capacity.

`calc_model<-lmer(calc.umol.mgAFDW.hr~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
calc_model<-lmer(calc.umol.mgAFDW.hr~timepoint*species*site+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(calc_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-154-1.png)<!-- -->

    ## 218 170 
    ## 212 166

``` r
#hist(residuals(antiox_model))
```

Residuals are not normally distributed. Attempt with log transformation.

`calc_model<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master)`

``` r
calc_model<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master)
qqPlot(residuals(calc_model))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-155-1.png)<!-- -->

    ## 170 218 
    ## 166 212

``` r
#hist(residuals(antiox_model))
```

We are going to need to find a stronger transformation for this data.

Generate a Type III Anova of model.

``` r
anova(calc_model, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint              1.90201 0.63400     3 320.51 90.2951 < 2.2e-16 ***
    ## site                   0.06561 0.03280     2 138.04  4.6719   0.01088 *  
    ## species                0.46485 0.23242     2 129.82 33.1017 2.417e-12 ***
    ## timepoint:site         0.11331 0.01889     6 318.39  2.6897   0.01464 *  
    ## timepoint:species      0.55302 0.09217     6 313.80 13.1269 2.924e-13 ***
    ## site:species           0.03745 0.00936     4 127.51  1.3336   0.26104    
    ## timepoint:site:species 0.15927 0.01327    12 312.08  1.8903   0.03487 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(calc_model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + calc.umol.mgAFDW.hr) ~ timepoint * site * species + (1 |  
    ##     colony_id_corr)
    ##    Data: master
    ## 
    ## REML criterion at convergence: -729.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6304 -0.2086 -0.0090  0.0785  7.1059 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0002821 0.01680 
    ##  Residual                   0.0070215 0.08379 
    ## Number of obs: 428, groups:  colony_id_corr, 140
    ## 
    ## Fixed effects:
    ##                                                     Estimate Std. Error
    ## (Intercept)                                         0.037731   0.023701
    ## timepointtimepoint2                                 0.162731   0.036590
    ## timepointtimepoint3                                 0.204826   0.064498
    ## timepointtimepoint4                                -0.030196   0.041642
    ## siteMahana                                         -0.022467   0.033514
    ## siteManava                                         -0.033458   0.032381
    ## speciesPocillopora                                 -0.020215   0.032382
    ## speciesPorites                                     -0.035678   0.032913
    ## timepointtimepoint2:siteMahana                     -0.001184   0.055435
    ## timepointtimepoint3:siteMahana                      0.266359   0.074909
    ## timepointtimepoint4:siteMahana                      0.022217   0.053685
    ## timepointtimepoint2:siteManava                     -0.033590   0.052002
    ## timepointtimepoint3:siteManava                      0.106507   0.073647
    ## timepointtimepoint4:siteManava                      0.029312   0.053386
    ## timepointtimepoint2:speciesPocillopora             -0.126400   0.048059
    ## timepointtimepoint3:speciesPocillopora             -0.153445   0.073066
    ## timepointtimepoint4:speciesPocillopora              0.013907   0.051675
    ## timepointtimepoint2:speciesPorites                 -0.049504   0.049311
    ## timepointtimepoint3:speciesPorites                 -0.140446   0.071630
    ## timepointtimepoint4:speciesPorites                  0.030239   0.052725
    ## siteMahana:speciesPocillopora                       0.013999   0.046602
    ## siteManava:speciesPocillopora                       0.017849   0.044970
    ## siteMahana:speciesPorites                           0.027303   0.046972
    ## siteManava:speciesPorites                           0.033236   0.045354
    ## timepointtimepoint2:siteMahana:speciesPocillopora   0.035345   0.071097
    ## timepointtimepoint3:siteMahana:speciesPocillopora  -0.249418   0.089306
    ## timepointtimepoint4:siteMahana:speciesPocillopora  -0.013381   0.070345
    ## timepointtimepoint2:siteManava:speciesPocillopora   0.113227   0.067906
    ## timepointtimepoint3:siteManava:speciesPocillopora  -0.109835   0.087835
    ## timepointtimepoint4:siteManava:speciesPocillopora  -0.014221   0.069261
    ## timepointtimepoint2:siteMahana:speciesPorites       0.049405   0.071948
    ## timepointtimepoint3:siteMahana:speciesPorites      -0.178236   0.088970
    ## timepointtimepoint4:siteMahana:speciesPorites      -0.025094   0.071120
    ## timepointtimepoint2:siteManava:speciesPorites       0.004716   0.069049
    ## timepointtimepoint3:siteManava:speciesPorites      -0.034257   0.086055
    ## timepointtimepoint4:siteManava:speciesPorites      -0.028824   0.070380
    ##                                                           df t value Pr(>|t|)
    ## (Intercept)                                       390.905963   1.592 0.112198
    ## timepointtimepoint2                               336.840998   4.447 1.18e-05
    ## timepointtimepoint3                               375.492879   3.176 0.001618
    ## timepointtimepoint4                               334.482902  -0.725 0.468876
    ## siteMahana                                        391.161649  -0.670 0.503017
    ## siteManava                                        390.946395  -1.033 0.302121
    ## speciesPocillopora                                390.781826  -0.624 0.532824
    ## speciesPorites                                    390.965627  -1.084 0.279027
    ## timepointtimepoint2:siteMahana                    335.655328  -0.021 0.982971
    ## timepointtimepoint3:siteMahana                    374.164609   3.556 0.000425
    ## timepointtimepoint4:siteMahana                    345.106923   0.414 0.679244
    ## timepointtimepoint2:siteManava                    336.721180  -0.646 0.518764
    ## timepointtimepoint3:siteManava                    366.619658   1.446 0.148977
    ## timepointtimepoint4:siteManava                    326.013793   0.549 0.583337
    ## timepointtimepoint2:speciesPocillopora            313.751532  -2.630 0.008958
    ## timepointtimepoint3:speciesPocillopora            363.289913  -2.100 0.036410
    ## timepointtimepoint4:speciesPocillopora            314.277018   0.269 0.788008
    ## timepointtimepoint2:speciesPorites                320.018093  -1.004 0.316172
    ## timepointtimepoint3:speciesPorites                361.742521  -1.961 0.050678
    ## timepointtimepoint4:speciesPorites                318.925087   0.574 0.566695
    ## siteMahana:speciesPocillopora                     391.057251   0.300 0.764029
    ## siteManava:speciesPocillopora                     390.799357   0.397 0.691651
    ## siteMahana:speciesPorites                         391.133628   0.581 0.561394
    ## siteManava:speciesPorites                         390.897006   0.733 0.464113
    ## timepointtimepoint2:siteMahana:speciesPocillopora 315.818343   0.497 0.619438
    ## timepointtimepoint3:siteMahana:speciesPocillopora 356.772198  -2.793 0.005506
    ## timepointtimepoint4:siteMahana:speciesPocillopora 321.151490  -0.190 0.849251
    ## timepointtimepoint2:siteManava:speciesPocillopora 313.025073   1.667 0.096434
    ## timepointtimepoint3:siteManava:speciesPocillopora 349.905383  -1.250 0.211964
    ## timepointtimepoint4:siteManava:speciesPocillopora 308.120241  -0.205 0.837460
    ## timepointtimepoint2:siteMahana:speciesPorites     318.713536   0.687 0.492787
    ## timepointtimepoint3:siteMahana:speciesPorites     356.819004  -2.003 0.045898
    ## timepointtimepoint4:siteMahana:speciesPorites     323.642507  -0.353 0.724440
    ## timepointtimepoint2:siteManava:speciesPorites     317.015621   0.068 0.945596
    ## timepointtimepoint3:siteManava:speciesPorites     347.484996  -0.398 0.690817
    ## timepointtimepoint4:siteManava:speciesPorites     311.939107  -0.410 0.682417
    ##                                                      
    ## (Intercept)                                          
    ## timepointtimepoint2                               ***
    ## timepointtimepoint3                               ** 
    ## timepointtimepoint4                                  
    ## siteMahana                                           
    ## siteManava                                           
    ## speciesPocillopora                                   
    ## speciesPorites                                       
    ## timepointtimepoint2:siteMahana                       
    ## timepointtimepoint3:siteMahana                    ***
    ## timepointtimepoint4:siteMahana                       
    ## timepointtimepoint2:siteManava                       
    ## timepointtimepoint3:siteManava                       
    ## timepointtimepoint4:siteManava                       
    ## timepointtimepoint2:speciesPocillopora            ** 
    ## timepointtimepoint3:speciesPocillopora            *  
    ## timepointtimepoint4:speciesPocillopora               
    ## timepointtimepoint2:speciesPorites                   
    ## timepointtimepoint3:speciesPorites                .  
    ## timepointtimepoint4:speciesPorites                   
    ## siteMahana:speciesPocillopora                        
    ## siteManava:speciesPocillopora                        
    ## siteMahana:speciesPorites                            
    ## siteManava:speciesPorites                            
    ## timepointtimepoint2:siteMahana:speciesPocillopora    
    ## timepointtimepoint3:siteMahana:speciesPocillopora ** 
    ## timepointtimepoint4:siteMahana:speciesPocillopora    
    ## timepointtimepoint2:siteManava:speciesPocillopora .  
    ## timepointtimepoint3:siteManava:speciesPocillopora    
    ## timepointtimepoint4:siteManava:speciesPocillopora    
    ## timepointtimepoint2:siteMahana:speciesPorites        
    ## timepointtimepoint3:siteMahana:speciesPorites     *  
    ## timepointtimepoint4:siteMahana:speciesPorites        
    ## timepointtimepoint2:siteManava:speciesPorites        
    ## timepointtimepoint3:siteManava:speciesPorites        
    ## timepointtimepoint4:siteManava:speciesPorites        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Since species is significant, run a model for site and time effects
within each species.

*Acropora*:

`calc_model_acr<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))`

``` r
calc_model_acr<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Acropora"))
qqPlot(residuals(calc_model_acr))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-157-1.png)<!-- -->

    ## 218 121 
    ##  64  47

Generate a Type III Anova of model.

``` r
anova(calc_model_acr, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## timepoint      1.31467 0.43822     3 88.838 29.2576 3.002e-13 ***
    ## site           0.06002 0.03001     2 53.378  2.0036    0.1449    
    ## timepoint:site 0.13953 0.02326     6 87.887  1.5526    0.1706    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(calc_model_acr)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + calc.umol.mgAFDW.hr) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Acropora")
    ## 
    ## REML criterion at convergence: -106.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4743 -0.2295 -0.0187  0.0827  4.3222 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0006834 0.02614 
    ##  Residual                   0.0149781 0.12239 
    ## Number of obs: 112, groups:  colony_id_corr, 50
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)                     0.037870   0.034705 99.861056   1.091  0.27781
    ## timepointtimepoint2             0.162513   0.053490 86.252312   3.038  0.00315
    ## timepointtimepoint3             0.205209   0.094353 95.355211   2.175  0.03211
    ## timepointtimepoint4            -0.030063   0.060873 85.595399  -0.494  0.62267
    ## siteMahana                     -0.022689   0.049075 99.908430  -0.462  0.64485
    ## siteManava                     -0.033521   0.047416 99.868425  -0.707  0.48124
    ## timepointtimepoint2:siteMahana -0.001063   0.081036 85.920816  -0.013  0.98957
    ## timepointtimepoint3:siteMahana  0.266182   0.109581 95.072027   2.429  0.01702
    ## timepointtimepoint4:siteMahana  0.022213   0.078490 88.186638   0.283  0.77783
    ## timepointtimepoint2:siteManava -0.033700   0.076020 86.156307  -0.443  0.65866
    ## timepointtimepoint3:siteManava  0.106013   0.107716 93.196233   0.984  0.32757
    ## timepointtimepoint4:siteManava  0.028936   0.078029 83.569842   0.371  0.71170
    ##                                  
    ## (Intercept)                      
    ## timepointtimepoint2            **
    ## timepointtimepoint3            * 
    ## timepointtimepoint4              
    ## siteMahana                       
    ## siteManava                       
    ## timepointtimepoint2:siteMahana   
    ## timepointtimepoint3:siteMahana * 
    ## timepointtimepoint4:siteMahana   
    ## timepointtimepoint2:siteManava   
    ## timepointtimepoint3:siteManava   
    ## timepointtimepoint4:siteManava   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.627                                                    
    ## tmpnttmpnt3      -0.352  0.235                                             
    ## tmpnttmpnt4      -0.545  0.357  0.201                                      
    ## siteMahana       -0.707  0.443  0.249  0.386                               
    ## siteManava       -0.732  0.459  0.258  0.399  0.518                        
    ## tmpnttmpnt2:stMh  0.414 -0.660 -0.155 -0.235 -0.582 -0.303                 
    ## tmpnttmpnt3:stMh  0.303 -0.202 -0.861 -0.173 -0.433 -0.222  0.267          
    ## tmpnttmpnt4:stMh  0.423 -0.277 -0.156 -0.776 -0.605 -0.310  0.366          
    ## tmpnttmpnt2:stMn  0.441 -0.704 -0.165 -0.251 -0.312 -0.601  0.464          
    ## tmpnttmpnt3:stMn  0.308 -0.206 -0.876 -0.176 -0.218 -0.422  0.136          
    ## tmpnttmpnt4:stMn  0.425 -0.278 -0.157 -0.780 -0.301 -0.582  0.184          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.274                                            
    ## tmpnttmpnt2:stMn  0.142            0.195                           
    ## tmpnttmpnt3:stMn  0.754            0.137            0.270          
    ## tmpnttmpnt4:stMn  0.135            0.605            0.366          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.258

*Pocillopora*:

`calc_model_poc<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))`

``` r
calc_model_poc<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Pocillopora"))
qqPlot(residuals(calc_model_poc))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-159-1.png)<!-- -->

    ## 170 172 
    ##  75  77

``` r
#hist(residuals(calc_model_poc))
```

Generate a Type III Anova of model.

``` r
anova(calc_model_poc, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## timepoint      0.203312 0.067771     3 115.568 13.0054 2.229e-07 ***
    ## site           0.001336 0.000668     2  44.859  0.1282    0.8800    
    ## timepoint:site 0.034105 0.005684     6 115.489  1.0908    0.3721    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(calc_model_poc)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + calc.umol.mgAFDW.hr) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Pocillopora")
    ## 
    ## REML criterion at convergence: -316.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5375 -0.2376 -0.0203  0.0845  8.0937 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev.
    ##  colony_id_corr (Intercept) 0.0003266 0.01807 
    ##  Residual                   0.0052110 0.07219 
    ## Number of obs: 159, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      0.017516   0.019214 145.822760   0.912
    ## timepointtimepoint2              0.036226   0.026850 111.844931   1.349
    ## timepointtimepoint3              0.051291   0.029626 119.650346   1.731
    ## timepointtimepoint4             -0.016289   0.026359 110.186023  -0.618
    ## siteMahana                      -0.008405   0.028192 146.059494  -0.298
    ## siteManava                      -0.015608   0.027173 145.822760  -0.574
    ## timepointtimepoint2:siteMahana   0.034203   0.038369 112.945815   0.891
    ## timepointtimepoint3:siteMahana   0.016640   0.041953 119.018077   0.397
    ## timepointtimepoint4:siteMahana   0.008677   0.039183 113.473039   0.221
    ## timepointtimepoint2:siteManava   0.079742   0.037626 111.030917   2.119
    ## timepointtimepoint3:siteManava  -0.003571   0.041299 118.767565  -0.086
    ## timepointtimepoint4:siteManava   0.014957   0.038027 112.163131   0.393
    ##                                Pr(>|t|)  
    ## (Intercept)                      0.3635  
    ## timepointtimepoint2              0.1800  
    ## timepointtimepoint3              0.0860 .
    ## timepointtimepoint4              0.5379  
    ## siteMahana                       0.7660  
    ## siteManava                       0.5666  
    ## timepointtimepoint2:siteMahana   0.3746  
    ## timepointtimepoint3:siteMahana   0.6924  
    ## timepointtimepoint4:siteMahana   0.8251  
    ## timepointtimepoint2:siteManava   0.0363 *
    ## timepointtimepoint3:siteManava   0.9312  
    ## timepointtimepoint4:siteManava   0.6948  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.673                                                    
    ## tmpnttmpnt3      -0.610  0.436                                             
    ## tmpnttmpnt4      -0.686  0.491  0.445                                      
    ## siteMahana       -0.682  0.459  0.416  0.467                               
    ## siteManava       -0.707  0.476  0.432  0.485  0.482                        
    ## tmpnttmpnt2:stMh  0.471 -0.700 -0.305 -0.343 -0.694 -0.333                 
    ## tmpnttmpnt3:stMh  0.431 -0.308 -0.706 -0.314 -0.634 -0.305  0.466          
    ## tmpnttmpnt4:stMh  0.461 -0.330 -0.299 -0.673 -0.679 -0.326  0.499          
    ## tmpnttmpnt2:stMn  0.481 -0.714 -0.311 -0.350 -0.328 -0.680  0.499          
    ## tmpnttmpnt3:stMn  0.438 -0.313 -0.717 -0.319 -0.298 -0.619  0.219          
    ## tmpnttmpnt4:stMn  0.475 -0.340 -0.308 -0.693 -0.324 -0.672  0.238          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.457                                            
    ## tmpnttmpnt2:stMn  0.220            0.236                           
    ## tmpnttmpnt3:stMn  0.507            0.215            0.447          
    ## tmpnttmpnt4:stMn  0.218            0.466            0.486          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.444

*Porites*:

`calc_model_por<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site*species+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))`

``` r
calc_model_por<-lmer(log(1+calc.umol.mgAFDW.hr)~timepoint*site+(1|colony_id_corr), na.action=na.omit, data=master, subset=c(species=="Porites"))
qqPlot(residuals(calc_model_por))
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-161-1.png)<!-- -->

    ## 210 433 
    ##  66 144

``` r
#hist(residuals(calc_model_por))
```

Generate a Type III Anova of model.

``` r
anova(calc_model_por, type="III")
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
    ## timepoint      0.55065 0.183548     3   145 55.3477 < 2.2e-16 ***
    ## site           0.03942 0.019712     2   145  5.9439  0.003303 ** 
    ## timepoint:site 0.06985 0.011641     6   145  3.5103  0.002847 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(calc_model_por)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: log(1 + calc.umol.mgAFDW.hr) ~ timepoint * site + (1 | colony_id_corr)
    ##    Data: master
    ##  Subset: c(species == "Porites")
    ## 
    ## REML criterion at convergence: -385.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9232 -0.1342 -0.0035  0.0983  3.6242 
    ## 
    ## Random effects:
    ##  Groups         Name        Variance  Std.Dev. 
    ##  colony_id_corr (Intercept) 1.002e-22 1.001e-11
    ##  Residual                   3.316e-03 5.759e-02
    ## Number of obs: 157, groups:  colony_id_corr, 45
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     1.572e-03  1.539e-02  1.450e+02   0.102
    ## timepointtimepoint2             1.139e-01  2.265e-02  1.450e+02   5.026
    ## timepointtimepoint3             6.486e-02  2.140e-02  1.450e+02   3.031
    ## timepointtimepoint4             5.807e-04  2.218e-02  1.450e+02   0.026
    ## siteMahana                      6.155e-03  2.218e-02  1.450e+02   0.277
    ## siteManava                      2.587e-04  2.140e-02  1.450e+02   0.012
    ## timepointtimepoint2:siteMahana  4.676e-02  3.145e-02  1.450e+02   1.486
    ## timepointtimepoint3:siteMahana  8.758e-02  3.289e-02  1.450e+02   2.663
    ## timepointtimepoint4:siteMahana -4.412e-03  3.199e-02  1.450e+02  -0.138
    ## timepointtimepoint2:siteManava -2.957e-02  3.116e-02  1.450e+02  -0.949
    ## timepointtimepoint3:siteManava  7.183e-02  3.056e-02  1.450e+02   2.350
    ## timepointtimepoint4:siteManava -2.635e-04  3.145e-02  1.450e+02  -0.008
    ##                                Pr(>|t|)    
    ## (Intercept)                     0.91878    
    ## timepointtimepoint2            1.46e-06 ***
    ## timepointtimepoint3             0.00289 ** 
    ## timepointtimepoint4             0.97915    
    ## siteMahana                      0.78180    
    ## siteManava                      0.99037    
    ## timepointtimepoint2:siteMahana  0.13934    
    ## timepointtimepoint3:siteMahana  0.00862 ** 
    ## timepointtimepoint4:siteMahana  0.89051    
    ## timepointtimepoint2:siteManava  0.34422    
    ## timepointtimepoint3:siteManava  0.02012 *  
    ## timepointtimepoint4:siteManava  0.99333    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##                  (Intr) tmpnt2 tmpnt3 tmpnt4 sitMhn sitMnv tmpnttmpnt2:stMh
    ## tmpnttmpnt2      -0.679                                                    
    ## tmpnttmpnt3      -0.719  0.489                                             
    ## tmpnttmpnt4      -0.694  0.471  0.499                                      
    ## siteMahana       -0.694  0.471  0.499  0.481                               
    ## siteManava       -0.719  0.489  0.517  0.499  0.499                        
    ## tmpnttmpnt2:stMh  0.489 -0.720 -0.352 -0.340 -0.705 -0.352                 
    ## tmpnttmpnt3:stMh  0.468 -0.318 -0.651 -0.325 -0.674 -0.337  0.476          
    ## tmpnttmpnt4:stMh  0.481 -0.327 -0.346 -0.693 -0.693 -0.346  0.489          
    ## tmpnttmpnt2:stMn  0.494 -0.727 -0.355 -0.343 -0.343 -0.687  0.524          
    ## tmpnttmpnt3:stMn  0.504 -0.342 -0.700 -0.349 -0.349 -0.700  0.246          
    ## tmpnttmpnt4:stMn  0.489 -0.332 -0.352 -0.705 -0.340 -0.680  0.239          
    ##                  tmpnttmpnt3:stMh tmpnttmpnt4:stMh tmpnttmpnt2:stMn
    ## tmpnttmpnt2                                                        
    ## tmpnttmpnt3                                                        
    ## tmpnttmpnt4                                                        
    ## siteMahana                                                         
    ## siteManava                                                         
    ## tmpnttmpnt2:stMh                                                   
    ## tmpnttmpnt3:stMh                                                   
    ## tmpnttmpnt4:stMh  0.468                                            
    ## tmpnttmpnt2:stMn  0.231            0.238                           
    ## tmpnttmpnt3:stMn  0.456            0.242            0.481          
    ## tmpnttmpnt4:stMn  0.229            0.489            0.467          
    ##                  tmpnttmpnt3:stMn
    ## tmpnttmpnt2                      
    ## tmpnttmpnt3                      
    ## tmpnttmpnt4                      
    ## siteMahana                       
    ## siteManava                       
    ## tmpnttmpnt2:stMh                 
    ## tmpnttmpnt3:stMh                 
    ## tmpnttmpnt4:stMh                 
    ## tmpnttmpnt2:stMn                 
    ## tmpnttmpnt3:stMn                 
    ## tmpnttmpnt4:stMn  0.476          
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

## Generating univariate response panels

First, specify symbiont and holobiont responses for plotting in separate
panels.

``` r
symbiont_responses<-c("cells.mgAFDW", "Sym_AFDW.mg.cm2", "Total_Chl", "Total_Chl_cell", "Am", "AQY", "Ratio_AFDW.mg.cm2") 
holobiont_responses<-c("cre.umol.mgafdw", "Host_AFDW.mg.cm2", "prot_mg.mgafdw", "Rd", "calc.umol.mgAFDW.hr")
```

Generate a panel of responses for symbiont metrics.

``` r
#build a plot for the legend
legend_plot<-chlplot_total_cell+theme(legend.position="bottom")

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  legend_plot + theme(legend.box.margin = margin(1,1,1,1))
)

symbiont_univ<-plot_grid(symbplot_afdw, sAFDWplot, am_plot, aqy_plot, chlplot_total, chlplot_total_cell, ncol=3, nrow=2, rel_widths = c(1,1), align="h")

symbiont_univ2<-plot_grid(symbiont_univ, legend, ncol=1, nrow=2, rel_heights = c(4, .2))

ggsave(filename="Figures/Symbiont_Univariate_Panel.png", plot=symbiont_univ2, dpi=150, width=30, height=11, units="in")
```

Generate a panel of responses for host/holobiont metrics.

``` r
#build a plot for the legend
legend_plot<-calcplot_afdw+theme(legend.position="bottom")
rd_plot<-rd_plot+theme(legend.position="none")
calcplot_afdw<-calcplot_afdw+theme(legend.position="none") 

holobiont_univ<-plot_grid(tacplot_afdw, hAFDWplot, ratioAFDWplot, protplot_afdw, rd_plot, calcplot_afdw, ncol=3, nrow=2, rel_widths = c(1, 1), align="h")

holobiont_univ2<-plot_grid(holobiont_univ, legend, ncol=1, nrow=2, rel_heights = c(4, .2))

ggsave(filename="Figures/Holobiont_Univariate_Panel.png", plot=holobiont_univ2, dpi=150, width=30, height=11, units="in")
```

# Correlations

Generate a correlation matrix to examine relationships between responses
of interest. Use normalizations by afdw for all responses.

Generating network plot with corrr package.

``` r
master%>%
  select(cre.umol.mgafdw, prot_mg.mgafdw, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chlc2.ug.mgAFDW, chla.ug.mgAFDW, cells.mgAFDW, chla.ug.cell, chlc2.ug.cell, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  corrr::correlate(., method="spearman", diagonal=1)%>% #generate correlation matrix
  rearrange(.)%>%
  #shave(.)%>% #remove upper triangle
  network_plot(min_cor=.2)
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-166-1.png)<!-- -->

Generating matrix with ggcorrplot for all colonies of all species.

``` r
p<-master%>%
  select(species, cre.umol.mgafdw, prot_mg.mgafdw, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chlc2.ug.mgAFDW, chla.ug.mgAFDW, cells.mgAFDW, chla.ug.cell, chlc2.ug.cell, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(!is.na(species))%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(prot_mg.mgafdw<1.5)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(calc.umol.mgAFDW.hr>-.06)%>%
  ggpairs(., columns=2:15, ggplot2::aes(colour=species)) #try to bold significant stats

#change colors
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("orange", "gray", "purple")) +
        scale_color_manual(values=c("orange", "gray", "purple")) +
        theme_classic()+
        theme(text=element_text(size=14, face="bold"))+
        theme(strip.text.x=element_text(face="bold", size=14))
  }
}

p
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-167-1.png)<!-- -->

``` r
ggsave(plot=p, "Figures/Correlation_Matrix.pdf", width=20, height=12, unit="in", dpi=300)
```

Acropora:

``` r
p<-master%>%
  filter(species=="Acropora")%>%
  select(site_code, cre.umol.mgafdw, prot_mg.mgafdw, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chlc2.ug.mgAFDW, chla.ug.mgAFDW, cells.mgAFDW, chla.ug.cell, chlc2.ug.cell, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(calc.umol.mgAFDW.hr>-.06)%>%
  ggpairs(., columns=2:15, ggplot2::aes(colour=site_code)) 

#change colors
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
        scale_color_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
        theme_classic()+
        theme(text=element_text(size=14, face="bold"))+
        theme(strip.text.x=element_text(face="bold", size=14))
  }
}

p
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-168-1.png)<!-- -->

``` r
ggsave(plot=p, "Figures/Correlation_Matrix_Acropora.pdf", width=20, height=12, unit="in", dpi=300)
```

Pocillopora

``` r
p<-master%>%
  filter(species=="Pocillopora")%>%
  select(site_code, cre.umol.mgafdw, prot_mg.mgafdw, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chlc2.ug.mgAFDW, chla.ug.mgAFDW, cells.mgAFDW, chla.ug.cell, chlc2.ug.cell, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(calc.umol.mgAFDW.hr>-.06)%>%
  ggpairs(., columns=2:15, ggplot2::aes(colour=site_code)) 

#change colors
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
        scale_color_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
        theme_classic()+
        theme(text=element_text(size=14, face="bold"))+
        theme(strip.text.x=element_text(face="bold", size=14))
  }
}

p
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-169-1.png)<!-- -->

``` r
ggsave(plot=p, "Figures/Correlation_Matrix_Pocillopora.pdf", width=20, height=12, unit="in", dpi=300)
```

Porites

``` r
p<-master%>%
  filter(species=="Porites")%>%
  select(site_code, cre.umol.mgafdw, prot_mg.mgafdw, Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Ratio_AFDW.mg.cm2, chlc2.ug.mgAFDW, chla.ug.mgAFDW, cells.mgAFDW, chla.ug.cell, chlc2.ug.cell, Am, AQY, Rd, calc.umol.mgAFDW.hr)%>%
  filter(Host_AFDW.mg.cm2>0)%>%
  filter(Sym_AFDW.mg.cm2>0)%>%
  filter(Ratio_AFDW.mg.cm2>0)%>%
  filter(prot_mg.mgafdw<1.5)%>%
  filter(Sym_AFDW.mg.cm2<20)%>%
  filter(calc.umol.mgAFDW.hr>-.06)%>%
  ggpairs(., columns=2:15, ggplot2::aes(colour=site_code)) 

#change colors
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
        scale_color_manual(values=c("#374d7c", "#00cccc", "#ff6633")) +
        theme_classic()+
        theme(text=element_text(size=14, face="bold"))+
        theme(strip.text.x=element_text(face="bold", size=14))
  }
}

p
```

![](2_univariate_analyses_files/figure-gfm/unnamed-chunk-170-1.png)<!-- -->

``` r
ggsave(plot=p, "Figures/Correlation_Matrix_Porites.pdf", width=20, height=12, unit="in", dpi=300)
```
