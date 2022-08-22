Modeling analysis of E5 time series biological data
================
Ariana S Huffmyer, E5 RoL Team
08/15/2022

NOTE: If you want to run this script, do not Knit, instead “run all
chunks below” and look at console output.

IDEA: Host phenotype = calcification (host fitness); could also do
symbiont fitness (perhaps max photosynthesis AM)?

Next script: Models of 1) best predictors of host fitness
(calcification); 2) best predictors of symbiont fitness (Am first, then
could look into Ik or E) or holobiont fitness (P:R ratios); 3) Lag times
between symbiont and host responses

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
if (!require("Hmisc")) install.packages("Hmisc")

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
library(Hmisc)
```

# Load dataframe

Load in master dataframe generated from 1_assemble_data.Rmd.

``` r
master<-read.csv("Output/master_timeseries.csv")
```

# Run models

## All species

Model the influence of physiological metrics on calcification.

``` r
df<-master%>%
  select(Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Am, AQY, Rd, calc.umol.mgAFDW.hr, cells.mgAFDW, prot_mg.mgafdw, cre.umol.mgafdw, Ratio_AFDW.mg.cm2, Total_Chl_cell, Total_Chl)

model1<-lm(calc.umol.mgAFDW.hr~ ., data=df)

summary(model1)
```

    ## 
    ## Call:
    ## lm(formula = calc.umol.mgAFDW.hr ~ ., data = df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36664 -0.07920 -0.02816  0.03310  1.05322 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -7.680e-02  4.093e-02  -1.876  0.06142 .  
    ## Host_AFDW.mg.cm2  -1.163e-02  6.334e-03  -1.836  0.06709 .  
    ## Sym_AFDW.mg.cm2   -4.857e-03  2.930e-03  -1.658  0.09827 .  
    ## Am                 6.754e-02  2.597e-02   2.601  0.00968 ** 
    ## AQY                2.413e+00  4.171e+00   0.578  0.56328    
    ## Rd                -1.752e-01  6.365e-02  -2.752  0.00621 ** 
    ## cells.mgAFDW      -5.882e-08  5.774e-08  -1.019  0.30900    
    ## prot_mg.mgafdw     4.651e-01  7.835e-02   5.935 6.76e-09 ***
    ## cre.umol.mgafdw   -2.094e-01  1.054e-01  -1.986  0.04778 *  
    ## Ratio_AFDW.mg.cm2  2.268e-01  1.026e-01   2.212  0.02760 *  
    ## Total_Chl_cell    -1.787e+02  1.397e+03  -0.128  0.89824    
    ## Total_Chl         -1.091e-02  6.049e-03  -1.804  0.07201 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1523 on 369 degrees of freedom
    ##   (67 observations deleted due to missingness)
    ## Multiple R-squared:  0.2042, Adjusted R-squared:  0.1805 
    ## F-statistic:  8.61 on 11 and 369 DF,  p-value: 1.375e-13

From this first look, Am, Rd, and protein.AFDW are most significant
(\<0.01) and antioxidant (cre.AFDW) and ratio host:symbiont AFDW are
also significant (\<0.05) as predictors of calcification.

Show effect plot

``` r
plot(allEffects(model1)) 
```

![](5_modeling_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Next, look for covariates and simplify the model by removing highly
correlated variables.

``` r
cor_all<-rcorr(as.matrix(df))
#combine all values 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

cor_all_df<-flattenCorrMatrix(cor_all$r, cor_all$P)
```

List correlations with absolute value r>0.9.

``` r
cor_all_df%>%
  filter(abs(cor)>0.85)
```

    ##   row column       cor p
    ## 1  Am     Rd 0.8657313 0

Respiration and Am are \> 0.85.

List correlations with high significance p\<0.01.

``` r
cor_all_df%>%
  filter(p<0.01)
```

    ##                    row            column        cor            p
    ## 1     Host_AFDW.mg.cm2   Sym_AFDW.mg.cm2  0.5268488 0.000000e+00
    ## 2     Host_AFDW.mg.cm2                Am  0.7715069 0.000000e+00
    ## 3      Sym_AFDW.mg.cm2                Am  0.4425720 0.000000e+00
    ## 4     Host_AFDW.mg.cm2               AQY  0.7000738 0.000000e+00
    ## 5      Sym_AFDW.mg.cm2               AQY  0.3856795 0.000000e+00
    ## 6                   Am               AQY  0.7220939 0.000000e+00
    ## 7     Host_AFDW.mg.cm2                Rd  0.7143785 0.000000e+00
    ## 8      Sym_AFDW.mg.cm2                Rd  0.4343864 0.000000e+00
    ## 9                   Am                Rd  0.8657313 0.000000e+00
    ## 10                 AQY                Rd  0.7339445 0.000000e+00
    ## 11    Host_AFDW.mg.cm2      cells.mgAFDW -0.1480500 2.009304e-03
    ## 12 calc.umol.mgAFDW.hr      cells.mgAFDW  0.1708633 4.094371e-04
    ## 13    Host_AFDW.mg.cm2    prot_mg.mgafdw  0.1368544 4.376829e-03
    ## 14     Sym_AFDW.mg.cm2    prot_mg.mgafdw  0.2345334 1.042934e-06
    ## 15                  Am    prot_mg.mgafdw  0.2089195 1.225315e-05
    ## 16                 AQY    prot_mg.mgafdw  0.1931586 5.425264e-05
    ## 17                  Rd    prot_mg.mgafdw  0.3015053 1.652924e-10
    ## 18 calc.umol.mgAFDW.hr    prot_mg.mgafdw  0.3413909 5.253575e-13
    ## 19        cells.mgAFDW    prot_mg.mgafdw  0.4763724 0.000000e+00
    ## 20    Host_AFDW.mg.cm2   cre.umol.mgafdw  0.6010372 0.000000e+00
    ## 21     Sym_AFDW.mg.cm2   cre.umol.mgafdw  0.4190958 0.000000e+00
    ## 22                  Am   cre.umol.mgafdw  0.5638934 0.000000e+00
    ## 23                 AQY   cre.umol.mgafdw  0.5042028 0.000000e+00
    ## 24                  Rd   cre.umol.mgafdw  0.5190912 0.000000e+00
    ## 25      prot_mg.mgafdw   cre.umol.mgafdw  0.5353145 0.000000e+00
    ## 26     Sym_AFDW.mg.cm2 Ratio_AFDW.mg.cm2  0.3972511 0.000000e+00
    ## 27                  Rd Ratio_AFDW.mg.cm2  0.1974731 3.888818e-05
    ## 28 calc.umol.mgAFDW.hr Ratio_AFDW.mg.cm2  0.2381165 7.713314e-07
    ## 29        cells.mgAFDW Ratio_AFDW.mg.cm2  0.4576379 0.000000e+00
    ## 30      prot_mg.mgafdw Ratio_AFDW.mg.cm2  0.5983127 0.000000e+00
    ## 31     cre.umol.mgafdw Ratio_AFDW.mg.cm2  0.3760915 1.021405e-14
    ## 32    Host_AFDW.mg.cm2    Total_Chl_cell -0.1965170 4.073329e-05
    ## 33     Sym_AFDW.mg.cm2    Total_Chl_cell -0.1513614 1.604603e-03
    ## 34                  Am    Total_Chl_cell -0.2558461 5.448323e-08
    ## 35                 AQY    Total_Chl_cell -0.1329581 5.266937e-03
    ## 36                  Rd    Total_Chl_cell -0.1671075 4.379996e-04
    ## 37        cells.mgAFDW    Total_Chl_cell -0.4282379 0.000000e+00
    ## 38     cre.umol.mgafdw    Total_Chl_cell -0.2188503 1.081184e-05
    ## 39    Host_AFDW.mg.cm2         Total_Chl -0.1645926 5.763184e-04
    ## 40 calc.umol.mgAFDW.hr         Total_Chl  0.1907257 7.598882e-05
    ## 41        cells.mgAFDW         Total_Chl  0.4944664 0.000000e+00
    ## 42      prot_mg.mgafdw         Total_Chl  0.6206999 0.000000e+00
    ## 43   Ratio_AFDW.mg.cm2         Total_Chl  0.4098812 0.000000e+00
    ## 44      Total_Chl_cell         Total_Chl  0.2525815 1.104011e-07

Many correlations are highly significant.

View the correlations between calcification and other physiological
parameters.

``` r
cor_all_df%>%
  filter(p<0.01)%>%
  filter(row=="calc.umol.mgAFDW.hr")
```

    ##                   row            column       cor            p
    ## 1 calc.umol.mgAFDW.hr      cells.mgAFDW 0.1708633 4.094371e-04
    ## 2 calc.umol.mgAFDW.hr    prot_mg.mgafdw 0.3413909 5.253575e-13
    ## 3 calc.umol.mgAFDW.hr Ratio_AFDW.mg.cm2 0.2381165 7.713314e-07
    ## 4 calc.umol.mgAFDW.hr         Total_Chl 0.1907257 7.598882e-05

There are significant, but weak positive correlations between
calcification & symbiont cell density, protein, symbiont:host biomass
ratio, and total chlorophyll.

## Acropora

Model the influence of physiological metrics on calcification.

``` r
df_acr<-master%>%
  filter(species=="Acropora")%>%
  select(Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Am, AQY, Rd, calc.umol.mgAFDW.hr, cells.mgAFDW, prot_mg.mgafdw, cre.umol.mgafdw, Ratio_AFDW.mg.cm2, Total_Chl_cell, Total_Chl)

model1_acr<-lm(calc.umol.mgAFDW.hr~ ., data=df_acr)

summary(model1_acr)
```

    ## 
    ## Call:
    ## lm(formula = calc.umol.mgAFDW.hr ~ ., data = df_acr)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40893 -0.11812 -0.03252  0.07664  0.72753 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -8.126e-01  3.275e-01  -2.482 0.014892 *  
    ## Host_AFDW.mg.cm2   2.490e-01  1.529e-01   1.628 0.106909    
    ## Sym_AFDW.mg.cm2   -3.360e-01  2.079e-01  -1.616 0.109542    
    ## Am                -3.902e-02  1.043e-01  -0.374 0.709223    
    ## AQY                1.275e+01  2.274e+01   0.561 0.576491    
    ## Rd                -3.237e-01  2.055e-01  -1.575 0.118652    
    ## cells.mgAFDW       2.156e-07  1.517e-07   1.421 0.158713    
    ## prot_mg.mgafdw     9.342e-01  2.000e-01   4.672 1.02e-05 ***
    ## cre.umol.mgafdw   -8.876e-01  5.230e-01  -1.697 0.093060 .  
    ## Ratio_AFDW.mg.cm2  1.504e+00  7.805e-01   1.927 0.057078 .  
    ## Total_Chl_cell     3.082e+04  8.675e+03   3.553 0.000604 ***
    ## Total_Chl         -7.268e-02  2.007e-02  -3.622 0.000479 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1979 on 92 degrees of freedom
    ##   (12 observations deleted due to missingness)
    ## Multiple R-squared:  0.5011, Adjusted R-squared:  0.4414 
    ## F-statistic: 8.399 on 11 and 92 DF,  p-value: 4.768e-10

From this first look, protein.AFDW, total chlorophyll.AFDW and total
chlorophyll per cell are significant (\<0.01).

Show effect plot

``` r
plot(allEffects(model1_acr)) 
```

![](5_modeling_analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Next, look for covariates and simplify the model by removing highly
correlated variables.

``` r
cor_acr<-rcorr(as.matrix(df_acr))
#combine all values 

cor_acr_df<-flattenCorrMatrix(cor_acr$r, cor_acr$P)
```

List correlations with absolute value r>0.9.

``` r
cor_acr_df%>%
  filter(abs(cor)>0.85)
```

    ## [1] row    column cor    p     
    ## <0 rows> (or 0-length row.names)

No correlations \> 0.85.

List correlations with high significance p\<0.01.

``` r
cor_acr_df%>%
  filter(p<0.01)
```

    ##                    row              column        cor            p
    ## 1     Host_AFDW.mg.cm2     Sym_AFDW.mg.cm2  0.3796031 3.108031e-05
    ## 2      Sym_AFDW.mg.cm2                  Am  0.3030589 9.939881e-04
    ## 3                   Am                 AQY  0.5874937 4.138245e-12
    ## 4      Sym_AFDW.mg.cm2                  Rd  0.3888123 1.753111e-05
    ## 5                   Am                  Rd  0.7689683 0.000000e+00
    ## 6                  AQY                  Rd  0.5679211 2.959921e-11
    ## 7     Host_AFDW.mg.cm2 calc.umol.mgAFDW.hr -0.3641716 7.906243e-05
    ## 8     Host_AFDW.mg.cm2        cells.mgAFDW -0.2841272 2.186799e-03
    ## 9                   Am        cells.mgAFDW  0.3475344 1.515258e-04
    ## 10                  Rd        cells.mgAFDW  0.2653496 4.327177e-03
    ## 11    Host_AFDW.mg.cm2      prot_mg.mgafdw -0.4998883 1.730008e-08
    ## 12 calc.umol.mgAFDW.hr      prot_mg.mgafdw  0.5499740 4.049636e-10
    ## 13        cells.mgAFDW      prot_mg.mgafdw  0.3557194 1.103793e-04
    ## 14                  Am     cre.umol.mgafdw -0.3131977 1.020758e-03
    ## 15 calc.umol.mgAFDW.hr     cre.umol.mgafdw  0.4333572 3.880995e-06
    ## 16      prot_mg.mgafdw     cre.umol.mgafdw  0.6715549 2.442491e-15
    ## 17    Host_AFDW.mg.cm2   Ratio_AFDW.mg.cm2 -0.4923268 2.642463e-08
    ## 18     Sym_AFDW.mg.cm2   Ratio_AFDW.mg.cm2  0.5642786 6.227996e-11
    ## 19                  Rd   Ratio_AFDW.mg.cm2  0.2547079 6.241356e-03
    ## 20 calc.umol.mgAFDW.hr   Ratio_AFDW.mg.cm2  0.3938683 1.734755e-05
    ## 21        cells.mgAFDW   Ratio_AFDW.mg.cm2  0.2774822 2.799149e-03
    ## 22      prot_mg.mgafdw   Ratio_AFDW.mg.cm2  0.5821236 1.355427e-11
    ## 23     cre.umol.mgafdw   Ratio_AFDW.mg.cm2  0.3816093 5.000188e-05
    ## 24 calc.umol.mgAFDW.hr      Total_Chl_cell  0.3573722 1.178826e-04
    ## 25        cells.mgAFDW      Total_Chl_cell -0.3737654 4.546085e-05
    ## 26      prot_mg.mgafdw      Total_Chl_cell  0.3629012 8.409070e-05
    ## 27     cre.umol.mgafdw      Total_Chl_cell  0.4769909 2.349843e-07
    ## 28   Ratio_AFDW.mg.cm2      Total_Chl_cell  0.3529941 1.256317e-04
    ## 29    Host_AFDW.mg.cm2           Total_Chl -0.3164609 6.382305e-04
    ## 30 calc.umol.mgAFDW.hr           Total_Chl  0.2994134 1.411070e-03
    ## 31        cells.mgAFDW           Total_Chl  0.4944600 2.602027e-08
    ## 32      prot_mg.mgafdw           Total_Chl  0.7418580 0.000000e+00
    ## 33     cre.umol.mgafdw           Total_Chl  0.4567787 8.581520e-07
    ## 34   Ratio_AFDW.mg.cm2           Total_Chl  0.5206205 3.402928e-09
    ## 35      Total_Chl_cell           Total_Chl  0.4942467 2.643705e-08

Many correlations are highly significant.

View the correlations between calcification and other physiological
parameters.

``` r
cor_acr_df%>%
  filter(p<0.01)%>%
  filter(row=="calc.umol.mgAFDW.hr")
```

    ##                   row            column       cor            p
    ## 1 calc.umol.mgAFDW.hr    prot_mg.mgafdw 0.5499740 4.049636e-10
    ## 2 calc.umol.mgAFDW.hr   cre.umol.mgafdw 0.4333572 3.880995e-06
    ## 3 calc.umol.mgAFDW.hr Ratio_AFDW.mg.cm2 0.3938683 1.734755e-05
    ## 4 calc.umol.mgAFDW.hr    Total_Chl_cell 0.3573722 1.178826e-04
    ## 5 calc.umol.mgAFDW.hr         Total_Chl 0.2994134 1.411070e-03

There are significant positive correlations between calcification &
protein, antioxidant capacity, host:symbiont biomass ratio, total
chlorophyll, and total chlorophyll per symbiont cell.

## Pocillopora

Model the influence of physiological metrics on calcification.

``` r
df_poc<-master%>%
  filter(species=="Pocillopora")%>%
  select(Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Am, AQY, Rd, calc.umol.mgAFDW.hr, cells.mgAFDW, prot_mg.mgafdw, cre.umol.mgafdw, Ratio_AFDW.mg.cm2, Total_Chl_cell, Total_Chl)

model1_poc<-lm(calc.umol.mgAFDW.hr~ ., data=df_poc)

summary(model1_poc)
```

    ## 
    ## Call:
    ## lm(formula = calc.umol.mgAFDW.hr ~ ., data = df_poc)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11233 -0.05256 -0.02525  0.02278  0.97889 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        1.480e-02  7.038e-02   0.210  0.83374   
    ## Host_AFDW.mg.cm2  -1.347e-03  1.584e-02  -0.085  0.93239   
    ## Sym_AFDW.mg.cm2   -1.964e-03  5.294e-02  -0.037  0.97046   
    ## Am                -2.599e-02  6.005e-02  -0.433  0.66591   
    ## AQY                8.819e+00  1.049e+01   0.841  0.40202   
    ## Rd                -6.179e-03  9.936e-02  -0.062  0.95051   
    ## cells.mgAFDW       2.531e-07  3.592e-07   0.704  0.48243   
    ## prot_mg.mgafdw     1.792e-02  1.966e-01   0.091  0.92753   
    ## cre.umol.mgafdw   -9.677e-01  3.175e-01  -3.048  0.00281 **
    ## Ratio_AFDW.mg.cm2  1.707e-01  2.136e-01   0.799  0.42567   
    ## Total_Chl_cell     3.861e+02  2.368e+03   0.163  0.87074   
    ## Total_Chl         -1.222e-02  2.076e-02  -0.588  0.55727   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1125 on 125 degrees of freedom
    ##   (28 observations deleted due to missingness)
    ## Multiple R-squared:  0.1103, Adjusted R-squared:  0.03196 
    ## F-statistic: 1.408 on 11 and 125 DF,  p-value: 0.1771

From this first look, antioxidant.AFDW is the only effect that is
significant (\<0.01).

Show effect plot

``` r
plot(allEffects(model1_poc)) 
```

![](5_modeling_analysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Next, look for covariates and simplify the model by removing highly
correlated variables.

``` r
cor_poc<-rcorr(as.matrix(df_poc))
#combine all values 

cor_poc_df<-flattenCorrMatrix(cor_poc$r, cor_poc$P)
```

List correlations with absolute value r>0.9.

``` r
cor_poc_df%>%
  filter(abs(cor)>0.85)
```

    ## [1] row    column cor    p     
    ## <0 rows> (or 0-length row.names)

No correlations \> 0.85.

List correlations with high significance p\<0.01.

``` r
cor_poc_df%>%
  filter(p<0.01)
```

    ##                    row            column        cor            p
    ## 1     Host_AFDW.mg.cm2   Sym_AFDW.mg.cm2  0.2580371 1.103197e-03
    ## 2      Sym_AFDW.mg.cm2               AQY  0.2271514 4.223731e-03
    ## 3                   Am               AQY  0.3257801 2.073490e-05
    ## 4      Sym_AFDW.mg.cm2                Rd  0.2408498 2.376938e-03
    ## 5                   Am                Rd  0.4251360 1.392652e-08
    ## 6                  AQY                Rd  0.3895323 2.530929e-07
    ## 7                   Am      cells.mgAFDW  0.2421793 2.100882e-03
    ## 8         cells.mgAFDW    prot_mg.mgafdw  0.5134015 5.279555e-12
    ## 9      Sym_AFDW.mg.cm2   cre.umol.mgafdw  0.3577725 9.272655e-06
    ## 10                  Am   cre.umol.mgafdw -0.2117951 8.808071e-03
    ## 11 calc.umol.mgAFDW.hr   cre.umol.mgafdw -0.2383081 3.537521e-03
    ## 12     Sym_AFDW.mg.cm2 Ratio_AFDW.mg.cm2  0.7324980 0.000000e+00
    ## 13        cells.mgAFDW Ratio_AFDW.mg.cm2  0.2842650 3.530600e-04
    ## 14      prot_mg.mgafdw Ratio_AFDW.mg.cm2  0.5869928 8.881784e-16
    ## 15     cre.umol.mgafdw Ratio_AFDW.mg.cm2  0.3225727 7.162062e-05
    ## 16        cells.mgAFDW    Total_Chl_cell -0.4318897 1.461302e-08
    ## 17                  Am         Total_Chl  0.2357263 2.610206e-03
    ## 18                 AQY         Total_Chl  0.2709837 5.072614e-04
    ## 19                  Rd         Total_Chl  0.2191367 5.223835e-03
    ## 20        cells.mgAFDW         Total_Chl  0.5785560 1.776357e-15
    ## 21      prot_mg.mgafdw         Total_Chl  0.4557949 1.391595e-09
    ## 22   Ratio_AFDW.mg.cm2         Total_Chl  0.4214824 4.719727e-08
    ## 23      Total_Chl_cell         Total_Chl  0.3710492 1.591107e-06

Many correlations are highly significant.

View the correlations between calcification and other physiological
parameters.

``` r
cor_poc_df%>%
  filter(p<0.01)%>%
  filter(row=="calc.umol.mgAFDW.hr")
```

    ##                   row          column        cor           p
    ## 1 calc.umol.mgAFDW.hr cre.umol.mgafdw -0.2383081 0.003537521

There is one significant negative correlation between calcification and
antioxidant capacity.

## Porites

Model the influence of physiological metrics on calcification.

``` r
df_por<-master%>%
  filter(species=="Porites")%>%
  select(Host_AFDW.mg.cm2, Sym_AFDW.mg.cm2, Am, AQY, Rd, calc.umol.mgAFDW.hr, cells.mgAFDW, prot_mg.mgafdw, cre.umol.mgafdw, Ratio_AFDW.mg.cm2, Total_Chl_cell, Total_Chl)

model1_por<-lm(calc.umol.mgAFDW.hr~ ., data=df_por)

summary(model1_por)
```

    ## 
    ## Call:
    ## lm(formula = calc.umol.mgAFDW.hr ~ ., data = df_por)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.19032 -0.05706 -0.01740  0.04014  0.26332 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        9.550e-02  5.823e-02   1.640   0.1034    
    ## Host_AFDW.mg.cm2  -1.917e-02  4.755e-03  -4.031 9.47e-05 ***
    ## Sym_AFDW.mg.cm2    1.154e-03  2.246e-03   0.514   0.6082    
    ## Am                 4.632e-02  1.847e-02   2.507   0.0134 *  
    ## AQY                2.287e+00  2.684e+00   0.852   0.3957    
    ## Rd                -3.380e-03  5.068e-02  -0.067   0.9469    
    ## cells.mgAFDW      -7.379e-09  1.156e-07  -0.064   0.9492    
    ## prot_mg.mgafdw    -7.288e-02  7.607e-02  -0.958   0.3399    
    ## cre.umol.mgafdw    1.580e-01  9.268e-02   1.705   0.0906 .  
    ## Ratio_AFDW.mg.cm2 -7.176e-02  1.391e-01  -0.516   0.6067    
    ## Total_Chl_cell    -1.586e+03  2.382e+03  -0.666   0.5065    
    ## Total_Chl          2.508e-03  6.410e-03   0.391   0.6963    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08591 on 128 degrees of freedom
    ##   (27 observations deleted due to missingness)
    ## Multiple R-squared:  0.1949, Adjusted R-squared:  0.1257 
    ## F-statistic: 2.817 on 11 and 128 DF,  p-value: 0.002525

From this first look, Host AFDW is the only effect that is higly
significant (\<0.01) and Am is also significant (\<0.05).

Show effect plot

``` r
plot(allEffects(model1_poc)) 
```

![](5_modeling_analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Next, look for covariates and simplify the model by removing highly
correlated variables.

``` r
cor_por<-rcorr(as.matrix(df_por))
#combine all values 

cor_por_df<-flattenCorrMatrix(cor_por$r, cor_por$P)
```

List correlations with absolute value r>0.9.

``` r
cor_por_df%>%
  filter(abs(cor)>0.85)
```

    ## [1] row    column cor    p     
    ## <0 rows> (or 0-length row.names)

No correlations \> 0.85.

List correlations with high significance p\<0.01.

``` r
cor_por_df%>%
  filter(p<0.01)
```

    ##                  row              column        cor            p
    ## 1   Host_AFDW.mg.cm2                  Am  0.3733510 1.253007e-06
    ## 2   Host_AFDW.mg.cm2                 AQY  0.2682030 6.306686e-04
    ## 3                 Am                 AQY  0.3775860 4.896892e-07
    ## 4   Host_AFDW.mg.cm2                  Rd  0.3567959 3.904731e-06
    ## 5                 Am                  Rd  0.7550583 0.000000e+00
    ## 6                AQY                  Rd  0.4785109 6.131184e-11
    ## 7   Host_AFDW.mg.cm2 calc.umol.mgAFDW.hr -0.2246752 4.670295e-03
    ## 8                 Am calc.umol.mgAFDW.hr  0.2480916 1.731083e-03
    ## 9                 Rd        cells.mgAFDW  0.2512130 1.402861e-03
    ## 10                Rd      prot_mg.mgafdw  0.2363604 2.879972e-03
    ## 11      cells.mgAFDW      prot_mg.mgafdw  0.7016420 0.000000e+00
    ## 12      cells.mgAFDW     cre.umol.mgafdw  0.4565224 1.000925e-08
    ## 13    prot_mg.mgafdw     cre.umol.mgafdw  0.5424470 2.601031e-12
    ## 14   Sym_AFDW.mg.cm2   Ratio_AFDW.mg.cm2  0.6400073 0.000000e+00
    ## 15      cells.mgAFDW   Ratio_AFDW.mg.cm2  0.5507244 6.528111e-14
    ## 16    prot_mg.mgafdw   Ratio_AFDW.mg.cm2  0.6068954 0.000000e+00
    ## 17   cre.umol.mgafdw   Ratio_AFDW.mg.cm2  0.5140594 6.022471e-11
    ## 18      cells.mgAFDW      Total_Chl_cell -0.2554101 1.157089e-03
    ## 19   cre.umol.mgafdw      Total_Chl_cell -0.2369545 4.380302e-03
    ## 20      cells.mgAFDW           Total_Chl  0.4598976 1.070705e-09
    ## 21    prot_mg.mgafdw           Total_Chl  0.6160563 0.000000e+00
    ## 22 Ratio_AFDW.mg.cm2           Total_Chl  0.2979827 1.431866e-04
    ## 23    Total_Chl_cell           Total_Chl  0.5830184 8.881784e-16

Many correlations are highly significant.

Based on correlation analysis there are no strongly correlated variables
and we can leave all responses in the model.

View the correlations between calcification and other physiological
parameters.

``` r
cor_por_df%>%
  filter(p<0.01)%>%
  filter(row=="calc.umol.mgAFDW.hr")
```

    ## [1] row    column cor    p     
    ## <0 rows> (or 0-length row.names)

There is nothing that correlates with calcification in Porites.
