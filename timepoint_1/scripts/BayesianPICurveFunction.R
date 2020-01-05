### Creating a function to run PI curves by individual in a bayesian framework

### Created by Nyssa Silbiger and Hollie Putnam
## Created on 10/26/19

###################################################################
rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("phytotools" %in% rownames(installed.packages()) == 'FALSE') install.packages('phytotools') 
if ("brms" %in% rownames(installed.packages()) == 'FALSE') install.packages('brms') 
if ("broom" %in% rownames(installed.packages()) == 'FALSE') install_github('broom') 
if ("tidybayes" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidybayes') 
if ("purrr" %in% rownames(installed.packages()) == 'FALSE') install.packages('purrr') 
if ("modelr" %in% rownames(installed.packages()) == 'FALSE') install.packages('modelr') 
if ("bayesplot" %in% rownames(installed.packages()) == 'FALSE') install.packages('bayesplot') 
if ("cowplot" %in% rownames(installed.packages()) == 'FALSE') install.packages('cowplot') 

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library(devtools)
library("tidyverse")
library("segmented")
library("plotrix")
library("gridExtra")
library("lubridate")
library("chron")
library('phytotools')
library(brms)
library(broom)
library(tidybayes)
library(purrr)
library(modelr)
library(bayesplot)
library(cowplot)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

################ Read in the Data ##################

Data<-read.csv(file = 'RAnalysis/Data/All_PI_Curve_rates.csv')
# add a column for site
Data<-separate(data = Data, col = Fragment.ID, into = "Site", sep =  "-", remove = TRUE)
# function to plot PI curves using Baysian analysis
PICurve_ind<-function(Data = Data, IndividualID, PlotDiagnostics = TRUE, PlotResults = TRUE, n_cores = 4,
                      PAR_col_name = "Light_Value", rate_col_name = "micromol.cm2.h", ind_col_name = "Sample.ID"){
  
  ## Data = Dataframe
  ## IndividualID = ID to sample from
  ## PlotDiagnostics = do you want to plot the diagnostic plotss?
  ## Plot Results = Do you want to plot the individual results?
  ## n_cores = how many cores do you want to use to run the bayesian analysis?
  ## PAR_col_name = what is the name of the column with light in quotes
  ## rate_col_name = what is the name of the column with rates
  ## ind_col_name = what is the name of the column with the individual IDs?
  
Individual.selected<-Data # subset the data by individual ID
# rename the Light and rate columns from the dataframe
colnames(Data)[colnames(Data)==PAR_col_name]<-"Light_Level"
colnames(Data)[colnames(Data)==rate_col_name]<-"micromol.cm2.h"
colnames(Data)[colnames(Data)==ind_col_name]<-"Sample.ID"
# pull out the names of the individual
spec.name<-unique(Individual.selected$Sample.ID)
y<-Individual.selected$micromol.cm2.h

options("scipen"=100,digits=12) # stan doesnt like scientific notation. This fixes that

#fit a model using a Bayesian model of a non-rectangular hyperbola (Marshall & Biscoe, 1980)

# this is for posterior predictove checks
stanvars <- stanvar(scode = "vector[N] y_new; real Ic;
  vector[N] nlp_theta = X_theta * b_theta;
  vector[N] nlp_AQY = X_AQY * b_AQY;
  vector[N] nlp_Am = X_Am * b_Am;
  vector[N] nlp_Rd = X_Rd * b_Rd;
  Ic = nlp_Rd[1]/nlp_AQY[1];

    for (n in 1:N) {  
        y_new[n] = normal_rng((1 / (2 * nlp_theta[n])) * (nlp_AQY[n] * C_1[n] + nlp_Am[n] - sqrt((nlp_AQY[n] * C_1[n] + nlp_Am[n]) ^ 2 - 4 * nlp_AQY[n] * nlp_theta[n] * nlp_Am[n] * C_1[n])) - nlp_Rd[n], sigma);};",
        block = "genquant")

# model fit
fit1<-brm(
  bf(micromol.cm2.h ~ (1/(2*theta))*(AQY*Light_Value+Am-sqrt((AQY*Light_Value+Am)^2-4*AQY*theta*Am*Light_Value))-Rd,
     AQY~1, Am~1, theta~1, Rd~1,
          nl = TRUE),
          data = Individual.selected,
          family = gaussian(),
  prior = c(
    prior(normal(0,10), nlpar = "AQY", lb=0),  # set the priors
    prior(normal(0,10), nlpar = "Am", lb=0),
    prior(normal(0, 1), nlpar = "theta"),
    prior(normal(0,1), nlpar = "Rd")
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 20), # force stan to take smaller steps to reduce divergent errors
cores = n_cores, chains = 4, seed = 126, iter = 3000, warmup = 2000, silent = TRUE, save_model = TRUE, stanvars = stanvars) 
  
# plot diagnostics if it is true
if (PlotDiagnostics==TRUE){
  ## assess the fits
  posterior <-  posterior_samples(fit1)

    # convergence diagnostics -------------------------------------------------
  # plot the traceplots
  pdf(paste0('RAnalysis/Output/BayesPICurves/traceplots_',spec.name,'.pdf'))
  plot(fit1, ask = FALSE, newpage = FALSE) # all of them
  dev.off()
  
  ## posterior predictive checks
  y_rep <- as.matrix(fit1, pars = "y_new")
  dim(y_rep)
  check1<-ppc_dens_overlay(y, y_rep[1:200, ]) # comparing density of y with densities of y over 200 posterior draws.
  check2<-ppc_stat(y = y, yrep = y_rep, stat = "mean") # compare estimates of summary statistics
  check3<-ppc_scatter_avg(y = y, yrep = y_rep) # observed vs predicted with 1:1 line
 # put them in one plot  
  ppcheckfig<-plot_grid(check1, check2, check3, labels = 'AUTO', nrow = 1 )
  # add title
  title <- ggdraw() + draw_label(spec.name, fontface='italic')
  ppcheckf<-plot_grid(title, ppcheckfig,nrow = 2 , rel_heights=c(0.1, 1))
  
  ggsave(filename = paste0('RAnalysis/Output/BayesPICurves/ppchecks_',spec.name,'.pdf'),plot = ppcheckf, width = 9, height = 6 )
}

if (PlotResults == TRUE){
  
pt1<-marginal_effects(fit1)

#population level 
pt2<-ggplot(as.data.frame(pt1$Light_Value))+ # pull pur the fitted data
  geom_line(aes(Light_Value, estimate__))+ # 
  geom_ribbon(aes(Light_Value, ymin = lower__, ymax = upper__), alpha=0.3)+
  geom_point(data = Individual.selected, aes(x = Light_Value, y = micromol.cm2.h)) +# add the raw data
  theme_bw()+
  xlab(expression(paste('PAR (', mu, "mol photons m"^-2, 'hr'^-1,")")))+
  ylab(expression(paste('Photosynthetic rate (', mu, "mol cm"^-2, 's'^-1,")")))+
  theme(text = element_text(size=18), title = element_text(face="italic"))

ggsave(filename = paste0('RAnalysis/Output/BayesPICurves/PICurve_',spec.name,'.pdf'),plot = pt2, width = 6, height = 6 )

}

# return the parameters
params<-fit1 %>%
  gather_draws(b_AQY_Intercept, b_Am_Intercept, b_theta_Intercept, b_Rd_Intercept, sigma, Ic) %>%
  median_qi()

return(params)

# print when each species is done
print(paste(spec.name, 'is done'))
rm(fit1) # clear the fit

} #end function

p<-as.data.frame(x = matrix(NA,nrow = 1, ncol = 7))
colnames(p)<-c('.variable','.value','.lower', '.upper', '.width', '.point','.interval')

#Same functions, but skips over species with an error
PICurve_ind_noerror <- possibly(PICurve_ind, otherwise = p)

#Run Bayesian PI for all Individuals, export plots, and save the 95% CI of each parameter
Param.output<-Data %>%
  group_by(Sample.ID, Species, Site) %>%
  do(PICurve_ind_noerror(Data = ., IndividualID = unique(.$Sample.ID)))

# add a column for plotting names
prettynames<-data.frame(.variable = unique(Param.output$.variable))
prettynames$varnames<-c("Am","AQY","Rd","Theta", "Sigma", "Ic")
Param.output<-left_join(Param.output,prettynames)# make the names easier for plotting
write.csv(file  = 'RAnalysis/Output/BayesPICurves/parameters.csv', x = Param.output)


Param.output <- read.csv(file='RAnalysis/Output/BayesPICurves/parameters.csv')
Param.output <- subset(Param.output, varnames!="Theta" & varnames!="Sigma")
Param.output$group <- paste0() 

ymin <- c(0,0,0)
ymax <- c(4,0.45,2)

## Make a plot of the means
Param.output%>%
  group_by(Species, Site, varnames)%>%
  summarise(mean.value = mean(.value), se = std.error(.value)) %>%
  ggplot(aes(x = Site, y = mean.value, group = Species, color = Species))+
  geom_point(size = 3)+
  geom_errorbar(aes(x = Site, ymin = mean.value-se, ymax = mean.value+se), width = 0.5)+
  facet_wrap(~varnames*Species, scales = "free_y", ncol = 4) +
  coord_cartesian(ylim = c(ymin, ymax)) 
  #coord_cartesian(ylim = c(0, 0.45)) +
 # coord_cartesian(ylim = c(0, 2.5))
  



