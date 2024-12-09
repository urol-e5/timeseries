![GitHub commit activity](https://img.shields.io/github/commit-activity/m/urol-e5/timeseries) ![GitHub last commit](https://img.shields.io/github/last-commit/urol-e5/timeseries) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/urol-e5/timeseries) ![GitHub issues](https://img.shields.io/github/issues/urol-e5/timeseries) ![GitHub closed issues](https://img.shields.io/github/issues-closed/urol-e5/timeseries) ![GitHub contributors](https://img.shields.io/github/contributors/urol-e5/timeseries)

# E5 URoL physiology timeseries project 

### Abstract  

**Holobiont composition and seasonal plasticity in physiology are key to understanding reef building coral performance** 

The symbiotic relationship between corals and Symbiodiniaceae generates holobiont sensitivity to fluctuations in temperature, light, and nutrients, which change seasonally in tropical oceans. However, prior work examining coral stress response has occurred primarily without consideration of the current physiological state due to seasonal timing. Here, we tested the effect of seasonal variation on three dominant and ecologically important genera, *Acropora, Pocillopora*, and *Porites* (n=15 tagged colonies per site) in Mo’orea French Polynesia across three lagoon sites on the North shore in January, March, September, and December of 2020. We quantified physiological variables for the coral and Symbiodiniaceae, as well as molecular identification of the Symbiodiniaceae community and host species. Physiology of the three genera differed significantly. Variance partitioning analyses identified seasonal timepoint as the dominant explanatory factor for both host and symbiont responses, with a small influence of site. Acropora colonies were composed of *A. pulchra*, but we found cryptic haplotypes in *Pocillopora* (*P. meandrina, P. tuahiniensis*) and *Porites* (*P. evermanni, P. lobata/lutea*). Although physiology did not vary strongly between *Pocillopora* haplotypes, Porites haplotypes exhibited strong differences in physiological responses and were associated with distinct Symbiodiniaceae communities. Overall, across functional groups of dominant coral genera, seasonality is the strongest driver of physiological state and particular attention should be paid to cryptic species variation. These data provide an essential picture of the need for considering physiological state and its seasonality, particularly in light of the plethora of climate change related stress test assays taking place throughout the year. 

### Experimental Design 

![Sites](https://github.com/urol-e5/timeseries/blob/master/metadata/E5_Sites.png?raw=true)

## Navigating this repository 

This repository is organized such that each sampling time point (1 = January 2020, 2 = March 2020, 3 = September 2020, 4 = November 2020) has a dedicated folder with all metadata, data files, and scripts necessary for calculating physiological responses. All data output in these individual time points is then read in, assembled, and analyzed in the `time_series_analysis` directory. 

### Contents 

- [`metadata`](https://github.com/urol-e5/timeseries/tree/master/metadata)
	- This folder contains metadata files that were generated during the study. If you are looking for final and corrected colony ID's or metadata, DO NOT use these files. See the "Looking for colony metadata" section below. 
	
- [`timepoint_1`](https://github.com/urol-e5/timeseries/tree/master/timepoint_1)
	- This folder contains `scripts`, `data`, and `output` for physiological responses collected in January 2020. The scripts in this folder read in raw data files and conduct normalization and calculations for each response. Calculated responses are output in the `output` folder and analyzed as described in the `time_series_analysis` folder. Open this project with the `timepoint_1.Rproj` file. 
	
- [`timepoint_2`](https://github.com/urol-e5/timeseries/tree/master/timepoint_2)
	- This folder contains `scripts`, `data`, and `output` for physiological responses collected in March 2020. The scripts in this folder read in raw data files and conduct normalization and calculations for each response. Calculated responses are output in the `output` folder and analyzed as described in the `time_series_analysis` folder. Open this project with the `timepoint_2.Rproj` file. 
	
- [`timepoint_3`](https://github.com/urol-e5/timeseries/tree/master/timepoint_3)
	- This folder contains `scripts`, `data`, and `output` for physiological responses collected in September 2020. The scripts in this folder read in raw data files and conduct normalization and calculations for each response. Calculated responses are output in the `output` folder and analyzed as described in the `time_series_analysis` folder. Open this project with the `timepoint_3.Rproj` file. 

- [`timepoint_4`](https://github.com/urol-e5/timeseries/tree/master/timepoint_4)
	- This folder contains `scripts`, `data`, and `output` for physiological responses collected in November 2020. The scripts in this folder read in raw data files and conduct normalization and calculations for each response. Calculated responses are output in the `output` folder and analyzed as described in the `time_series_analysis` folder. Open this project with the `timepoint_4.Rproj` file. 

- [`ITS2_sequencing`](https://github.com/urol-e5/timeseries/tree/master/ITS2_sequencing)
	- This folder contains raw data files from Symportal analysis of ITS2 sequences. Data are analyzed in the `time_series_analysis` folder described below. 
	
- [`time_series_analysis`](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis)
	- Open this project with the `time_series_analysis.Rproj` file. This directory contains scripts, a `Figures` directory, and a `Output` directory and subfolders within these for each data type.  
	
	- [`its2_analysis.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/its2_analysis.Rmd) 
		- This script analyzes ITS2 sequence data produced through the SymPortal workflow. Output figures are stored in the [Figures/ITS2 folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/ITS2).  
		
	- [`environmental_data_analyses.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/environmental_data_analyses.Rmd)
		- This script analyzes available temperature, light, pH, and salinity data. Figures and summary data are stored in the  [Figures/Environmental folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/Environmental) and [Output folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Output). 

	- [`1_assemble_data.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/1_assemble_data.Rmd)
		- This script reads in all data files for each physiological response from each individual time point and assembles the data into one master dataframe. This script corrects colony ID errors, and QC's data, calculates different normalizations for each metric, and subsets data for other analyses. Master data frame is stored in the [Output directory](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Output). 

	- [`2_univariate_analyses.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/2_univariate_analyses.Rmd)
		- This script plots univariate response metrics for all physiological variables across site and time and species. Each response is analyzed using linear mixed effect models accounting for repeated measures. Figures are stored in the [Figures/Univariate folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/Univariate) and summaries of data are stored in the [Output directory](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Output). 

	- [`3_species_pca.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/3_species_pca.Rmd)
		- This script generates PCA plots and conducts PERMANOVA analyses of multivariate physiology by species. PCA's were also generated for each species for responses normalized to each possible normalizer (surface area, protein, and biomass). Figures are stored in the [Figures/NormalizerPCA folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/NormalizerPCA)

	- [`4_multivariate_analysis.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/4_multivariate_analysis.Rmd)
		- This script generates PCA plots with trajectory arrows to track change in group centroid between sites and across seasons. This script generates biplot loadings for these plots and conducts PERMANOVA analyses. Figures are stored in the [Figures/Multivariate folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/Multivariate). 

	- [`5_plasticity_analysis_travel_distance.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/5_plasticity_analysis_travel_distances.Rmd)
		- This scripts calculates a "plasticity" metric that is the total "distance" that species trajectories travel across time from plots generated in script 4. Figures are stored in the [Figures/Plasticity folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/Plasticity). Note that this analysis is kept in the repository but is not used in further analysis or presented in the publication. We instead pursued variance partitioning analyses as described below. 

	- [`6_modeling_analysis.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/6_modeling_analysis.Rmd)
		- This script conducts preliminary linear model selection analyses to analyze the relationships between physiological metrics and a phenotype-level response (calcification). 

	- [`7_variance_partitioning.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/7_variance_partitioning.Rmd)
		- This script conducts variance partitioning and redundancy analyses (RDA) on physiological metrics. Figures are stored in the [Figures/Multivariate folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/Multivariate). 

	- [`8_plasticity_analysis_centroid_distances.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/8_plasticity_analysis_centroid_distances.Rmd)
		- This script also conducts an analysis to generate a "plasticity" metric calculated as the distance between group centroids. Figures are stored in the [Figures/Multivariate folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/Multivariate). Note that this analysis is not used in further analyses or presented in the publication. We instead used variance partitioning analyses. 

	- [`9_multivariate_analysis_betadispersion.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/9_multivariate_analysis_betadispersion.Rmd)
		- This scripts conducts multivariate analyses of betadispersion. 

	- [`10_rda_analysis.Rmd`](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/10_rda_analysis.Rmd)
		- This script conducts distance-based redundancy analyses (dbRDA) to investigate the relationships between symbiont community composition and physiological responses. Figures are stored in the [Figures/ITS2 folder](https://github.com/urol-e5/timeseries/tree/master/time_series_analysis/Figures/ITS2). 

### Responses measured 

We measured the following physiological responses:  

<img width="554" alt="responses" src="https://github.com/urol-e5/timeseries/assets/32178010/733de652-9b2d-49c0-974d-30778f81fbb6">

### Running scripts 

**To run .Rmd scripts, set Knit directory to "Project Directory" prior to running full scripts. Restart your R environment prior to running scripts.** 

*To run all data from start to finish or re run all multivariate and univariate analyses after revising individual time point analyses:* 

1. Open the R project for the time points individually (e.g., `timepoint_1.Rproj` within each time point folder). Run all scripts in the `scripts` folders in individual time point folders. First run the `surface_area.Rmd` and then the `protein.Rmd` to generate output for normalizers needed in other data sets. Note that output of calculated metrics will be stored in the `output` folder within each time opint with the timepoint number in the file name: `1_surface_area.csv`. Run these scripts for all four time points in turn (timepoint 1, 2, 3, and 4). 

2. Open the `time_series_analysis.Rproj` project in the time series analysis folder. This project contains all scripts to read in data from each time point, conduct QC and metadata merging, and perform univariate and multivariate analyses. Run the scripts in this project in the numerical order they are named (e.g., `1_assemble_data.Rmd` then `2_univariate_analyses.Rmd`.   

*If you are looking to run individual steps of the analysis without changing any inputs, you can open and run scripts at any point in the workflow.*  

### Looking for colony metadata? 

QC'd colony metadata is available in the [time series analysis repository here](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/Output/master_timeseries.csv). It is critical to use this file for colony metadata. Colony metadata in individual time point folders includes errors and typos that are corrected in the time series analysis. 

Use the `colony_id_corr` column for the corrected colony_id. 

#### Data standards

- 'site'  *(site where sample was collected)*
	- 'site1' = Manava
	- 'site2' = Mahana
	- 'site3' = Hilton

- 'lon' *(longitude coordinate)*
	- Use decimal degrees: *dd.ddddd°*

- 'lat' *(latitude coordinate)*
	- Use decimal degrees: *dd.ddddd°*

- 'date'
	- Use format: *YYYYMMDD*

- 'time'
	- Use format: *hh:mm* with 24-hour clock

- 'species'. *(species of coral sampled)*
	- 'POC' = *Pocillopora*
	- 'POR' = *Porites*
	- 'ACR' = *Acropora*

- 'colony_id' 	*(coral tag number)*
	- Always include species and number separated by a dash, e.g. 'POC-201'

# Contact 

If you have any questions, contact Ariana Huffmyer at ashuffmyer (at) uri.edu. 
	

	
