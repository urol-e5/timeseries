  library(tidyverse)
  
  # Which samples are represented in which datasets?
  coral_metadata <- read_csv("../metadata/coral_metadata.csv") %>% 
    filter(!colony_id %in% c("POC-215", "POC-232", "ACR-360", "ACR-370", "ACR-376", "ACR-178", "ACR-347")) %>%   # know we didn't collect!
    pull(colony_id)
  
  homog_vols <- read_csv("data/1_homogenate_vols.csv") %>% 
    pull(colony_id)
  
  surface_area <- read_csv("data/1_Wax_dipping.csv") %>% 
    filter(!grepl("Standard", Sample)) %>%
    pull(colony_id)
  
  chl <- list.files("data/1_Chl/", "platemap", full.names = TRUE) %>%
    map_df(read_csv) %>%
    drop_na() %>%
    filter(!colony_id == "BK") %>%
    pull(colony_id)
  
  afdw <- read_csv("data/1_Biomass.csv") %>%
    drop_na() %>%
    pull(colony_id)
  
  molec <- read_csv("data/1_DNARNAShield_sample_metadata.csv") %>%
    distinct(colony_id) %>%
    pull()
  
  picurve <- list.files("data/1_PICurves/", "*.csv") %>%
    gsub("_.*.csv", "", .) %>%
    subset(!grepl("BK", .))
  
  sym_counts <- read_csv("data/1_Sym_Counts_Data.csv") %>%
    pull(colony_id)
  
  
  #####################################
  
  diffs <- function(x, y) {
    list(
      setdiff(x, y),
      setdiff(y, x)
    )
  }
  
  # Missing from DNA/RNA shield tubes
  diffs(coral_metadata, molec)
  
  # Missing from homogenate volumes
  diffs(coral_metadata, homog_vols)
  
  # Missing from surface area
  diffs(coral_metadata, surface_area)
  
  # Missing from chlorophyll
  diffs(coral_metadata, chl)
  
  # Missing from ash free dry weight
  diffs(coral_metadata, afdw)
  
  # Missing from PI curves
  diffs(coral_metadata, picurve)
  
  # Missing from sym_counts
   diffs(coral_metadata, sym_counts)
  

  
  
