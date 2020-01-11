library(tidyverse)


# Function to read in chl data
read_chl <- function(file) {
  chl_data <- read_csv(file, skip = 24, n_max = 24) %>%
    select(-1) %>%
    magrittr::set_colnames(c("row", 1:12, "wavelength")) %>%
    fill(row) %>%
    gather("col", "absorbance", -wavelength, -row) %>%
    unite("well", c(row, col), sep = "")
}

# List chlorophyll data files
chl_path = "data/1_Chl/"                                                 # Path to chl data directory
all_chl_files <- list.files(path = chl_path, pattern = "*.csv")          # List all files in directory
chl_platemaps <- list.files(path = chl_path, pattern = "platemap")       # List platemap files
chl_data_files <- setdiff(all_chl_files, chl_platemaps)                  # List absorbance data files

# Read in all files into tibble
df <- tibble(file = chl_data_files) %>%
  #filter(chl_data_files == "1_20200106_Chl.csv") %>%
  mutate(platemap = map(file, ~ read_csv(paste0(chl_path, tools::file_path_sans_ext(.), "_platemap.csv"))),
         chl_data = map(file, ~ read_chl(paste0(chl_path, .))))

# Merge platemap and data for each plate
df <- df %>%
  mutate(merged = map2(platemap, chl_data, ~ right_join(.x, .y)))


####

# average all technical replicates for each colony/wavelength, including all acetone blanks together
df <- df %>%
  unnest(merged) %>%
  filter(!is.na(colony_id)) %>%                         # remove empty wells (colony_id is NA)
  group_by(file, colony_id, wavelength) %>%
  summarise(n = n(), mean_abs = mean(absorbance)) %>%
  spread(wavelength, mean_abs)

# get the acetone blank 750 absorbace for each file (i.e., plate), and subtract from 630 and 663 values for each sample
df <- df %>%
  group_by(file) %>%
  mutate(blank750 = `750`[colony_id == "BK"]) %>%
  ungroup() %>%
  mutate(adj630 = `630` - blank750,
         adj663 = `663` - blank750)

# calculate chla and chlc2 values based on equations from Jeffrey and Humphrey 1975
#units µg/ml
df <- df %>%
  mutate(chla = 11.43 * adj663 - 0.64 * adj630,
        chlc2 = 27.09 * adj630 - 3.63 * adj663)

#Multiply by the homogenate volume
#Load homogenate volume
homog.vol <- read.csv("data/1_homogenate_vols.csv", header=TRUE)
Meta <- read.csv("data/1_Sample_Info.csv")
Data <-merge(df, Meta, by="colony_id")
Data <- merge(Data, homog.vol, by="colony_id")

#calculate surface area standard curve
wax.data <- read.csv("data/1_Wax_dipping.csv", header=TRUE)
wax.data$delta.mass.g <- wax.data$weight2.g-wax.data$weight1.g
stnds <- subset(wax.data, Sample=="Standard")
stnds <- stnds[-1,]
stnds$rad <- stnds$Diameter/2
stnds$surface.area.cm2 <- 4*pi*(stnds$rad)^2
stnd.curve <- lm(surface.area.cm2~delta.mass.g, data=stnds)
plot(surface.area.cm2~delta.mass.g, data=stnds)

stnd.curve$coefficients
summary(stnd.curve)$r.squared

#Calculate surface area
smpls <- subset(wax.data, Sample=="Coral")
smpls$surface.area.cm2 <- stnd.curve$coefficients[2] * smpls$delta.mass.g + stnd.curve$coefficients[1]
plot(smpls$surface.area.cm2)
range(smpls$surface.area.cm2)
range(stnds$surface.area.cm2)

Data <- merge(Data, smpls, by="colony_id")

Data$chla.cm2 <-  (Data$chla * Data$homog_vol_ml)/Data$surface.area.cm2
Data$chlc2.cm2 <-  (Data$chlc2 * Data$homog_vol_ml)/Data$surface.area.cm2

Data%>%
  group_by(Species, Site)%>%
  summarise(mean.value = mean(chla.cm2), se = std.error(chla.cm2)) %>%
  ggplot(aes(x = Site, y = mean.value, group = Species, color = Species))+
  ylab("CHL-a µg cm-2")+
  geom_point(size = 3)+
  geom_errorbar(aes(x = Site, ymin = mean.value-se, ymax = mean.value+se), width = 0.5)+
  facet_grid(~Species, scales = "free_y")

Data%>%
  group_by(Species, Site)%>%
  summarise(mean.value = mean(chlc2.cm2 ), se = std.error(chlc2.cm2 )) %>%
  ggplot(aes(x = Site, y = mean.value, group = Species, color = Species))+
  ylab("CHL-c2 µg cm-2")+
  geom_point(size = 3)+
  geom_errorbar(aes(x = Site, ymin = mean.value-se, ymax = mean.value+se), width = 0.5)+
  facet_grid(~Species, scales = "free_y")
