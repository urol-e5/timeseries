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

# List files with cholorophyll data
files <- list.files(path = "data/1_Chl", pattern = "*.csv", full.names = TRUE)

# Read in all files into tibble
df <- tibble(file = files) %>%
  mutate(chl_data = map(file, ~ read_chl(.)))

  
df$chl_data[[1]]
  

