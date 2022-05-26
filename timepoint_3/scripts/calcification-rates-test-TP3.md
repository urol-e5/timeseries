Calcification-test.Rmd
================
Emma Strand
3/10/2022

# Time point 3 - September

# Calcification test

We want to avoid manual entry when possible. So this script tests
loading in:  
- `raw_files` folder with raw data file (example =
`TA_Output_20220303_PutnamTitrations_PutnamLab`) -
`3_DeltaTA_metadata`  
- *finish this list*

### Load all packages needed

``` r
## install packages if you dont already have them in your library
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse') 
if ("broom" %in% rownames(installed.packages()) == 'FALSE') install.packages('broom') 
if ("purrr" %in% rownames(installed.packages()) == 'FALSE') install.packages('purrr') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("nlstools" %in% rownames(installed.packages()) == 'FALSE') install.packages('nlstools')
if ("stringr" %in% rownames(installed.packages()) == 'FALSE') install.packages('stringr') 

#Read in required libraries

library(broom)
library(purrr)
library(lubridate)
library(tidyverse)
library(nlstools)
library(stringr)
```

### Import raw data

Putnam_Sample_TA file (output from CSUN titrations) is only used for
salinity. We are pulling TA values straight from raw output files in
‘data/3_calcification/raw_files’ folder.  
I’m taking out initial 1 run 3 from CSUN for now (done at both CSUN and
URI) – circle back to which one to keep for analysis.

``` r
#bring in calcification data file with TA and chamber pH, temp, salinity measurements
raw.data <- list.files(path = 'data/3_calcification/raw_files', pattern = ".csv", full.names = TRUE) %>% 
  set_names(.) %>% 
  map_dfr(read.table, .id = "titration.run", header=TRUE, sep=",") %>%
  filter(!SampleID == "JUNK 1") %>% filter(!SampleID == "JUNK 2") %>% filter(!SampleID == "CRM") %>% 
  select(-Sample.Index, -TA_evap) %>%
  rename(Salinity.lab = Salinity) 

# load salinity values from CSUN
Putnam_Sample_TA <- read.csv("data/3_calcification/Putnam_Sample_TA.csv") %>% select(SampleID, Salinity) %>% filter(!SampleID == "Initial1_20200912_3")  

# merging CSUN and URI raw data 
raw.data <- full_join(raw.data, Putnam_Sample_TA, by = "SampleID") %>% 
  gather("salinity.origin", "Salinity.lab", 5:6) %>% # this will combine both salinity columns into one
  filter(!is.na(Salinity.lab)) %>% select(-salinity.origin) %>%
  separate(SampleID, c("colony_id", "Date", "Run.Number"), sep = "_", convert = FALSE) # keeps new column as old string (character)

# load field data for calcification runs 
deltaTA <- read.csv("data/3_calcification/3_DeltaTA_metadata.csv", header = TRUE, sep = ",",
                    colClasses=c(rep('character', 11), rep('numeric', 5), rep('character', 3))) %>%
           rename(Salinity.chamber = Salinity) %>% select(-Surface.Area)

# load surface area values 
SA <- read.csv("output/3_surface_area.csv") 
SA$surface.area.cm2 <- as.numeric(SA$surface.area.cm2)

# join all df together
data <- full_join(deltaTA, raw.data, by = c("colony_id", "Date", "Run.Number")) %>% 
  mutate(sample.type = case_when(
    startsWith(colony_id, "Initial") ~ "Initial",
    startsWith(colony_id, "BK") ~ "Blank",
    startsWith(colony_id, "P") ~ "Sample", #this covers POC, POR
    startsWith(colony_id, "A") ~ "Sample")) # this covers ACR

## fix mis-labeled coral fragments 
## ACR-143 to ACR-145
## POC-365 to POR-365 
## POC-367 to POR-367
## POR-387 to POR-385 
## POR-373 to POC-373
data$colony_id <- gsub("ACR-143", "ACR-145", data$colony_id)
data$colony_id <- gsub("POC-365", "POR-365", data$colony_id)
data$colony_id <- gsub("POC-367", "POR-367", data$colony_id)
data$colony_id <- gsub("POR-387", "POR-385", data$colony_id)
data$colony_id <- gsub("POR-373", "POC-373", data$colony_id)

data <- full_join(data, SA, by = "colony_id")

## read in corrected salinity values 
corr.salinity <- read.csv("data/3_calcification/3_corrected_salinity.csv", header=TRUE) %>%
  separate(colony_id, c("colony_id", "Date", "Run.Number"), sep = "_", convert = FALSE)

# merging correct salinity and data df
data <- left_join(data, corr.salinity)

# if there is an NA in the corrected salinity value, fill in with previous salinity value
data <- data %>%
  mutate(corrected_salinity = if_else(is.na(corrected_salinity), Salinity.lab, corrected_salinity),
         corrected_salinity2 = if_else(corrected_salinity > 37.5, 36.5, corrected_salinity)) ## circle back to make sure this is correct

data %>% subset(Run.Number == "3")
```

    ##         Date Run.Number Start.Time Stop.Time Chamber Channel Position
    ## 21  20200912          3      11:25     12:58       1       1        4
    ## 22  20200912          3      11:24     13:02       2       2        5
    ## 23  20200912          3      11:25     13:04       3       3        6
    ## 24  20200912          3      11:25     13:07       4       4        1
    ## 25  20200912          3      11:23     13:10       5       5        3
    ## 26  20200912          3      11:23     13:13       6       6        9
    ## 27  20200912          3      11:23     13:16       7       7        8
    ## 28  20200912          3      11:25     13:19       8       8        7
    ## 29  20200912          3      11:23     13:21       9       9       12
    ## 30  20200912          3      11:26     13:24      10      10       10
    ## 111 20200912          3       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 112 20200912          3       <NA>      <NA>    <NA>    <NA>     <NA>
    ##         Species colony_id Plug.ID Site Vol.ml pH.mV Temp.C Salinity.chamber
    ## 21      Porites   POR-354     354    3    595 -60.6   27.4            36.53
    ## 22      Porites   POR-340     340    3    580 -60.8   27.6            36.35
    ## 23      Porites   POR-385     387    3    610 -60.6   27.7            36.37
    ## 24  Pocillopora   POC-391     391    3    600 -60.2   27.6            36.33
    ## 25     Acropora   ACR-364     364    3    590 -60.5   27.8            36.28
    ## 26     Acropora   ACR-343     343    3    590 -60.2   27.8            36.45
    ## 27      Porites   POR-362     362    3    595 -60.5   27.5            36.27
    ## 28      Porites   POR-381     381    3    580 -58.6   27.7            36.30
    ## 29      Porites   POC-373     373    3    595 -61.9   27.9            36.67
    ## 30         BK-3      BK-3    BK-3    3    605 -59.9   27.5            36.24
    ## 111        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 112        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ##     TA.Start.Time TA.Stop.Time Notes
    ## 21          11:25        12:58      
    ## 22          11:24        13:02      
    ## 23          11:25        13:04      
    ## 24          11:25        13:07      
    ## 25          11:23        13:10      
    ## 26          11:23        13:13      
    ## 27          11:23        13:16      
    ## 28          11:25        13:19      
    ## 29          11:23        13:21      
    ## 30          11:26        13:24      
    ## 111          <NA>         <NA>  <NA>
    ## 112          <NA>         <NA>  <NA>
    ##                                                                             titration.run
    ## 21  data/3_calcification/raw_files/TA_Output_20220308_Run1_PutnamTitrations_PutnamLab.csv
    ## 22       data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 23       data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 24  data/3_calcification/raw_files/TA_Output_20220308_Run1_PutnamTitrations_PutnamLab.csv
    ## 25                                 data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv
    ## 26       data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 27       data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 28       data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 29       data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 30                                 data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv
    ## 111      data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ## 112      data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ##           TA    Mass Salinity.lab sample.type surface.area.cm2  timepoint
    ## 21  2267.029 59.9790        39.94      Sample         16.36341 timepoint3
    ## 22  2218.778 59.6260        40.16      Sample         20.61319 timepoint3
    ## 23  2295.439 60.4780        40.35      Sample         12.00889 timepoint3
    ## 24  2319.235 60.1190        39.73      Sample         43.91962 timepoint3
    ## 25  2231.662 59.9025        36.18      Sample         49.98752 timepoint3
    ## 26  2333.810 60.0810        40.35      Sample         26.30325 timepoint3
    ## 27  2265.813 60.0310        40.29      Sample         13.04889 timepoint3
    ## 28  2145.788 59.9340        40.30      Sample         24.77692 timepoint3
    ## 29  2286.196 60.1730        40.26      Sample         82.71005 timepoint3
    ## 30  2346.754 59.9889        36.28       Blank               NA       <NA>
    ## 111 2337.256 60.3790        40.34     Initial               NA       <NA>
    ## 112 2348.791 60.4460        40.41     Initial               NA       <NA>
    ##     corrected_salinity corrected_salinity2
    ## 21               38.12               36.50
    ## 22               38.32               36.50
    ## 23               38.25               36.50
    ## 24               38.34               36.50
    ## 25               36.18               36.18
    ## 26               37.79               36.50
    ## 27               36.57               36.57
    ## 28               37.35               37.35
    ## 29               38.30               36.50
    ## 30               36.28               36.28
    ## 111              38.02               36.50
    ## 112              38.22               36.50

*Volume of blanks — variations.* Blank 3 and Blank 4 are both from
chamber 10 but have different volumes – this shouldn’t be the case.

## Normalize TA values to salinity and separate blanks from samples.

``` r
# Normalize data to the salinity measurement taken at the same time as titration
data <- data %>% mutate(Ta.norm = TA * corrected_salinity/36,
                        Ta.norm_test = TA *corrected_salinity2/36)

## select for Initial 1 Run 3 output (comment out the filter function further above in this chunk prior to running this)
data %>% subset(colony_id == "Initial1") %>% subset(Run.Number == "3")
```

    ##         Date Run.Number Start.Time Stop.Time Chamber Channel Position Species
    ## 111 20200912          3       <NA>      <NA>    <NA>    <NA>     <NA>    <NA>
    ##     colony_id Plug.ID Site Vol.ml pH.mV Temp.C Salinity.chamber TA.Start.Time
    ## 111  Initial1    <NA> <NA>     NA    NA     NA               NA          <NA>
    ##     TA.Stop.Time Notes
    ## 111         <NA>  <NA>
    ##                                                                        titration.run
    ## 111 data/3_calcification/raw_files/TA_Output_20220303_PutnamTitrations_PutnamLab.csv
    ##           TA   Mass Salinity.lab sample.type surface.area.cm2 timepoint
    ## 111 2337.256 60.379        40.34     Initial               NA      <NA>
    ##     corrected_salinity corrected_salinity2  Ta.norm Ta.norm_test
    ## 111              38.02                36.5 2468.402     2369.718

``` r
# Calculate initial TA norm averages for each run
initial <- data %>% subset(sample.type == "Initial") %>% dplyr::group_by(Run.Number) %>%
  mutate(Ta.norm_initial_avg = mean(Ta.norm)) %>% select(Run.Number, Ta.norm_initial_avg) %>% distinct()

# Rejoining initial data back to the full df   
calc.data <- full_join(data, initial, by = "Run.Number") %>% subset(!sample.type == "Initial")

# Calculating delta TA for the blank samples
blanks <- calc.data %>% subset(sample.type == "Blank") %>%
  mutate(delta.TA.blank = Ta.norm_initial_avg - Ta.norm) %>% select(Run.Number, delta.TA.blank)

# Rejoining blank data back to the full df 
calc.data <- full_join(calc.data, blanks, by = "Run.Number") %>% subset(!sample.type == "Blank")   
```

#### Looking at variation in samples across blanks, initials, and titration data source (CSUN vs URI)

``` r
### those that are under 37 for the "problem runs": 1, 2, and 3 - could use these values?
data %>% subset(Run.Number == "1" | Run.Number == "2" | Run.Number == "3") %>% group_by(Run.Number) %>%
  filter(corrected_salinity < 37)
```

    ## # A tibble: 5 × 29
    ## # Groups:   Run.Number [3]
    ##   Date     Run.Number Start.Time Stop.Time Chamber Channel Position Species    
    ##   <chr>    <chr>      <chr>      <chr>     <chr>   <chr>   <chr>    <chr>      
    ## 1 20200909 1          11:57      13:49     6       6       9        Pocillopora
    ## 2 20200909 2          17:14      18:59     5       5       3        Pocillopora
    ## 3 20200912 3          11:23      13:10     5       5       3        Acropora   
    ## 4 20200912 3          11:23      13:16     7       7       8        Porites    
    ## 5 20200912 3          11:26      13:24     10      10      10       BK-3       
    ## # … with 21 more variables: colony_id <chr>, Plug.ID <chr>, Site <chr>,
    ## #   Vol.ml <dbl>, pH.mV <dbl>, Temp.C <dbl>, Salinity.chamber <dbl>,
    ## #   TA.Start.Time <chr>, TA.Stop.Time <chr>, Notes <chr>, titration.run <chr>,
    ## #   TA <dbl>, Mass <dbl>, Salinity.lab <dbl>, sample.type <chr>,
    ## #   surface.area.cm2 <dbl>, timepoint <chr>, corrected_salinity <dbl>,
    ## #   corrected_salinity2 <dbl>, Ta.norm <dbl>, Ta.norm_test <dbl>

``` r
### average salinity value from CSUN = 36.164
data %>% mutate(titration.location = case_when(
    endsWith(titration.run, "CSUN.csv") ~ "CSUN",
    endsWith(titration.run, "Lab.csv") ~ "URI")) %>%
  subset(titration.location == "CSUN") %>%
  mutate(mean=mean(corrected_salinity))
```

    ##         Date Run.Number Start.Time Stop.Time Chamber Channel Position
    ## 6   20200909          1      11:57     13:49       6       6        9
    ## 25  20200912          3      11:23     13:10       5       5        3
    ## 30  20200912          3      11:26     13:24      10      10       10
    ## 31  20200912          4      16:20     17:52       1       1        4
    ## 32  20200912          4      16:20     17:55       2       2        5
    ## 33  20200912          4      16:20     17:58       3       3        6
    ## 34  20200912          4      16:20     18:01       4       4        1
    ## 35  20200912          4      16:19     18:03       5       5        3
    ## 36  20200912          4      16:19     18:06       6       6        9
    ## 37  20200912          4      16:19     18:09       7       7        8
    ## 38  20200912          4      16:19     18:11       8       8        7
    ## 39  20200912          4      16:20     18:14       9       9       12
    ## 40  20200912          4      16:21     18:16      10      10       10
    ## 41  20200913          5      11:46     13:21       1       1        4
    ## 42  20200913          5      11:46     13:24       2       2        5
    ## 43  20200913          5      11:46     13:27       3       3        6
    ## 45  20200913          5      11:44     13:32       5       5        3
    ## 46  20200913          5      11:44     13:35       6       6        9
    ## 47  20200913          5      11:44     13:37       7       7        8
    ## 48  20200913          5      11:44     13:39       8       8        7
    ## 49  20200913          5      11:46     13:42       9       9       12
    ## 50  20200913          5      11:46     13:44      10      10       10
    ## 51  20200913          6      16:31     18:02       1       1        4
    ## 52  20200913          6      16:31     18:04       2       2        5
    ## 53  20200913          6      16:32     18:07       3       3        6
    ## 54  20200913          6      16:30     18:10       4       4        1
    ## 55  20200913          6      16:31     18:12       5       5        3
    ## 56  20200913          6      16:30     18:15       6       6        9
    ## 57  20200913          6      16:30     18:17       7       7        8
    ## 58  20200913          6      16:30     18:20       8       8        7
    ## 59  20200913          6      16:30     18:23       9       9       12
    ## 60  20200913          6      16:31     18:26      10      10       10
    ## 61  20200914          7      11:12     12:45       1       1        4
    ## 62  20200914          7      11:12     12:47       2       2        5
    ## 63  20200914          7      11:11     12:50       3       3        6
    ## 64  20200914          7      11:11     12:52       4       4        1
    ## 65  20200914          7      11:10     12:55       5       5        3
    ## 66  20200914          7      11:12     12:57       6       6        9
    ## 67  20200914          7      11:12     12:59       7       7        8
    ## 68  20200914          7      11:10     13:02       8       8        7
    ## 69  20200914          7      11:10     13:05       9       9       12
    ## 70  20200914          7      11:10     13:07      10      10       10
    ## 71  20200914          8      16:00     17:35       1       1        4
    ## 72  20200914          8      15:59     17:37       2       2        5
    ## 73  20200914          8      16:00     17:40       3       3        6
    ## 74  20200914          8      15:59     17:42       4       4        1
    ## 75  20200914          8      16:00     17:44       5       5        3
    ## 76  20200914          8      16:01     17:46       6       6        9
    ## 77  20200914          8      15:59     17:48       7       7        8
    ## 78  20200914          8      15:59     17:51       8       8        7
    ## 79  20200914          8      16:01     17:54       9       9       12
    ## 80  20200914          8      16:01     17:57      10      10       10
    ## 81  20200915          9      11:12     12:47       1       1        4
    ## 82  20200915          9      11:11     12:49       2       2        5
    ## 83  20200915          9      11:11     12:52       3       3        6
    ## 84  20200915          9      11:13     12:54       4       4        1
    ## 85  20200915          9      11:13     12:57       5       5        3
    ## 86  20200915          9      11:12     12:59       6       6        9
    ## 87  20200915          9      11:12     13:02       7       7        8
    ## 88  20200915          9      11:11     13:04       8       8        7
    ## 89  20200915          9      11:11     13:07       9       9       12
    ## 90  20200915          9      11:12     13:09      10      10       10
    ## 91  20200915         10      16:01     17:34       1       1        4
    ## 92  20200915         10      16:01     17:37       2       2        5
    ## 93  20200915         10      16:00     17:39       3       3        6
    ## 94  20200915         10      16:00     17:41       4       4        1
    ## 95  20200915         10      16:02     17:43       5       5        3
    ## 96  20200915         10      16:03     17:45       6       6        9
    ## 97  20200915         10      16:00     17:48       7       7        8
    ## 98  20200915         10      16:00     17:50       8       8        7
    ## 99  20200915         10      16:02     17:52       9       9       12
    ## 100 20200915         10      16:02     17:54      10      10       10
    ## 101 20200916         11      10:54     12:28       1       1        4
    ## 102 20200916         11      10:53     12:30       2       2        5
    ## 103 20200916         11      10:53     12:32       3       3        6
    ## 104 20200916         11      10:54     12:35       4       4        1
    ## 105 20200916         11      10:54     12:37       5       5        3
    ## 106 20200916         11      10:54     12:39       6       6        9
    ## 107 20200916         11      10:54     12:41       7       7        8
    ## 108 20200916         11      10:53     12:43       8       8        7
    ## 109 20200916         11      10:53     12:45       9       9       12
    ## 110 20200916         11      10:54     12:47      10      10       10
    ## 116 20200915          9       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 117 20200914          8       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 118 20200914          7       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 119 20200915         10       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 120 20200914          7       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 121 20200913          6       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 122 20200913          6       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 123 20200913          5       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 124 20200916         11       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 125 20200914          8       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 126 20200913          5       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 127 20200916         11       <NA>      <NA>    <NA>    <NA>     <NA>
    ## 128 20200915          9       <NA>      <NA>    <NA>    <NA>     <NA>
    ##         Species colony_id Plug.ID Site Vol.ml pH.mV Temp.C Salinity.chamber
    ## 6   Pocillopora   POC-371     371    3    590 -60.9 28.500            36.05
    ## 25     Acropora   ACR-364     364    3    590 -60.5 27.800            36.28
    ## 30         BK-3      BK-3    BK-3    3    605 -59.9 27.500            36.24
    ## 31     Acropora   ACR-229     229    2    595 -62.7 27.100            36.73
    ## 32  Pocillopora   POC-255     255    2    595 -61.8 27.100            36.49
    ## 33  Pocillopora   POC-219     219    2    590 -60.6 27.100            36.28
    ## 34  Pocillopora   POC-377     377    3    595 -60.6 27.000            36.34
    ## 35     Acropora   ACR-220     220    2    590 -60.1 27.000            36.30
    ## 36      Porites   POR-214     214    2    585 -61.3 27.000            36.38
    ## 37      Porites   POR-209     209    2    595 -60.9 27.000            36.46
    ## 38  Pocillopora   POC-375     375    3    595 -62.5 26.900            36.55
    ## 39  Pocillopora   POC-366     366    3    600 -63.3 27.100            36.42
    ## 40         BK-4      BK-4    BK-4    3    615 -60.7 26.800            36.32
    ## 41         BK-5      BK-5    BK-5    3    610 -61.3 27.629            36.47
    ## 42  Pocillopora   POC-200     200    2    600 -61.5 27.800            36.56
    ## 43     Acropora   ACR-210     210    2    605 -58.5 27.700            36.27
    ## 45     Acropora   ACR-241     241    2    610 -58.7 27.500            36.28
    ## 46  Pocillopora   POC-222     222    2    595 -63.0 27.600            36.35
    ## 47  Pocillopora   POC-201     201    2    590 -60.8 27.700            36.47
    ## 48      Porites   POR-245     245    2    595 -57.4 27.500            36.47
    ## 49     Acropora   ACR-244     244    2    605 -58.0 27.700            36.50
    ## 50      Porites   POR-242     242    2    580 -54.4 27.600            36.48
    ## 51         BK-6      BK-6    BK-6    2    615 -61.5 26.700            36.51
    ## 52      Porites   POR-251     251    2    600 -61.1 26.800            36.49
    ## 53      Porites   POR-262     262    2    610 -59.9 26.800            36.27
    ## 54      Porites   POR-253     253    2    600 -58.6 26.900            36.28
    ## 55  Pocillopora   POC-205     205    2    605 -59.2 26.900            36.29
    ## 56     Acropora   ACR-265     265    2    595 -59.5 26.700            36.51
    ## 57  Pocillopora   POC-239     239    2    605 -60.6 26.800            36.49
    ## 58     Acropora   ACR-213     213    2    605 -57.7 26.600            36.55
    ## 59      Porites   POR-221     221    2    605 -58.8 26.700            36.54
    ## 60  Pocillopora   POC-248     248    2    600 -58.1 26.700            36.40
    ## 61         BK-7      BK-7    BK-7    2    610 -59.2 28.400            36.63
    ## 62  Pocillopora   POC-238     238    2    605 -58.2 28.700            36.41
    ## 63      Porites   POR-224     224    2    605 -58.7 28.500            36.38
    ## 64      Porites   POR-216     216    2    600 -57.1 28.600            36.16
    ## 65  Pocillopora   POC-259     259    2    600 -60.6 28.700            36.42
    ## 66     Acropora   ACR-237     237    2    605 -58.4 28.800            36.51
    ## 67     Acropora   ACR-225     225    2    610 -57.3 28.700            36.41
    ## 68      Porites   POR-260     260    2    600 -58.3 28.700            36.50
    ## 69  Pocillopora   POC-257     257    2    605 -61.0 28.800            36.49
    ## 70     Acropora   ACR-218     210    2    590 -55.3 28.500            36.46
    ## 71         BK-8      BK-8    BK-8    2    615 -61.3 27.300            36.56
    ## 72     Acropora   ACR-175     175    1    610 -60.2 27.400            36.34
    ## 73     Acropora   ACR-150     150    1    615 -59.7 27.500            36.92
    ## 74  Pocillopora   POC-207     207    2    605 -60.2 27.500            36.39
    ## 75  Pocillopora    POC-57      57    1    600 -59.2 27.600            36.45
    ## 76  Pocillopora    POC-45      45    1    595 -59.3 27.500            36.54
    ## 77      Porites    POR-80      80    1    600 -61.2 27.400            36.56
    ## 78  Pocillopora   POC-254     254    2    600 -59.0 27.400            36.27
    ## 79      Porites    POR-81      81    1    590 -59.2 27.400            36.27
    ## 80      Porites    POR-76      76    1    590 -60.4 27.300            36.51
    ## 81         BK-9      BK-9    BK-9    1    615 -59.4 27.800            36.39
    ## 82  Pocillopora    POC-48      48    1    590 -59.9 28.000            36.68
    ## 83  Pocillopora    POC-52      52    1    595 -60.7 27.900            36.64
    ## 84  Pocillopora    POC-55      55    1    595 -59.1 27.900            36.77
    ## 85  Pocillopora    POC-40      40    1    595 -59.0 27.900            36.49
    ## 86  Pocillopora    POC-43      43    1    590 -59.5 27.900            36.50
    ## 87  Pocillopora    POC-44      44    1    600 -58.0 27.900            36.40
    ## 88  Pocillopora    POC-68      68    1    600 -57.4 27.900            36.50
    ## 89  Pocillopora    POC-53      53    1    590 -56.5 28.000            36.27
    ## 90  Pocillopora    POC-42      42    1    600 -58.9 27.900            36.62
    ## 91        BK-10     BK-10   BK-10    1    615 -61.9 27.100            36.39
    ## 92     Acropora   ACR-186     186    1    605 -59.9 27.200            36.76
    ## 93     Acropora   ACR-178     178    1    600 -59.2 27.300            36.53
    ## 94     Acropora   ACR-190     190    1    610 -59.7 27.100            36.31
    ## 95     Acropora   ACR-139     139    1    605 -59.1 27.200            36.66
    ## 96      Porites    POR-75      75    1    600 -59.9 27.100            36.57
    ## 97      Porites    POR-72      72    1    605 -60.2 27.100            36.49
    ## 98      Porites    POR-83      83    1    600 -59.1 27.100            36.74
    ## 99      Porites    POR-79      79    1    595 -59.6 27.100            36.28
    ## 100 Pocillopora    POC-47      47    1    595 -59.7 27.100            36.51
    ## 101     Porites    POR-69      69    1    595 -58.9 28.000            36.57
    ## 102     Porites    POR-73      73    1    600 -59.7 28.200            36.34
    ## 103    Acropora   ACR-173     173    1    605 -58.4 28.000            36.20
    ## 104       BK-11     BK-11   BK-11    1    620 -60.0 27.900            36.21
    ## 105     Porites    POR-82      82    1    600 -58.5 28.000            36.50
    ## 106     Porites    POR-74      74    1    595 -56.4 27.900            36.48
    ## 107     Porites    POR-70      70    1    590 -58.8 28.100            36.43
    ## 108     Porites    POR-71      71    1    600 -55.9 27.900            36.37
    ## 109    Acropora   ACR-140     140    1    605 -56.4 27.900            36.37
    ## 110    Acropora   ACR-145     143    1    600 -55.8 27.900            36.52
    ## 116        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 117        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 118        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 119        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 120        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ## 121        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ## 122        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 123        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ## 124        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ## 125        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ## 126        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 127        <NA>  Initial1    <NA> <NA>     NA    NA     NA               NA
    ## 128        <NA>  Initial2    <NA> <NA>     NA    NA     NA               NA
    ##     TA.Start.Time TA.Stop.Time Notes
    ## 6           11:57        13:49      
    ## 25          11:23        13:10      
    ## 30          11:26        13:24      
    ## 31          16:20        17:52      
    ## 32          16:20        17:55      
    ## 33          16:20        17:58      
    ## 34          16:20        18:01      
    ## 35          16:19        18:03      
    ## 36          16:19        18:06      
    ## 37          16:19        18:09      
    ## 38          16:19        18:11      
    ## 39          16:20        18:14      
    ## 40          16:21        18:16      
    ## 41          11:46        13:21      
    ## 42          11:46        13:24      
    ## 43          11:46        13:27      
    ## 45          11:44        13:32      
    ## 46          11:44        13:35      
    ## 47          11:44        13:37      
    ## 48          11:44        13:39      
    ## 49          11:46        13:42      
    ## 50          11:46        13:44      
    ## 51          16:31        18:02      
    ## 52          16:31        18:04      
    ## 53          16:32        18:07      
    ## 54          16:30        18:10      
    ## 55          16:31        18:12      
    ## 56          16:30        18:15      
    ## 57          16:30        18:17      
    ## 58          16:30        18:20      
    ## 59          16:30        18:23      
    ## 60          16:31        18:26      
    ## 61          11:12        12:45      
    ## 62          11:12        12:47      
    ## 63          11:11        12:50      
    ## 64          11:11        12:52      
    ## 65          11:10        12:55      
    ## 66          11:12        12:57      
    ## 67          11:12        12:59      
    ## 68          11:10        13:02      
    ## 69          11:10        13:05      
    ## 70          11:10        13:07      
    ## 71          16:00        17:35      
    ## 72          15:59        17:37      
    ## 73          16:00        17:40      
    ## 74          15:59        17:42      
    ## 75          16:00        17:44      
    ## 76          16:01        17:46      
    ## 77          15:59        17:48      
    ## 78          15:59        17:51      
    ## 79          16:01        17:54      
    ## 80          16:01        17:57      
    ## 81          11:12        12:47      
    ## 82          11:11        12:49      
    ## 83          11:11        12:52      
    ## 84          11:13        12:54      
    ## 85          11:13        12:57      
    ## 86          11:12        12:59      
    ## 87          11:12        13:02      
    ## 88          11:11        13:04      
    ## 89          11:11        13:07      
    ## 90          11:12        13:09      
    ## 91          16:01        17:34      
    ## 92          16:01        17:37      
    ## 93          16:00        17:39      
    ## 94          16:00        17:41      
    ## 95          16:02        17:43      
    ## 96          16:03        17:45      
    ## 97          16:00        17:48      
    ## 98          16:00        17:50      
    ## 99          16:02        17:52      
    ## 100         16:02        17:54      
    ## 101         10:54        12:28      
    ## 102         10:53        12:30      
    ## 103         10:53        12:32      
    ## 104         10:54        12:35      
    ## 105         10:54        12:37      
    ## 106         10:54        12:39      
    ## 107         10:54        12:41      
    ## 108         10:53        12:43      
    ## 109         10:53        12:45      
    ## 110         10:54        12:47      
    ## 116          <NA>         <NA>  <NA>
    ## 117          <NA>         <NA>  <NA>
    ## 118          <NA>         <NA>  <NA>
    ## 119          <NA>         <NA>  <NA>
    ## 120          <NA>         <NA>  <NA>
    ## 121          <NA>         <NA>  <NA>
    ## 122          <NA>         <NA>  <NA>
    ## 123          <NA>         <NA>  <NA>
    ## 124          <NA>         <NA>  <NA>
    ## 125          <NA>         <NA>  <NA>
    ## 126          <NA>         <NA>  <NA>
    ## 127          <NA>         <NA>  <NA>
    ## 128          <NA>         <NA>  <NA>
    ##                                              titration.run       TA    Mass
    ## 6   data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2272.142 59.9766
    ## 25  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2231.662 59.9025
    ## 30  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2346.754 59.9889
    ## 31  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2330.218 60.2157
    ## 32  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2328.536 59.9063
    ## 33  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2313.089 60.3524
    ## 34  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2329.505 59.8302
    ## 35  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2263.332 59.8119
    ## 36  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2253.585 59.8515
    ## 37  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2243.001 59.9955
    ## 38  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2294.772 59.6028
    ## 39  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2339.664 60.4108
    ## 40  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2353.593 59.8394
    ## 41  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2370.052 60.2492
    ## 42  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2316.281 60.5494
    ## 43  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2266.667 59.7542
    ## 45  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2317.167 60.2639
    ## 46  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2291.329 60.0259
    ## 47  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2280.257 59.9221
    ## 48  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2187.665 60.0428
    ## 49  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2305.212 60.1405
    ## 50  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2063.667 59.8281
    ## 51  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2355.166 59.8130
    ## 52  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2215.544 59.5779
    ## 53  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2309.144 60.1321
    ## 54  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2288.573 59.8592
    ## 55  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2347.324 59.7517
    ## 56  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2311.926 60.3983
    ## 57  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2331.014 59.6093
    ## 58  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2319.038 59.7731
    ## 59  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2280.655 60.2386
    ## 60  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2337.568 59.9144
    ## 61  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2340.889 59.7729
    ## 62  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2326.911 60.0789
    ## 63  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2292.691 60.4517
    ## 64  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2249.790 60.0695
    ## 65  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2335.429 59.7525
    ## 66  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2263.221 60.2248
    ## 67  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2315.557 59.8001
    ## 68  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2216.004 59.9607
    ## 69  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2266.521 60.0345
    ## 70  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2263.035 60.2002
    ## 71  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2323.287 60.0365
    ## 72  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2324.283 59.5282
    ## 73  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2302.141 60.1609
    ## 74  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2334.991 59.5757
    ## 75  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2342.115 60.3579
    ## 76  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2313.635 59.5629
    ## 77  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2300.964 59.7662
    ## 78  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2287.649 59.7061
    ## 79  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2264.843 60.2914
    ## 80  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2226.915 60.4601
    ## 81  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2346.377 60.2229
    ## 82  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2328.735 59.6577
    ## 83  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2302.834 59.8449
    ## 84  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2330.806 59.8678
    ## 85  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2304.699 59.6919
    ## 86  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2338.907 59.7265
    ## 87  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2345.213 59.5399
    ## 88  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2316.764 60.1294
    ## 89  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2300.858 60.0902
    ## 90  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2351.202 59.6550
    ## 91  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2354.515 60.3043
    ## 92  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2302.768 59.9428
    ## 93  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2282.768 59.5835
    ## 94  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2318.415 59.5506
    ## 95  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2300.485 60.1274
    ## 96  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2266.898 59.9771
    ## 97  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2292.240 60.0200
    ## 98  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2259.310 59.7957
    ## 99  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2259.508 59.8509
    ## 100 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2347.516 60.3074
    ## 101 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2201.071 59.8883
    ## 102 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2273.030 60.1943
    ## 103 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2287.242 59.8397
    ## 104 data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2354.485 60.4868
    ## 105 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2242.790 60.3408
    ## 106 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2243.620 60.0087
    ## 107 data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2293.703 60.0661
    ## 108 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2271.170 59.9535
    ## 109 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2344.319 59.8587
    ## 110 data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2309.475 59.7765
    ## 116 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2358.128 59.7355
    ## 117 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2345.932 60.1714
    ## 118 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2346.717 60.0851
    ## 119 data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2357.366 59.8400
    ## 120 data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2320.098 60.4887
    ## 121 data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2332.678 60.2730
    ## 122 data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2330.260 60.2292
    ## 123 data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2356.138 59.9884
    ## 124 data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2346.250 60.3172
    ## 125 data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2346.726 59.8898
    ## 126 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2363.889 60.0473
    ## 127 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2350.892 59.8138
    ## 128 data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2366.407 60.3133
    ##     Salinity.lab sample.type surface.area.cm2  timepoint corrected_salinity
    ## 6          36.21      Sample         47.10321 timepoint3              36.21
    ## 25         36.18      Sample         49.98752 timepoint3              36.18
    ## 30         36.28       Blank               NA       <NA>              36.28
    ## 31         36.26      Sample         18.48082 timepoint3              36.26
    ## 32         36.13      Sample         39.55013 timepoint3              36.13
    ## 33         36.16      Sample         40.98667 timepoint3              36.16
    ## 34         36.24      Sample         35.15071 timepoint3              36.24
    ## 35         36.21      Sample         35.82409 timepoint3              36.21
    ## 36         36.16      Sample         18.25262 timepoint3              36.16
    ## 37         36.11      Sample         16.27363 timepoint3              36.11
    ## 38         36.27      Sample         63.52620 timepoint3              36.27
    ## 39         36.20      Sample         60.76534 timepoint3              36.20
    ## 40         36.27       Blank               NA       <NA>              36.27
    ## 41         36.20       Blank               NA       <NA>              36.20
    ## 42         36.20      Sample         58.64045 timepoint3              36.20
    ## 43         36.16      Sample         25.77951 timepoint3              36.16
    ## 45         36.14      Sample         21.95995 timepoint3              36.14
    ## 46         36.09      Sample         53.32075 timepoint3              36.09
    ## 47         36.13      Sample         63.19699 timepoint3              36.13
    ## 48         36.19      Sample         28.27476 timepoint3              36.19
    ## 49         36.14      Sample         18.99708 timepoint3              36.14
    ## 50         36.29      Sample         28.31591 timepoint3              36.29
    ## 51         36.21       Blank               NA       <NA>              36.21
    ## 52         36.20      Sample         28.98555 timepoint3              36.20
    ## 53         36.13      Sample         10.60976 timepoint3              36.13
    ## 54         36.26      Sample         19.25147 timepoint3              36.26
    ## 55         36.20      Sample         24.24570 timepoint3              36.20
    ## 56         35.96      Sample         14.68370 timepoint3              35.96
    ## 57         36.22      Sample         49.71817 timepoint3              36.22
    ## 58         36.30      Sample         19.54326 timepoint3              36.30
    ## 59         36.13      Sample         20.49348 timepoint3              36.13
    ## 60         36.19      Sample         34.34640 timepoint3              36.19
    ## 61         36.06       Blank               NA       <NA>              36.06
    ## 62         36.10      Sample         33.41489 timepoint3              36.10
    ## 63         36.03      Sample         12.16975 timepoint3              36.03
    ## 64         36.11      Sample         15.74989 timepoint3              36.11
    ## 65         36.09      Sample         32.92108 timepoint3              36.09
    ## 66         36.08      Sample         22.15448 timepoint3              36.08
    ## 67         36.09      Sample         12.91047 timepoint3              36.09
    ## 68         35.93      Sample         27.39936 timepoint3              35.93
    ## 69         36.05      Sample         50.88910 timepoint3              36.05
    ## 70         36.12      Sample         28.89202 timepoint3              36.12
    ## 71         36.16       Blank               NA       <NA>              36.16
    ## 72         36.19      Sample         18.70902 timepoint3              36.19
    ## 73         36.18      Sample         32.49835 timepoint3              36.18
    ## 74         36.19      Sample         33.11187 timepoint3              36.19
    ## 75         36.19      Sample         42.25487 timepoint3              36.19
    ## 76         36.16      Sample         65.45281 timepoint3              36.16
    ## 77         36.11      Sample         26.90929 timepoint3              36.11
    ## 78         36.13      Sample         53.68737 timepoint3              36.13
    ## 79         36.25      Sample         16.96946 timepoint3              36.25
    ## 80         36.14      Sample         13.62500 timepoint3              36.14
    ## 81         36.25       Blank               NA       <NA>              36.25
    ## 82         36.17      Sample         65.91670 timepoint3              36.17
    ## 83         36.14      Sample         54.91816 timepoint3              36.14
    ## 84         36.17      Sample         41.33085 timepoint3              36.17
    ## 85         36.17      Sample         79.04387 timepoint3              36.17
    ## 86         36.25      Sample         44.46955 timepoint3              36.25
    ## 87         36.25      Sample         38.38294 timepoint3              36.25
    ## 88         36.23      Sample         36.66208 timepoint3              36.23
    ## 89         36.09      Sample         89.05852 timepoint3              36.09
    ## 90         36.13      Sample         28.91821 timepoint3              36.13
    ## 91         36.15       Blank               NA       <NA>              36.15
    ## 92         36.20      Sample         45.35616 timepoint3              36.20
    ## 93         36.21      Sample         38.80567 timepoint3              36.21
    ## 94         36.27      Sample         23.13088 timepoint3              36.27
    ## 95         36.13      Sample         34.25287 timepoint3              36.13
    ## 96         36.17      Sample         17.02557 timepoint3              36.17
    ## 97         36.15      Sample         16.02298 timepoint3              36.15
    ## 98         36.27      Sample         20.79276 timepoint3              36.27
    ## 99         36.04      Sample         17.64658 timepoint3              36.04
    ## 100        36.19      Sample         37.33920 timepoint3              36.19
    ## 101        36.20      Sample         21.39880 timepoint3              36.20
    ## 102        36.22      Sample         17.41089 timepoint3              36.22
    ## 103        36.07      Sample         41.44682 timepoint3              36.07
    ## 104        36.14       Blank               NA       <NA>              36.14
    ## 105        36.18      Sample         15.46183 timepoint3              36.18
    ## 106        36.16      Sample         15.86960 timepoint3              36.16
    ## 107        36.13      Sample         23.22067 timepoint3              36.13
    ## 108        36.03      Sample         16.94327 timepoint3              36.03
    ## 109        36.12      Sample         24.37289 timepoint3              36.12
    ## 110        36.16      Sample         11.88544 timepoint3              36.16
    ## 116        36.13     Initial               NA       <NA>              36.13
    ## 117        36.23     Initial               NA       <NA>              36.23
    ## 118        36.17     Initial               NA       <NA>              36.17
    ## 119        36.28     Initial               NA       <NA>              36.28
    ## 120        36.05     Initial               NA       <NA>              36.05
    ## 121        36.21     Initial               NA       <NA>              36.21
    ## 122        36.20     Initial               NA       <NA>              36.20
    ## 123        36.11     Initial               NA       <NA>              36.11
    ## 124        36.21     Initial               NA       <NA>              36.21
    ## 125        36.19     Initial               NA       <NA>              36.19
    ## 126        36.13     Initial               NA       <NA>              36.13
    ## 127        36.15     Initial               NA       <NA>              36.15
    ## 128        36.10     Initial               NA       <NA>              36.10
    ##     corrected_salinity2  Ta.norm Ta.norm_test titration.location   mean
    ## 6                 36.21 2285.396     2285.396               CSUN 36.164
    ## 25                36.18 2242.820     2242.820               CSUN 36.164
    ## 30                36.28 2365.006     2365.006               CSUN 36.164
    ## 31                36.26 2347.047     2347.047               CSUN 36.164
    ## 32                36.13 2336.944     2336.944               CSUN 36.164
    ## 33                36.16 2323.369     2323.369               CSUN 36.164
    ## 34                36.24 2345.035     2345.035               CSUN 36.164
    ## 35                36.21 2276.535     2276.535               CSUN 36.164
    ## 36                36.16 2263.601     2263.601               CSUN 36.164
    ## 37                36.11 2249.855     2249.855               CSUN 36.164
    ## 38                36.27 2311.983     2311.983               CSUN 36.164
    ## 39                36.20 2352.663     2352.663               CSUN 36.164
    ## 40                36.27 2371.245     2371.245               CSUN 36.164
    ## 41                36.20 2383.219     2383.219               CSUN 36.164
    ## 42                36.20 2329.149     2329.149               CSUN 36.164
    ## 43                36.16 2276.741     2276.741               CSUN 36.164
    ## 45                36.14 2326.178     2326.178               CSUN 36.164
    ## 46                36.09 2297.057     2297.057               CSUN 36.164
    ## 47                36.13 2288.491     2288.491               CSUN 36.164
    ## 48                36.19 2199.211     2199.211               CSUN 36.164
    ## 49                36.14 2314.177     2314.177               CSUN 36.164
    ## 50                36.29 2080.291     2080.291               CSUN 36.164
    ## 51                36.21 2368.904     2368.904               CSUN 36.164
    ## 52                36.20 2227.852     2227.852               CSUN 36.164
    ## 53                36.13 2317.482     2317.482               CSUN 36.164
    ## 54                36.26 2305.102     2305.102               CSUN 36.164
    ## 55                36.20 2360.365     2360.365               CSUN 36.164
    ## 56                35.96 2309.357     2309.357               CSUN 36.164
    ## 57                36.22 2345.259     2345.259               CSUN 36.164
    ## 58                36.30 2338.363     2338.363               CSUN 36.164
    ## 59                36.13 2288.891     2288.891               CSUN 36.164
    ## 60                36.19 2349.905     2349.905               CSUN 36.164
    ## 61                36.06 2344.790     2344.790               CSUN 36.164
    ## 62                36.10 2333.374     2333.374               CSUN 36.164
    ## 63                36.03 2294.601     2294.601               CSUN 36.164
    ## 64                36.11 2256.664     2256.664               CSUN 36.164
    ## 65                36.09 2341.267     2341.267               CSUN 36.164
    ## 66                36.08 2268.251     2268.251               CSUN 36.164
    ## 67                36.09 2321.345     2321.345               CSUN 36.164
    ## 68                35.93 2211.695     2211.695               CSUN 36.164
    ## 69                36.05 2269.669     2269.669               CSUN 36.164
    ## 70                36.12 2270.579     2270.579               CSUN 36.164
    ## 71                36.16 2333.613     2333.613               CSUN 36.164
    ## 72                36.19 2336.550     2336.550               CSUN 36.164
    ## 73                36.18 2313.652     2313.652               CSUN 36.164
    ## 74                36.19 2347.314     2347.314               CSUN 36.164
    ## 75                36.19 2354.476     2354.476               CSUN 36.164
    ## 76                36.16 2323.918     2323.918               CSUN 36.164
    ## 77                36.11 2307.995     2307.995               CSUN 36.164
    ## 78                36.13 2295.910     2295.910               CSUN 36.164
    ## 79                36.25 2280.571     2280.571               CSUN 36.164
    ## 80                36.14 2235.576     2235.576               CSUN 36.164
    ## 81                36.25 2362.672     2362.672               CSUN 36.164
    ## 82                36.17 2339.732     2339.732               CSUN 36.164
    ## 83                36.14 2311.790     2311.790               CSUN 36.164
    ## 84                36.17 2341.813     2341.813               CSUN 36.164
    ## 85                36.17 2315.582     2315.582               CSUN 36.164
    ## 86                36.25 2355.150     2355.150               CSUN 36.164
    ## 87                36.25 2361.500     2361.500               CSUN 36.164
    ## 88                36.23 2331.566     2331.566               CSUN 36.164
    ## 89                36.09 2306.611     2306.611               CSUN 36.164
    ## 90                36.13 2359.693     2359.693               CSUN 36.164
    ## 91                36.15 2364.326     2364.326               CSUN 36.164
    ## 92                36.20 2315.561     2315.561               CSUN 36.164
    ## 93                36.21 2296.084     2296.084               CSUN 36.164
    ## 94                36.27 2335.803     2335.803               CSUN 36.164
    ## 95                36.13 2308.793     2308.793               CSUN 36.164
    ## 96                36.17 2277.602     2277.602               CSUN 36.164
    ## 97                36.15 2301.791     2301.791               CSUN 36.164
    ## 98                36.27 2276.255     2276.255               CSUN 36.164
    ## 99                36.04 2262.019     2262.019               CSUN 36.164
    ## 100               36.19 2359.905     2359.905               CSUN 36.164
    ## 101               36.20 2213.299     2213.299               CSUN 36.164
    ## 102               36.22 2286.921     2286.921               CSUN 36.164
    ## 103               36.07 2291.690     2291.690               CSUN 36.164
    ## 104               36.14 2363.641     2363.641               CSUN 36.164
    ## 105               36.18 2254.004     2254.004               CSUN 36.164
    ## 106               36.16 2253.592     2253.592               CSUN 36.164
    ## 107               36.13 2301.985     2301.985               CSUN 36.164
    ## 108               36.03 2273.063     2273.063               CSUN 36.164
    ## 109               36.12 2352.133     2352.133               CSUN 36.164
    ## 110               36.16 2319.739     2319.739               CSUN 36.164
    ## 116               36.13 2366.643     2366.643               CSUN 36.164
    ## 117               36.23 2360.920     2360.920               CSUN 36.164
    ## 118               36.17 2357.799     2357.799               CSUN 36.164
    ## 119               36.28 2375.701     2375.701               CSUN 36.164
    ## 120               36.05 2323.320     2323.320               CSUN 36.164
    ## 121               36.21 2346.286     2346.286               CSUN 36.164
    ## 122               36.20 2343.206     2343.206               CSUN 36.164
    ## 123               36.11 2363.337     2363.337               CSUN 36.164
    ## 124               36.21 2359.936     2359.936               CSUN 36.164
    ## 125               36.19 2359.112     2359.112               CSUN 36.164
    ## 126               36.13 2372.425     2372.425               CSUN 36.164
    ## 127               36.15 2360.687     2360.687               CSUN 36.164
    ## 128               36.10 2372.980     2372.980               CSUN 36.164

``` r
data %>% filter(!is.na(Ta.norm)) %>% gather("TA_test", "value", 28:29) %>%
  ggplot(aes(x=Run.Number, y=value, group = sample.type, color = sample.type)) + facet_wrap(~TA_test) + theme_bw() + geom_point()
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
data %>% filter(!is.na(Ta.norm)) %>%
  ggplot(aes(x=Run.Number, y=Ta.norm, group = sample.type, color = sample.type)) + theme_bw() + geom_point()
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
data %>% subset(!sample.type == "Sample") %>%
  ggplot(aes(x=Run.Number, y=Ta.norm, group = sample.type, color = sample.type)) + theme_bw() + geom_point()
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
data %>% filter(!is.na(Ta.norm)) %>% subset(sample.type=="Sample") %>% mutate(titration.location = case_when(
    endsWith(titration.run, "CSUN.csv") ~ "CSUN",
    endsWith(titration.run, "Lab.csv") ~ "URI")) %>%
  ggplot(aes(x=sample.type, y=Ta.norm, color = titration.location))  + theme_bw() + geom_point() + facet_grid(~Species)
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
data %>% filter(!is.na(Ta.norm)) %>%
  ggplot(aes(x=Mass, y=Ta.norm, group = sample.type, color = Run.Number)) + theme_bw() + geom_point() + facet_wrap(~sample.type)
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

## Calculate Net Ecosystem Calcification (NEC) from total alkalinity method

1.  Convert time from character to time for both start and stop time.  
2.  Calculate deltaTA (for each coral) as the sample TA.norm subtracted
    from the initial from that run, then subtracting the blank from that
    run.  
3.  Calculate the time difference between start and stop times for each
    chamber(coral).  
4.  Calculate calcification rate according the Net Ecosystem
    Calcification rates equation, normalizing to the surface area of the
    coral.  
5.  Make any negative values equal to zero.  
6.  Add the time point for the exported df (to then bring into a larger
    df with all time points in a different script).

``` r
NEC <- calc.data %>%
  mutate(start.time = strptime(as.character(TA.Start.Time), "%H:%M")) %>% #convert time from character to time
  mutate(stop.time = strptime(as.character(TA.Stop.Time), "%H:%M")) %>%
  mutate(deltaTA = (Ta.norm_initial_avg - Ta.norm) - delta.TA.blank) %>% #calculate difference in TA corrected to blanks
  mutate(timediff = as.numeric(stop.time - start.time)) %>% #calculate time difference
  mutate(umol.cm2.hr = (deltaTA/2)*(1.023)*(Vol.ml/surface.area.cm2)*(1/timediff)*(1/1000)) %>% # calculate net ecosystem calcification rates
  mutate(umol.cm2.hr=if_else(umol.cm2.hr<0, 0, umol.cm2.hr)) %>% #make any negative values equal to 0
  mutate(timepoint="timepoint3") 

write.csv(NEC, 'output/3_calcification_rates.csv') 

NEC %>% filter(is.na(umol.cm2.hr)) # missing calc rates 
```

    ##        Date Run.Number Start.Time Stop.Time Chamber Channel Position
    ## 1  20200912          4      16:20     17:52       1       1        4
    ## 2  20200912          4      16:20     17:55       2       2        5
    ## 3  20200912          4      16:20     17:58       3       3        6
    ## 4  20200912          4      16:20     18:01       4       4        1
    ## 5  20200912          4      16:19     18:03       5       5        3
    ## 6  20200912          4      16:19     18:06       6       6        9
    ## 7  20200912          4      16:19     18:09       7       7        8
    ## 8  20200912          4      16:19     18:11       8       8        7
    ## 9  20200912          4      16:20     18:14       9       9       12
    ## 10 20200913          5      11:47     13:29       4       4        1
    ##        Species colony_id Plug.ID Site Vol.ml pH.mV Temp.C Salinity.chamber
    ## 1     Acropora   ACR-229     229    2    595 -62.7   27.1            36.73
    ## 2  Pocillopora   POC-255     255    2    595 -61.8   27.1            36.49
    ## 3  Pocillopora   POC-219     219    2    590 -60.6   27.1            36.28
    ## 4  Pocillopora   POC-377     377    3    595 -60.6   27.0            36.34
    ## 5     Acropora   ACR-220     220    2    590 -60.1   27.0            36.30
    ## 6      Porites   POR-214     214    2    585 -61.3   27.0            36.38
    ## 7      Porites   POR-209     209    2    595 -60.9   27.0            36.46
    ## 8  Pocillopora   POC-375     375    3    595 -62.5   26.9            36.55
    ## 9  Pocillopora   POC-366     366    3    600 -63.3   27.1            36.42
    ## 10     Porites   POR-240     240    2    590 -59.8   27.6            36.28
    ##    TA.Start.Time TA.Stop.Time Notes
    ## 1          16:20        17:52      
    ## 2          16:20        17:55      
    ## 3          16:20        17:58      
    ## 4          16:20        18:01      
    ## 5          16:19        18:03      
    ## 6          16:19        18:06      
    ## 7          16:19        18:09      
    ## 8          16:19        18:11      
    ## 9          16:20        18:14      
    ## 10         11:47        13:29      
    ##                                             titration.run       TA    Mass
    ## 1  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2330.218 60.2157
    ## 2  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2328.536 59.9063
    ## 3  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2313.089 60.3524
    ## 4  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2329.505 59.8302
    ## 5  data/3_calcification/raw_files/TA2021-10-20 - CSUN.csv 2263.332 59.8119
    ## 6  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2253.585 59.8515
    ## 7  data/3_calcification/raw_files/TA2021-10-15 - CSUN.csv 2243.001 59.9955
    ## 8  data/3_calcification/raw_files/TA2021-10-08 - CSUN.csv 2294.772 59.6028
    ## 9  data/3_calcification/raw_files/TA2021-10-21 - CSUN.csv 2339.664 60.4108
    ## 10                                                   <NA>       NA      NA
    ##    Salinity.lab sample.type surface.area.cm2  timepoint corrected_salinity
    ## 1         36.26      Sample         18.48082 timepoint3              36.26
    ## 2         36.13      Sample         39.55013 timepoint3              36.13
    ## 3         36.16      Sample         40.98667 timepoint3              36.16
    ## 4         36.24      Sample         35.15071 timepoint3              36.24
    ## 5         36.21      Sample         35.82409 timepoint3              36.21
    ## 6         36.16      Sample         18.25262 timepoint3              36.16
    ## 7         36.11      Sample         16.27363 timepoint3              36.11
    ## 8         36.27      Sample         63.52620 timepoint3              36.27
    ## 9         36.20      Sample         60.76534 timepoint3              36.20
    ## 10           NA      Sample         30.46324 timepoint3                 NA
    ##    corrected_salinity2  Ta.norm Ta.norm_test Ta.norm_initial_avg delta.TA.blank
    ## 1                36.26 2347.047     2347.047                  NA             NA
    ## 2                36.13 2336.944     2336.944                  NA             NA
    ## 3                36.16 2323.369     2323.369                  NA             NA
    ## 4                36.24 2345.035     2345.035                  NA             NA
    ## 5                36.21 2276.535     2276.535                  NA             NA
    ## 6                36.16 2263.601     2263.601                  NA             NA
    ## 7                36.11 2249.855     2249.855                  NA             NA
    ## 8                36.27 2311.983     2311.983                  NA             NA
    ## 9                36.20 2352.663     2352.663                  NA             NA
    ## 10                  NA       NA           NA            2367.881      -15.33736
    ##             start.time           stop.time deltaTA timediff umol.cm2.hr
    ## 1  2022-05-26 16:20:00 2022-05-26 17:52:00      NA 1.533333          NA
    ## 2  2022-05-26 16:20:00 2022-05-26 17:55:00      NA 1.583333          NA
    ## 3  2022-05-26 16:20:00 2022-05-26 17:58:00      NA 1.633333          NA
    ## 4  2022-05-26 16:20:00 2022-05-26 18:01:00      NA 1.683333          NA
    ## 5  2022-05-26 16:19:00 2022-05-26 18:03:00      NA 1.733333          NA
    ## 6  2022-05-26 16:19:00 2022-05-26 18:06:00      NA 1.783333          NA
    ## 7  2022-05-26 16:19:00 2022-05-26 18:09:00      NA 1.833333          NA
    ## 8  2022-05-26 16:19:00 2022-05-26 18:11:00      NA 1.866667          NA
    ## 9  2022-05-26 16:20:00 2022-05-26 18:14:00      NA 1.900000          NA
    ## 10 2022-05-26 11:47:00 2022-05-26 13:29:00      NA 1.700000          NA

#### Check sample variation one more time post normalizing to surface area

#### TA as a function of chamber temperature, salinity, pH, initial measurements temperature and pH

Do we need to do any correction for this?

``` r
NEC %>% filter(!is.na(umol.cm2.hr)) %>%
  ggplot(aes(x=Run.Number, y=umol.cm2.hr, color = Species)) + theme_bw() + geom_point() + facet_grid(~Species)
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
NEC %>% filter(!is.na(umol.cm2.hr)) %>%
  ggplot(aes(x=Temp.C, y=umol.cm2.hr, color = Species)) + theme_bw() + geom_point() 
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
NEC %>% filter(!is.na(umol.cm2.hr)) %>%
  ggplot(aes(x=pH.mV, y=umol.cm2.hr, color = Species)) + theme_bw() + geom_point() 
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
NEC %>% filter(!is.na(umol.cm2.hr)) %>%
  ggplot(aes(x=Salinity.chamber, y=umol.cm2.hr, color = Species)) + theme_bw() + geom_point() 
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
NEC %>% filter(!is.na(umol.cm2.hr)) %>%
  ggplot(aes(x=Chamber, y=umol.cm2.hr, color = Species)) + theme_bw() + geom_point() # this accounts for chamber and channel # 
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->

``` r
NEC %>% filter(!is.na(umol.cm2.hr)) %>%
  ggplot(aes(x=titration.run, y=umol.cm2.hr, color = Species)) + theme_bw() + geom_point()  
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->

## Visuals

``` r
NEC %>%
  filter(!Run.Number=="NA")%>%
  filter(!Species=="NA")%>%
  filter(!Site=="NA")%>%
  ggplot(aes(x = as.factor(Site), y = umol.cm2.hr, color = Species), position=position_dodge(0.5)) +
  geom_boxplot() +
  geom_point(size=3, position=position_jitterdodge(0.6))+
  facet_wrap(~Species) +
  labs(x = "Site", y = expression(paste(mu,"mol Ca", CO[3], " cm"^-2, " hr"^-1))) +                                       # Plot all points
  theme_classic() 
```

![](calcification-rates-test-TP3_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Write data to file for time series analysis.

``` r
# write calcification data to file for use in time series output
NEC %>%
  select(colony_id, umol.cm2.hr) %>%
  mutate(timepoint="timepoint3")%>%
  filter(!is.na(umol.cm2.hr))%>%
  write_csv(path = "output/3_calcification_output.csv")
```
