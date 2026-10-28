## N uptake estimate filtering and corrections
# JAK and JRB

## Import packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse"), require, character.only=T)

## Import data
dat <- read.csv("../om_corrected_data/N_Incubation_Uptake_Rates_20261019.csv")
colnames(dat)

#### FIT MM MODELS WITH FILTERING ####

## Now, we are going to manually filter outliers after looking at plots
# This was decided because the previous filtering was biased to removing
# the background concentrations since they had >20% deviation
# Filtering was done be JB and JK visually inspecting data
dat$flag <- (
  (dat$Site == "BW0.5m" & dat$Inc_month == "July" & dat$Type == "biofilm" & dat$Analyte == "NH3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "BW0.5m" & dat$Inc_month == "July" & dat$Type == "biofilm" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "BW0.5m" & dat$Inc_month == "July" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0 & dat$Spike_µg_L != 0) |
    (dat$Site == "BW0.5m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NH3" & dat$net_delta_Conc_ugNLhr_OM > 200) |
    (dat$Site == "BW0.5m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM > 200) |
    (dat$Site == "BW10m" & dat$Inc_month == "June" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$Spike_µg_L == 800) |
    (dat$Site == "BW3m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "GB0.5m" & dat$Inc_month == "May" & dat$Type == "biofilm" & dat$Analyte == "NO3" & dat$Spike_µg_L == 800) |
    (dat$Site == "GB3m" & dat$Inc_month == "June" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < -50) |
    (dat$Site == "GB3m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "GB10m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "GB10m" & dat$Inc_month == "June" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM > 100) |
    (dat$Site == "SS3m" & dat$Inc_month == "June" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM > 40) |
    (dat$Site == "SS3m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "SH3m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$Spike_µg_L == 800)
)

# count flags
sum(dat$flag, na.rm = TRUE)

# visualize
dat %>%
  filter(flag) %>%
  count(Site, Inc_month, Analyte, Type, name = "Filtered_Count")


### Now that we have visualized the flagged outliers
# Create a new dataframe excluding the flagged rows
flagged_rows <- dat %>% filter(flag == TRUE)
dat_filtered <- dat %>% filter(flag == FALSE)

# create vertical line for "true" spike conc 
# read the p val and ver

# Preprocess data to ensure consistency in 'Site' names
dat_filtered <- dat_filtered %>%
  mutate(Site = gsub("_", "", Site))

## Save dat_filtered csv
write.csv(dat_filtered, "../om_corrected_data/Filtered_NUptake_2026_10_28.csv")
