## N uptake estimate filtering and corrections
# JAK and JRB

## Import packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse"), require, character.only=T)

## Import data
Incubation_all <- readRDS("../data/N_Incubation_Uptake_Rates_20240924.rds")

#### CORRECT FOR AFDW ####
# Multiply by the volume of lake water in L
# Divide by the grams of AFDW 
# Make uptake rates mostly positive to visualize.
# Final units are in ugN/g AFDW-hr
dat <- Incubation_all %>%
  mutate(net_delta_Conc_ugNLhr_OM = (((-1 * net_delta_Conc_ugNLhr * Volume_lake_water_L)/AFDM_g))) # correct OM here

# dat$Inc_month <- as.character(dat$Inc_month)

# visualize raw uptake net_delta_Conc_ugNLhr, with correction for Sample_dryweight_g, and net_delta_Conc_ugNLhr_OM
# Raw uptake (converted to positive)
hist(dat$net_delta_Conc_ugNLhr * -1,
     main = "Raw Uptake", xlab = "µg N/L/hr")

# Normalized by dry mass (converted to positive)
hist((dat$net_delta_Conc_ugNLhr * -1) / dat$Weight_g,
     main = "Per gram", xlab = "µg N/L/hr/g")

# Normalized by OM (already positive)
hist(dat$net_delta_Conc_ugNLhr_OM,
     main = "Per OM", xlab = "µg N/L/hr/g OM")


# plot together
dat_uptake <- dat %>%
  mutate(
    raw_uptake = net_delta_Conc_ugNLhr * -1,
    weight_corrected = net_delta_Conc_ugNLhr / Weight_g * -1,
    afdm_corrected = net_delta_Conc_ugNLhr_OM
  ) %>%
  pivot_longer(
    cols = c(raw_uptake, weight_corrected, afdm_corrected),
    names_to = "Metric", values_to = "Uptake"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c(
      "raw_uptake", 
      "weight_corrected", 
      "afdm_corrected"
    ))
  )

correction_bio_sed <- ggplot(dat_uptake, aes(x = Type, y = Uptake, fill = Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(x = "Sample Type", y = "Uptake (µg N/L/hr)") +
  theme_bw() +
  theme(legend.position = "none")
correction_bio_sed

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


### Now that we have visualized the flagged outliers let's rerun models
# Create a new dataframe excluding the flagged rows
flagged_rows <- dat %>% filter(flag == TRUE)
dat_filtered <- dat %>% filter(flag == FALSE)

# create vertical line for "true" spike conc 
# read the p val and ver

# Preprocess data to ensure consistency in 'Site' names
dat_filtered <- dat_filtered %>%
  mutate(Site = gsub("_", "", Site))

## Save dat_filtered csv
write.csv(dat_filtered, "../data/Filtered_NUptake_2026_04_20.csv")







