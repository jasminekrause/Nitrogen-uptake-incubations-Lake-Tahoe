#### Uptake Rate Calculation Workflow
#### Script created by: Heili Lowman, KAL, JAK

#1. Data is adjusted by nutrient concentration and time
#If the NO3 or NH4 value is below detection we cut the value in half
#If there is no recorded time we use the averaged time across experiments ~6:46:00
#Filter out values less than 0
#East shore June NH4 background was anomalously high - replaced with mean of May & July
#Spike concentrations of 0 are set to 0 (not mean water conc)
#Missing BW3m July 100ug/L spike conc hardcoded to 117.5

#2. Calculate change in incubation concentrations
#Filter out all incubated samples and group them by triplicates, calc the mean
#Next, we filter the spiked mean values, which represent the baseline prior to incubation
#SH samples assigned to GB water corrections, SS samples assigned to BW water corrections

#3. Calc uptake rates
#We will subtract each mean spike conc from each incubated samples and divide by
#time incubated
#Subtract the uptake that occurred in the water column
#And finally corrected for organic matter (OM) in sediments and biofilm
#OM correction: multiply by volume lake water (L), divide by grams OM
#Units: ugN / g OM·hr, sign flipped so uptake is positive

#4. Tidy data and visualize
#Select key columns for export

### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(here)
library(drc)
library(patchwork)
library(hms)

# And creating function here to turn time-formatted hours into decimal hours.
# https://stackoverflow.com/questions/21781311/how-to-convert-time-to-decimal
decimateTime <- function(time) {
  time_num = as.numeric(unlist(strsplit(time, ":")))
  time_new = time_num[1] + time_num[2]/60 
  return(time_new)
}

setwd("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening/Step_2/")
getwd()

# Load raw datasets.
dat_raw <- read.csv("Step2_NH4corrected_Ninc_data_20260515.csv", header = T)

# Check data structure.
str(dat_raw)

# set theme for plots
custom_theme <- theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none"
  )

dat_raw1 <- dat_raw %>%
  mutate(
    real_icu_time = format(as.POSIXct(Actual_incubation_hrs, format = "%H:%M"), format = "%H:%M:%S")
  )
str(dat_raw1)

# Calculate the average incubation hours
dat_raw$decimal_hours <- sapply(dat_raw$Actual_incubation_hrs, decimateTime)

# Calculate the average of these decimal hours
average_incubation_hours <- mean(dat_raw$decimal_hours, na.rm = TRUE)
# Print the average incubation hours
print(average_incubation_hours)

#### Data quality check ####

# Before looking at any of the triplicate incubations, we must calculate
# the mean starting concentration for each spike (NH4 & NO3) in both
# kinds of lake water to account for any pre-existing concentrations of either.

# # if they are below the limit of detection, set them at half of the limit of detection.
# # Ammonia - 0.002 mg N/L
# # Nitrate - 0.003 mg N/L


dat_raw2 <- dat_raw1 %>%
  mutate(
    Conc_mgNL_ed = pmin(NH4_mgNL, NO2_NO3_mgNL, na.rm = TRUE),  # Populate with non-NA values from NH4_mgNL and NO2_NO3_mgNL
    Conc_mgNL_ed = case_when(
      Conc_mgNL_ed < 0.002 ~ ifelse(Analyte == "NH3", 0.001, Conc_mgNL_ed),  # If NH4_mgNL < 0.002, set to 0.001
      Conc_mgNL_ed < 0.003 & Analyte == "NO3" ~ 0.0015,  # If NO2_NO3_mgNL < 0.003 and Analyte is NO3, set to 0.0015
      TRUE ~ Conc_mgNL_ed  # Otherwise, keep the original value
    ),
    Conc_ugNLK = 1000 * Conc_mgNL_ed,
    real_icu_time = ifelse(is.na(real_icu_time), "06:46:00", 
      real_icu_time)
  ) %>%
  ungroup()

#filter out values less than 0
dat_raw_nozero <- dat_raw2 %>%
  filter(Conc_ugNLK > 0)


Chemcheck <- ggplot(dat_raw2, aes(x = Vial_num, y = NH4_mgNL, color=Location)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

Chemcheck2 <- ggplot(dat_raw2, aes(x = Vial_num, y = NO2_NO3_mgNL, color=Location)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

Chemcheck3 <- ggplot(dat_raw2, aes(x = Vial_num, y = Conc_ugNLK, color=Location, 
                                  shape=Analyte)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)


#### Calculate change in incubation concentrations ####

# Next, we take the difference between final measured concentrations and
# original concentrations to examine the change over time. We also calculate
# this as a rate based on the amount of time everything was incubated.

# In order to do this, we first append the dataset created above to the
# larger dataset, creating essentially a new column with our mean "before"
# concentrations.

# Calculate mean of all  triplicates in dataset
dat_avg_before <- dat_raw_nozero %>%
  group_by(Site, Location, Inc_month, Type, Analyte, Spike_µg_L) %>%
  summarise(
    Mean_Conc = mean(Conc_ugNLK, na.rm = TRUE),
    SD_Conc = sd(Conc_ugNLK, na.rm = TRUE),
    .groups = "drop"
  )

# we need to make a new column to facilitate joining.
dat_raw_nozeroa <- dat_raw_nozero %>%
  left_join(dat_avg_before, by = c("Site", "Location", "Type","Inc_month", "Analyte", "Spike_µg_L"))

##SPIKES
# filter out mean spike concentrations and join  them back to dataset
mean_spike_conc_df <- dat_raw_nozeroa %>%
  filter(Type == "spike") %>%
  group_by(Shore, Inc_month, Analyte, Spike_µg_L) %>%
  summarise(Mean_Spike_Conc = mean(Conc_ugNLK, na.rm = TRUE))


### visualize
# background concentrations only (0 spike additions)
background_df <- mean_spike_conc_df %>%
  filter(Spike_µg_L == 0) %>%
  mutate(
    Analyte = factor(Analyte, levels = c("NH3", "NO3")),
    Inc_month = factor(Inc_month, levels = c("May", "June", "July"))
  )

# plot
ggplot(background_df,
       aes(x = Inc_month,
           y = Mean_Spike_Conc,
           color = Shore,
           group = Shore)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~Analyte, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Incubation month",
    y = expression("Background concentration ("*mu*"g/L)"),
    color = "Shore"
  )

### east shore june is very high. take the average of may and july
# calculate mean spike conc for GB June NH3 Spike=0 using May and July
GB_june_NH3_fix <- mean_spike_conc_df %>%
  filter(Shore == "E",
         Analyte == "NH3",
         Spike_µg_L == 0,
         Inc_month %in% c("May", "July")) %>%
  pull(Mean_Spike_Conc) %>%
  mean(na.rm = TRUE)

GB_june_NH3_fix

# join as normal then overwrite the problematic value
dat_raw_nozerob <- dat_raw_nozeroa %>%
  left_join(mean_spike_conc_df, by = c("Shore", "Inc_month", "Analyte", "Spike_µg_L")) %>%
  mutate(
    Mean_Spike_Conc = case_when(
      Shore == "E" & Inc_month == "June" & 
        Analyte == "NH3" & Spike_µg_L == 0 ~ GB_june_NH3_fix,
      TRUE ~ Mean_Spike_Conc
    )
  )

# verify fix
dat_raw_nozerob %>%
  filter(Shore == "E", Inc_month == "June", 
         Analyte == "NH3", Spike_µg_L == 0) %>%
  dplyr::select(Shore, Inc_month, Analyte, Spike_µg_L, Mean_Spike_Conc) %>%
  distinct()

# The original dataset had the 0 spike concentration as the average water concentration
# which is incorrect. We need to change this column.
dat_raw_nozerob$Mean_Spike_Conc[dat_raw_nozerob$Spike_µg_L == 0] <- 0

# There was some missing spike concentrations for 100 ug/L for BW3m July
dat_raw_nozerob$Mean_Spike_Conc[dat_raw_nozerob$Spike_µg_L == 100 & is.na(dat_raw_nozerob$Mean_Spike_Conc)] <- 117.5





#### Calculate net uptake rates ####

# Finally, we subtract mean rates of change in lakewater alone from the
# triplicate mean rates of change of subtrate incubations to determine the
# true net uptake rates of the substrate (i.e., sediment or biofilm) alone.

# Correct the calculation of change in concentration
dat_raw_nozeroc <- dat_raw_nozerob %>%
  mutate(
    delta_Conc_ugNL = Conc_ugNLK - Mean_Spike_Conc, # Subtract mean spike concentration
    real_icu_hr = hour(as_hms(real_icu_time)) + (minute(as_hms(real_icu_time)) / 60), # Calculate total hours
    real_icu_hrc = ifelse(is.na(real_icu_time) | real_icu_hr > 7 + 41/60, 6 + 46/60, real_icu_hr), # Correct ICU hours
    delta_Conc_ugNLhr = delta_Conc_ugNL / real_icu_hrc # Calculate delta Conc_?gNL per hour
  )

# Calculate mean uptake rates for water using Assigned_Location for grouping
dat_avg_water <- dat_raw_nozeroc %>%
  filter(Type == "water") %>%
  group_by(Inc_month, Location, Analyte, Mean_Spike_Conc) %>%
  summarize(water_delta_Conc_ugNLhr = mean(delta_Conc_ugNLhr, na.rm = TRUE), .groups = 'drop')

# We need to make sure SH and SS are included
dat_raw_nozeroc_adjusted <- dat_raw_nozeroc %>%
  mutate(Assigned_Location = case_when(
    Location == "SH" ~ "GB",  # This assumes you want to use GB water data for SH
    Location == "SS" ~ "BW",  # This assumes you want to use BW water data for SS
    TRUE ~ Location
  ))

# Join with the larger dataset using Assigned_Location
dat_raw_nozerod <- left_join(dat_raw_nozeroc_adjusted, dat_avg_water,
                             by = c("Assigned_Location" = "Location", "Inc_month", "Analyte", "Mean_Spike_Conc"))

# Calculate net uptake rates including uptake in the water column using the original Location
dat_raw_nozero <- dat_raw_nozerod %>%
  mutate(net_delta_Conc_ugNLhr = delta_Conc_ugNLhr - water_delta_Conc_ugNLhr) # Use mean spike conc

### visualize water uptake
dat_raw_nozeroc %>%
  filter(Type == "water", 
         Spike_µg_L == 0,
         !is.na(Analyte)) %>%
  mutate(Inc_month = factor(Inc_month, levels = c("May", "June", "July"))) %>%
  ggplot(aes(x = Location, y = delta_Conc_ugNLhr, fill = Location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_grid(Analyte ~ Inc_month, scales = "free_y") +
  labs(x = "Location",
       y = expression(N~uptake(µg~N~L^-1~h^-1)),
       fill = "Location") +
  scale_fill_manual(values = c(BW = "#3283a8", GB = "#a67d17")) +
  custom_theme +
  theme(legend.position = "none")


#### correct spike values for organic matter (percent) by multiplying by dry weight
# Multiply by the volume of lake water in L
# Divide by the grams of OM 
# Make uptake rates mostly positive to visualize.
# Final units are in ugN/g OM-hr
# read in OM
om <- read_csv("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening/Raw_dat/Step2_clean_ninc_om.csv")

# standardize site names in om to match dat_raw_nozero
om_clean <- om %>%
  mutate(Site = case_when(
    Site %in% c("BW0.5m_a", "BW0.5m") ~ "BW0.5m",
    Site %in% c("BWNS1", "BWNS2", "BWNS3") ~ "BW_3m",
    Site == "BW10m"  ~ "BW_10m",
    Site %in% c("SSNS1", "SSNS2", "SSNS3", "SSNS") ~ "SS_3m",
    Site == "GB0.5m" ~ "GB0.5m",
    Site %in% c("GBNS1", "GBNS2", "GBNS3") ~ "GB_3m",
    Site == "GB10m"  ~ "GB_10m",
    Site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH_3m",
    TRUE ~ Site
  )) %>%
  group_by(Site, Date) %>%
  summarise(pct_OM = mean(pct_OM, na.rm = TRUE), .groups = "drop")

# join and calculate OM-corrected uptake
dat_raw_nozero_om <- dat_raw_nozero %>%
  mutate(Sample_Date = as.Date(Sample_Date)) %>%
  left_join(om_clean, by = c("Site", "Sample_Date" = "Date")) %>%
  mutate(
    OM_g = Sample_dryweight_g * pct_OM,
    net_delta_Conc_ugNLhr_OM = ((-1 * net_delta_Conc_ugNLhr * Volume_lake_water_L) / OM_g)
  )

#
dat_raw_nozero_om %>%
filter(Type %in% c("sediment", "biofilm"),
       !is.na(Analyte),
       !is.na(net_delta_Conc_ugNLhr_OM)) %>%
  mutate(Inc_month = factor(Inc_month, levels = c("May", "June", "July"))) %>%
  pivot_longer(cols = c(net_delta_Conc_ugNLhr, net_delta_Conc_ugNLhr_OM),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = Location, y = value, fill = Location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_grid(metric ~ Inc_month + Analyte, scales = "free_y") +
  scale_fill_manual(values = c(BW = "#3283a8", GB = "#a67d17",
                               SS = "#136F63", SH = "#c76640")) +
  labs(x = "Location", y = NULL) +
  custom_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


#### Export data. ####
names(dat_raw_nozero_om)

# Let's tidy the dataset a bit so that we can better see the steps
dat_tidy <- dat_raw_nozero_om %>%
  filter(Type %in% c("biofilm", "sediment")) %>%
  dplyr::select(Inc_month, Site, Location, Depth, Analyte, Spike_µg_L, Mean_Spike_Conc,
                Volume_lake_water_L, Type, NO2_NO3_mgNL, NH4_mgNL,
                NH4_mgNL,Conc_ugNLK,pH,Weight_g,
                OM_g, real_icu_hrc, water_delta_Conc_ugNLhr, 
                net_delta_Conc_ugNLhr, net_delta_Conc_ugNLhr_OM)

# Visualize
dat_tidy <- dat_tidy %>%
  mutate(Inc_month = factor(Inc_month, levels = c("May", "June", "July")))

check_dat_tidy <- ggplot(dat_tidy, aes(x = real_icu_hrc,
                                       y = net_delta_Conc_ugNLhr, 
                                       color = Location, shape = Type)) +
  geom_point() + 
  theme_bw() + 
  facet_grid(. ~ Inc_month)  # Months will now appear in the correct order

check_dat_tidy

#ggsave("check_dat_tidy_plot.png", check_dat_tidy, width = 10, height = 6, dpi = 300)

# Export for use in Michaelis-Menten calculation script.
# saveRDS(dat_tidy, "Step2_N_Incubation_Uptake_Rates_20250615.rds")
# write_csv(dat_tidy, "N_Incubation_Uptake_Rates_20240924.csv")


# End of script.