#Ammonium Calculation
#Adopted from Meredith Brehob 2021
#Krause and Loria 2023 edits
#This code will calculate ammonium values from aqueous ammonia concentrations
#We are pulling in raw data read from the YSI during incubations
#The data is uploaded from individual google sheets from each incubation 
#and merged into one dataframe

##### Read in data #####

# Load packages.
library(tidyverse)
library(lubridate)
library(here)
library(drc)
library(patchwork)
library(readxl)
library(dplyr)
library(readxl)

setwd("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening/Step_1/")
getwd()

# This actually reads in with all incubation.
May_inc <- read.csv("N_inc_May_20240220.csv") [, -c(32)]
June_inc <- read.csv("N_inc_June_20240220.csv")
July_inc <- read.csv("N_inc_July_20240220.csv")

colnames(May_inc)
colnames(June_inc)
colnames(July_inc)

merged_data <- rbind(May_inc,June_inc,July_inc)
summary(merged_data)
unique(merged_data$Inc_month) # still some data missing from the .cvs

merged_data$Sample_Date <- as.Date((merged_data$Sample_Date), format = "%Y-%m-%d")


# read in ammonium conversion data for later
# pH and temperature reading from the outlet and taken with YSI on day of incubation
df <- read_excel("Lake_data_ammonium_conversion.xlsx")
# df <- read.csv("Raw/Lake_data_ammonium_conversion.csv")

# merge data
ammonium_df_sub <- merged_data %>%
  #filter(NH3_mgNL >= 0) %>% # Isolate just NH3, let's skip for now
  mutate(Shore = ifelse(grepl("^GB|^SH", Location), "E",
                        ifelse(grepl("^DI_blank", Location), "L",
                               ifelse(grepl("^BW|^SS", Location), "W", 
                                      NA_character_))))

str(ammonium_df_sub)

# merge data sets
ammonium_df <-merge(x=ammonium_df_sub, y=df[,c("Shore","Inc_month","Temp_degC","pH")],
                    by = c("Shore","Inc_month"), all.x=TRUE)
str(ammonium_df)

# arrange and filter empty rows
ammonium_df <- ammonium_df %>%
  arrange(Inc_month, Vial_num) %>%
  slice(71:n())




##### Assign variables #####
ammonium_df$Temp_degC= as.numeric(ammonium_df$Temp_degC)
ammonium_df$pH=as.numeric(ammonium_df$pH)
ammonium_df$ammon_cal <- ifelse(ammonium_df$Analyte == "NH3", ammonium_df$Conc_mgNL, NA)

##### Calculate NH4 #####
#Calculate pKa
ammonium_df$pKa = 0.09018 + 2727.92/(ammonium_df$Temp_degC+273.15)
#Calculate fraction of NH3
ammonium_df$f = 1/(10^(ammonium_df$pKa-ammonium_df$pH)+1)
#Calculate concentration of NH4
ammonium_df$NH4_mgNL <- (1 - ammonium_df$f) * ammonium_df$ammon_cal

hist(ammonium_df$NH3_mgNL)
hist(ammonium_df$NH4_mgNL)

Amoncheck <- ggplot(ammonium_df, aes(x = NH3_mgNL, y = NH4_mgNL, color=Location)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

Amoncheck <- ggplot(ammonium_df, aes(x = NH3_mgNL, y = NH4_mgNL, color=Type)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

Chemcheck <- ggplot(ammonium_df, aes(x = NO3_spike_µg_L, y = NO2_NO3_mgNL, color=Location)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

Chemcheck2 <- ggplot(ammonium_df, aes(x = NH3_spike_µg_L, y = NH4_mgNL, color=Location)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

Chemcheck3 <- ggplot(ammonium_df, aes(x = NH3_spike_µg_L, y = NH3_mgNL, color=Location)) +
  geom_point() + theme_bw() + facet_grid(.~Inc_month)

##### Write data #####
# clean up df
ammonium_save <- ammonium_df %>%
  rename(OM_percent = OM_.) %>%
  select(-ammon_cal, -pKa, -f, -NH3_mgNL)
# write.csv(ammonium_save, file = "NH4corrected_Ninc_data_20240220.csv", row.names = FALSE)

# ggsave("Amoncheck.png", plot = Amoncheck, width = 15, height = 5, units = "in")