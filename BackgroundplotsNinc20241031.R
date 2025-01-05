## Main figures and supp figures for N manuscript
## The following datasets are included: 
###### 1. Metabolism data from Loria et al. This is all 3m sites filtered for incubation periods
###### 2. chla/pheo data 
###### 3. Discharge and weather station data pulled from SNOTEL 
###### 4. Temperature data from 3m sensors
###### 5. Water chemistry data collected from stream and lake

# Last updated 2025-01-05

setwd("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening")

library(snotelr)
library(dplyr)
library(dataRetrieval)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(scales)

#### 1. METAB ####
# Script adopted by KAL
se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

getwd()
SFS_datQ <- readRDS("./SFS24_analysis_dat.rds")
summary(SFS_datQ)
str(SFS_datQ)

SFS_datQ<- SFS_datQ%>%
  drop_na(site)

unique(SFS_datQ$shore)
unique(SFS_datQ$site)

## useful variable creation
SFS_datQ1 <- SFS_datQ %>%
  mutate(
    position = case_when(
      site %in% c("BWNS1","GBNS1", "SSNS1", "SHNS1") ~ "north",
      site %in% c("BWNS2","GBNS2", "SSNS2", "SHNS2" ) ~ "center",
      site %in% c("BWNS3","GBNS3", "SHNS3") ~ "south",
      TRUE ~ as.character(site)
    )
  )

## create column for year and doy 
SFS_datQ1$year <- year(SFS_datQ1$date)
SFS_datQ1$yday <- yday(SFS_datQ1$date)
SFS_datQ1$week <- week(SFS_datQ1$date)
SFS_datQ1$month <- month(SFS_datQ1$date)


# Quality filter to remove high days 
# This removed 38 datapoints
SFS_datQ2 <- SFS_datQ1 %>%
  filter(middle_GPP<30 &  middle_ER>-30)
summary(SFS_datQ2)

# slim the dataset for things you want: 
SFS_datQ3 <-  SFS_datQ2 %>%
  select(site, shore, position, year, month, yday, date, 
         middle_GPP, middle_ER, lake_tempC, lake_DO, Kd_fill, ppt_mm, windsp_mean)

# Average by shore instead of sensor: 
SFS_shore <- SFS_datQ3 %>%
  arrange(shore, yday, date, month, year) %>%
  group_by(shore, yday, date, month, year) %>%
  summarise(
    GPP_m = mean(middle_GPP, na.rm=T),
    ER_m = mean(c(middle_ER * -1), na.rm=T),
    temp_m = mean(lake_tempC, na.rm=T),
    LakeDO_m = mean(lake_DO, na.rm=T),
    Kd_m = mean(Kd_fill, na.rm=T),
    ppt_m = mean(ppt_mm, na.rm=T),
    windsp_m = mean(windsp_mean, na.rm=T)
  ) 

## Jasmine adding on for N incubation plots
# Filter data for April to August 2023
SFS_datQ_filtered <- SFS_shore %>%
  filter(year == 2023 & month %in% c(5, 6, 7))

# Create a factor for the months to ensure they appear in the correct order in the plots
SFS_datQ_filtered$month <- factor(SFS_datQ_filtered$month, levels = c(5, 6, 7), labels = c("May", "June", "July"))

# Boxplot for GPP by site over each month
gpp_site_boxplot <- ggplot(SFS_datQ_filtered, aes(x = shore, y = GPP_m, fill = shore)) +
  geom_boxplot() +
  facet_wrap(~month, scales = "free_x") +
  labs(y = expression(GPP~(mmol~O[2]~m^-3~d^-1)),
       x = "Site") +
  theme_minimal() +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(legend.position = "none")  # remove legend from this plot


# Boxplot for ER by site over each month
er_site_boxplot <- ggplot(SFS_datQ_filtered, aes(x = shore, y = ER_m, fill = shore)) +
  geom_boxplot() +
  facet_wrap(~month, scales = "free_x") +
  labs(y = expression(ER~(mmol~O[2]~m^-3~d^-1)),
       x = "Site") +
  theme_minimal() +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640"))  +
  theme(legend.position = "none")  # remove legend from this plot

# NEP calculation
SFS_datQ_filtered <- SFS_datQ_filtered %>%
  mutate(NEP = GPP_m - ER_m)

# Boxplot for NEP by site over each month
nep_site_boxplot <- ggplot(SFS_datQ_filtered, aes(x = shore, y = NEP, fill = shore)) +
  geom_boxplot() +
  facet_wrap(~month, scales = "free_x") +
  labs(y = expression(NEP~(mmol~O[2]~m^-3~d^-1)),
       x = "Site") +
  theme_minimal() +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) 
#theme(legend.position = "bottom")  # keep legend for this plot

nep_site_boxplot

## plot all together
# Increase the font size in the plots
custom_theme <- theme_minimal() +
  theme(
    text = element_text(size = 16),               # General text size
    axis.title = element_text(size = 18),         # Axis title size
    axis.text = element_text(size = 14),          # Axis text size
    strip.text = element_text(size = 18),         # Facet label size
    legend.title = element_text(size = 16),       # Legend title size
    legend.text = element_text(size = 14)         # Legend text size
  )

# Apply this theme to each plot
gpp_site_boxplot <- gpp_site_boxplot + custom_theme
er_site_boxplot <- er_site_boxplot + custom_theme
nep_site_boxplot <- nep_site_boxplot + custom_theme


#ggsave("GPP.png", plot = gpp_site_boxplot, width = 12, height = 8, dpi = 300)
#ggsave("ER.png", plot = er_site_boxplot, width = 12, height = 8, dpi = 300)
#ggsave("NEP.png", plot = nep_site_boxplot, width = 12, height = 8, dpi = 300)


# Plot kd
kd_time_series_plot <- ggplot(SFS_datQ_filtered, aes(x = date, y = Kd_m, color = shore)) +
  geom_line() +  # Line plot
  geom_point() +  # Add points for better visibility
  labs(y = expression(Kd~(m^-1)),
       x = "Month") +
  theme_minimal() +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme(legend.position = "none")  # Remove legend

kd_time_series_plot <- kd_time_series_plot + custom_theme


#ggsave("Kd.png", plot = kd_time_series_plot, width = 12, height = 8, dpi = 300)


#### 2. Chl-a data ####
NS_W_ChlaPheo <- read.csv("NS_W_ChlaPheo_processed.csv")
NS_RS_ChlaPheo <- read.csv("NS_RS_ChlaPheo_processed.csv")
# Combine the datasets
NS_ChlaPheo_combined$date <- as.Date(NS_ChlaPheo_combined$date)
head(NS_ChlaPheo_combined$date)

# Extract the month from the 'date' column
NS_ChlaPheo_combined$month <- month(NS_ChlaPheo_combined$date)

# Convert the month to a factor for ordered plotting
NS_ChlaPheo_combined$month <- factor(NS_ChlaPheo_combined$month, 
                                     levels = c(5, 6, 7), 
                                     labels = c("May", "June", "July"))

# Filter out sites containing "10m", "15m", or "20m"
NS_ChlaPheo_filtered <- NS_ChlaPheo_combined %>%
  filter(!grepl("15m|20m", site))

# Group the sites into categories
NS_ChlaPheo_filtered <- NS_ChlaPheo_filtered %>%
  mutate(site_group = case_when(
    grepl("GB", site) ~ "GB",
    grepl("BW", site) ~ "BW",
    grepl("SS", site) ~ "SS",
    grepl("SH", site) ~ "SH",
    TRUE ~ site
  ))

# Convert site_group to factor for ordered plotting
NS_ChlaPheo_filtered$site_group <- factor(NS_ChlaPheo_filtered$site_group, levels = c("GB", "BW", "SS", "SH"))

# Filter out rows where the month is NA
NS_ChlaPheo_filtered <- NS_ChlaPheo_filtered %>%
  filter(!is.na(month))

# Create the boxplot with RS on top and W on bottom
chla_boxplot_grouped <- ggplot(NS_ChlaPheo_filtered, aes(x = site_group, y = Chla_ugL, fill = Sample_Type)) +
  geom_boxplot() +
  facet_grid(Sample_Type ~ month, scales = "free_y") +
  labs(y = expression(Chla~(µg~L^-1)),
       x = "Site Group",
       title = "Chlorophyll-a (Chla) by Site Group and Month") +
  theme_minimal() +
  scale_fill_manual(values = c("W" = "#3283a8", "RS" = "#136F63"))

# Display the plot
chla_boxplot_grouped


### PHEO
# Calculate the ratio of Chla_ugL to Pheo_ugL
NS_ChlaPheo_filtered <- NS_ChlaPheo_filtered %>%
  mutate(Chla_Pheo_Ratio = ifelse(Pheo_ugL > 0, Chla_ugL / Pheo_ugL, NA)) %>% 
  filter(Chla_Pheo_Ratio >= 0)

# Create the boxplot with RS on top and W on bottom, for the Chla/Pheo ratio
ratio_boxplot_grouped <- ggplot(NS_ChlaPheo_filtered, aes(x = site_group, y = Chla_Pheo_Ratio, fill = Sample_Type)) +
  geom_boxplot() +
  facet_grid(Sample_Type ~ month, scales = "free_y") +
  labs(y = expression(Chla/Pheo~Ratio),
       x = "Site Group",
       title = "Chla to Pheo Ratio by Site Group and Month") +
  theme_minimal() +
  scale_fill_manual(values = c("W" = "#3283a8", "RS" = "#136F63"))

# Display the plot
ratio_boxplot_grouped



#### 3. Flow and climate ####
# Site and parameter setup
siteNo_GB <- "10336730"
siteNo_BW <- "10336660"
pCode_flow <- "00060"
pCode_stage <- "00065"
start.date <- "2023-04-01"
end.date <- "2023-09-01"

# Retrieve and convert flow data from CFS to CMS
flow_data_GB <- readNWISdata(siteNumbers = siteNo_GB, parameterCd = pCode_flow, startDate = start.date, endDate = end.date) %>%
  mutate(Site = "GB", label = "east") %>%
  dplyr::rename(date = dateTime, dischargeCFS = X_00060_00003) %>%
  mutate(
    dischargeCMS = dischargeCFS * 0.0283168,
    scale_Q = ((dischargeCFS * 0.0283168) / 10.64485)) %>%
  select(date, dischargeCFS, Site, label, dischargeCMS, scale_Q)

flow_data_BW <- readNWISdata(siteNumbers = siteNo_BW, parameterCd = pCode_flow, startDate = start.date, endDate = end.date) %>%
  mutate(Site = "BW", label = "west") %>%
  rename(date = dateTime, dischargeCFS = X_00060_00003) %>%
  mutate(
    dischargeCMS = dischargeCFS * 0.0283168,
    scale_Q = ((dischargeCFS * 0.0283168) / 29.00787)
  ) %>%
  select(date, dischargeCFS, Site, label, dischargeCMS, scale_Q)

# Combine flow data
flow_data <- rbind(flow_data_GB, flow_data_BW)

flow_data <- flow_data %>%
  mutate(date = as.Date(date))


### Retrieve SWE data from SNOTEL
target_snow_data <- snotelr::snotel_download(site_id = c(848, 615), internal = TRUE) %>%
  mutate(date = as.Date(date),
         Site = ifelse(site_id == 615, "GB", "BW")) %>%
  filter(date >= as.Date(start.date) & date <= as.Date(end.date)) %>%
  arrange(Site, date) %>%
  group_by(Site) %>%
  mutate(
    daily_delta_SWE = snow_water_equivalent - lag(snow_water_equivalent, default = 0),
    snow_water_equivalent_cumulative = cumsum(daily_delta_SWE),
    precipitation_quality = ifelse(daily_delta_SWE > 0, "Accumulate", "Melt"),
    rain = ifelse(daily_delta_SWE < 0, 1, 0)
  ) %>%
  ungroup()

# Plot daily delta SWE
delta_swe_plot <- ggplot(target_snow_data, aes(x = date, y = daily_delta_SWE, color = Site)) +
  geom_line() +
  scale_x_date(expand = c(0, 0), breaks = seq(as.Date("2023-01-01"), as.Date("2023-09-01"), "2 months"), date_labels = "%b-%y") +
  labs(y = "Delta SWE (mm)", x = "Date") +
  scale_color_manual(values = alpha(c("#3283A8", "#A67D17"), 0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        plot.subtitle = element_text(size = 12))
  

#### 4. Temperature ####
# Read in raw miniDOT data
BW_temp <- read.csv("Raw_dat/NS_BW_DO_dat_24.csv")
GB_temp <- read.csv("Raw_dat/NS_GB_DO_dat_24.csv")
SS_temp <- read.csv("Raw_dat/NS_SS_DO_dat_24.csv")
SH_temp <- read.csv("Raw_dat/NS_SH_DO_dat_24.csv")

# Rename the column 'location.type' to 'replicate' in SS_temp
colnames(SS_temp)[colnames(SS_temp) == "location.type"] <- "replicate"

# Combine all the datasets into one
df_temp <- rbind(BW_temp, GB_temp, SS_temp, SH_temp)

# Correct time
df_temp$Pacific_Standard_Time <- as.Date(df_temp$Pacific_Standard_Time)

# Filter out flagged data
df_temp <- subset(df_temp, Flag1 != "YES")

# Filter date range
df_temp <- df_temp %>%
  filter(Pacific_Standard_Time >= as.Date(start.date) & Pacific_Standard_Time <= as.Date(end.date))

# Aggregate to daily data
daily_data <- df_temp %>%
  group_by(Pacific_Standard_Time, shore) %>%
  summarise(Temperature_deg_C = mean(Temperature_deg_C, na.rm = TRUE))



#### PLOT OF DISCHARGE SWE AND TEMP
# Define color palette for GB and BW
temp_colors <- c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640") 

# Define highlight dates
highlight_dates <- as.Date(c("2023-05-24", "2023-06-26", "2023-07-21"))

# Custom theme across plots
theme_custom <- theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(t = 20, b = 10)  # Adjust plot margins for better spacing
  )

# Common scale_x_date for all plots 
scale_x_custom <- scale_x_date(
  expand = c(0, 0), 
  breaks = seq(as.Date("2023-04-01"), as.Date("2023-09-01"), by = "1 month"),  # Adjusted for monthly breaks
  date_labels = "%b"  # Show only the month abbreviation (e.g., Apr, May, Jun)
)

# Flow plot
flow_plot <- ggplot(flow_data, aes(x = date, color = as.factor(Site))) +
  geom_line(aes(y = dischargeCMS), linetype = "solid", linewidth = 1) +
  geom_point(data = flow_data %>% filter(date %in% highlight_dates), 
             aes(y = dischargeCMS), color = "white", shape = 21, fill = "black", size = 2) +
  scale_x_custom +
  labs(y = expression('Discharge ('*m^3~s^-1~km^-2*')'), x = 'Month', color = "Site") +
  scale_color_manual(values = temp_colors) +
  theme_custom +
  theme(legend.position = "none", plot.margin = margin(t = 10, b = 20))

# SWE plot
swe_plot <- ggplot(target_snow_data, aes(x = date, y = snow_water_equivalent_cumulative, color = as.factor(Site))) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(data = target_snow_data %>% filter(date %in% highlight_dates), 
             aes(y = snow_water_equivalent_cumulative), color = "white", shape = 21, fill = "black", size = 2) +
  scale_x_custom +
  labs(y = "Cumulative SWE (mm)", x = 'Month', color = "Site") +
  scale_color_manual(values = temp_colors) +
  theme_custom +
  theme(legend.position = "none")

# Temperature month plot
temp_month_plot <- ggplot(daily_data, aes(x = Pacific_Standard_Time, y = Temperature_deg_C, color = shore)) +
  geom_line(linewidth = 1, linetype = "solid") +
  scale_x_date(expand = c(0, 0), 
               breaks = seq(as.Date("2023-04-01"), as.Date("2023-09-01"), by = "1 month"),
               date_labels = "%b") +  # Show only the month abbreviation
  labs(y = "Water Temperature (°C)", x = "Month", color = "Shore") +
  scale_color_manual(values = temp_colors) +
  theme_custom

# Combine the plots 
combined_plot <- ((flow_plot / swe_plot) | temp_month_plot) + 
  plot_layout(widths = c(2, 1), heights = c(1, 2)) &  # Increase height for the temperature plot
  theme(plot.margin = margin(10, 10, 10, 10))

combined_plot

# Save the combined plot
#ggsave("combined_plot.png", combined_plot, width = 12, height = 8, dpi = 300)




#### 5. Read in Chemistry Data ####
# Load the dataset for chemistry data
ns_chem_data <- read.csv("Raw_dat/NS_chem_dat_24.csv")
head(ns_chem_data)

# Filter data for the Water Year 2023
data_long <- ns_chem_data %>%
  filter(WaterYear == 2023)

# Display column names for reference
colnames(data_long)

## Clean Data ####
data_sum <- data_long %>%
  # Convert PO4 from µg/L to mg/L (division by 1000)
  mutate(PO4_mgL_dl = PO4_ugL_dl / 1000) %>%
  
  # Exclude specific stream sites not involved in incubation ("GBL" and "BWL")
  filter(!site %in% c("GBL", "BWL")) %>%
  
  # Add 'position' column based on site depth categories
  mutate(position = case_when(
    site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH(NS)", 
    site %in% c("SSNS1", "SSNS2", "SSNS3") ~ "SS(NS)",    
    site %in% c("GBNS1", "GBNS2", "GBNS3") ~ "GB(NS)",   
    site %in% c("BWNS1", "BWNS2", "BWNS3") ~ "BW(NS)",    
    site == "GB0.5m" ~ "GB(0.5m)",                
    site == "BW0.5m" ~ "BW(0.5m)",             
    site == "GBO" ~ "GB(Outlet)",            
    site == "BWO" ~ "BW(Outlet)",            
    TRUE ~ NA_character_            
  )) %>%
  
  # Convert 'date' to Date format
  mutate(date = as.Date(date, format = "%m/%d/%y"))

# Display the cleaned data structure
head(data_sum)

# Reshape the data to long format for plotting
data_long_clean <- data_sum %>%
  pivot_longer(
    cols = c("PO4_mgL_dl", "NO3_mgL_dl", "NH4_mgL_dl", "DOC_mgL_dl"),
    names_to = "Variable", 
    values_to = "Concentration"
  )

# Remove rows with any NA values in the concentration columns
data_long_clean <- data_long_clean %>%
  filter(
    !is.na(Concentration)
  )

# Add custom y-axis limits based on the variable
data_long_clean <- data_long_clean %>%
  mutate(
    y_limit = case_when(
      Variable == "DOC_mgL_dl" ~ 5,  # Set y-limit for DOC to 5
      TRUE ~ 0.125                  # Set y-limit for NH4, NO3, PO4 to 0.125
    )
  )

# Ensure 'position' is a factor
data_long_clean <- data_long_clean %>%
  mutate(
    position = factor(position, levels = c(
      "GB(Outlet)", "GB(0.5m)", "GB(NS)", 
      "BW(Outlet)", "BW(0.5m)", "BW(NS)",
      "SS(NS)", "SH(NS)"
    ))
  )



#### CHEM FIGURE 
p <- ggplot(data_long_clean %>% 
              filter(!is.na(Concentration) & position != "NA"), aes(x = position, y = Concentration, fill = position)) +
  # Boxplot with outliers
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +  # Show outliers
  facet_wrap(
    ~ Variable, 
    scales = "free_y", 
    labeller = labeller(Variable = c(
      "DOC_mgL_dl" = "DOC",
      "NH4_mgL_dl" = "NH4",
      "NO3_mgL_dl" = "NO3",
      "PO4_mgL_dl" = "PO4"
    ))
  ) +
  labs(
    y = "Concentration (mg/L)", 
    x = "Position"
  ) +
  scale_fill_manual(
    values = c(
      "GB(Outlet)" = "#A67D17",  
      "GB(0.5m)" = "#A67D17",     
      "GB(NS)" = "#A67D17",  
      "BW(Outlet)" = "#3283A8",    
      "BW(0.5m)" = "#3283A8",     
      "BW(NS)" = "#3283A8",      
      "SH(NS)" = "#c76640",      
      "SS(NS)" = "#136F63"
    ),
    na.value = "gray"  
  ) +
  theme_minimal() + 
  theme(
    strip.background = element_rect(fill = "gray80", color = NA), 
    strip.text = element_text(size = 12, face = "bold"),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14), 
    legend.position = "none", 
    panel.grid.major = element_line(color = "gray85", size = 0.25),  
    panel.grid.minor = element_blank(),  
    plot.title = element_blank()
  )

# Save
# ggsave("nut_plot_20241216.png", plot = p, width = 12, height = 8, dpi = 300)




#### END OF SCRIPT ####






