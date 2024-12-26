## Pull data for figures for N incubation manuscript
# 2024-10-31

library(snotelr)
library(dplyr)
library(dataRetrieval)
library(ggplot2)
library(patchwork)
library(tidyr)

#### FLOW FIG ####
# Site and parameter setup
siteNo_GB <- "10336730"
siteNo_BW <- "10336660"
pCode_flow <- "00060"
pCode_stage <- "00065"
start.date <- "2023-01-01"
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
  





#### Read in Chemistry Data ####
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
    site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH(NS)",     # Shallow Littoral
    site %in% c("SSNS1", "SSNS2", "SSNS3") ~ "SS(NS)",     # Shallow Littoral
    site %in% c("GBNS1", "GBNS2", "GBNS3") ~ "GB(NS)",     # Gravel Bars
    site %in% c("BWNS1", "BWNS2", "BWNS3") ~ "BW(NS)",     # Benthos
    site == "GB0.5m" ~ "GB(0.5m)",                         # Shallow Gravel Bar
    site == "BW0.5m" ~ "BW(0.5m)",                         # Shallow Benthos
    site == "GBO" ~ "GB(Outlet)",                          # Outlet Gravel Bar
    site == "BWO" ~ "BW(Outlet)",                          # Outlet Benthos
    TRUE ~ NA_character_                                   # Optional: keep original if needed
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

# Ensure 'position' is a factor with the expected levels
data_long_clean <- data_long_clean %>%
  mutate(
    position = factor(position, levels = c(
      "GB(Outlet)", "GB(0.5m)", "GB(NS)", 
      "BW(Outlet)", "BW(0.5m)", "BW(NS)",
      "SS(NS)", "SH(NS)"
    ))
  )



#### CHEM FIGURE ####
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
      "GB(Outlet)" = "#A67D17",    # Gold for GB
      "GB(0.5m)" = "#A67D17",      # Gold for GB
      "GB(NS)" = "#A67D17",        # Gold for GB
      "BW(Outlet)" = "#3283A8",    # Blue for BW
      "BW(0.5m)" = "#3283A8",      # Blue for BW
      "BW(NS)" = "#3283A8",        # Blue for BW
      "SH(NS)" = "#c76640",        # Red for SS
      "SS(NS)" = "#136F63"         # Green for SH
    ),
    na.value = "gray"  # Optional: specify color for NA if any
  ) +
  theme_minimal() +  # Minimal theme for clarity
  theme(
    strip.background = element_rect(fill = "gray80", color = NA),  # Facet strip background
    strip.text = element_text(size = 12, face = "bold"),  # Facet labels styling
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    legend.position = "none",  # Remove the legend
    panel.grid.major = element_line(color = "gray85", size = 0.25),  # Reduce gridlines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_blank()  # Remove plot title
  )

# Save
# ggsave("nut_plot_20241216.png", plot = p, width = 12, height = 8, dpi = 300)








#### Read in temperature ####
BW_temp <- read.csv("Raw_dat/NS_BW_DO_dat_24.csv")
GB_temp <- read.csv("Raw_dat/NS_GB_DO_dat_24.csv")
SS_temp <- read.csv("Raw_dat/NS_SS_DO_dat_24.csv")
SH_temp <- read.csv("Raw_dat/NS_SH_DO_dat_24.csv")


    
    
## Make plots
## Could add on Kd

#### PLOT OF DISCHARGE SWE AND TEMP
# Define color palette for GB and BW
temp_colors <- c("GB" = "#A67D17", "BW" = "#3283A8")  # Yellow for GB, Blue for BW

# Define highlight dates
highlight_dates <- as.Date(c("2023-05-24", "2023-06-26", "2023-07-21"))

# Modify flow_plot for GB and BW categories with updated colors
flow_plot <- ggplot(flow_data, aes(x = date, color = as.factor(Site))) +
  geom_line(aes(y = dischargeCMS), linetype = "solid", linewidth = 1) +
  geom_point(data = flow_data %>% filter(date %in% highlight_dates), 
             aes(y = dischargeCMS), color = "white", shape = 21, fill = "black", size = 2) +
  scale_x_date(expand = c(0, 0), breaks = seq(as.Date("2015-01-01"), as.Date("2023-10-01"), by = "4 months"),
               date_labels = "%b-%y") +
  labs(y = expression('Discharge ('*m^3~s^-1~km^-2*')'), x = 'Date', color = "Site") +
  scale_color_manual(values = temp_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(t = 10, b = 20)
  )

# Modify swe_plot for GB and BW categories
swe_plot <- ggplot(target_snow_data, aes(x = date, y = snow_water_equivalent_cumulative, color = as.factor(Site))) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(data = target_snow_data %>% filter(date %in% highlight_dates), 
             aes(y = snow_water_equivalent_cumulative), color = "white", shape = 21, fill = "black", size = 2) +
  scale_x_date(expand = c(0, 0), breaks = seq(as.Date("2015-01-01"), as.Date("2023-10-01"), by = "4 months"),
               date_labels = "%b-%y") +
  labs(y = "Cumulative SWE (mm)", x = 'Date', color = "Site") +
  scale_color_manual(values = temp_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(t = 20, b = 10) 
  )

# Assign positions and colors for temperature plot
data_long <- data_long %>%
  mutate(
    position = case_when(
      site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH(NS)",
      site %in% c("SSNS1", "SSNS2", "SSNS3") ~ "SS(NS)",
      site %in% c("GBNS1", "GBNS2", "GBNS3") ~ "GB(NS)",
      site %in% c("BWNS1", "BWNS2", "BWNS3") ~ "BW(NS)",
      site == "GB0.5m" ~ "GB(0.5m)",
      site == "BW0.5m" ~ "BW(0.5m)",
      site == "GBO" ~ "GB(Outlet)",
      site == "BWO" ~ "BW(Outlet)",
      TRUE ~ NA_character_  # Ensures unmatched sites are marked as NA
    ),
    color_group = ifelse(grepl("GB", position), "GB", "BW"),  # Assign color group
    month = factor(format(as.Date(date), "%b"), levels = month.abb[1:10])
  )

# Filter data for January to October and remove rows with NA in 'position'
data_long_filtered <- data_long %>% 
  filter(month %in% month.abb[1:9], !is.na(position))

# Reorder 'position' based on desired levels and assign shapes for each
data_long_filtered$position <- factor(data_long_filtered$position, levels = c("SH(NS)", "SS(NS)", "GB(NS)", "BW(NS)",
                                                                              "GB(0.5m)", "BW(0.5m)", "GB(Outlet)", "BW(Outlet)"))

# Plot temperature over month with single color per group and different shapes for positions
temp_month_plot <- ggplot(data_long_filtered, aes(x = month, y = temp_C, color = color_group, shape = position)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(aes(group = color_group, color = color_group), method = "loess", se = FALSE, linetype = "dashed", linewidth = 0.8) +
  labs(y = "Water Temperature (°C)", x = "Month", color = "Group", shape = "Position Group") +
  scale_color_manual(values = temp_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 1)) +  # Custom shapes for each position
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom"
  )

# Combine the plots: stack discharge and SWE on the left, temperature plot on the right
combined_plot <- ((flow_plot / swe_plot) | temp_month_plot) + 
  plot_layout(widths = c(2, 1), heights = c(1, 1)) &
  theme(plot.margin = margin(10, 10, 10, 10))

# Display the combined plot
combined_plot

# Save the combined plot
ggsave("combined_plot.png", combined_plot, width = 12, height = 8, dpi = 300)




