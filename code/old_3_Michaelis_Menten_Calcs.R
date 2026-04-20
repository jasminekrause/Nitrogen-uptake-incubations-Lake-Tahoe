#### MM Calc Workflow
#### March 14, 2024
#### Script created by: HL, KL, JK 

# This script is designed to use calculated uptake rates to estimate parameters
# for Michaelis-Menten model of uptake rates, based on this workflow:
# https://rpubs.com/RomanL/6752

#First we standardize the volume of lake water in each vial and correct for OM
#Then fit the mm model and save plots to a new folder
#Extract model outputs in model_info_df
#Compare GPP and ER 

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(here)
library(drc)
library(patchwork)
library(ggplot2)

getwd()
setwd("/Users/jasminekrause/Desktop/Lab Tech Position UNR/Nearshore Greening Tahoe/RProj Nearshore Greening/Step_3/")

# Load data
Incubation_all <- readRDS("N_Incubation_Uptake_Rates_20240924.rds")
print(names(Incubation_all))

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

# ggsave("correction_bio_sed.png",plot=correction_bio_sed)



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

# Ensure the plots directory exists
plots_directory <- "MM_plots_20250225_filtered"
dir.create(plots_directory, showWarnings = FALSE, recursive = TRUE)
print(paste("Directory created: ", plots_directory))

# Define unique combinations from the dataset
unique_combinations <- unique(dat_filtered[, c("Site", "Inc_month", "Type", "Analyte")])
unique_combinations <- data.frame(lapply(unique_combinations, as.character), stringsAsFactors = FALSE)
print(head(unique_combinations))


# Function to fit and plot MM models with confidence intervals and p-values
fit_and_plot_MM_model_annotate <- function(site, inc_month, type, analyte, data) {
  
  # Filter data for the specific combination
  specific_data <- data %>%
    filter(Site == site,
           Inc_month == inc_month,
           Type == type,
           Analyte == analyte)
  
  # Print the number of data points for the model
  print(paste("Number of data points:", nrow(specific_data)))
  
  # Return NULL otherwise
  if (nrow(specific_data) == 0) {
    print(paste("No data for:", site, inc_month, type, analyte))
    return(NULL)
  }
  
  # Identify unique spike concentrations for visualization
  unique_spikes <- unique(specific_data$Spike_ug_L)
  
  # Attempt to fit MM model
  tryCatch({
    # Fit the MM model
    mm_model <- drm(net_delta_Conc_ugNLhr_OM ~ Mean_Spike_Conc, 
                    data = specific_data, 
                    fct = MM.2())
    
    # Extract coefficient summary
    coef_table <- summary(mm_model)$coefficients
    params <- coef(mm_model)
    
    # Extract Km and Vmax values
    Km <- params[1]   
    Vmax <- params[2]
    
    # Extract standard errors
    Km_se <- coef_table[1, 2]
    Vmax_se <- coef_table[2, 2]
    
    # Calculate 95% confidence intervals
    Km_CI <- c(Km - 1.96 * Km_se, Km + 1.96 * Km_se)
    Vmax_CI <- c(Vmax - 1.96 * Vmax_se, Vmax + 1.96 * Vmax_se)
    
    # Extract p-values
    Km_pval <- coef_table[1, 4]
    Vmax_pval <- coef_table[2, 4]
    
    # Generate predicted values for the plot
    pred_data <- data.frame(
      Mean_Spike_Conc = seq(0, max(specific_data$Mean_Spike_Conc, na.rm = TRUE), length.out = 100)
    )
    pred_data$net_delta_Conc_ugNLhr_OM <- predict(mm_model, newdata = pred_data)
    
    # Create the plot with data points, model predictions, and annotations
    plot <- ggplot() +
      geom_point(data = specific_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "blue") +
      geom_line(data = pred_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "red") +
      geom_vline(xintercept = unique_spikes, linetype = "dashed", color = "green", size = 1) +  # Add dashed lines for spike concentrations
      theme_minimal() +
      xlab("Concentration [ug/L]") +
      ylab("Uptake Rate [ugN/gAFDW-hr]") +
      ggtitle(paste(site, inc_month, type, analyte, sep = " - ")) +
      annotate("text", x = Inf, y = Inf, 
               label = paste("Vmax:", round(Vmax, 3), "\nCI:", round(Vmax_CI[1], 3), "-", round(Vmax_CI[2], 3),
                             "\np-value:", format(Vmax_pval, digits = 3, scientific = TRUE)), 
               vjust = 2, hjust = 2, size = 4, color = "black", parse = FALSE) +
      annotate("text", x = Inf, y = Inf, 
               label = paste("Km:", round(Km, 3), "\nCI:", round(Km_CI[1], 3), "-", round(Km_CI[2], 3),
                             "\np-value:", format(Km_pval, digits = 3, scientific = TRUE)), 
               vjust = 4, hjust = 2, size = 4, color = "black", parse = FALSE)
    
    # Return both the fitted model and the generated plot
    return(list(model = mm_model, plot = plot))
    
  }, error = function(e) {
    # Print errors
    print(paste("Error fitting model for:", site, inc_month, type, analyte, "Error message:", e$message))
    
    # Create a fallback plot with only the data points
    plot <- ggplot() +
      geom_point(data = specific_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "blue") +
      theme_minimal() +
      xlab("Concentration [?g/L]") +
      ylab("Uptake Rate [?gN/gAFDW-hr]") +
      ggtitle(paste(site, inc_month, type, analyte, " - Data Only"))
    
    # Return NULL for the model and the fallback plot
    return(list(model = NULL, plot = plot))
  })
}

# Initialize a list to store model results
model_results <- list()

# Loop through each unique combination to fit models and generate plots
for(i in 1:nrow(unique_combinations)) {
  comb <- unique_combinations[i, ]
  
  # Fit model and generate plot
  result <- fit_and_plot_MM_model_annotate(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, dat_filtered)
  
  if(!is.null(result)) {
    # Save the model
    model_key <- paste(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, sep = "_")
    model_results[[model_key]] <- result$model
    
    # Save the plot
    filename <- paste(plots_directory, sprintf("MM_plot_%s_%s_%s_%s.png", comb$Site, comb$Inc_month, comb$Type, comb$Analyte), sep = "/")
    print(paste("Saving plot to: ", filename))
    ggsave(filename, result$plot, width = 10, height = 8)
  } else {
    print(paste("Plotting failed for: ", paste(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, sep = " - ")))
  }
}


# Function to extract coefficients, standard errors, confidence intervals, and p-values
extract_model_info <- function(model_results) {
  
  # Check if model_results is empty or NULL
  if (length(model_results) == 0) {
    warning("No models found in model_results.")
    return(NULL)
  }
  
  # Iterate over all model names in the list
  model_info <- do.call(rbind, lapply(names(model_results), function(model_name) {
    model <- model_results[[model_name]]
    
    # Check if the model is NULL or failed
    if (is.null(model) || class(model) != "drc") {
      return(data.frame(model_name = model_name, 
                        Vmax = NA, Km = NA, 
                        Std_Error_Km = NA, p_value_Km = NA, 
                        Std_Error_Vmax = NA, p_value_Vmax = NA,
                        CI_Low_Km = NA, CI_High_Km = NA,
                        CI_Low_Vmax = NA, CI_High_Vmax = NA,
                        stringsAsFactors = FALSE))
    }
    
    # Extract summary information
    model_summary <- summary(model)
    coef_table <- coef(model_summary)
    
    # Ensure the coefficient table has expected structure
    if (nrow(coef_table) < 2) {
      warning(paste("Unexpected coefficient structure in model:", model_name))
      return(data.frame(model_name = model_name, 
                        Vmax = NA, Km = NA, 
                        Std_Error_Km = NA, p_value_Km = NA, 
                        Std_Error_Vmax = NA, p_value_Vmax = NA,
                        CI_Low_Km = NA, CI_High_Km = NA,
                        CI_Low_Vmax = NA, CI_High_Vmax = NA,
                        stringsAsFactors = FALSE))
    }
    
    # Extract estimates
    Vmax <- coef_table[1, 1]
    Km <- coef_table[2, 1]
    
    # Extract standard errors
    Std_Error_Vmax <- coef_table[1, 2]
    Std_Error_Km <- coef_table[2, 2]
    
    # Extract p-values
    p_value_Vmax <- coef_table[1, 4]
    p_value_Km <- coef_table[2, 4]
    
    # Calculate 95% Confidence Intervals (CI)
    CI_Low_Km <- Km - 1.96 * Std_Error_Km
    CI_High_Km <- Km + 1.96 * Std_Error_Km
    
    CI_Low_Vmax <- Vmax - 1.96 * Std_Error_Vmax
    CI_High_Vmax <- Vmax + 1.96 * Std_Error_Vmax
    
    # Return structured dataframe
    return(data.frame(model_name = model_name, 
                      Vmax = Vmax, Km = Km, 
                      Std_Error_Km = Std_Error_Km, p_value_Km = p_value_Km, 
                      Std_Error_Vmax = Std_Error_Vmax, p_value_Vmax = p_value_Vmax,
                      CI_Low_Km = CI_Low_Km, CI_High_Km = CI_High_Km,
                      CI_Low_Vmax = CI_Low_Vmax, CI_High_Vmax = CI_High_Vmax,
                      stringsAsFactors = FALSE))
  }))
  
  return(model_info)
}

# Extract model information
model_MM_filtered_20250225 <- extract_model_info(model_results)

# Save the model_info_df dataframe to a CSV file
output_file <- "model_MM_filtered_20250225.csv"
#write.csv(model_MM_filtered_20250225, file = output_file, row.names = FALSE)









#### FIT LINEAR MODELS ####
## fit linear models to above to poor MM fits


##########
# Preprocess data to ensure consistency in 'Site' names
dat_filtered <- dat_filtered %>%
  mutate(Site = gsub("_", "", Site))

# Ensure the plots directory exists
plots_directory <- "Linear_plots_20250225_filtered"
dir.create(plots_directory, showWarnings = FALSE, recursive = TRUE)
print(paste("Directory created: ", plots_directory))

# Define unique combinations from the dataset
unique_combinations <- unique(dat_filtered[, c("Site", "Inc_month", "Type", "Analyte")])
unique_combinations <- data.frame(lapply(unique_combinations, as.character), stringsAsFactors = FALSE)
print(head(unique_combinations))

# Function to fit and plot linear models
fit_and_plot_linear_model_annotate <- function(site, inc_month, type, analyte, data) {
  
  # Filter data for the specific combination
  specific_data <- data %>%
    filter(Site == site,
           Inc_month == inc_month,
           Type == type,
           Analyte == analyte)
  
  # Print the number of data points for the model
  print(paste("Number of data points:", nrow(specific_data)))
  
  # Return NULL if no data points are available for the combination
  if (nrow(specific_data) == 0) {
    print(paste("No data for:", site, inc_month, type, analyte))
    return(NULL)
  }
  
  # Attempt to fit the linear model
  tryCatch({
    # Fit the linear model
    linear_model <- lm(net_delta_Conc_ugNLhr_OM ~ Mean_Spike_Conc, data = specific_data)
    
    # Extract coefficients and calculate uptake
    coef_table <- summary(linear_model)$coefficients
    intercept <- coef_table[1, 1]       # Intercept
    slope <- coef_table[2, 1]           # Slope
    p_value <- coef_table[2, 4]         # p-value for slope
    uptake <- as.numeric(slope) * -1    # Calculate uptake (negative of slope)
    
    # Generate predicted values for the plot
    pred_data <- data.frame(
      Mean_Spike_Conc = seq(0, max(specific_data$Mean_Spike_Conc, na.rm = TRUE), length.out = 100)
    )
    pred_data$net_delta_Conc_ugNLhr_OM <- predict(linear_model, newdata = pred_data)
    
    # Create the plot with data points, model predictions, and annotations
    plot <- ggplot() +
      geom_point(data = specific_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "blue") +
      geom_line(data = pred_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "red") +
      theme_minimal() +
      xlab("Concentration [ug/L]") +
      ylab("Uptake Rate [ugN/gAFDW-hr]") +
      ggtitle(paste(site, inc_month, type, analyte, sep = " - ")) +
      annotate("text", x = Inf, y = Inf, 
               label = paste("p-value:", format(p_value, digits = 6, scientific = FALSE)), 
               vjust = 2, hjust = 2, size = 4, color = "black", parse = FALSE) +
      annotate("text", x = Inf, y = Inf, 
               label = paste("Intercept:", format(intercept, digits = 6, scientific = FALSE)), 
               vjust = 4, hjust = 2, size = 4, color = "black", parse = FALSE)
    
    # Return both the fitted model and the generated plot
    return(list(model = linear_model, plot = plot))
    
  }, error = function(e) {
    # Handle errors during model fitting
    print(paste("Error fitting model for:", site, inc_month, type, analyte, "Error message:", e$message))
    
    # Create a fallback plot with only the data points
    plot <- ggplot() +
      geom_point(data = specific_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "blue") +
      theme_minimal() +
      xlab("Concentration [ug/L]") +
      ylab("Uptake Rate [ugN/gAFDW-hr]") +
      ggtitle(paste(site, inc_month, type, analyte, " - Data Only"))
    
    # Return NULL for the model and the fallback plot
    return(list(model = NULL, plot = plot))
  })
}

# Initialize a list to store model results
model_results <- list()

# Loop through each unique combination to fit models and generate plots
for(i in 1:nrow(unique_combinations)) {
  comb <- unique_combinations[i, ]
  
  # Fit model and generate plot
  result <- fit_and_plot_linear_model_annotate(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, dat_filtered)
  
  if(!is.null(result)) {
    # Save the model in the results list
    model_key <- paste(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, sep = "_")
    model_results[[model_key]] <- result$model
    
    # Save the plot
    filename <- paste(plots_directory, sprintf("Linear_plot_%s_%s_%s_%s.png", comb$Site, comb$Inc_month, comb$Type, comb$Analyte), sep = "/")
    print(paste("Saving plot to: ", filename))
    ggsave(filename, result$plot, width = 10, height = 8)
  } else {
    print(paste("Plotting failed for: ", paste(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, sep = " - ")))
  }
}

# Function to extract coefficients, standard errors, and p-values from model summaries
extract_model_info <- function(model_results) {
  model_info <- do.call(rbind, lapply(names(model_results), function(model_name) {
    model_summary <- summary(model_results[[model_name]])
    coef_table <- coef(model_summary)
    
    data.frame(model_name = model_name,
               intercept = coef_table[1, 1],
               slope = coef_table[2, 1],
               Std_Error_slope = coef_table[2, 2],
               p_value_slope = coef_table[2, 4],
               Std_Error_intercept = coef_table[1, 2],
               p_value_intercept = coef_table[1, 4],
               stringsAsFactors = FALSE)
  }))
  
  return(model_info)
}

# Extract model information
model_linear_filtered_20250225 <- extract_model_info(model_results)

# Save the model_info_df dataframe to a CSV file
output_file <- "model_info_df_filtered_linear_20250225.csv"
# write.csv(model_linear_filtered_20250225, file = output_file, row.names = FALSE)





#### analysis and visualization
## change site names to keep consistency 
## 0.5m --> O (outlet)
## 3m --> NS (nearshore)
## 10m --> SL (shallow littoral)

# Use dat_filtered which are outliers excluded from analysis
# Transform depths into deployment locations
dat_filtered <- dat_filtered %>%
  mutate(Location = case_when(
    Depth == 0.5  ~ "INT",
    Depth == 10.0 ~ "SL",
    Depth == 3.0  ~ "NS",
    TRUE ~ NA_character_  # Handles unexpected values
  ))

# Extract site and create new label based on location
dat_filtered <- dat_filtered %>%
  mutate(
    Site_prefix = str_extract(Site, "^[A-Za-z]+"),  # Extracts the first alphabetic part of Site
    New_Site = paste0(Site_prefix, Location)  # Merges with Location
  )

# Define site order
site_order <- c("BWINT", "GBINT", "BWNS", "GBNS", "SSNS", "SHNS", "BWSL", "GBSL")

#### background figures ####
# visualize OM
# Create boxplots of AFDM_g across Inc_month and Analyte
ggplot(dat_filtered, aes(x = New_Site, y = AFDM_g)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 2, outlier.color = "black") +
  facet_wrap(~ Inc_month, scales = "fixed", ncol = 4) +  # Fix scales across facets
  labs(title = "AFDM_g Boxplots Across Inc_month and Analyte",
       x = "Analyte", y = "AFDM (g)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity



## Uptake plots 
# updated 2025-06-24
## standardize the uptake across sediment and biofilms for the spike concentration
dat_filtered <- dat_filtered %>%
  mutate(
    is_control = ifelse(Mean_Spike_Conc == 0, "control", "spiked"),
    standardized_uptake = ifelse(Mean_Spike_Conc > 0, net_delta_Conc_ugNLhr_OM / Mean_Spike_Conc, NA)
  )

# raw uptake
ggplot(dat_filtered %>% filter(Type %in% c("sediment", "biofilm")),
       aes(x = Type, y = net_delta_Conc_ugNLhr_OM, fill = Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  labs(
    title = "Raw Uptake Rates",
    x = "Sample Type",
    y = expression(paste("Uptake Rate (?gN gAFDM"^-1, " hr"^-1, ")"))
  ) +
  scale_fill_manual(values = c("sediment" = "#d95f02", "biofilm" = "#1b9e77")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# standardized uptake
ggplot(dat_filtered %>% filter(Type %in% c("sediment", "biofilm")),
       aes(x = Type, y = standardized_uptake, fill = Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  labs(
    title = "Standardized Uptake Rates",
    x = "Sample Type",
    y = expression(paste("Standardized Uptake (?gN gAFDM"^-1, " hr"^-1, " per ?gN L"^-1, ")"))
  ) +
  scale_fill_manual(values = c("sediment" = "#d95f02", "biofilm" = "#1b9e77")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# plot uptake concentration vs uptake rate
uptake_spike <- ggplot(dat_filtered %>% filter(Type %in% c("sediment", "biofilm")),
                       aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM, color = Type)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Uptake Rate vs. Spike Concentration",
    x = "Mean Spike Concentration (?gN/L)",
    y = expression(paste("Uptake Rate (?gN gAFDM"^-1, " hr"^-1, ")")),
    color = "Sample Type"
  ) +
  theme_bw(base_size = 14)

# ggsave("Manuscript_figures/uptake_spike.png",plot=uptake_spike, width = 10, height = 8, dpi = 300, units = "in")








#### boxplot for sediment
# Filter for sediment only
dat_sediment <- dat_filtered %>% filter(Type == "sediment")

## AFDM
# Determine the global y-axis limits for AFDM_g
y_min_afdm <- min(dat_sediment$AFDM_g, na.rm = TRUE)
y_max_afdm <- max(dat_sediment$AFDM_g, na.rm = TRUE)

# Create the boxplot for AFDM_g
afdm_boxplot <- ggplot(dat_sediment, aes(x = factor(New_Site, levels = site_order), y = AFDM_g)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 2, outlier.color = "black") +  
  facet_wrap(~ Inc_month, scales = "fixed", ncol = 4) +  # Fix scales across facets
  labs(
    x = "Site", 
    y = expression(paste("AFDM (g)")),
    title = "Sediment AFDM"
  ) +
  ylim(y_min_afdm, y_max_afdm) +  # Set fixed y-axis limits
  theme_bw(base_size = 14) +  
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"  # Hide legend since there's only one type
  )

# Print the plot
print(afdm_boxplot)

# ggsave("Manuscript_figures/afdm_boxplot.png",plot=afdm_boxplot, width = 10, height = 8, dpi = 300, units = "in")


# Filter sediment data to only include spiked treatments
dat_sediment_spiked <- dat_sediment %>%
  filter(is_control != "control")  # or is_control == "spiked"

# Determine the global y-axis limits for sediment
y_min_sed <- min(dat_sediment_spiked$net_delta_Conc_ugNLhr_OM, na.rm = TRUE)
y_max_sed <- max(dat_sediment_spiked$net_delta_Conc_ugNLhr_OM, na.rm = TRUE)

# Boxplot for sediment
sediment_boxplot <- ggplot(dat_sediment_spiked, aes(x = factor(New_Site, levels = site_order), y = net_delta_Conc_ugNLhr_OM, fill = Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 2, outlier.color = "black") +  
  facet_wrap(~ Inc_month + Analyte, scales = "fixed", ncol = 4) +
  scale_fill_manual(values = c("sediment" = "#d95f02")) +  
  labs(
    x = "Site", 
    y = expression(paste("Uptake Rate (?gN gAFDM"^-1, " hr"^-1, ")")),
    title = "Sediment Uptake Rates",
    fill = "Sample Type"
  ) +
  ylim(y_min_sed, y_max_sed) +
  theme_bw(base_size = 14) +  
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

# Plot
print(sediment_boxplot)

# ggsave("Manuscript_figures/boxplot_sed.png",plot=sediment_boxplot, width = 10, height = 8, dpi = 300, units = "in")



## biofilm boxplot
## AFDM
# Filter for biofilm only
dat_biofilm <- dat_filtered %>% filter(Type == "biofilm")

# Determine the global y-axis limits for biofilm
y_min_bio_afdm <- min(dat_biofilm$AFDM_g, na.rm = TRUE)
y_max_bio_afdm <- max(dat_biofilm$AFDM_g, na.rm = TRUE)


# Create the boxplot for AFDM_g
afdm_boxplot_bio <- ggplot(dat_biofilm, aes(x = factor(New_Site, levels = site_order), y = AFDM_g)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 2, outlier.color = "black") +  
  facet_wrap(~ Inc_month, scales = "fixed", ncol = 4) +  # Fix scales across facets
  labs(
    x = "Site", 
    y = expression(paste("AFDM (g)")),
    title = "Biofilm AFDM"
  ) +
  ylim(y_min_bio_afdm, y_max_bio_afdm) +  # Set fixed y-axis limits
  theme_bw(base_size = 14) +  
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"  # Hide legend since there's only one type
  )

# Print the plot
print(afdm_boxplot_bio)

# ggsave("Manuscript_figures/afdm_boxplot_bio.png",plot=afdm_boxplot_bio, width = 10, height = 8, dpi = 300, units = "in")



## Uptake
# Filter biofilm data to only include spiked treatments
dat_biofilm_spiked <- dat_biofilm %>%
  filter(is_control != "control")  # or is_control == "spiked"

# Determine the global y-axis limits for biofilm
y_min_bio <- min(dat_biofilm_spiked$net_delta_Conc_ugNLhr_OM, na.rm = TRUE)
y_max_bio <- max(dat_biofilm_spiked$net_delta_Conc_ugNLhr_OM, na.rm = TRUE)

# boxplot for biofilm
biofilm_boxplot <- ggplot(dat_biofilm_spiked, aes(x = factor(New_Site, levels = site_order), y = net_delta_Conc_ugNLhr_OM, fill = Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 2, outlier.color = "black") +  
  facet_wrap(~ Inc_month + Analyte, scales = "fixed", ncol = 4) +
  scale_fill_manual(values = c("biofilm" = "#1b9e77")) +  
  labs(
    x = "Site", 
    y = expression(paste("Uptake Rate (?gN gAFDM"^-1, " hr"^-1, ")")),
    title = "Biofilm Uptake Rates",
    fill = "Sample Type"
  ) +
  ylim(y_min_bio, y_max_bio) +
  theme_bw(base_size = 14) +  
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

# plot
print(biofilm_boxplot)

# ggsave("Manuscript_figures/boxplot_bio.png",plot=biofilm_boxplot, width = 10, height = 8, dpi = 300, units = "in")



#### End of script