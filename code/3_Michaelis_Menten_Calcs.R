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

setwd("R:/Blaszczak_Lab/Ongoing Projects/Tahoe Data/Raw Data/N_inc/Raw/Uptake_Calculations/")
getwd()

# Load data
Incubation_all <- readRDS("N_Incubation_Uptake_Rates_nofilter_20240924.rds")
print(names(Incubation_all))

#### CORRECT FOR AFDW ####
# Multiply by the volume of lake water in L
# Divide by the grams of AFDW 
# Make uptake rates mostly positive to visualize.
# Final units are in ugN/g AFDW-hr
dat <- Incubation_all %>%
  mutate(net_delta_Conc_ugNLhr_OM = (((-1 * net_delta_Conc_�gNLhr * Volume_lake_water_L)/AFDM_g))) # correct OM here

# dat$Inc_month <- as.character(dat$Inc_month)


#### FIT MM MODELS ####
# Preprocess data to ensure consistency in 'Site' names
dat$Site <- gsub("_", "", dat$Site)
# Ensure the plots directory exists
plots_directory <- "MM_plots_20241121_nofilter"
dir.create(plots_directory, showWarnings = FALSE, recursive = TRUE)
print(paste("Directory created: ", plots_directory))

# Define unique combinations from the dataset
unique_combinations <- unique(dat[, c("Site", "Inc_month", "Type", "Analyte")])
unique_combinations <- data.frame(lapply(unique_combinations, as.character), stringsAsFactors = FALSE)
print(head(unique_combinations))


# Function to fit and plot MM models
fit_and_plot_MM_model <- function(site, inc_month, type, analyte, data) {
  
  # Filter data for the specific combination
  specific_data <- data %>%
    filter(Site == site,
           Inc_month == inc_month,
           Type == type,
           Analyte == analyte)
  
  # to debug any missing data points
  print(paste("Number of data points:", nrow(specific_data)))
  
  if (nrow(specific_data) == 0) {
    print(paste("No data for:", site, inc_month, type, analyte))
    return(NULL)
  }
  
  # Fit the MM model
  tryCatch({
    mm_model <- drm(net_delta_Conc_ugNLhr_OM ~ Mean_Spike_Conc, data = specific_data, fct = MM.2())
    
    # Predict values for plotting
    pred_data <- data.frame(Mean_Spike_Conc = seq(0, max(specific_data$Mean_Spike_Conc, na.rm = TRUE), length.out = 100))
    pred_data$net_delta_Conc_ugNLhr_OM <- predict(mm_model, newdata = pred_data)
    
    # Generate the plot
    plot <- ggplot() +
      geom_point(data = specific_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "blue") +
      geom_line(data = pred_data, aes(x = Mean_Spike_Conc, y = net_delta_Conc_ugNLhr_OM), colour = "red") +
      theme_minimal() +
      xlab("Concentration [ug/L]") +
      ylab("Uptake Rate [ugN/gAFDW-hr]") +
      ggtitle(paste(site, inc_month, type, analyte, sep = " - "))
    
    # Return both the model and the plot
    return(list(model = NULL, plot = plot))
  }, error = function(e) {
    print(paste("Error fitting model for:", site, inc_month, type, analyte, "Error message:", e$message))
    return(NULL)
  })
}

# Initialize a list to store model results
model_results <- list()

# Loop through each unique combination to fit models and generate plots
for(i in 1:nrow(unique_combinations)) {
  comb <- unique_combinations[i, ]
  
  # Fit model and generate plot
  result <- fit_and_plot_MM_model(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, dat)
  
  if(!is.null(result)) {
    # Save the model in the results list
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

# Function to extract coefficients, standard errors, and p-values from model summaries
extract_model_info <- function(model_results) {
  model_info <- do.call(rbind, lapply(names(model_results), function(model_name) {
    model_summary <- summary(model_results[[model_name]])
    coef_table <- coef(model_summary)
    
    data.frame(model_name = model_name,
               d_Intercept = coef_table[1, 1],
               e_Intercept = coef_table[2, 1],
               Std_Error_e = coef_table[2, 2],
               p_value_e = coef_table[2, 4],
               Std_Error_d = coef_table[1, 2],
               p_value_d = coef_table[1, 4],
               stringsAsFactors = FALSE)
  }))
  
  # Calculate uptake from the d_Intercept
  model_info$uptake <- as.numeric(model_info$d_Intercept) * -1
  return(model_info)
}

# Extract model information
model_info_df <- extract_model_info(model_results)

# Save the model_info_df dataframe to a CSV file
output_file <- "model_info_df.csv"
# write.csv(model_info_df, file = output_file, row.names = FALSE)









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
    (dat$Site == "GB10m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "SS3m" & dat$Inc_month == "June" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM > 40) |
    (dat$Site == "SS3m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$net_delta_Conc_ugNLhr_OM < 0) |
    (dat$Site == "SH3m" & dat$Inc_month == "May" & dat$Type == "sediment" & dat$Analyte == "NO3" & dat$Spike_µg_L == 800)
)


### Now that we have visualized the flagged outliers let's rerun models
# Create a new dataframe excluding the flagged rows
flagged_rows <- dat %>% filter(flag == TRUE)
dat_filtered <- dat %>% filter(flag == FALSE)

# create vertical line for "true" spike conc 
# read the p val and ver

#### FIT MM MODELS ####
# Preprocess data to ensure consistency in 'Site' names
dat_filtered <- dat_filtered %>%
  mutate(Site = gsub("_", "", Site))

# Ensure the plots directory exists
plots_directory <- "MM_plots_20241121_filtered"
dir.create(plots_directory, showWarnings = FALSE, recursive = TRUE)
print(paste("Directory created: ", plots_directory))

# Define unique combinations from the dataset
unique_combinations <- unique(dat_filtered[, c("Site", "Inc_month", "Type", "Analyte")])
unique_combinations <- data.frame(lapply(unique_combinations, as.character), stringsAsFactors = FALSE)
print(head(unique_combinations))


# Function to fit and plot MM models
fit_and_plot_MM_model_annotate <- function(site, inc_month, type, analyte, data) {
  
  # Filter data for the specific combination
  specific_data <- data %>%
    filter(Site == site,
           Inc_month == inc_month,
           Type == type,
           Analyte == analyte)
  
  # Debugging: Print the number of data points for the current combination
  print(paste("Number of data points:", nrow(specific_data)))
  
  # Return NULL if no data points are available for the combination
  if (nrow(specific_data) == 0) {
    print(paste("No data for:", site, inc_month, type, analyte))
    return(NULL)
  }
  
  # Identify unique spike concentrations for visualization
  unique_spikes <- unique(specific_data$Spike_µg_L)
  
  # Attempt to fit the Michaelis-Menten (MM) model
  tryCatch({
    # Fit the MM model
    mm_model <- drm(net_delta_Conc_ugNLhr_OM ~ Mean_Spike_Conc, 
                    data = specific_data, 
                    fct = MM.2())
    
    # Extract coefficients and calculate uptake
    coef_table <- summary(mm_model)$coefficients
    d_Intercept <- coef_table[1, 1]       # Parameter d (intercept)
    p_value_d <- coef_table[1, 4]         # p-value for parameter d
    uptake <- as.numeric(d_Intercept) * -1  # Calculate uptake (negative of d)
    
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
               label = paste("p-value:", format(p_value_d, digits = 6, scientific = FALSE)), 
               vjust = 2, hjust = 2, size = 4, color = "black", parse = FALSE) +
      annotate("text", x = Inf, y = Inf, 
               label = paste("Uptake:", format(uptake, digits = 6, scientific = FALSE)), 
               vjust = 4, hjust = 2, size = 4, color = "black", parse = FALSE)
    
    # Return both the fitted model and the generated plot
    return(list(model = mm_model, plot = plot))
    
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
  result <- fit_and_plot_MM_model_annotate(comb$Site, comb$Inc_month, comb$Type, comb$Analyte, dat_filtered)
  
  if(!is.null(result)) {
    # Save the model in the results list
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

# Function to extract coefficients, standard errors, and p-values from model summaries
extract_model_info <- function(model_results) {
  model_info <- do.call(rbind, lapply(names(model_results), function(model_name) {
    model_summary <- summary(model_results[[model_name]])
    coef_table <- coef(model_summary)
    
    data.frame(model_name = model_name,
               d_Intercept = coef_table[1, 1],
               e_Intercept = coef_table[2, 1],
               Std_Error_e = coef_table[2, 2],
               p_value_e = coef_table[2, 4],
               Std_Error_d = coef_table[1, 2],
               p_value_d = coef_table[1, 4],
               stringsAsFactors = FALSE)
  }))
  
  # Calculate uptake from the d_Intercept
  model_info$uptake <- as.numeric(model_info$d_Intercept) * -1
  return(model_info)
}

# Extract model information
model_info_df_filtered_20241121 <- extract_model_info(model_results)

# Save the model_info_df dataframe to a CSV file
output_file <- "model_info_df_filtered_20241121.csv"
#write.csv(model_info_df_filtered_20241121, file = output_file, row.names = FALSE)


## End of script