# Nitrogen-uptake-incubations-Lake-Tahoe

## Overview
This project calculates ammonium values, estimates nitrogen uptake rates, and applies Michaelis-Menten models to analyze nitrogen uptake dynamics in aquatic systems. The analysis is based on raw data collected from sensors and incubation experiments of water, sediment, and biofilms from Lake Tahoe. 

## Workflows

### 1. **Ammonium Calculation**
- **Purpose**: Calculate ammonium values from aqueous ammonia concentrations.
- **Data Source**: Raw data pulled from individual Google Sheets for each incubation and merged into a single dataframe.

### 2. **Uptake Rate Calculation**
- **Purpose**: Edit raw nitrogen incubation data to calculate uptake rates.
- **Steps**:
  - Adjust data for nutrient concentrations and time.
  - Handle below-detection values and missing time data.
  - Calculate change in concentrations and uptake rates, subtracting water column uptake.

### 3. **Michaelis-Menten (MM) Model Calculation**
- **Purpose**: Use calculated uptake rates to estimate parameters for the Michaelis-Menten model of nutrient uptake.
- **Steps**:
  - Standardize water volume in each vial and correct for organic matter.
  - Fit the Michaelis-Menten model, generate plots, and extract model parameters.

### 4. **Analysis of MM Models**
- Analyze the results from the fitted Michaelis-Menten models.
- Compare Gross Primary Production (GPP) and Ecosystem Respiration (ER).

### 5. **Figures and Plotting**
- Generate and save plots related to the analysis, stored in `BackgroundplotsNinc20241031`.

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/jasminekrause/Nitrogen-uptake-incubations-Lake-Tahoe.git
