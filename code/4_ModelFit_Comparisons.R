## Compare model fits per site based on output from JAK (need that other code)
## JRB

# Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2"),require, character.only=T)

# Import model fits
fit_MM <- read.csv("../data/model_outputs_MM_filtered.csv", header = T) #59
fit_lin <- read.csv("../data/model_outputs_linear_filtered.csv", header = T) #60

# Merge
colnames(fit_MM); colnames(fit_lin)
fits <- left_join(fit_MM, fit_lin, by = "model_name")

##########################################
## Prep Data
##########################################

# Remove messy incubations based on visual assessments on Nov. 19th
inc_to_exclude <- c("BW0.5m_June_sediment_NO3",
                    "BW0.5m_May_sediment_NO3",
                    "BW3m_June_sediment_NO3",
                    "GB0.5m_July_sediment_NO3",
                    "GB0.5m_June_sediment_NO3",
                    "GB0.5m_May_biofilm_NH3",
                    "GB0.5m_May_biofilm_NO3",
                    "GB0.5m_May_sediment_NH3",
                    "GB0.5m_May_sediment_NO3",
                    "SS3m_May_sediment_NO3")
fits <- fits[-which(fits$model_name %in% inc_to_exclude),]

## Split model name
fits <- fits %>%
  separate(model_name, 
           into = c("Site", "Month", "Substrate", "N_Type"), 
           sep = "_", 
           remove = FALSE)

## Rename Site column
levels(as.factor(fits$Site))
fits$SiteID <- revalue(fits$Site, replace = c("BW0.5m" = "BW_INT",
                                              "BW10m" = "BW_SL",
                                              "BW3m" = "BW_NS",
                                              "GB0.5m" = "GB_INT",
                                              "GB10m" = "GB_SL",
                                              "GB3m" = "GB_NS",
                                              "SH3m" = "SH_NS",
                                              "SS3m" = "SS_NS"))
## Split SiteID
fits <- fits %>%
  separate(SiteID, 
           into = c("Stream", "Depth_Loc"), 
           sep = "_", 
           remove = FALSE)


## Summarize total acceptable incubations by different groupings
fits_total_incub <- fits %>%
  group_by(Stream, Depth_Loc, Substrate, N_Type) %>%
  count()
#write.csv(fits_total_incub, "../data/final_incub_by_SiteLocSubstrateNType.csv")

##########################################################
## Filter which sites had MM models that fit
##########################################################
# Vmax is the rate of uptake at saturating levels of S
# Km which represents the half-saturation constant
## What proportion of incubations had MM models that fit Vmax and Km?

# add category for Vmax and Km p values
fits$Vmax_sig_0.05 <- ifelse(fits$p_value_Vmax < 0.05, yes = "yes", no = "no")
fits$Km_sig_0.05 <- ifelse(fits$p_value_Km < 0.05, yes = "yes", no = "no")
## which have both significant






















