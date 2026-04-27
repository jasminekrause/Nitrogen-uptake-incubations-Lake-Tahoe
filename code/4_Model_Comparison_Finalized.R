## Model comparison code for Lake Tahoe N uptake measurements
## Comparison of mean, linear, and Michaelis-Menten models
## for each location and date
# JRB

## Import packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse","shinystan"), require, character.only=T)
library(rstan)
library(bayesplot)
library(loo)
library(posterior)

## Import data
dat_filtered <- read.csv("../data/Filtered_NUptake_2026_04_20.csv")

#####################################################################
## Create lists of site-month-type-analyte data frames
#####################################################################
head(dat_filtered)

# Create a unique ID row
dat_filtered$Inc_ID <- paste(dat_filtered$Site,"_",
                             dat_filtered$Inc_month,"_",
                             dat_filtered$Type,"_",
                             dat_filtered$Analyte, sep = "")

# Subset to columns for model
df <- dat_filtered[, c("Inc_ID",
                       "Mean_Spike_Conc",
                       "net_delta_Conc_ugNLhr_OM")]

# Rename to more simple column names
colnames(df) <- c("Inc_ID","Conc","Uptake")

# Create a list
l <- split(df, df$Inc_ID)

# Visualize
lapply(l[50:60], function(x) ggplot(x, aes(Conc, Uptake))+
         geom_point(size=2)+theme_bw(base_size = 14)+
         labs(title = x$Inc_ID[1]))

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

colnames(l$BW0.5m_July_biofilm_NH3)

stan_data_compile <- function(data){
  mean_model_data <- list(N =length(data$Inc_ID), y = data$Uptake)
  model_data <- list(N =length(data$Inc_ID), y = data$Uptake, x = data$Conc)
  
  list <- list("Mean_Model_Data" = mean_model_data,
               "Model_Data" = model_data)
  return(list)
}

stan_data_l <- lapply(l, function(x) stan_data_compile(x))

#######################################
## Initial tests to fit each model
#######################################

## Initial tests
#Mean Model
test_mean <- stan("Mean_Model_Nuptake.stan",
                  data=stan_data_l$BW0.5m_June_sediment_NH3$Mean_Model_Data,
                  chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_mean)
summarize_draws(as_draws_df(test_mean))

#1stOrder Model
test_1stOrder <- stan("1stOrder_Model_Nuptake.stan",
                      data=stan_data_l$BW0.5m_June_sediment_NH3$Model_Data,
                      chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_1stOrder)
summarize_draws(as_draws_df(test_1stOrder))


#MM Model
test_MM <- stan("MM_Model_Nuptake.stan",
                data=stan_data_l$BW0.5m_June_sediment_NH3$Model_Data,
                chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_MM)
summarize_draws(as_draws_df(test_MM))

#MM Model
test_MM <- stan("MM_log_Model_Nuptake.stan",
                data=stan_data_l$BW0.5m_June_sediment_NH3$Model_Data,
                chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_MM)
summarize_draws(as_draws_df(test_MM))


#Hill Model
test_Hill <- stan("Hill_Model_Nuptake.stan",
                  data=stan_data_l$BW0.5m_June_sediment_NH3$Model_Data,
                  chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_Hill)
summarize_draws(as_draws_df(test_Hill))





#######################################
## Fit each model to each incubation
#######################################

all_model_fit <- function(x){
  
  dat_mean <- x$Mean_Model_Data
  dat <- x$Model_Data
  
  #Mean Model
  fit_mean <- stan("Mean_Model_Nuptake.stan",
                   data=dat_mean,
                   chains=4,iter=5000,
                   control=list(max_treedepth=12))
  
  #1stOrder Model
  fit_1stOrder <- stan("1stOrder_Model_Nuptake.stan",
                       data=dat,
                       chains=4,iter=5000,
                       control=list(max_treedepth=12))

  #MM Model
  fit_MM <- stan("MM_Model_Nuptake.stan",
                 data=dat,
                 chains=4,iter=5000,
                 control=list(max_treedepth=12))
  

  
  
  ## Extract parameters
  parsout_mean<-extract(fit_mean,pars=c('mu','sigma'))
  parsout_1stOrder<-extract(fit_1stOrder,pars=c('a','sigma','ymu'))
  parsout_MM<-extract(fit_MM,pars=c('Vmax','K','sigma','ymu'))


  
  ## Extract model fit diagnostics and compile
  diag_mean <- summarize_draws(as_draws_df(fit_mean))[1:3,]
  diag_linear <- summarize_draws(as_draws_df(fit_linear))[1:3,]
  diag_MM <- summarize_draws(as_draws_df(fit_MM))[1:3,]
  diag_Hill <- summarize_draws(as_draws_df(fit_Hill))[1:4,]
  diag_1stOrder <- summarize_draws(as_draws_df(fit_1stOrder))[1:3,]
  
  diag_mean$model <- "mean"
  diag_linear$model <- "linear"
  diag_MM$model <- "MM"
  diag_Hill$model <- "Hill"
  diag_1stOrder$model <- "1stOrder"
  diagnostics <- rbind(diag_mean, diag_linear, diag_MM, diag_Hill, diag_1stOrder)
  
  parsout_list <- list("Pars_Mean" = parsout_mean,
                       "Pars_Linear" = parsout_linear,
                       "Pars_MM" = parsout_MM,
                       "Pars_Hill" = parsout_Hill,
                       "Pars_1stOrder" = parsout_1stOrder,
                       "Diagnostics" = diagnostics)
  
  return(parsout_list)
  
}

all_modelfits_pars_diag <- lapply(stan_data_l, function(x) all_model_fit(x))
# Save output locally (but in folder under gitignore)
saveRDS(all_modelfits_pars_diag, "../rds files/modelfits_iter5k_chains4_2026_04_27.rds")







