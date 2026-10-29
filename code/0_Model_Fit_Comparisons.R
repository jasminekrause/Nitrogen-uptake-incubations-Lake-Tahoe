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
dat_filtered <- read.csv("../om_corrected_data/Filtered_NUptake_2026_10_28.csv")

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

# Visualize subset
lapply(l[30:40], function(x) ggplot(x, aes(Conc, Uptake))+
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

#Linear Model
test_linear <- stan("Linear_Model_Nuptake.stan",
                  data=stan_data_l$BW0.5m_June_sediment_NH3$Model_Data,
                  chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_linear)
summarize_draws(as_draws_df(test_linear))

#MM Model
test_MM <- stan("MM_Model_Nuptake.stan",
                  data=stan_data_l$BW0.5m_June_sediment_NH3$Model_Data,
                  chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_MM)
summarize_draws(as_draws_df(test_MM))


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
  
  #Linear Model
  fit_linear <- stan("Linear_Model_Nuptake.stan",
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
  parsout_linear<-extract(fit_linear,pars=c('b0','b1','sigma','ymu'))
  parsout_MM<-extract(fit_MM,pars=c('Vmax','K','sigma','ymu'))

  ## Extract model fit diagnostics and compile
  diag_mean <- summarize_draws(as_draws_df(fit_mean))[1:2,]
  diag_1stOrder <- summarize_draws(as_draws_df(fit_1stOrder))[1:2,]
  diag_linear <- summarize_draws(as_draws_df(fit_linear))[1:3,]
  diag_MM <- summarize_draws(as_draws_df(fit_MM))[1:3,]

  
  diag_mean$model <- "mean"
  diag_1stOrder$model <- "1stOrder"
  diag_linear$model <- "linear"
  diag_MM$model <- "MM"

  diagnostics <- rbind(diag_mean, diag_1stOrder, diag_linear, diag_MM)
  
  parsout_list <- list("Pars_Mean" = parsout_mean,
                       "Pars_1stOrder" = parsout_1stOrder,
                       "Pars_Linear" = parsout_linear,
                       "Pars_MM" = parsout_MM,
                       "Diagnostics" = diagnostics)
  
  return(parsout_list)
  
}

all_modelfits_pars_diag <- lapply(stan_data_l, function(x) all_model_fit(x))
# Save output locally (but in folder under gitignore)
saveRDS(all_modelfits_pars_diag, "../rds files/modelfits_iter5k_chains4_2026_10_28.rds")

#all_modelfits_pars_diag <- readRDS("../rds files/modelfits_iter5k_chains4_2026_10_28.rds")

##############################################
## Extract and combine model diagnostics
##############################################
extract_diag <- function(y){
  df <- y$Diagnostics
  return(df)
}

## Extract
allpars_allmodels_diag <- ldply(lapply(all_modelfits_pars_diag, function(x) extract_diag(x)), data.frame)

## Save
write.csv(allpars_allmodels_diag, "../data/allpars_diag_iter5k_chains4_2026_10_28.csv")


## Explore
a <- allpars_allmodels_diag

# Mean (null)
View(a[which(a$variable == "mu" & a$rhat < 1.1),])
# 1st order
View(a[which(a$variable == "a" & a$rhat < 1.1),])
# Linear
View(a[which(a$variable == "b0" & a$rhat < 1.1),])
View(a[which(a$variable == "b1" & a$rhat < 1.1),])
# M-M
View(a[which(a$variable == "Vmax" & a$rhat < 1.1),])
View(a[which(a$variable == "K" & a$rhat < 1.1),])


ggplot(a[which(a$variable == "a"),], aes(mean))+
  geom_histogram()+
  scale_x_continuous(trans = "log")


##############################################
## Log-likelihood function
##############################################

## First merge lists by Inc_ID
combined <- c(stan_data_l, all_modelfits_pars_diag)
merged_list <- tapply(combined, names(combined), function(x) unlist(x, recursive = FALSE))

listnames <- c("Mean_Model_Data","Model_Data",
               "Pars_Mean","Pars_1stOrder",
               "Pars_Linear","Pars_MM",
               "Diagnostics")

renaming <- function(x) {
  names(x) <- listnames
  return(x)
}

merged_list <- lapply(merged_list, function(y) renaming(y))


## Log likelihood function
log_lik_fn<-function(data,pred,sigma,iter) {
  
  log_lik_out<-matrix("NA",iter,length(data))
  
  for(i in 1:iter){
    log_lik_out[i,]<-dnorm(data, pred[i,], sigma[i], log=T)
  }
  
  return(log_lik_out)
}


## Initial test - figure out how to run comparison
test <- merged_list[1]
View(test[[1]]$Diagnostics)
ll_test <- apply(log_lik_fn(test$BW0.5m_July_biofilm_NH3$Model_Data$y,
           test$BW0.5m_July_biofilm_NH3$Pars_Linear$ymu,
           test$BW0.5m_July_biofilm_NH3$Pars_Linear$sigma, 10000),
           2, as.numeric)

ll_waic <- loo::waic.matrix(ll_test) #lower the value the better
#waic does not appear to work well

loo(ll_test)


########################################################
## Function to run log likelihood function across models
########################################################

ll_models <- function(x, iter){

  ## Run Log likelihood function (data,pred,sigma,iter)
  ll_mean <- apply(log_lik_fn(x$Model_Data$y,
                              matrix(x$Pars_Mean$mu,iter,x$Model_Data$N),
                              matrix(x$Pars_Mean$sigma,iter,x$Model_Data$N), iter),2,as.numeric)
  ll_1storder <- apply(log_lik_fn(x$Model_Data$y,
                            x$Pars_1stOrder$ymu,
                            x$Pars_1stOrder$sigma, iter),2,as.numeric)
  ll_linear <- apply(log_lik_fn(x$Model_Data$y,
                          x$Pars_Linear$ymu,
                          x$Pars_Linear$sigma, iter),2,as.numeric)
  ll_MM <- apply(log_lik_fn(x$Model_Data$y,
                      x$Pars_MM$ymu,
                      x$Pars_MM$sigma, iter),2,as.numeric)
  
  ll_list <- list(ll_mean, ll_1storder, ll_linear, ll_MM)
  names(ll_list) <- c("ll_mean","ll_1storder","ll_linear","ll_MM")
  
  return(ll_list)
  
}

ll_output <- lapply(merged_list, function(z) ll_models(z, 10000))
saveRDS(ll_output,"../rds files/ll_output_2026_10_28.rds")

##################
## Loo comparison
##################

ll_comparison_loo <- function(x){
  
  #run a loo test on each ll output
  loo_mean <- loo(x$ll_mean)
  loo_1storder <- loo(x$ll_1storder)
  loo_linear <- loo(x$ll_linear)
  loo_MM <- loo(x$ll_MM)
  
  comparison <- loo_compare(loo_mean, #model 1
                            loo_1storder, #model 2
                            loo_linear, #model 3
                            loo_MM) #model 4
  return(comparison)
  
}


output_loo_comparisons <- lapply(ll_output, function(y) ll_comparison_loo(y))
output_loo_df <- ldply(output_loo_comparisons, data.frame)

## save both
saveRDS(output_loo_comparisons, "../rds files/loo_comparisons_2026_10_29.rds")
write.csv(output_loo_df, "../data/loo_comparisons_df_2026_10_29.csv")


