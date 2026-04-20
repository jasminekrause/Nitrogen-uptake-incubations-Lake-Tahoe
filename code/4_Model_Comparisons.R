## Model comparison code for Lake Tahoe N uptake measurements
## Comparison of mean, linear, and Michaelis-Menten models
## for each location and date
# JAK and JRB

## Import packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse"), require, character.only=T)
library(rstan)
library(bayesplot)
library(loo)

## Import data
dat_filtered <- read.csv("../data/Filtered_NUptake_2026_04_20.csv")

#####################################################################
## Create lists of site-month-type-analyte data frames
#####################################################################
head(dat_filtered)
# Define unique combinations from the dataset
unique_combinations <- unique(dat_filtered[, c("Site", "Inc_month", "Type", "Analyte")])
unique_combinations <- data.frame(lapply(unique_combinations, as.character), stringsAsFactors = FALSE)
print(head(unique_combinations))

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
ggplot(l$BW0.5m_June_sediment_NH3, aes(Conc, Uptake))+
  geom_point(size=2)+theme_bw(base_size = 14)

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

colnames(l$BW0.5m_July_biofilm_NH3)

stan_data_compile <- function(data){
  mean_model_data <- list(N =length(data$Inc_ID), y = data$Uptake)
  linear_model_data <- list(N =length(data$Inc_ID), y = data$Uptake, x = data$Conc)
  MM_model_data <- list(N =length(data$Inc_ID), y = data$Uptake, x = data$Conc)
  
  list <- list("Mean_Model" = mean_model_data,
               "Linear_Model" = linear_model_data,
               "MM_Model" = MM_model_data)
  return(list)
}

stan_data_l <- lapply(l, function(x) stan_data_compile(x))


###########################
## Fit each model
###########################

## Initial tests
#Mean Model
test_mean <- stan("Mean_Model_Nuptake.stan",
                data=stan_data_l$BW0.5m_June_sediment_NH3$Mean_Model,
                chains=3,iter=3000, control=list(max_treedepth=12))
launch_shinystan(test_mean)
#Linear Model
test_linear <- stan("Mean_Model_Nuptake.stan",
                  data=stan_data_l$BW0.5m_June_sediment_NH3$Mean_Model,
                  chains=3,iter=3000, control=list(max_treedepth=12))
launch_shinystan(test_mean)
#MM Model
test_mean <- stan("Mean_Model_Nuptake.stan",
                  data=stan_data_l$BW0.5m_June_sediment_NH3$Mean_Model,
                  chains=3,iter=3000, control=list(max_treedepth=12))
launch_shinystan(test_mean)



fit1<-stan(file='GrowthCurves_WAIC.stan',data=modeldata, chains=3,iter=2000, control = list(adapt_delta = 0.9999, max_treedepth =12))
launch_shinystan(fit1)
parsout1<-extract(fit1,pars=c('ymu','sigma'))
fit1ll<-log_lik(Y,parsout1$ymu,parsout1$sigma,iter=2000)
loo::waic.matrix(fit1ll)





#######################
## Log-likelihood
#######################


log_lik<-function(data,pred,sigma,iter)
{
  log_lik_out<-matrix("NA",iter,length(data))

  for(i in 1:iter){
    log_lik_out[i,]<-dnorm(data,pred[i,],sigma[i],log=T)
    
      }
   
  
  return(log_lik_out)
}


gdata<-read.csv("avg_biovolume_data.csv", header=T)
nstrains<-6
ndays<-46
use<-matrix(NA,ndays,nstrains)

for (i in 1:nstrains){
  for (ii in 1:ndays){
    use[ii,i]<-length(which(gdata$sID==i & gdata$day==ii))
    
    
  }
  
  
}

Ydata<-array(-99,dim=c(3,ndays,nstrains))

for (i in 1:nstrains){
  for (ii in 1:ndays){
    if(length(which(gdata$sID==i & gdata$day==ii))>0){
    Ydata[1:use[ii,i],ii,i]<-log(gdata$bv[which(gdata$sID==i & gdata$day==ii)])
    }
    
  }
  
  
}



YdataNA<-Ydata
YdataNA[which(YdataNA==-99,arr.ind = T)]<-NA

####Model 1: All the same growth rates####
smatch<-c(1,1,1,1,1,1)

modeldata<-list("ndays"=ndays, "nstrains"=nstrains, 'nreps'=3, 'use'=use, "Yobs"=Ydata, 'smatch'=smatch, 'nrates'=length(unique(smatch)))

options(mc.cores = parallel::detectCores())
fit1<-stan(file='GrowthCurves_WAIC.stan',data=modeldata, chains=3,iter=2000, control = list(adapt_delta = 0.9999, max_treedepth =12))
launch_shinystan(fit1)
parsout1<-extract(fit1,pars=c('ymu','sigma'))
fit1ll<-log_lik(Y,parsout1$ymu,parsout1$sigma,iter=2000)
loo::waic.matrix(fit1ll)

####Model 2: All Different growth rates####
smatch<-c(1,2,3,4,5,6)

modeldata<-list("ndays"=ndays, "nstrains"=nstrains, 'nreps'=3, 'use'=use, "Yobs"=Ydata, 'smatch'=smatch, 'nrates'=length(unique(smatch)))

options(mc.cores = parallel::detectCores())
fit2<-stan(file='GrowthCurves_WAIC.stan',data=modeldata, chains=3,iter=10000, control = list(adapt_delta = 0.9999, max_treedepth =12))
launch_shinystan(fit2)

parsout2<-extract(fit2,pars=c('ymu','sigma_o','sigma_p'))

fit2ll<-log_lik(YdataNA,parsout2$ymu,sqrt((parsout2$sigma_o^2)+(parsout2$sigma_p^2)),iter=15000)
loo::waic.array(fit2ll)

########## Model 3: Different Growth Rates Between Toxic and Nontoxic ########
smatch<-c(1,1,2,2,2,2)

modeldata<-list("ndays"=ndays, "nstrains"=nstrains, 'nreps'=3, 'use'=use, "Yobs"=Ydata, 'smatch'=smatch, 'nrates'=length(unique(smatch)))

options(mc.cores = parallel::detectCores())
fit3<-stan(file='GrowthCurves_WAIC.stan',data=modeldata, chains=3,iter=10000, control = list(adapt_delta = 0.9999, max_treedepth =12))
launch_shinystan(fit3)
parsout1<-extract(fit3,pars=c('ymu','sigma_o','sigma_p'))
fit3ll<-log_lik(YdataNA,parsout3$ymu,sqrt((parsout3$sigma_o^2)+(parsout3$sigma_p^2)),iter=15000)
loo::waic.matrix(fit3ll)

####### Model 4: Different growth rates between two congeners of atx, one, and no toxin production #########
smatch<-c(1,1,2,2,3,3)

modeldata<-list("ndays"=ndays, "nstrains"=nstrains, 'nreps'=3, 'use'=use, "Yobs"=Ydata, 'smatch'=smatch, 'nrates'=length(unique(smatch)))

options(mc.cores = parallel::detectCores())
fit4<-stan(file='GrowthCurves_WAIC.stan',data=modeldata, chains=3,iter=10000, control = list(adapt_delta = 0.9999, max_treedepth =12))
launch_shinystan(fit4)
parsout4<-extract(fit4,pars=c('ymu','sigma_o','sigma_p'))
fit4ll<-log_lik(YdataNA,parsout4$ymu,sqrt((parsout4$sigma_o^2)+(parsout4$sigma_p^2)),iter=15000)
loo::waic.matrix(fit4ll)

#######Model 5: Different Growth Rates Based on Toxin Concentration #########
smatch<-c(1,1,2,3,2,3)

modeldata<-list("ndays"=ndays, "nstrains"=nstrains, 'nreps'=3, 'use'=use, "Yobs"=Ydata, 'smatch'=smatch, 'nrates'=length(unique(smatch)))

options(mc.cores = parallel::detectCores())
fit5<-stan(file='GrowthCurves_WAIC.stan',data=modeldata, chains=3,iter=10000, control = list(adapt_delta = 0.9999, max_treedepth =12))
launch_shinystan(fit5)
parsout5<-extract(fit5,pars=c('ymu','sigma_o','sigma_p'))
fit5ll<-log_lik(YdataNA,parsout5$ymu,sqrt((parsout5$sigma_o^2)+(parsout5$sigma_p^2)),iter=15000)
loo::waic.matrix(fit5ll)


#####Start here with workspace file###
pars<-extract(fit1,c("y",'a','b'))

for (i in 1:6){
  

s1mean<-apply(pars$y[,,i],2,mean)
s1upper<-apply(pars$y[,,i],2,quantile, prob=.975)
s1lower<-apply(pars$y[,,i],2,quantile, prob=.025)

if(i==1){plot(1:46,s1mean, type='l', lwd=2, ylab="ln (Biovolume)", xlab='Day', ylim=c(10,20), col=i)} 
else{lines(1:46,s1mean, type='l', lwd=2,col=i)}
#lines(1:46,s1upper, lty=2, lwd=2,col=i)
#lines(1:46,s1lower, lty=2, lwd=2,col=i)

}


for (i in 1:6){
  
  
  s1mean<-apply(apply(pars$y[,,i],1,diff),1,mean)
  s1upper<-apply(apply(pars$y[,,i],1,diff),1,quantile, prob=.975)
  s1lower<-apply(apply(pars$y[,,i],1,diff),1,quantile, prob=.025)
  
  if(i==1){plot(1:45,s1mean, type='l', lwd=2, ylab="Growth Rate", xlab='Day', ylim=c(0,1), col=i)} 
  else{lines(1:45,s1mean, type='l', lwd=2,col=i)}
  #lines(1:45,s1upper, lty=2, lwd=2,col=i)
  #lines(1:45,s1lower, lty=2, lwd=2,col=i)
  
}

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(fit5, regex_pars = "a\\[[1-6]\\]",
           prob = 0.8) + plot_title


#################ggplots#####################
samples<-as.data.frame(extract(fit5))
a5<-as.data.frame(fit5,c("a"))
write.csv(a5, "est_gr_5.csv")


gr.2<-read.csv("est_gr_2.csv", header=T)
cgr.2 <- gr.2 %>%
  group_by(strain) %>%
  summarise(mean = mean(gr))


qgr.2<-gr.2%>%
  group_by(strain)%>%
  summarize_at(vars(gr),
               list(Q1=~quantile(.,probs=0.05), Q3=~quantile(., probs = 0.9)))


ggplot(gr.2, mapping=aes(x=gr, fill=strain, alpha =0.9)) +
  scale_fill_manual(values=c("cadetblue", "steelblue4", "green3","gold","orange", "olivedrab4")) +
  geom_density(color="black") + 
  geom_segment(data=cgr.2, color="black", linetype = "dashed", aes(x=mean, xend=mean, y=0, yend=1)) +
  facet_grid(rows = vars(strain)) + theme_classic()+
  labs(title="Posterior Distribution of Maximum Growth Rates", x="Growth rate", y="Density") + 
  theme(legend.position = "none") 


fit5.gr<-read.csv("fit5.csv", header=T)
fit5.gr

gr5 <- fit5.gr %>%
  group_by(strain) %>%
  summarise(mean = mean(gr))



qgr5<-fit5.gr%>%
  group_by(strain)%>%
  summarize_at(vars(gr),
               list(Q1=~quantile(.,probs=0.05), Q3=~quantile(., probs = 0.9)))



ggplot(fit5.gr, mapping=aes(x=gr, fill=strain, alpha =0.9)) +
  scale_fill_manual(values=c("blue4", "gold2","green4")) +
  geom_density(color="black") + 
  geom_segment(data=gr5, color="black", linetype = "dashed", aes(x=mean, xend=mean, y=0, yend=1.75)) +
  facet_grid(rows = vars(strain)) + theme_classic()+
  labs(title="Posterior Distribution of Maximum Growth Rates", x="Growth rate", y="Density") + 
  theme(legend.position = "none") 

ggsave("Posteriordistributions5.png", height=4, width=6)
