
## Data needed --> symptom onset and last day of exposure

## MCMC, at each iteration: 
# Sample a date of infection for each confirmed case --> fixed to last day of exposure if known, otherwise based on travel history
# Compute incubation period for each confirmed case
# Sample value for distributional parameters
# Compute likelihood of observing computed incubation periods given these parameters
# Accept/reject parameters

load('data/contact_data.RData')
names(data.contact.clean); dim(data.contact.clean)

data <- data.contact.clean

##-----------------------------------------------------------------------
## Exposure windows

summary(as.Date(data$date)) # date of enquete

# Dates of last contact
data$contact1_lastcontact <- as.Date(data$date) - data$days_since_last_contact1
data$contact2_lastcontact <- as.Date(data$date) - data$days_since_last_contact2
data$contact3_lastcontact <- as.Date(data$date) - data$days_since_last_contact3
data$contact4_lastcontact <- as.Date(data$date) - data$days_since_last_contact4
summary(data$contact1_lastcontact)

# Exposure window defined as [SO - 21 days ; last contact in 21 days before SO]
data <- data[!is.na(data$contact1_lastcontact) | !is.na(data$contact2_lastcontact) | !is.na(data$contact3_lastcontact) | !is.na(data$contact4_lastcontact), ]
data <- data[!is.na(data$symptom.onset),]
dim(data)

max.int <- 18 # change for sensitivity analysis

# ## Using date of first symptoms
# library(dplyr)
# data <- data %>%
#   group_by(ID) %>%
#   mutate(
#     # exposure_lower = as.Date(symptom.onset - 21),
#     exposure_lower = min(contact1_lastcontact[(contact1_lastcontact < symptom.onset) & (contact1_lastcontact > (symptom.onset - max.int))], 
#                          contact2_lastcontact[(contact2_lastcontact < symptom.onset) & (contact2_lastcontact > (symptom.onset - max.int))], 
#                          contact3_lastcontact[(contact3_lastcontact < symptom.onset) & (contact3_lastcontact > (symptom.onset - max.int))], 
#                          contact4_lastcontact[(contact4_lastcontact < symptom.onset) & (contact4_lastcontact > (symptom.onset - max.int))], 
#                          na.rm = T),
#     exposure_upper = max(contact1_lastcontact[(contact1_lastcontact < symptom.onset) & (contact1_lastcontact > (symptom.onset - max.int))], 
#                          contact2_lastcontact[(contact2_lastcontact < symptom.onset) & (contact2_lastcontact > (symptom.onset - max.int))], 
#                          contact3_lastcontact[(contact3_lastcontact < symptom.onset) & (contact3_lastcontact > (symptom.onset - max.int))], 
#                          contact4_lastcontact[(contact4_lastcontact < symptom.onset) & (contact4_lastcontact > (symptom.onset - max.int))], 
#                          na.rm = T)
#   )
# data <- data[!is.infinite((data$exposure_upper)), ]
# dim(data); head(data[,c("ID","date","symptom.onset","exposure_lower","exposure_upper")])
# sum(data$exposure_upper >= data$symptom.onset, na.rm = T)
# data$exposure_lower <- ifelse(data$exposure_lower == data$exposure_upper, as.Date(data$symptom.onset - max.int), data$exposure_lower)
# data$exposure_lower <- as.Date(data$exposure_lower); summary(data$exposure_lower)
# 
# data$exposureDuration <- as.numeric(data$exposure_upper - data$exposure_lower)
# summary(data$exposureDuration)
# sum(data$exposureDuration == 0); sum(data$exposureDuration != 0)
# 
# data.incubation <- data[,c(1,2,5,9,10,11,13,14,15,17,18,19,100:104,111,112,128:130)]
# save(data.incubation, file = 'data/data_exposure_sens24d.RData')
# 
# # Dates in numeric format
# minDate <- min(data.incubation$symptom.onset, data.incubation$exposure_lower, data.incubation$exposure_upper)
# data.incubation$symptom.onset.num <- as.numeric(data.incubation$symptom.onset - minDate) + 1
# data.incubation$exposure.lower.num <- as.numeric(data.incubation$exposure_lower - minDate) + 1
# data.incubation$exposure.upper.num <- as.numeric(data.incubation$exposure_upper - minDate) + 1
# 
# summary(data.incubation$exposureDuration)
# sum(data.incubation$exposureDuration > 0)
# sum(data.incubation$symptom.onset < data.incubation$exposure.upper.num)
# sum(data.incubation$symptom.onset < data.incubation$exposure.lower.num)
# 
# data.stan <- data.incubation[data.incubation$exposureDuration > 0, ]
# dim(data.stan)

# ##-------------------------------------------------
# ## Using onset of rash
# data <- data[!is.na(data$rash_onset), ]
# dim(data)

# summary(as.numeric(as.Date(data$rash_onset) - as.Date(data$symptom.onset)))

library(dplyr)
data <- data %>%
  group_by(ID) %>%
  mutate(
    # exposure_lower = as.Date(symptom.onset - 21),
    exposure_lower = min(contact1_lastcontact[(contact1_lastcontact < symptom.onset) & (contact1_lastcontact > (symptom.onset - max.int))], 
                         contact2_lastcontact[(contact2_lastcontact < symptom.onset) & (contact2_lastcontact > (symptom.onset - max.int))], 
                         contact3_lastcontact[(contact3_lastcontact < symptom.onset) & (contact3_lastcontact > (symptom.onset - max.int))], 
                         contact4_lastcontact[(contact4_lastcontact < symptom.onset) & (contact4_lastcontact > (symptom.onset - max.int))], 
                         na.rm = T),
    exposure_upper = max(contact1_lastcontact[(contact1_lastcontact < symptom.onset) & (contact1_lastcontact > (symptom.onset - max.int))], 
                         contact2_lastcontact[(contact2_lastcontact < symptom.onset) & (contact2_lastcontact > (symptom.onset - max.int))], 
                         contact3_lastcontact[(contact3_lastcontact < symptom.onset) & (contact3_lastcontact > (symptom.onset - max.int))], 
                         contact4_lastcontact[(contact4_lastcontact < symptom.onset) & (contact4_lastcontact > (symptom.onset - max.int))], 
                         na.rm = T)
  )
data <- data[!is.infinite((data$exposure_upper)), ]
dim(data); head(data[,c("ID","date","symptom.onset","exposure_lower","exposure_upper")])
sum(data$exposure_upper >= data$symptom.onset, na.rm = T)
data$exposure_lower <- ifelse(data$exposure_lower == data$exposure_upper, as.Date(data$symptom.onset - max.int), data$exposure_lower)
data$exposure_lower <- as.Date(data$exposure_lower); summary(data$exposure_lower)

data$exposureDuration <- as.numeric(data$exposure_upper - data$exposure_lower)
summary(data$exposureDuration)
sum(data$exposureDuration == 0); sum(data$exposureDuration != 0)

# data.incubation <- data[,c(1,2,5,9,10,11,13,14,15,17,18,19,100:104,111,112,124:133)]
data.incubation <- data[,c(1,22,23,26,27,114:128)]
# save(data.incubation, file = 'data/data_exposure_rash.RData')

# Dates in numeric format
minDate <- min(data.incubation$symptom.onset, data.incubation$exposure_lower, data.incubation$exposure_upper)
data.incubation$symptom.onset.num <- as.numeric(data.incubation$symptom.onset - minDate) + 1
data.incubation$exposure.lower.num <- as.numeric(data.incubation$exposure_lower - minDate) + 1
data.incubation$exposure.upper.num <- as.numeric(data.incubation$exposure_upper - minDate) + 1

summary(data.incubation$exposureDuration)
sum(data.incubation$exposureDuration > 0)
sum(data.incubation$exposureDuration == 23)
sum(data.incubation$symptom.onset < data.incubation$exposure.upper.num)
sum(data.incubation$symptom.onset < data.incubation$exposure.lower.num)

data.stan <- data.incubation[data.incubation$exposureDuration > 0, ]
dim(data.stan)

##--------------------------------------------------------------------------
## Estimate incubation period using MCMC
## https://github.com/fmiura/MpxInc_2022/blob/main/src/1_StanModel_MPXinc.R

library(patchwork)
library(rstan)
library(loo)
library(tidyverse)#added
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

# Input for Stan model
input_data <- list(N = length(data.stan$exposure.lower.num),
                   tStartExposure = data.stan$exposure.lower.num,
                   tEndExposure = data.stan$exposure.upper.num,
                   tSymptomOnset = data.stan$symptom.onset.num
)

# Stan model for different distributions
distributions <- c("weibull", "gamma", "lognormal")
code <- sprintf("
  data{
    int<lower=1> N;
    vector[N] tStartExposure;
    vector[N] tEndExposure;
    vector[N] tSymptomOnset;
  }
  parameters{
    real<lower=0> par[2];
    vector<lower=0, upper=1>[N] uE;	// Uniform value for sampling between start and end exposure
  }
  transformed parameters{
    vector[N] tE; 	// infection moment
    tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
  }
  model{
    // Contribution to likelihood of incubation period
    target += %s_lpdf(tSymptomOnset -  tE  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tSymptomOnset[i] -  tE[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code = code)
# source('code/stanmodIncub.RData')

# Fit stan models
fit <- mapply(sampling, models, list(input_data), iter = 10000, warmup = 3000, chain = 4)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY = FALSE)

# Summary of model fits
means <- cbind(pos$weibull[,2]*gamma(1+1/pos$weibull[,1]),
               pos$gamma[,1] / pos$gamma[,2],
               exp(pos$lognormal[,1]+pos$lognormal[,2]^2/2))
a_percentile <- c(0.025, 0.5, 0.975)
res <- apply(means, 2, quantile, a_percentile)
ll <- mapply(function(z) loo(extract_log_lik(z))$looic, fit)
waic <- mapply(function(z) waic(extract_log_lik(z))$waic, fit)
rbind(res, looIC=ll, WAIC=waic)

variances <- cbind((pos$weibull[,2]^2)*((gamma(1 + (2/pos$weibull[,1])) - (gamma(1 + (1/pos$weibull[,1])))^2)),
                   pos$gamma[,1] / (pos$gamma[,2]^2),
                   (exp(pos$lognormal[,2]^2) - 1) * (exp(2*pos$lognormal[,1] + pos$lognormal[,2]^2))
)
res2 <- apply(variances, 2, quantile, a_percentile)
sqrt(res2)

cens_w_percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qweibull(p = p, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025, 0.5, 0.975)))
colnames(cens_w_percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
cens_g_percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qgamma(p = p, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025, 0.5, 0.975)))
colnames(cens_g_percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
cens_ln_percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qlnorm(p = p, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025, 0.5, 0.975)))
colnames(cens_ln_percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)

# Plots
#make data frames for visualization
df <- data.frame(
  #Take mean values to draw emprical CDF
  inc_day = ((input_data$tSymptomOnset-input_data$tEndExposure)+(input_data$tSymptomOnset-input_data$tStartExposure))/2
)
x_plot <- seq(0,30,by=0.1)
Gam_plot <- as.data.frame(list(dose= x_plot, 
                               pred= sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.5))),
                               low = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025))),
                               upp = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.975)))
))
Wei_plot <- as.data.frame(list(dose= x_plot, 
                               pred= sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.5))),
                               low = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025))),
                               upp = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.975)))
))
ln_plot <- as.data.frame(list(dose= x_plot, 
                              pred= sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.5))),
                              low = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025))),
                              upp = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.975)))
))
library(ggplot2)
gamma_ggplot <- ggplot(df, aes(x=inc_day)) +
  stat_ecdf(geom = "step")+ 
  xlim(c(0, 30))+
  geom_line(data=Gam_plot, aes(x=x_plot, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Gam_plot, aes(x=x_plot,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 24)+
  labs(x="Incubation period (days)", y = "Proportion")+
  ggtitle("Gamma")
weibul_ggplot <- ggplot(df, aes(x=inc_day)) +
  stat_ecdf(geom = "step")+ 
  xlim(c(0, 30))+
  geom_line(data=Wei_plot, aes(x=x_plot, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Wei_plot, aes(x=x_plot,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 24)+
  labs(x="Incubation period (days)", y = "Proportion")+
  ggtitle("Weibull")
lognorm_ggplot <- ggplot(df, aes(x=inc_day)) +
  stat_ecdf(geom = "step")+ 
  xlim(c(0, 30))+
  geom_line(data=ln_plot, aes(x=x_plot, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=ln_plot, aes(x=x_plot,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 24)+
  labs(x="Incubation period (days)", y = "Proportion")+
  ggtitle("Lognormal")
# (lognorm_ggplot|gamma_ggplot|weibul_ggplot) + plot_annotation(tag_levels = 'A') #16inchi x 6 inchi
library(ggpubr)
fInc <- ggarrange(weibul_ggplot, gamma_ggplot, lognorm_ggplot, ncol = 3)
ggsave('fitIncubation.jpeg', plot = fInc, width = 50, height = 30, units = 'cm', dpi = 300)

##--------------------------------------------------------------
## EpiLPS semi-parametric estimation

library(EpiLPS)
?estimIncub
# https://statsandr.com/blog/epilps-for-estimation-of-incubation-times/
# Model computes a semi-parametric fit to the data and compares it with classic parametric fits (lognormal, weibull, gamma)
# candidate with the lowest BIC is selected

data.stan$tL <- data.stan$symptom.onset.num - data.stan$exposure.upper.num
data.stan$tR <- data.stan$symptom.onset.num - data.stan$exposure.lower.num
# data.stan$tL <- data.stan$symptom.onset.num - data.stan$exposure.upper.num
# data.stan$tR <- data.stan$symptom.onset.num - data.stan$exposure.lower.num

dataIncub <- data.frame(tL = data.stan$tL, tR = data.stan$tR)
head(dataIncub)
summary(dataIncub$tR - dataIncub$tL)

set.seed(222)
out <- estimIncub(dataIncub, verbose = TRUE)
out$stats
library(ggplot2)
library(gridExtra)
grid.arrange(plot(out, typ = "incubwin"), plot(out, type = "pdf"), nrow = 1)
