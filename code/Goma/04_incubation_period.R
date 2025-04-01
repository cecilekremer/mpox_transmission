# 
# load('data/contact_data_100325.RData')
# dataKamituga <- data.contact.clean
# 
# load('data/contact_data_Goma_100325.RData')
# dataGoma <- data.contact.clean
# 
# names(dataKamituga)[!(names(dataKamituga)%in%names(dataGoma))]
# names(dataGoma)[!(names(dataGoma)%in%names(dataKamituga))]
# 
# dataKamituga <- dataKamituga[,c(1:5,18:23,32:104,114,116,117,118,119,120,121)]
# dim(dataKamituga); names(dataKamituga)
# dataKamituga$ID_code <- NA
# dataKamituga <- dataKamituga[,c(1:86,92,87:91)]
# 
# dataGoma <- dataGoma[,c(1:5,17:22,31:103, 112, 114, 115, 116, 117, 118, 119, 120)]
# dim(dataGoma); names(dataGoma)
# 
# all(names(dataGoma) == names(dataKamituga))
# # which(names(dataGoma) != names(dataKamituga))
# 
# dataGoma$loc <- 'Goma'
# dataKamituga$loc <- 'Kamituga'
# 
# data.all <- rbind(dataGoma, dataKamituga)      
# dim(data.all)
# 
# # Create ID combined with location
# data.all$ID2 <- ifelse(!is.na(data.all$ID_code), paste0(data.all$ID_code, data.all$ID),
#                        paste0('KM', data.all$ID))
# table(nchar(data.all$ID2))
# 
# table(data.all$vagina_lesion)
# data.all$genital.lesion <- ifelse(data.all$n.lesion.genital>0, 1, 0)
# data.all$anal.lesion <- ifelse(data.all$n.lesion.anal>0, 1, 0)
# 
# table(data.all$genital.lesion, data.all$vagina_lesion)
# table(data.all$genital.lesion, data.all$penis_lesion)
# 
# table(data.all$vagina_lesion)/sum(table(data.all$vagina_lesion))
# table(data.all$penis_lesion)/sum(table(data.all$penis_lesion))
# 
# save(data.all, file = "data/contact_data_ALL_100325.RData")

load('data/contact_data_ALL_100325.RData')

##-----------------------------------------------------------------------
## Exposure windows

data <- data.all

summary(as.Date(data$date)) # date of enquete
summary(data$days_since_last_contact1)

# Dates of last contact
data$contact1_lastcontact <- as.Date(data$date) - data$days_since_last_contact1
data$contact2_lastcontact <- as.Date(data$date) - data$days_since_last_contact2
data$contact3_lastcontact <- as.Date(data$date) - data$days_since_last_contact3
data$contact4_lastcontact <- as.Date(data$date) - data$days_since_last_contact4
summary(data$contact1_lastcontact)

# Exposure window defined as [SO - 21 days ; last contact in 21 days before SO]
data <- data[!is.na(data$contact1_lastcontact) | !is.na(data$contact2_lastcontact) | !is.na(data$contact3_lastcontact) | !is.na(data$contact4_lastcontact), ]
# individuals reporting only one contact
# data <- data[!is.na(data$contact1_lastcontact) & (is.na(data$contact2_lastcontact) & is.na(data$contact3_lastcontact) & is.na(data$contact4_lastcontact)), ]
data <- data[!is.na(data$symptom.onset),]
dim(data)

length(unique(data$ID2))
sum(is.na(data$contact1_lastcontact) & is.na(data$contact2_lastcontact) & is.na(data$contact3_lastcontact) & is.na(data$contact4_lastcontact))
sum(is.na(data$symptom.onset))

table(data$transm_sexual, data$contact1_sexual)
table(data$contact1_sexual) # 1 = yes, 2 = no, 3 = don't remember/no response
table(data$transm_sexual)

sum(data$contact1_lastcontact > data$symptom.onset, na.rm = T)

## Remove case if all last known exposures are after symptom onset
data <- data[(data$contact1_lastcontact <= data$symptom.onset) | (data$contact2_lastcontact <= data$symptom.onset) |
       (data$contact3_lastcontact <= data$symptom.onset) | (data$contact3_lastcontact <= data$symptom.onset), ]
data <- data[!is.na(data$ID2), ]
dim(data); length(unique(data$ID2))


table(data$contact1_type) # 1 = single contact, 2 = multiple occasions / ongoing


### Set exposure limits

max.int <- 60 # change for sensitivity analysis

library(dplyr)
data <- data %>%
  group_by(ID2) %>%
mutate(
  # exposure_lower = first date of contact before symptom onset and after the max incubation period
  exposure_lower = min(contact1_lastcontact[(contact1_lastcontact <= symptom.onset) & (contact1_lastcontact >= (symptom.onset - max.int))],
                       contact2_lastcontact[(contact2_lastcontact <= symptom.onset) & (contact2_lastcontact >= (symptom.onset - max.int))],
                       contact3_lastcontact[(contact3_lastcontact <= symptom.onset) & (contact3_lastcontact >= (symptom.onset - max.int))],
                       contact4_lastcontact[(contact4_lastcontact <= symptom.onset) & (contact4_lastcontact >= (symptom.onset - max.int))],
                       na.rm = T),
  # exposure_upper = last date of contact before symptom onset and after the max incubation period
  exposure_upper = max(contact1_lastcontact[(contact1_lastcontact <= symptom.onset) & (contact1_lastcontact >= (symptom.onset - max.int))],
                       contact2_lastcontact[(contact2_lastcontact <= symptom.onset) & (contact2_lastcontact >= (symptom.onset - max.int))],
                       contact3_lastcontact[(contact3_lastcontact <= symptom.onset) & (contact3_lastcontact >= (symptom.onset - max.int))],
                       contact4_lastcontact[(contact4_lastcontact <= symptom.onset) & (contact4_lastcontact >= (symptom.onset - max.int))],
                       na.rm = T),
  # exposure_lower2 = first date of contact before symptom onset
  exposure_lower2 = min(contact1_lastcontact[(contact1_lastcontact <= symptom.onset)],
                        contact2_lastcontact[(contact2_lastcontact <= symptom.onset)],
                        contact2_lastcontact[(contact2_lastcontact <= symptom.onset)],
                        contact2_lastcontact[(contact2_lastcontact <= symptom.onset)],
                        na.rm = T),
  # exposure_upper2 = last date of contact before symptom onset
  exposure_upper2 = max(contact1_lastcontact[(contact1_lastcontact <= symptom.onset)],
                        contact2_lastcontact[(contact2_lastcontact <= symptom.onset)],
                        contact2_lastcontact[(contact2_lastcontact <= symptom.onset)],
                        contact2_lastcontact[(contact2_lastcontact <= symptom.onset)],
                        na.rm = T)
)
# # reporting only one contact
# mutate(
#   exposure_upper = ifelse(length(contact1_lastcontact[(contact1_lastcontact < symptom.onset) & (contact1_lastcontact > (symptom.onset - max.int))])>0,
#                           contact1_lastcontact[(contact1_lastcontact < symptom.onset) & (contact1_lastcontact > (symptom.onset - max.int))],
#                           symptom.onset
#   )
# )
sum(!is.infinite(data$exposure_lower))

## Fix exposure for those reporting single contact and one-time occasion

## Reporting multiple exposure: set lower bound to two months before last exposure (Perez-Saez et al)

summary(as.numeric(data$exposure_upper-data$symptom.onset))

data <- data[!is.infinite((data$exposure_upper)), ]

dim(data); head(data[,c("ID","date","symptom.onset","exposure_lower","exposure_upper")])
# dim(data); head(data[,c("ID","date","symptom.onset","exposure_upper")])
data$exposure_upper <- as.Date(data$exposure_upper)

sum(data$exposure_upper > data$symptom.onset, na.rm = T)
sum(data$exposure_upper < data$symptom.onset, na.rm = T)
sum(data$exposure_upper == data$symptom.onset, na.rm = T)
sum(data$exposure_lower == data$exposure_upper, na.rm = T)
sum(data$exposure_lower < data$exposure_upper, na.rm = T)

data$exposure_lower <- ifelse(data$exposure_lower == data$exposure_upper, as.Date(data$symptom.onset - max.int), data$exposure_lower)
data$exposure_lower <- as.Date(data$exposure_lower); summary(data$exposure_lower)
# data$exposure_lower <- as.Date(data$symptom.onset - max.int)
summary(data$exposure_lower)

data$exposureDuration <- as.numeric(data$exposure_upper - data$exposure_lower)
summary(data$exposureDuration)
sum(data$exposureDuration == 0); sum(data$exposureDuration != 0)

table(data$transm_sexual)

summary(as.numeric(data$symptom.onset - data$exposure_lower))
summary(as.numeric(data$symptom.onset - data$exposure_upper))

# data.incubation <- data[,c(1,2,5,9,10,11,13,14,15,17,18,19,100:104,111,112,124:133)]
# data.incubation <- data[,c(1,22,23,26,27,114:128)]
data.incubation <- data
# data.incubation <- data[data$transm_sexual == 0, ]
dim(data.incubation)
# save(data.incubation, file = 'data/data_exposure_rash.RData')

# Dates in numeric format
minDate <- min(data.incubation$symptom.onset, data.incubation$exposure_lower, data.incubation$exposure_upper)
data.incubation$symptom.onset.num <- as.numeric(data.incubation$symptom.onset - minDate) + 1
data.incubation$exposure.lower.num <- as.numeric(data.incubation$exposure_lower - minDate) + 1
data.incubation$exposure.upper.num <- as.numeric(data.incubation$exposure_upper - minDate) + 1

summary(data.incubation$exposureDuration)
sum(data.incubation$exposureDuration > 0)
# sum(data.incubation$exposureDuration == 23)
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
ggsave('fitIncubationAll.jpeg', plot = fInc, width = 50, height = 30, units = 'cm', dpi = 300)

##--------------------------------------------------------------
## EpiLPS semi-parametric estimation

library(EpiLPS)
# ?estimIncub
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
