# 
# ## Combine data from Kamituga and Goma
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


## Data containing one row for each case x contact x exposure
load('./Final code/data/contact_data.RData')
data <- data.all
dim(data); length(unique(data$ID2))

sum(is.na(data$ID))
sum(is.na(data$ID2))

## Those reporting symptom onset and at least one contact, incl. last exposure date
data <- data[!is.na(data$symptom.onset), ]
# Dates of last contact
data$contact1_lastcontact <- as.Date(data$date) - data$days_since_last_contact1
data$contact2_lastcontact <- as.Date(data$date) - data$days_since_last_contact2
data$contact3_lastcontact <- as.Date(data$date) - data$days_since_last_contact3
data$contact4_lastcontact <- as.Date(data$date) - data$days_since_last_contact4
summary(data$contact1_lastcontact)

data <- data[!is.na(data$contact1_lastcontact) | !is.na(data$contact2_lastcontact) | !is.na(data$contact3_lastcontact) | !is.na(data$contact4_lastcontact), ]
## Remove case if all last known exposures are after symptom onset
data <- data[(data$contact1_lastcontact <= data$symptom.onset) | (data$contact2_lastcontact <= data$symptom.onset) |
               (data$contact3_lastcontact <= data$symptom.onset) | (data$contact4_lastcontact <= data$symptom.onset), ]
# Remove ???
data <- data[!is.na(data$ID2), ]
dim(data); length(unique(data$ID2))
table(data$ID_code)

## Format data: one line per contact
library(tidyr)
library(dplyr)

df <- data %>%
  mutate(across(matches("^contact\\d+_"), as.character))

df_long <- df %>%
  pivot_longer(cols = matches("^contact\\d+_"),
               names_to = c("contact_num", ".value"),
               names_pattern = "(contact\\d+)_(.*)")
df_long <- df_long[!is.na(df_long$type), ]
df_long <- df_long[!is.na(df_long$lastcontact), ]
length(unique(df_long$ID2))
dim(df_long)

table(df_long$num_contact_mpox)

# ## Raw incubation periods
# raw_incubation_periods <- df_long %>%
#   pivot_longer(cols = c('symptom.onset', 'rash_onset'),
#                names_to = 'to_what',
#                values_to = 'to_date') %>%
#   mutate(to_date = as.Date(to_date)) %>%
#   pivot_longer(cols = c('lastcontact'),
#                names_to = "from_what",
#                values_to = 'from_date') %>%
#   mutate(time_diff = difftime(to_date, from_date, units = 'days') %>% as.numeric()) 
# raw_incubation_periods <- raw_incubation_periods[raw_incubation_periods$time_diff >= 0, ]
# library(ggplot2)
# raw_incubation_periods %>%
#   ggplot(aes(x = time_diff, fill = sexual)) +
#   geom_histogram() +
#   facet_grid(~to_what) +
#   # ggthemes::scale_fill_few("Dark") +
#   theme_bw() +
#   labs(x = 'time difference (days)', title = 'Distribution of times from last contact to symptom/rash onset')
# summary(raw_incubation_periods$time_diff)
# 
# data %>%
#   mutate(time_diff = difftime(rash_onset, symptom.onset, units = 'days') %>%
#            as.numeric()) %>%
#   ggplot(aes(x = time_diff, fill = as.factor(agecat))) +
#   geom_histogram() +
#   facet_wrap(~transm_sexual) +
#   # ggthemes::scale_fill_few('Dark') +
#   theme_bw() +
#   labs(x = 'time (days) between symptom and rash onset')

###------------------------------------------------------------
### Data for Stan model

df_long <- df_long[df_long$symptom.onset >= df_long$lastcontact, ]
dim(df_long); length(unique(df_long$ID2))
sum(df_long$rash_onset >= df_long$lastcontact, na.rm = T)

# df.unique <- df_long %>%
#   distinct(ID2, .keep_all = TRUE)
# summary(as.numeric(df.unique$rash_onset - df.unique$symptom.onset, na.rm = T))

# 114 individuals that report only one contact
df_long <- df_long[df_long$num_contact_mpox == 1, ]
dim(df_long); length(unique(df_long$ID2))
table(df_long$contact_num)

table(df_long$ID_code)
sum(is.na(df_long$ID_code))

table(df_long$type) # 1 = single exposure, 2 = multiple/ongoing exposure

# # Individuals with only single exposure
# table(df_long$type)
# df_long <- df_long[df_long$type == 1, ]
# dim(df_long); length(unique(df_long$ID2))

# df_long <- df_long[df_long$type == 1, ]
table(df_long$sexual) # 1 = yes, 2 = no

table(df_long$included) # 1 = contact included in study

# Scenario 1: set to sexual if missing and hypothesis is sexual (if not family)
# Scenario 2: set to non-sexual if missing
# Scenario 3: set to sexual if anal/genital lesions present (if not family)
df_long$sexual1 <- df_long$sexual
df_long$sexual2 <- df_long$sexual
df_long$sexual3 <- df_long$sexual
for(i in 1:dim(df_long)[1]){
  if(is.na(df_long$sexual[i])){
    if(df_long$rel[i] %in% c(1,2)){
      df_long$sexual1[i] <- 2 # not sexual if family, assuming spouses have indicated sexual contact anyways
      df_long$sexual2[i] <- 2
      df_long$sexual3[i] <- 2
    }else{
      df_long$sexual1[i] <- 1 # sexual if missing and not family
      df_long$sexual2[i] <- 2 # non-sexual if missing
      if(!is.na(df_long$anal.lesion[i]) | !is.na(df_long$genital.lesion[i])){
        if(df_long$anal.lesion[i] == 1 | df_long$genital.lesion[i] == 1){
          df_long$sexual3[i] <- 1 # sexual
        }else{
          df_long$sexual3[i] <- 2 # non-sexual
        }
      }else{
        df_long$sexual3[i] <- 2 # non-sexual
      }
    }
  }
}

## Set exposure windows
mindate <- min(as.Date(df_long$lastcontact)) - 60
# exposure
df_long$a_minus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$symptom.onset) - mindate - 21), # multiple exposure (type = 2)
                          as.numeric(as.Date(df_long$lastcontact) - mindate))
df_long$a_plus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$lastcontact) - mindate + 1),
                         as.numeric(as.Date(df_long$lastcontact) - mindate + 1))
# symptom onset
df_long$b_minus <- as.numeric(as.Date(df_long$symptom.onset) - mindate)
df_long$b_plus <- as.numeric(as.Date(df_long$symptom.onset) - mindate + 1)

summary(df_long$a_plus - df_long$a_minus)
sum(df_long$a_plus-df_long$a_minus < 0)
summary(df_long$b_plus - df_long$a_plus)

summary(df_long$a_plus); summary(df_long$a_minus)

df_longg <- df_long[df_long$a_plus-df_long$a_minus >= 0, ]

# df_longg <- df_longg[df_longg$type == 2, ]
# df_longg <- df_longg[is.na(df_longg$ID_code), ] # Kamituga = ID code NA

dim(df_longg)


##-------------------------------------------------------------------------------
## Weibull posterior estimates

data_stan <- list(
  N = nrow(df_longg),
  a_minus = df_longg$a_minus,
  a_plus = df_longg$a_plus,
  b_minus = df_longg$b_minus,
  b_plus = df_longg$b_plus,
  # incubation = 0,
  upper_bound = as.numeric(max(data.all$date) - mindate), # to account for truncation
  incl_cov = 2, # 1 = overall (no covariate), 2 = include two-level covariate
  # sexual = ifelse(df_longg$sexual1 == 2, 0, 1) # 1 = sexual, 0 = non-sexual
  # sexual = ifelse(df_longg$agecat == 3, 0, 1) # 0 = child, 1 = adult
  # sexual = ifelse(df_longg$agenum < 15, 0, 1)
  sexual = ifelse(df_longg$gender == 1, 0, 1) # 0 = male, 1 = female
)

# mod_estIncub <- rstan::stan_model('./Final code/stan_incub_ward.stan', model_name = 'mod_estIncub')
# save(mod_estIncub, file = 'mod_est_incub.RData')

load('./mod_est_incub.RData')

nrchains = 4
options(mc.cores=parallel::detectCores())
fit_Weibull <- rstan::sampling(mod_estIncub, data = (data_stan), 
                               verbose = TRUE,
                               chains = 4, iter = 10000)

## Results
rstan::traceplot(fit_Weibull, pars = c("mean_","alpha"), inc_warmup = TRUE)
quantile(as.matrix(fit_Weibull)[,'mean_[2]'], probs = c(0.025,0.5,0.975)) # mean sexual
quantile(as.matrix(fit_Weibull)[,'sd_[2]'], probs = c(0.025,0.5,0.975)) # mean sexual
quantile(as.matrix(fit_Weibull)[,'mean_[1]'], probs = c(0.025,0.5,0.975)) # mean non-sexual
quantile(as.matrix(fit_Weibull)[,'sd_[1]'], probs = c(0.025,0.5,0.975)) # mean non-sexual



##-------------------------------------------------------------------------------
## Bayesian regression model

df_longg$agecat2 <- ifelse(df_longg$agenum < 15, 0, 1) # 0 = child, 1 = adult

data_stan <- list(
  N = nrow(df_longg),
  a_minus = df_longg$a_minus,
  a_plus = df_longg$a_plus,
  b_minus = df_longg$b_minus,
  b_plus = df_longg$b_plus,
  # incubation = 0,
  upper_bound = as.numeric(max(data.all$date) - mindate), # to account for truncation
  sexual = ifelse(df_longg$sexual1 == 2, 0, 1),
  # age = df_longg$agenum,
  # age = ifelse(df_longg$agecat == 3, 0, 1), # 0 = child, 1 = adult
  age = df_longg$agecat2, # child cutoff at 15y
  # kamituga = ifelse(is.na(df_longg$ID_code), 1, 0),
  kamituga = rep(0, nrow(df_longg)),
  n_steps_prior = df_longg$a_plus - df_longg$a_minus
)
table(data_stan$sexual)
table(data_stan$age)
table(data_stan$kamituga)

# mod_Weibull_reg <- rstan::stan_model('./Final code/stan_incub_reg.stan', model_name = 'mod_Weibull_reg')
# save(mod_Weibull_reg, file = 'mod_Weibull_reg.RData')

# table(df_longg$type) # 1 = single exposure

load('./mod_Weibull_reg.RData')

nrchains = 4
options(mc.cores=parallel::detectCores())
fit_Weibull <- rstan::sampling(mod_Weibull_reg, data = (data_stan), 
                               verbose = TRUE,
                               chains = 4, iter = 10000)
# rstan::pairs.stanfit(fit_Weibull)

## Results
rstan::traceplot(fit_Weibull, pars = c("theta[1]","theta[2]","theta[3]","alpha"), inc_warmup = TRUE)
quantile(as.matrix(fit_Weibull)[,'theta[1]'], probs = c(0.025,0.5,0.975)) # intercept
quantile(as.matrix(fit_Weibull)[,'theta[2]'], probs = c(0.025,0.5,0.975)) # sexual
quantile(as.matrix(fit_Weibull)[,'theta[3]'], probs = c(0.025,0.5,0.975)) # age
# quantile(as.matrix(fit_Weibull)[,'theta[4]'], probs = c(0.025,0.5,0.975)) # kamituga
# quantile(as.matrix(fit_Weibull)[,'theta[4]'], probs = c(0.025,0.5,0.975)) # age x sexual
quantile(as.matrix(fit_Weibull)[,'alpha'], probs = c(0.025,0.5,0.975)) # Weibull shape alpha

### Individual incubation periods?
posterior_samples <- rstan::extract(fit_Weibull)
a_samples <- posterior_samples$a
b_samples <- posterior_samples$b
incubation_samples <- b_samples - a_samples # matrix with 20000 draws (rows) for each individual (columns)
summary(apply(incubation_samples, 2, mean)) # posterior mean incubation for each individual
hist(apply(incubation_samples, 2, mean)) # posterior mean incubation for each individual
alpha_samples <- posterior_samples$alpha
beta_samples <- posterior_samples$beta
# pooled incubation period across posterior draws and individuals
incubation_all <- as.vector(incubation_samples)
# estimate posterior mean Weibull parameters
alpha_mean <- mean(as.vector(alpha_samples))
beta_mean <- mean(as.vector(beta_samples))
# plot
library(ggplot2)
df <- data.frame(incubation = incubation_all)
ggplot(df, aes(x = incubation)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "white", alpha = 0.7) +
  stat_function(fun = dweibull, args = list(shape = alpha_mean, scale = beta_mean),
                color = "red", size = 1.2, linetype = 'dashed') +
  labs(title = "Posterior incubation periods with fitted Weibull",
       x = "Incubation period (days)", y = "Density") +
  theme_minimal()
# # posterior predictive overlay
# ggplot(df, aes(x = incubation)) +
#   geom_density(fill = "skyblue", alpha = 0.5) +
#   lapply(samples, function(i) {
#     stat_function(fun = dweibull,
#                   args = list(shape = alpha_samples[i], scale = beta_samples[i]),
#                   color = "red", alpha = 0.05)
#   }) +
#   labs(title = "Incubation periods with posterior Weibull overlays",
#        x = "Days", y = "Density") +
#   theme_minimal()

# # hist(as.matrix(fit_Weibull)[,'limit_val_[1]']); summary(as.matrix(fit_Weibull)[,'limit_val_[1]'])
# library(posterior)
# ?rhat
# rhat(fit_Weibull)

# ###---------------------------------------------------------------------------------------
# ### EpiLPS semi-parametric estimation
# 
# library(EpiLPS)
# # ?estimIncub
# # https://statsandr.com/blog/epilps-for-estimation-of-incubation-times/
# # Model computes a semi-parametric fit to the data and compares it with classic parametric fits (lognormal, weibull, gamma)
# # candidate with the lowest BIC is selected
# 
# mindate <- min(as.Date(df_long$lastcontact)) - 60
# # df_long$a_minus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$symptom.onset) - mindate - 25),
# #                           as.numeric(as.Date(df_long$lastcontact) - mindate))
# # df_long$a_plus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$lastcontact) - mindate + 1),
# #                          as.numeric(as.Date(df_long$lastcontact) - mindate + 1))
# # df_long$b_minus <- as.numeric(as.Date(df_long$symptom.onset) - mindate)
# # df_long$b_plus <- as.numeric(as.Date(df_long$symptom.onset) - mindate + 1)
# 
# data.epi <- df_longg # excluded cases with upper < lower bound
# data.epi$symptom.onset.num <- as.numeric(as.Date(data.epi$symptom.onset) - mindate); summary(data.epi$symptom.onset.num)
# data.epi$lastcontact.num <- as.numeric(as.Date(data.epi$lastcontact) - mindate); summary(data.epi$lastcontact.num)
# # data.epi$exposure.upper.num <- as.numeric(as.Date(data.epi$lastcontact) - mindate); summary(data.epi$exposure.upper.num)
# # data.epi$exposure.lower.num <- ifelse(df_longg$type == 2, as.numeric(as.Date(data.epi$symptom.onset)-mindate-21),
# #                                       as.numeric(as.Date(data.epi$lastcontact)-mindate))
# data.epi$exposure.lower.num <- data.epi$a_minus
# data.epi$exposure.upper.num <- data.epi$a_plus #- 1 # because EpiLPS inherently accounts for censoring
# summary(as.numeric(data.epi$exposure.upper.num - data.epi$exposure.lower.num))
# summary(data.epi$a_plus - data.epi$a_minus)
# 
# # incubation period window
# data.epi$tL <- data.epi$symptom.onset.num - data.epi$exposure.upper.num
# # data.epi$tR <- data.epi$symptom.onset.num - data.epi$exposure.lower.num + 1
# data.epi$tR <- data.epi$symptom.onset.num - data.epi$exposure.lower.num + 1
# data.epi$duration <- (data.epi$tR - data.epi$tL)
# hist(data.epi$duration)
# summary(data.epi$duration)
# 
# data.epi$tL <- ifelse(data.epi$tL == -1, 0, data.epi$tL)
# 
# # data.epi <- data.epi[!is.na(data.epi$sexual), ]
# # data.epi <- data.epi[data.epi$duration > 5, ]
# 
# dataIncub <- data.frame(tL = data.epi$tL, tR = data.epi$tR, sexual = data.epi$sexual)
# head(dataIncub); dim(dataIncub)
# summary(dataIncub$tR - dataIncub$tL)
# 
# set.seed(222)
# # out <- estimIncub(dataIncub[dataIncub$sexual == 2, c(1:2)], verbose = TRUE)
# out <- estimIncub(dataIncub[,c(1:2)], K = 10, verbose = TRUE)
# gridExtra::grid.arrange(plot(out, typ = 'incubwin'), plot(out, type = "pdf"), nrow = 1)
# out$stats
# library(ggplot2)
# library(gridExtra)
# grid.arrange(plot(out, typ = "incubwin"), plot(out, type = "pdf"), nrow = 1)
