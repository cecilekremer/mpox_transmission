
## Data containing one row for each case x contact x exposure
load('data/contact_data_ALL_100325.RData')
data <- data.all
dim(data); length(unique(data$ID2))


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
               (data$contact3_lastcontact <= data$symptom.onset) | (data$contact3_lastcontact <= data$symptom.onset), ]
data <- data[!is.na(data$ID2), ]
dim(data); length(unique(data$ID2))

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

## Raw incubation periods
raw_incubation_periods <- df_long %>%
  pivot_longer(cols = c('symptom.onset', 'rash_onset'),
               names_to = 'to_what',
               values_to = 'to_date') %>%
  mutate(to_date = as.Date(to_date)) %>%
  pivot_longer(cols = c('lastcontact'),
               names_to = "from_what",
               values_to = 'from_date') %>%
  mutate(time_diff = difftime(to_date, from_date, units = 'days') %>% as.numeric()) 
raw_incubation_periods <- raw_incubation_periods[raw_incubation_periods$time_diff >= 0, ]
library(ggplot2)
raw_incubation_periods %>%
  ggplot(aes(x = time_diff, fill = sexual)) +
  geom_histogram() +
  facet_grid(~to_what) +
  ggthemes::scale_fill_few("Dark") +
  theme_bw() +
  labs(x = 'time difference (days)', title = 'Distribution of times from last contact to symptom/rash onset')
summary(raw_incubation_periods$time_diff)

data %>%
  mutate(time_diff = difftime(rash_onset, symptom.onset, units = 'days') %>%
           as.numeric()) %>%
  ggplot(aes(x = time_diff, fill = as.factor(agecat))) +
  geom_histogram() +
  facet_wrap(~transm_sexual) +
  ggthemes::scale_fill_few('Dark') +
  theme_bw() +
  labs(x = 'time (days) between symptom and rash onset')

###------------------------------------------------------------
### Data for Stan model

df_long <- df_long[df_long$symptom.onset >= df_long$lastcontact, ]
dim(df_long); length(unique(df_long$ID2))
sum(df_long$rash_onset >= df_long$lastcontact, na.rm = T)

df.unique <- df_long %>%
  distinct(ID2, .keep_all = TRUE)
summary(as.numeric(df.unique$rash_onset - df.unique$symptom.onset, na.rm = T))

# df_long <- df_long[df_long$contact_num == 'contact1', ]

# 114 individuals that report only one contact
df_long <- df_long[df_long$num_contact_mpox == 1, ]
dim(df_long); length(unique(df_long$ID2))
table(df_long$contact_num)

table(df_long$ID_code)
sum(is.na(df_long$ID_code))

# # Individuals with only single exposure
# table(df_long$type)
# df_long <- df_long[df_long$type == 1, ]
# dim(df_long); length(unique(df_long$ID2))

# df_long <- df_long[df_long$type == 1, ]
table(df_long$sexual) # 1 = yes, 2 = no

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
      # if(!is.na(df_long$transm_sexual[i]) & df_long$transm_sexual[i] == 1){
      #   df_long$sexual2[i] <- 1 # otherwise sexual if this is the most likely hypothesis
      # }else{
      #   df_long$sexual2[i] <- 2
      # }
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
table(df_long$sexual)
table(df_long$sexual1)
table(df_long$sexual2)
table(df_long$sexual3)
table(df_long$agecat)

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
# df_longg <- df_longg[df_longg$sexual2 == 2, ]
# df_longg <- df_longg[df_longg$agecat == 1, ]

# data_stan <- list(
#   N = nrow(df_longg),
#   a_minus = df_longg$a_minus,
#   a_plus = df_longg$a_plus,
#   b_minus = df_longg$b_minus,
#   b_plus = df_longg$b_plus,
#   # incubation = 0,
#   upper_bound = as.numeric(max(data.all$date) - mindate), # to account for truncation
#   incl_cov = 2, # include covariate? 1 = no, 2 = yes
#   sexual = ifelse(df_longg$sexual1 == 2, 0, 1),
#   # sexual = ifelse(df_longg$agecat == 1, 0, 1), # 0 = adult, 1 = child
#   n_steps_prior = df_longg$a_plus - df_longg$a_minus
# )
# table(data_stan$sexual)
# df_longg <- df_longg[!is.na(df_longg$ID_code), ]
dim(df_longg)

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
  age = ifelse(df_longg$agecat == 3, 0, 1), # 0 = child, 1 = adult
  # kamituga = ifelse(is.na(df_longg$ID_code), 1, 0),
  kamituga = rep(0, nrow(df_longg)),
  n_steps_prior = df_longg$a_plus - df_longg$a_minus
)
table(data_stan$sexual)
table(data_stan$age)
table(data_stan$kamituga)
# table(data_stan$kamituga, data_stan$age)
# table(data_stan$kamituga, data_stan$sexual)

# mod_Weibull <- rstan::stan_model('./code/stan_incub_ward.stan', model_name = 'mod_Weibull')
# save(mod_Weibull, file = 'mod_Weibull.RData')
# mod_Weibull_noTrunc <- rstan::stan_model('./code/stan_incub_ward_noTrunc.stan', model_name = 'mod_Weibull_noTrunc')
# save(mod_Weibull_noTrunc, file = 'mod_Weibull_noTrunc.RData')
# mod_Weibull_prob <- rstan::stan_model('./code/stan_incub_ward_prob_v2.stan', model_name = 'mod_Weibull_prob')
# save(mod_Weibull_prob, file = 'mod_Weibull_prob.RData')

# mod_Weibull_reg <- rstan::stan_model('./code/incubation_final/stan_incub_ward_reg.stan', model_name = 'mod_Weibull_reg')
# save(mod_Weibull_reg, file = 'mod_Weibull_reg.RData')

table(df_longg$type) # 1 = single exposure

# load('mod_Weibull.RData')
# load('mod_Weibull_noTrunc.RData')
# load('mod_Weibull_prob.RData')
load('mod_Weibull_reg.RData')

# initf <- function(chain_id = 1){
#   meaninit = rep(5, data_stan$incl_cov)
#   alphainit = rep(1, data_stan$incl_cov)
#   dim(meaninit) = data_stan$incl_cov
#   dim(alphainit) = data_stan$incl_cov
#   
#   list(mean_ = meaninit,
#        alpha = alphainit
#   )
# }
nrchains = 4
# init_ll <- lapply(1:nrchains, function(id) initf(chain_id = id))

options(mc.cores=parallel::detectCores())
fit_Weibull <- rstan::sampling(mod_Weibull_reg, data = (data_stan), 
                               # init = init_ll,
                               # init = '0',
                               # init_r = 10,
                               verbose = TRUE,
                               chains = 4, iter = 10000)
# rstan::pairs.stanfit(fit_Weibull)

## Results

### Regression model
rstan::traceplot(fit_Weibull, pars = c("theta[1]","theta[2]","theta[3]","alpha"), inc_warmup = TRUE)
quantile(as.matrix(fit_Weibull)[,'theta[1]'], probs = c(0.025,0.5,0.975)) # intercept
quantile(as.matrix(fit_Weibull)[,'theta[2]'], probs = c(0.025,0.5,0.975)) # sexual
quantile(as.matrix(fit_Weibull)[,'theta[3]'], probs = c(0.025,0.5,0.975)) # age
# quantile(as.matrix(fit_Weibull)[,'theta[4]'], probs = c(0.025,0.5,0.975)) # kamituga
# quantile(as.matrix(fit_Weibull)[,'theta[4]'], probs = c(0.025,0.5,0.975)) # age x sexual
quantile(as.matrix(fit_Weibull)[,'alpha'], probs = c(0.025,0.5,0.975)) # Weibull shape alpha

### Incubation period stratified

# sumfit <- summary(fit_Weibull, pars = c("median_[1]","mean_[1]","sd_[1]"), probs = c(0.025, 0.5, 0.975))
# sumfit$summary
# rstan::traceplot(fit_Weibull, pars = c("median_[1]","mean_[1]","sd_[1]"), inc_warmup = TRUE) ## Non-sexual
# rstan::traceplot(fit_Weibull, pars = c("median_[2]","mean_[2]","sd_[2]"), inc_warmup = TRUE) ## Sexual

# # fit_Weibull
# # colnames(as.matrix(fit_Weibull))
# hist(as.matrix(fit_Weibull)[,'median_[1]']); summary(as.matrix(fit_Weibull)[,'median_[1]']) ## Non-sexual
# quantile(as.matrix(fit_Weibull)[,'median_[1]'], probs = c(0.025,0.5,0.975))
# hist(as.matrix(fit_Weibull)[,'median_[2]']); summary(as.matrix(fit_Weibull)[,'median_[2]']) ## Sexual
# quantile(as.matrix(fit_Weibull)[,'median_[2]'], probs = c(0.025,0.5,0.975))
# 
# hist(as.matrix(fit_Weibull)[,'mean_[1]']); summary(as.matrix(fit_Weibull)[,'mean_[1]']) ## Non-sexual
# quantile(as.matrix(fit_Weibull)[,'mean_[1]'], probs = c(0.025,0.5,0.975))
# hist(as.matrix(fit_Weibull)[,'mean_[2]']); summary(as.matrix(fit_Weibull)[,'mean_[2]']) ## Sexual
# quantile(as.matrix(fit_Weibull)[,'mean_[2]'], probs = c(0.025,0.5,0.975))
# 
# hist(as.matrix(fit_Weibull)[,'sd_[1]']); summary(as.matrix(fit_Weibull)[,'sd_[1]']) ## Non-sexual
# quantile(as.matrix(fit_Weibull)[,'sd_[1]'], probs = c(0.025,0.5,0.975))
# hist(as.matrix(fit_Weibull)[,'sd_[2]']); summary(as.matrix(fit_Weibull)[,'sd_[2]']) ## Sexual
# quantile(as.matrix(fit_Weibull)[,'sd_[2]'], probs = c(0.025,0.5,0.975))


# hist(as.matrix(fit_Weibull)[,'limit_val_[1]']); summary(as.matrix(fit_Weibull)[,'limit_val_[1]'])
library(posterior)
?rhat
rhat(fit_Weibull)

###---------------------------------------------------------------------------------------
### EpiLPS semi-parametric estimation

library(EpiLPS)
# ?estimIncub
# https://statsandr.com/blog/epilps-for-estimation-of-incubation-times/
# Model computes a semi-parametric fit to the data and compares it with classic parametric fits (lognormal, weibull, gamma)
# candidate with the lowest BIC is selected

mindate <- min(as.Date(df_long$lastcontact)) - 60
# df_long$a_minus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$symptom.onset) - mindate - 25),
#                           as.numeric(as.Date(df_long$lastcontact) - mindate))
# df_long$a_plus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$lastcontact) - mindate + 1),
#                          as.numeric(as.Date(df_long$lastcontact) - mindate + 1))
# df_long$b_minus <- as.numeric(as.Date(df_long$symptom.onset) - mindate)
# df_long$b_plus <- as.numeric(as.Date(df_long$symptom.onset) - mindate + 1)

data.epi <- df_longg # excluded cases with upper < lower bound
data.epi$symptom.onset.num <- as.numeric(as.Date(data.epi$symptom.onset) - mindate); summary(data.epi$symptom.onset.num)
data.epi$lastcontact.num <- as.numeric(as.Date(data.epi$lastcontact) - mindate); summary(data.epi$lastcontact.num)
# data.epi$exposure.upper.num <- as.numeric(as.Date(data.epi$lastcontact) - mindate); summary(data.epi$exposure.upper.num)
# data.epi$exposure.lower.num <- ifelse(df_longg$type == 2, as.numeric(as.Date(data.epi$symptom.onset)-mindate-21),
#                                       as.numeric(as.Date(data.epi$lastcontact)-mindate))
data.epi$exposure.lower.num <- data.epi$a_minus
data.epi$exposure.upper.num <- data.epi$a_plus #- 1 # because EpiLPS inherently accounts for censoring
summary(as.numeric(data.epi$exposure.upper.num - data.epi$exposure.lower.num))
summary(data.epi$a_plus - data.epi$a_minus)

# incubation period window
data.epi$tL <- data.epi$symptom.onset.num - data.epi$exposure.upper.num
# data.epi$tR <- data.epi$symptom.onset.num - data.epi$exposure.lower.num + 1
data.epi$tR <- data.epi$symptom.onset.num - data.epi$exposure.lower.num + 1
data.epi$duration <- (data.epi$tR - data.epi$tL)
hist(data.epi$duration)
summary(data.epi$duration)

data.epi$tL <- ifelse(data.epi$tL == -1, 0, data.epi$tL)

dataIncub <- data.frame(tL = data.epi$tL, tR = data.epi$tR)
head(dataIncub); dim(dataIncub)
summary(dataIncub$tR - dataIncub$tL)

set.seed(222)
out <- estimIncub(dataIncub, verbose = TRUE)
out$stats
library(ggplot2)
library(gridExtra)
grid.arrange(plot(out, typ = "incubwin"), plot(out, type = "pdf"), nrow = 1)
