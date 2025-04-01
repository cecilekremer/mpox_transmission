
library(rstan)
load('code/stanmod_adj.RData')

# ## Stan model for different distributions
# distributions <- c("weibull", "gamma", "lognormal")
# code <- sprintf("
#   data{
#     int<lower=1> N; // number of observations: individuals x visits x possible contacts
#   int<lower=1> M; // number of individuals x visits
#   
#   // Data
#   array[N] real<lower=0> t_s; // symptom onset
#   array[N] real<lower=0> t_e_r; // upper bound exposure
#   array[N] real<lower=0> t_e_l; // lower bound exposure
#   array[N] int<lower=0,upper=1> censored; // single (0) or multiple (1) exposure
#   
#   array[M] int<lower=1> n_contacts; // number of contacts for each case
#   
#   // prior on duration of lesions; assuming infectious period corresponds to lesion period
#   real<lower=0> log_mu_lesions; 
#   real<lower=0> sigma_lesions;
#   
#   // Auxiliary data
#   array[N] int<lower=0> n_steps_prior; // number of days over which to integrate
#   }
#   transformed data{
#   vector[N] log_prior_contacts; // prior on who is the infecting contact
#   array[N] vector[max(n_steps_prior)] log_prior_exposures; // prior on when exposure occurred
#   
#   for(m in 1:M){
#     // for(i in 1:n_contacts[m]){
#       log_prior_contacts[m] = log(1./n_contacts[m]); // uniform prior over contacts
#       // }
#   }
#   
#   for(i in 1:N){
#     int n = n_steps_prior[i];
#     array[n] real ts;
#     ts = linspaced_array(n, t_e_r[i] - n, t_e_r[i]);
#     
#     if(n > 0){
#       for(j in 1:n){
#         log_prior_exposures[i][j] = lognormal_lccdf(t_e_r[i] - ts[j] + 1 | log_mu_lesions, sigma_lesions);
#       }
#       // normalize
#       log_prior_exposures[i][1:n] = log(softmax(log_prior_exposures[i][1:n]));
#       // if(i == 1){
#       //   print(log_prior_exposures[i]);
#       // }
#     }
#   }
# }
#   parameters{
#       real<lower=0> par[2]; // distributional parametsrs
#   }
#   transformed parameters{
#     // vector[M] tE; 	// infection moment
#   vector[M] log_lik; // log-likelihood vector
#   
#   for(m in 1:M){
#     int n = n_contacts[m];
#     int cnt = 1;
#     vector[n] ll; // placeholder for logliks for each contact
#     
#     // integrate over all contacts
#     for(i in 1:n){
#       
#       real ll_c;
#       
#       if(censored[i] == 1){ // multiple occasions
#       
#         int k = n_steps_prior[i];
#         vector[k] ll_exp; // loglik of exposures prior to most recent contact
#         vector[k] ll_exp2 = rep_vector(-log(k), k); // ll denominator to account for truncation
#         array[k] real ts; // times of possible exposure to integrate over
#         
#         ts = linspaced_array(k, t_e_r[i] - k, t_e_r[i]);
#         
#         for(j in 1:k){
#           
#           real tE = t_s[i] - ts[j];
#           // real bound = t_e_r[i] - ts[j] + 21;
#           real this_ll;
#           
#           this_ll = %s_lpdf(tE | par[1], par[2]);
#           
#           ll_exp[j] = log_prior_exposures[i][j] + this_ll;
#           
#         }
#         
#         ll_c = log_sum_exp(ll_exp);
#         
#       }else{ // single occasion
#         
#         real tE = t_s[i] - t_e_l[i]; // exposure time = infection time
#         ll_c = %s_lpdf(tE | par[1], par[2]);
#         
#       }
#       
#       // this contact's loglik
#       ll[cnt] = log_prior_contacts[i] + ll_c;
#       cnt += 1;
#       
#     }
#     
#     // sum over contacts
#     log_lik[m] = log_sum_exp(ll);
#   }
#   }
#   model{
#     // Contribution to likelihood of incubation period
#   target += sum(log_lik);
#   }
# ", distributions, distributions)
# 
# names(code) <- distributions
# 
# models <- mapply(stan_model, model_code = code)

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

data %>%
  mutate(time_diff = difftime(rash_onset, symptom.onset, units = 'days') %>%
           as.numeric()) %>%
  ggplot(aes(x = time_diff, fill = as.factor(agecat))) +
  geom_histogram() +
  facet_wrap(~transm_sexual) +
  ggthemes::scale_fill_few('Dark') +
  theme_bw() +
  labs(x = 'time (days) between symptom and rash onset')

##--------------------------------------------------------------------------------------------------
## Data for Stan model
# https://github.com/HopkinsIDD/mpox_cladeIb_incubation_period/blob/main/analysis/01_prepare_data.R

df_long <- df_long[df_long$symptom.onset >= df_long$lastcontact, ]

# Set EL to one month before first reported case since we only have date of most recent contact
min_exp_date <- min(as.Date(df_long$lastcontact) - 30)

library(purrr)
make_map <- function(vec, u_ref) {
  map_dbl(vec, ~ which(u_ref == .))
}
make_starts <- function(map_to_i) {
  u_i <- unique(map_to_i)
  starts <- rep(0, length(u_i))
  for (i in 1:length(u_i)) {
    starts[i] <- which(map_to_i == i)[1]
  }
  starts
}
make_ends <- function(map_to_i) {
  u_i <- unique(map_to_i)
  ends <- rep(0, length(u_i))
  for (i in 1:length(u_i)) {
    ends[i] <- which(map_to_i == i) %>% last()
  }
  ends
}
u_ids <- unique(df_long$ID2)
map_to_id <- make_map(vec = df_long$ID2, u_ref = u_ids)
starts <- make_starts(map_to_id)
ends <- make_ends(map_to_id)
n_contacts <- ends - starts + 1

ref_time <- min(as.Date(df_long$symptom.onset), as.Date(df_long$lastcontact)) - 60
t_s <- as.numeric(difftime(df_long$symptom.onset, ref_time, units = 'days'))
t_e_r <- as.numeric(difftime(as.Date(df_long$lastcontact), ref_time, units = 'days'))
t_e_l <- as.numeric(difftime(as.Date(df_long$lastcontact), ref_time, units = 'days'))
O <- 60 # window for incubation period
Q <- 30 # window for reporting
all_times <- seq(min(t_s) - O - Q - 1, max(as.Date(df_long$date)))

t_e_l[df_long$type == 2] <- t_e_r[df_long$type == 2] - O

stan_data <- list(
  N = nrow(df_long),
  M = df_long %>% distinct(ID2) %>% nrow(),
  # n_contacts <- (df_long %>% group_by(ID2) %>% mutate(n_contacts = n()) %>% ungroup())$n_contacts
  t_s = t_s,
  t_e_r = t_e_r,
  t_e_l = t_e_l,
  censored = ifelse(df_long$type == 1, 0, 1),
  n_contacts = n_contacts,
  log_mu_lesions = log(14.1),
  sigma_lesions = 0.5,
  n_steps_prior = t_e_r - t_e_l
)

fit <- mapply(sampling, models, list(stan_data), iter = 10000, warmup = 3000, chain = 4)

