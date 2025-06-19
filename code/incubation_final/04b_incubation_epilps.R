
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

## Select data to include in estimation
df_long <- df_long[df_long$symptom.onset >= df_long$lastcontact, ]
dim(df_long); length(unique(df_long$ID2))
sum(df_long$rash_onset >= df_long$lastcontact, na.rm = T)

df.unique <- df_long %>%
  distinct(ID2, .keep_all = TRUE)
summary(as.numeric(df.unique$rash_onset - df.unique$symptom.onset, na.rm = T))

# 114 individuals that report only one contact
df_long <- df_long[df_long$num_contact_mpox == 1, ]
dim(df_long); length(unique(df_long$ID2))
table(df_long$contact_num)

table(df_long$ID_code)
sum(is.na(df_long$ID_code))

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
# exposure lower
df_long$a_minus <- ifelse(df_long$type == 2, as.numeric(as.Date(df_long$symptom.onset) - mindate - 21), # multiple exposure (type = 2)
                          as.numeric(as.Date(df_long$lastcontact) - mindate))
# exposure upper
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
summary(df_longg$a_plus - df_longg$a_minus)

##----------------------------------------------
## Non-parametric estimation using EpiLPS

library(EpiLPS)

data.epi <- df_longg
data.epi$symptom.onset.num <- as.numeric(as.Date(data.epi$symptom.onset) - mindate); summary(data.epi$symptom.onset.num)
data.epi$lastcontact.num <- as.numeric(as.Date(data.epi$lastcontact) - mindate); summary(data.epi$lastcontact.num)

summary(data.epi$symptom.onset.num - data.epi$lastcontact.num)
hist(data.epi$symptom.onset.num - data.epi$lastcontact.num)

mindate = mindate + 1

## lower and upper bound on exposure
data.epi$exposure.lower.num <- ifelse(data.epi$type == 2, as.numeric(as.Date(data.epi$symptom.onset) - mindate - 21), # multiple exposures (type = 2)
                                      as.numeric(as.Date(data.epi$lastcontact) - mindate) - 0.5)
data.epi$exposure.upper.num <- ifelse(data.epi$type == 2, as.numeric(as.Date(data.epi$lastcontact) - mindate),
                                      as.numeric(as.Date(data.epi$lastcontact) - mindate) + 0.5) # day after reported exposure (to avoid zero-width intervals)
summary(as.numeric(data.epi$exposure.upper.num - data.epi$exposure.lower.num))

# incubation period window
data.epi$tL <- data.epi$symptom.onset.num - data.epi$exposure.upper.num
# data.epi$tR <- data.epi$symptom.onset.num - data.epi$exposure.lower.num + 1
data.epi$tR <- data.epi$symptom.onset.num - data.epi$exposure.lower.num
data.epi$duration <- (data.epi$tR - data.epi$tL)
hist(data.epi$duration)
summary(data.epi$duration)

# data.epi <- data.epi[data.epi$duration < 20, ]
# data.epi <- data.epi[data.epi$type == 2, ]
data.epi <- data.epi[data.epi$sexual1 == 2, ]

# dataIncub <- data.frame(tL = data.epi$tL+1, tR = data.epi$tR+1)#, sexual = data.epi$sexual)
dataIncub <- data.frame(tL = data.epi$tL, tR = data.epi$tR)
head(dataIncub); dim(dataIncub)
hist(dataIncub$tR)
hist(dataIncub$tL)

sum(dataIncub$tL == dataIncub$tR)

set.seed(222)
# out <- estimIncub(dataIncub[dataIncub$sexual == 2, c(1:2)], verbose = TRUE)
out <- estimIncub(dataIncub, verbose = TRUE)
library(ggplot2)
library(gridExtra)
gridExtra::grid.arrange(plot(out, typ = 'incubwin'), plot(out, type = "pdf"), nrow = 1)
out$stats




