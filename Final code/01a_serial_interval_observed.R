
library(dplyr)
library(tidyr)

# setwd('./Final code/')


###------------------------------------------------------------
### GOMA

# load('data/contact_data_ALL_100325.RData')
load('./Final code/data/contact_data.RData')
dataGM <- data.all[!is.na(data.all$ID_code), ]
table(dataGM$contact_mpox) # 1 = yes
sum(dataGM$contacts != "")
sum(!(grepl(',', dataGM$contacts)) & dataGM$contacts != "")

# table(dataGM$contact1_sexual)

# symptom onset available
summary(data.all$symptom.onset)
sum(!is.na(data.all$symptom.onset))
sum(!is.na(dataGM$symptom.onset))

###------------------------------------------------------------
### KAMITUGA

load('./Final code/data/contact_data_all.RData')
data <- data.contact.clean

table(data$contact_mpox) # 1 = yes
sum(data$contacts != "") # number of participants for which at least one contact was included in the study
sum(!(grepl(',', data$contacts)) & data$contacts != "")

sum(!is.na(data$symptom.onset))

##---------------------------------------------------------------
## Individuals reporting only one contact 

data.si <- data[data$contacts != '', ]
excl <- which(grepl(',', data.si$contacts)) # exclude those reporting multiple contacts
data.si <- data.si[-excl, ]
all(nchar(data.si$contacts) <= 4 )
data.si$contacts <- as.numeric(data.si$contacts)
dim(data.si)
df <- data.si %>%
  mutate(across(matches("^contact\\d+_"), as.character))
df_long <- df %>%
  pivot_longer(cols = matches("^contact\\d+_"),
               names_to = c("contact_num", ".value"),
               names_pattern = "(contact\\d+)_(.*)")
df_long <- df_long[!is.na(df_long$symptom.onset), ]
# df_long <- df_long[df_long$num_contact_mpox == 1, ]
length(unique(df_long$ID))
df_long <- df_long[df_long$id != '', ]
table(df_long$contact_num[!is.na(df_long$type)]) # 1 = single, 2 = multiple/ongoing
# df_long <- df_long[!is.na(df_long$rel), ]
length(unique(df_long$ID))

data.si <- df_long

summary(data.si$symptom.onset)
data.si$symptom.onset.index <- NA
for(i in 1:dim(data.si)[1]){
  idcon <- data.si$contacts[i] # probable index case
  dt <- data$symptom.onset[data$ID == idcon] # take from full dataset
  data.si$symptom.onset.index[i] <- ifelse(is.na(dt), NA, as.Date(dt))
}
data.si$symptom.onset.index <- as.Date(data.si$symptom.onset.index)

data.si$serial.interval <- as.numeric(data.si$symptom.onset - data.si$symptom.onset.index)

sum(!is.na(data.si$serial.interval))
summary(data.si$serial.interval)
hist(data.si$serial.interval, breaks = 50)


## Exclude outlier (last contact was 45 days after index symptom onset)
data.si <- data.si[data.si$serial.interval <= 30, ]
data.si <- data.si[data.si$serial.interval >= -5, ]
data.si <- data.si[!is.na(data.si$serial.interval), ]
dim(data.si)
table(data.si$sexual)
summary(data.si$serial.interval)
sum(data.si$serial.interval < 0)

# Check whether contact != ID
sum(data.si$contacts == data.si$ID)
data.si <- data.si[data.si$contacts != data.si$ID, ]
dim(data.si)

# si.IDs <- data.si$ID
# save(si.IDs, file = 'IDs_serial_interval.RData')

# # Contacts reporting each other?
# ddd <- data.si[data.si$contacts %in% data.si$ID, ]
# ids <- c(ddd$contacts, ddd$ID)
# View(data.si[data.si$ID %in% ids,][,c("ID","contacts")])

# remove duplicate
data.si <- data.si[-which(data.si$ID==839)[1], ]
dim(data.si)

##-------------------------------------------------------------------------------------------
## Characteristics of the pairs with sexual transmission

summary(data.si$serial.interval)
dim(data.si)
sum(data.si$serial.interval < 0)

table(data.si$sexual) # 1 = sexual, 2 = non-sexual
table(data.si$sexual[data.si$serial.interval<0])

# data.si.id <- data.si[,c("ID","contacts")]
# save(data.si.id, file = './Final code/data/data_obs_SI.RData')

table(data.si$sexual, data.si$n.lesion.genital)
table(data.si$n.lesion.genital[is.na(data.si$sexual)])
table(data.si$n.lesion.anal[is.na(data.si$sexual)])

data.si$lesion_intimate <- ifelse(data.si$n.lesion.anal>0 | data.si$n.lesion.genital>0, 1, 0)
table(data.si$lesion_intimate[is.na(data.si$sexual)])
table(data.si$lesion_intimate[data.si$sexual == 2])

# HH members
table(data.si$rel[data.si$sexual == 2])
# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)
table(data.si$rel[is.na(data.si$sexual)])
table(data.si$famtype[is.na(data.si$sexual)]) # check if all rel=1 or rel=2 are parent/child/sibling

table(data.si$sexual)

## plot matrix?
table(data.si$occupation_prim[data.si$sexual == 1], data.si$rel[data.si$sexual == 1])
# xdat <- as.data.frame((table(data.si$occupation_prim[data.si$contact1_sexual == 1], data.si$contact1_rel[data.si$contact1_sexual == 1])))
# var1 <- xdat$Var1[xdat$Freq != 0]; var2 <- xdat$Var2[xdat$Freq != 0]
# xdat <- xdat[xdat$Var1 %in% unique(var1) & xdat$Var2 %in% unique(var2), ]
# pairs <- data.frame(xtabs(~Var1 + Var2, xdat))
# trans.mat <- matrix(xdat$Freq, ncol=5, nrow=7) # row = index, col = contact
# removeRow <- c()
# for(i in 1:dim(trans.mat)[1]){
#   if(sum(trans.mat[i, 1:5]) == 0){
#     removeRow <- c(removeRow, i)
#   }
# }
# trans.mat <- trans.mat[-removeRow, ]

# Presymptomatic transmission
sum(data.si$serial.interval < 0)
table(data.si$sexual[data.si$serial.interval < 0])
sum(data.si$serial.interval[is.na(data.si$sexual)] < 0)
sum(data.si$serial.interval[!is.na(data.si$sexual)] < 0)

##-------------------------------------------------------------------------------------------
##
## EpiLPS to obtain nonparametric estimate of the serial interval
##
##-------------------------------------------------------------------------------------------

# library(EpiLPS)
# source('code/serialinterval_final/estimSI_boot.R')
library(EpiDelays)

hist(data.si$serial.interval, breaks = 20)

xS <- data.frame(sL = data.si$serial.interval - 0.5, sR = data.si$serial.interval + 0.5)
set.seed(2022)
fitS <- estimSI(x = xS)
round(fitS$npestim, 2)

##---------------------------------------------------------------------------------------
## Sexual vs non-sexual transmission: assuming NAs are sexual (SCENARIO 1)

# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)
table(data.si$rel[is.na(data.si$sexual)])

data.sens1 <- data.si
# no sexual contact if child or sibling
data.sens1$sexual <- ifelse((is.na(data.sens1$sexual) & data.sens1$rel %in% c(1,2)) | data.sens1$sexual == 2, 2, 1)
data.sens1$sexual <- ifelse(is.na(data.sens1$sexual) | data.sens1$sexual == 1, 1, 2)
table(data.sens1$sexual)

## Non-sexual transmission
xNonSexual <- data.frame(sL = data.sens1$serial.interval[data.sens1$sexual == 2] - 0.5, 
                         sR = data.sens1$serial.interval[data.sens1$sexual == 2] + 0.5)
set.seed(2022)
fitNonSexual <- estimSI(x = xNonSexual)
round(fitNonSexual$npestim, 2)

## Sexual transmission
xSexual <- data.frame(sL = data.sens1$serial.interval[data.sens1$sexual == 1] - 0.5, 
                      sR = data.sens1$serial.interval[data.sens1$sexual == 1] + 0.5)
set.seed(2022)
fitSexual <- estimSI(x = xSexual)
round(fitSexual$npestim, 2)

##---------------------------------------------------------------------------------------
## Sexual vs non-sexual transmission: assuming NAs are non-sexual (SCENARIO 2)

data.sens2 <- data.si
data.sens2$sexual <- ifelse(is.na(data.sens2$sexual) | data.sens2$sexual == 2, 2, 1)
table(data.sens2$sexual)

## Non-sexual transmission
xNonSexual <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 2] - 0.5, 
                         sR = data.sens2$serial.interval[data.sens2$sexual == 2] + 0.5)
set.seed(2022)
fitNonSexual <- estimSI(x = xNonSexual)
round(fitNonSexual$npestim, 2)

## Sexual transmission
xSexual <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 1] - 0.5, 
                      sR = data.sens2$serial.interval[data.sens2$sexual == 1] + 0.5)
set.seed(2022)
fitSexual <- estimSI(x = xSexual)
round(fitSexual$npestim, 2)

##------------------------------------------------------------------------------------------------------------------
## Sexual vs non-sexual transmission: assuming presence of genital/anal lesions indicate sexual contact (SCENARIO 3)

data.sens3 <- data.si
for(i in 1:dim(data.sens3)[1]){
  if(is.na(data.sens3$sexual[i])){
    if(data.sens3$rel[i] %in% c(1,2)){
      data.sens3$sexual[i] <- 2
    }else if(is.na(data.sens3$n.lesion.anal[i]) & is.na(data.sens3$n.lesion.genital[i])){
      data.sens3$sexual[i] <- 2
    }else{ 
      if(data.sens3$n.lesion.anal[i] > 0 | data.sens3$n.lesion.genital[i] > 0){
        data.sens3$sexual[i] <- 1
      }else{
        data.sens3$sexual[i] <- 2
      }
    }
  }
}
table(data.sens3$sexual)

## Non-sexual transmission
xNonSexual <- data.frame(sL = data.sens3$serial.interval[data.sens3$sexual == 2] - 0.5, 
                         sR = data.sens3$serial.interval[data.sens3$sexual == 2] + 0.5)
set.seed(2022)
fitNonSexual <- estimSI(x = xNonSexual)
round(fitNonSexual$npestim, 2)

## Sexual transmission
xSexual <- data.frame(sL = data.sens3$serial.interval[data.sens3$sexual == 1] - 0.5, 
                      sR = data.sens3$serial.interval[data.sens3$sexual == 1] + 0.5)
set.seed(2022)
fitSexual <- estimSI(x = xSexual)
round(fitSexual$npestim, 2)

###------------------------------------------------------------
### Age? Non-sexual for adults vs children

table(data.sens2$agecat, data.sens2$sexual)
xChild <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat == 3] - 0.5, 
                      sR = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat == 3] + 0.5)
set.seed(2022)
fitChild <- estimSI(x = xChild)

xAdult <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat != 3] - 0.5, 
                     sR = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat != 3] + 0.5)
set.seed(2022)
fitAdult <- estimSI(x = xAdult)

round(fitAdult$npestim, 2)
round(fitChild$npestim, 2)

## If cutoff for children is at 15y
data.sens2$agecat2 <- ifelse(data.sens2$agenum < 15, 2, 1) # 1 = adult, 2 = child
table(data.sens2$agecat2, data.sens2$sexual)

xChild <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat2 == 2] - 0.5, 
                     sR = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat2 == 2] + 0.5)
set.seed(2022)
fitChild <- estimSI(x = xChild)

xAdult <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat2 != 2] - 0.5, 
                     sR = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat2 != 2] + 0.5)
set.seed(2022)
fitAdult <- estimSI(x = xAdult)

round(fitAdult$npestim, 2)
round(fitChild$npestim, 2)

###----------------------------------------------------------------
### Sexual transmission among adults

data.sens1$agecat2 <- ifelse(data.sens1$agenum < 15, 2, 1) # 1 = adult, 2 = child
data.sens3$agecat2 <- ifelse(data.sens3$agenum < 15, 2, 1) # 1 = adult, 2 = child

table(data.sens2$sexual[data.sens2$agecat == 1])

## Adult cutoff at 12y
xSexual <- data.frame(sL = data.sens3$serial.interval[data.sens3$sexual == 1 & data.sens3$agecat != 3] - 0.5,
                      sR = data.sens3$serial.interval[data.sens3$sexual == 1 & data.sens3$agecat != 3] + 0.5)
set.seed(2022)
fitSexual <- estimSI(x = xSexual)

xNonSexual <- data.frame(sL = data.sens3$serial.interval[data.sens3$sexual == 2 & data.sens3$agecat != 3] - 0.5,
                      sR = data.sens3$serial.interval[data.sens3$sexual == 2 & data.sens3$agecat != 3] + 0.5)
set.seed(2022)
fitNonSexual <- estimSI(x = xNonSexual)

round(fitSexual$npestim, 2)
round(fitNonSexual$npestim, 2)

## Adult cutoff at 15y
xSexual <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 1 & data.sens2$agecat2 != 2] - 0.5,
                      sR = data.sens2$serial.interval[data.sens2$sexual == 1 & data.sens2$agecat2 != 2] + 0.5)
set.seed(2022)
fitSexual <- estimSI(x = xSexual)

xNonSexual <- data.frame(sL = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat2 != 2] - 0.5,
                         sR = data.sens2$serial.interval[data.sens2$sexual == 2 & data.sens2$agecat2 != 2] + 0.5)
set.seed(2022)
fitNonSexual <- estimSI(x = xNonSexual)

round(fitSexual$npestim, 2)
round(fitNonSexual$npestim, 2)

###--------------------------------------------------
### Gender

xMale <- data.frame(sL = data.si$serial.interval[data.si$gender == 1] - 0.5,
                      sR = data.si$serial.interval[data.si$gender == 1] + 0.5)
set.seed(2022)
fitMale <- estimSI(x = xMale)

xFemale <- data.frame(sL = data.si$serial.interval[data.si$gender == 2] - 0.5,
                         sR = data.si$serial.interval[data.si$gender == 2] + 0.5)
set.seed(2022)
fitFemale <- estimSI(x = xFemale)

round(fitMale$npestim, 2)
round(fitFemale$npestim, 2)

###--------------------------------------------------
### Time between symptom onset and rash

summary(data.sens1$rash_onset)
summary(data.sens1$symptom.onset)
data.sens1$time.symptom.to.rash <- as.numeric(data.sens1$rash_onset - data.sens1$symptom.onset)
summary(data.sens1$time.symptom.to.rash)
aggregate(data.sens1$time.symptom.to.rash, by = list(data.sens1$sexual), FUN = sd, na.rm = T)
t.test(data.sens1$time.symptom.to.rash ~ data.sens1$sexual)

data.sens1$adult <- ifelse(data.sens1$agecat == 3, 0, 1)
aggregate(data.sens1$time.symptom.to.rash[data.sens1$sexual == 2], by = list(data.sens1$adult[data.sens1$sexual == 2]), FUN = sd, na.rm = T)
t.test(data.sens1$time.symptom.to.rash[data.sens1$sexual == 2] ~ data.sens1$adult[data.sens1$sexual == 2])


aggregate(data.sens1$time.symptom.to.rash[data.sens1$adult == 1], by = list(data.sens1$sexual[data.sens1$adult == 1]), FUN = sd, na.rm = T)
t.test(data.sens1$time.symptom.to.rash[data.sens1$adult == 1] ~ data.sens1$sexual[data.sens1$adult == 1])

