
load('data/contact_data.RData')
data <- data.contact.clean

##---------------------------------------------------------------
## Individuals reporting only one contact (N = 86)

data.si <- (data[data$contact1_id != '' & data$contact2_id == '' & data$contact3_id == '', ])
data.si <- data.si[data.si$contacts != '', ]
all(nchar(data.si$contacts) <= 3 )
data.si$contacts <- as.numeric(data.si$contacts)

summary(data.si$symptom.onset)
data.si$symptom.onset.secondary <- NA
for(i in 1:dim(data.si)[1]){
  idcon <- data.si$contacts[i]
  dt <- data$symptom.onset[data$ID == idcon] # take from full dataset
  data.si$symptom.onset.secondary[i] <- ifelse(is.na(dt), NA, as.Date(dt))
}
data.si$symptom.onset.secondary <- as.Date(data.si$symptom.onset.secondary)

data.si$serial.interval <- as.numeric(data.si$symptom.onset.secondary - data.si$symptom.onset)

summary(data.si$serial.interval)
hist(data.si$serial.interval[data.si$contact1_sexual == 1])
hist(data.si$serial.interval[data.si$contact1_sexual == 2])

##--------------------------------------------------------------
## Reporting multiple contacts

data.si <- data[data$contacts != '', ]

# Make sure all contacts are included as cases
all.contact.ids <- unique(c(as.numeric(gsub("SCREEN", "", data.si$contact1_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact2_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact3_id))
))
table(all.contact.ids %in% data.si$ID)
all.case.ids <- unique(c(all.contact.ids, unique(data$ID[data$contacts != ''])))

data.si <- data[data$ID %in% all.case.ids, ]
data.si <- data.si[,c(1:3,5,9:11,13:15,17,18,19,22,23,30:33,38:41,46:49,62:89,100)]

data.si$ID_orig <- data.si$ID # save original IDs
data.si$ID <- 1:dim(data.si)[1] # IDs need to be 1 to 171

table(data.si$contact1_rel)
# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)

## Initial infector-infectee matrix
NCases <- length(unique(data.si$ID))
Case <- data.si$ID
infectors <- list()
PossibleInfector <- matrix(nrow = NCases, ncol = 4)
for(i in 1:NCases){
  inf <- as.numeric(unlist(strsplit(data.si$contacts[i], ",")))
  infectors[[i]] <- data.si$ID[data.si$ID_orig %in% inf]
  if(length(infectors[[i]]) > 0){
    PossibleInfector[i, 1:length(infectors[[i]])] <- infectors[[i]]
  }else{
    PossibleInfector[i, 1] <- 0 # no contacts among included cases
  }
}

data.si$contact_list <- lapply(infectors, function(x){
  if(length(x) == 0){
    return(NA)
  }else{
    return(as.numeric(x))
  }
})

min.si <- 5 # absolute value of the max. allowed negative serial interval

rm(Case); rm(NCases); rm(i)
library(igraph)

# ## Sample networks
# source('code/functions/fun_network.R')
# 
# num.nets <- 1000
# trees <- matrix(NA, nrow = num.nets, ncol = length(unique(data.si$ID)) + 1)
# onsets <- matrix(NA, nrow = num.nets, ncol = length(unique(data.si$ID)) + 1)
# 
# # ptm <- proc.time()
# for(i in 1:num.nets){
#   net <- setup_network(case.ids = data.si$ID,
#                        cluster = NA,
#                        contact.list = data.si$contact_list,
#                        infector.mat = PossibleInfector,
#                        symptom.onset = as.Date(data.si$symptom.onset),
#                        helper.date = as.Date(data.si$date), # date of questionaire
#                        min.si = min.si, # absolute value of the max. allowed negative serial interval
#                        max.si = 30
#   )
# 
#   trees[i, ] <- c(i, net$network)
#   onsets[i, ] <- c(i, net$onset.times)
# 
#   if(i%%100 == 0){
#     print(i)
#   }
# 
# }
# # proc.time() - ptm
# save(trees, file = 'code/trees.RData')
# save(onsets, file = 'code/onsets.RData')

## MCMC to estimate serial interval using normal distribution
load('code/trees.RData')
load('code/onsets.RData')
source('code/functions/fun_mcmc_si.R')
nrun <- 2000000
burnin <- 0.4
thin <- 100
updatefr <- 2 # update network every other run
out <- estimate_si(case.ids = data.si$ID,
                   # cluster = NA,
                   # contact.list = data.si$contact_list,
                   # infector.mat = PossibleInfector,
                   # symptom.onset = as.Date(data.si$symptom.onset),
                   # helper.date = as.Date(data.si$date), # date of questionaire
                   # min.si = min.si, # absolute value of the max. allowed negative serial interval
                   # max.si = 200,
                   networks = trees,
                   onsets = onsets,
                   max.si = 30,
                   start.parms = c(1, 1), # mean = theta[1], sd = theta[2]
                   tuning.parms = c(0.2, 0.2),
                   mcmc.runs = nrun,
                   burnin = burnin,
                   thin = thin,
                   update.freq = updatefr
)

length(unique(out$networkIDs))

## Check MCMC convergence
library(coda)
# thin <- 1
chain1 <- coda::mcmc(out$parms[seq(1, (nrun-(burnin*nrun))/thin), ])
colnames(chain1) = c("mean SI","sd SI")
autocorr.plot(chain1[,c(1,2)]) # see if thinning has to be increased
plot(chain1[,c(1,2)])
summary(chain1)

par(mfrow=c(2,1))
plot(out$log_post[seq(1, (nrun-(burnin*nrun))/thin)], type="l", ylab="Posterior")
plot(out$log_lik[seq(1, (nrun-(burnin*nrun))/thin)], type="l", ylab="Likelihood")



