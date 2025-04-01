
# https://khayatrayen.github.io/MCMC.html

# https://medium.com/@tinonucera/bayesian-linear-regression-from-scratch-a-metropolis-hastings-implementation-63526857f191


load('/lustre1/scratch/326/vsc32693/simNet/contact_data_100325.RData')
data <- data.contact.clean

data.si <- data[data$contacts != '', ]
# table(data.si$contact1_sexual)
# table(data.si$transm_sexual) # not specified which contact was sexual

# Make sure all contacts are included as cases
all.contact.ids <- unique(c(as.numeric(gsub("SCREEN", "", data.si$contact1_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact2_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact3_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact4_id))
))
# table(all.contact.ids %in% data.si$ID)
all.case.ids <- unique(c(all.contact.ids, unique(data$ID[data$contacts != ''])))

data.si <- data[data$ID %in% all.case.ids, ]
# data.si <- data.si[,c(1:3,5,9:11,13:15,17,18,19,22,23,30:33,38:41,46:49,62:89,100)]
# dim(data.si)

data.si$ID_orig <- data.si$ID # save original IDs
data.si$ID <- 1:dim(data.si)[1] # IDs need to be 1 to N

# summary(data.si$agenum)
# 
# table(data.si$contact1_rel)
# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)

## Initial infector-infectee matrix + transmission routes
NCases <- length(unique(data.si$ID))
Case <- data.si$ID
infectors <- list()
routes <- list(); routes.hh <- list()
PossibleInfector <- matrix(nrow = NCases, ncol = 4)
TransmissionSexual <- matrix(NA, nrow = NCases, ncol = 4)
TransmissionHousehold <- matrix(NA, nrow = NCases, ncol = 4)

for(i in 1:NCases){
  inf <- as.numeric(unlist(strsplit(data.si$contacts[i], ",")))
  infectors[[i]] <- data.si$ID[data.si$ID_orig %in% inf]
  if(length(infectors[[i]]) > 0){
    PossibleInfector[i, 1:length(infectors[[i]])] <- infectors[[i]]
  }else{
    PossibleInfector[i, 1] <- 0 # no contacts among included cases
  }
  
  ## Sexual transmission
  sextrans <- numeric()
  # transmission routes: 1 = sexual, 2 = other (also if not specifically sexual reported) (data columns: 83,84,85,86)
  if(length(inf) > 0){
    for(c in 1:length(inf)){
      # sex <- ifelse(data.si[i, 82 + c] == 1 & !is.na(data.si[i, 82 + c]), 1, 2) # Scenario 2: NA = non-sexual
      # Scenario 1: NA = sexual if not family
      sex <- ifelse((is.na(data.si[i, 82 + c]) & data.si[i, 58 + c] %in% c(1,2)) | data.si[i, 82 + c] == 2, 2, 1)
      
      sextrans <- c(sextrans, sex)
    }
  }
  routes[[i]] <- sextrans
  if(length(routes[[i]]) > 0){
    TransmissionSexual[i, 1:length(routes[[i]])] <- routes[[i]]
  }
  
  ## Household transmission
  HHtrans <- numeric()
  # transmission routes: 1 = HH member, 2 = other (also if not specifically reported) (data columns: 59,60,61,62)
  if(length(inf) > 0){
    for(c in 1:length(inf)){
      HH <- ifelse(data.si[i, 58 + c] == 1 & !is.na(data.si[i, 58 + c]), 1, 2)
      HHtrans <- c(HHtrans, HH)
    }
  }
  
  routes.hh[[i]] <- HHtrans
  if(length(routes.hh[[i]]) > 0){
    TransmissionHousehold[i, 1:length(routes.hh[[i]])] <- routes.hh[[i]]
  }
  
}

data.si$contact_list <- lapply(infectors, function(x){
  if(length(x) == 0){
    return(NA)
  }else{
    return(as.numeric(x))
  }
})

data.si$transmission_list_sexual <- lapply(routes, function(x){
  if(length(x) == 0){
    return(NA)
  }else{
    return(as.numeric(x))
  }
})

data.si$transmission_list_household <- lapply(routes.hh, function(x){
  if(length(x) == 0){
    return(NA)
  }else{
    return(as.numeric(x))
  }
})


min.si <- 5 # absolute value of the max. allowed negative serial interval
max.si <- 30

rm(Case); rm(NCases); rm(i)
library(igraph)

# ## Sample networks --> VSC
# source('code/01_baseline/functions/fun_network.R')

## MCMC to estimate serial interval using normal distribution

# Change network when changing min and max si !!!
load('/lustre1/scratch/326/vsc32693/simNet/trees_100325_base.RData')
load('/lustre1/scratch/326/vsc32693/simNet/onsets_100325_base.RData')
# load('code/01_baseline/routes_221024_base.RData')

# # check if all networks are unique
# library(mgcv)
# dim(uniquecombs(trees[,c(-1)]))

trees <- trees[!is.na(trees[,1]),]
onsets <- onsets[trees[,1],]

# trees <- trees[1:1000, ]
# onsets <- onsets[1:1000, ]
# transroutes <- transroutes[1:1000, ]

trees[,1] <- 1:dim(trees)[1]
onsets[,1] <- 1:dim(onsets)[1]
# transroutes[,1] <- 1:dim(transroutes)[1]

source('/lustre1/scratch/326/vsc32693/simNet/fun_reg_si_age.R')

nrun <- 5000000
burnin <- 0.4
thin <- 200
updatefr <- 2 # update network every other run
out <- estimate_si_reg(case.ids = data.si$ID,
                   networks = trees,
                   onsets = onsets,
                   # routes = transroutes,
                   data = data.si,
                   max.si = max.si, ## Change network when changing min and max si !!!
                   start.parms = rep(1, 5),
                   tuning.parms = rep(0.05, 5),
                   mcmc.runs = nrun,
                   burnin = burnin,
                   thin = thin,
                   update.freq = updatefr
)


save.image(file = '/lustre1/scratch/326/vsc32693/simNet/SerialIntervalReg_agegroup_240325.RData')
