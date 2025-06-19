
load('/lustre1/scratch/326/vsc32693/simNet/contact_data.RData')
data <- data.contact.clean

##--------------------------------------------------------------
## Reporting multiple contacts

data.si <- data[data$contacts != '', ]
# table(data.si$contact1_sexual)
# table(data.si$transm_sexual) # not specified which contact was sexual

# Make sure all contacts are included as cases
all.contact.ids <- unique(c(as.numeric(gsub("SCREEN", "", data.si$contact1_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact2_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact3_id))
))
# table(all.contact.ids %in% data.si$ID)
all.case.ids <- unique(c(all.contact.ids, unique(data$ID[data$contacts != ''])))

data.si <- data[data$ID %in% all.case.ids, ]
data.si <- data.si[,c(1:3,5,9:11,13:15,17,18,19,22,23,30:33,38:41,46:49,62:89,100,111:112)]

data.si$ID_orig <- data.si$ID # save original IDs
data.si$ID <- 1:dim(data.si)[1] # IDs need to be 1 to 171

# table(data.si$contact1_rel)
# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)

## Set to only the sexual infector(s) if both variables indicate sexual transmission?
for(i in 1:dim(data.si)[1]){
  inf <- as.numeric(unlist(strsplit(data.si$contacts[i], ",")))
  if(length(inf) > 1){
    if(length(inf) == 2){
      routes <- c(data.si$contact1_sexual[i], data.si$contact2_sexual[i])
      if(1 %in% routes){
        if(data.si$transm_sexual[i] == 1){
          cnts <- which(routes == 1)
          data.si$contacts[i] <- paste(c(rep(NA, 2-length(inf[cnts])), inf[cnts]), collapse = ',')
        }
      }
    }else if(length(inf) == 3){
      routes <- c(data.si$contact1_sexual[i], data.si$contact2_sexual[i], data.si$contact3_sexual[i])
      if(1 %in% routes){
        if(data.si$transm_sexual[i] == 1){
          cnts <- which(routes == 1)
          data.si$contacts[i] <- paste(c(rep(NA, 3-length(inf[cnts])), inf[cnts]), collapse = ',')
        }
      }
    }
  }
}

## Initial infector-infectee matrix + transmission routes
NCases <- length(unique(data.si$ID))
Case <- data.si$ID
infectors <- list()
routes <- list()
PossibleInfector <- matrix(nrow = NCases, ncol = 4)
TransmissionRoutes <- matrix(NA, nrow = NCases, ncol = 4)


for(i in 1:NCases){
  
  inf <- as.numeric(unlist(strsplit(data.si$contacts[i], ",")))
  
  sextrans <- numeric()
  # transmission routes: 1 = sexual, 2 = other (also if not specifically sexual reported) (data columns: 36,37,38,39)
  if(length(inf) > 0){
    for(c in 1:length(inf)){
      sex <- ifelse(data.si[i, 35 + c] == 1 & !is.na(data.si[i, 35 + c]), 1, 2)
      sextrans <- c(sextrans, sex)
    }
  }
  sextrans <- sextrans[!is.na(inf)]
  
  infectors[[i]] <- data.si$ID[data.si$ID_orig %in% inf]
  if(length(infectors[[i]]) > 0){
    PossibleInfector[i, 1:length(infectors[[i]])] <- infectors[[i]]
  }else{
    PossibleInfector[i, 1] <- 0 # no contacts among included cases
  }
  
  routes[[i]] <- sextrans
  if(length(routes[[i]]) > 0){
    TransmissionRoutes[i, 1:length(routes[[i]])] <- routes[[i]]
  }
  
}

data.si$contact_list <- lapply(infectors, function(x){
  if(length(x) == 0){
    return(NA)
  }else{
    return(as.numeric(x))
  }
})

data.si$transmission_list <- lapply(routes, function(x){
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

## Sample networks
source('/lustre1/scratch/326/vsc32693/simNet/fun_network.R')

num.nets <- 10000

library(foreach)
library(doParallel)

cl <- makeCluster(96)
registerDoParallel(cl)

foreach(t = 1:num.nets, .packages = c('igraph'), .errorhandling = 'remove') %dopar% {
  
  net <- setup_network(case.ids = data.si$ID,
                       cluster = NA,
                       contact.list = data.si$contact_list,
                       transm.list = data.si$transmission_list,
                       infector.mat = PossibleInfector,
                       transm.route = TransmissionRoutes,
                       symptom.onset = as.Date(data.si$symptom.onset),
                       helper.date = as.Date(data.si$date), # date of questionaire
                       min.si = min.si, # absolute value of the max. allowed negative serial interval
                       max.si = max.si
  )
  
  out.trees <- paste0("/lustre1/scratch/326/vsc32693/simNet/out/network_", t, ".csv")
  out.times <- paste0("/lustre1/scratch/326/vsc32693/simNet/out/onset_", t, ".csv")
  out.trans <- paste0("/lustre1/scratch/326/vsc32693/simNet/out/transm_", t, ".csv")
  
  write.table(c(t, net$network), file = out.trees, append = F, col.names = F, row.names = F, sep = ",", quote = F)
  write.table(c(t, net$onset.times), file = out.times, append = F, col.names = F, row.names = F, sep = ",", quote = F)
  write.table(c(t, net$transm.route), file = out.trans, append = F, col.names = F, row.names = F, sep = ",", quote = F)
  
}

stopCluster(cl)

