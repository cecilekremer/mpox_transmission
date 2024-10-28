
load("G:/My Drive/UHasselt/Research/MPOX/mpoxDRC/results/WS_serialInterval_all.RData")

load('data/contact_data.RData')
data <- data.contact.clean

# Make sure all contacts are included as cases
all.contact.ids <- unique(c(as.numeric(gsub("SCREEN", "", data.si$contact1_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact2_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact3_id))
))
table(all.contact.ids %in% data.si$ID)
all.case.ids <- unique(c(all.contact.ids, unique(data$ID[data$contacts != ''])))

data.si <- data[data$ID %in% all.case.ids, ]
data.si <- data.si[,c(1:3,5,9:11,13:15,17,18,19,22,23,30:33,38:41,46:49,62:89,100,101:104,124:126)]

data.si$ID_orig <- data.si$ID # save original IDs
data.si$ID <- 1:dim(data.si)[1] # IDs need to be 1 to 171

# sexual transmission vs lesions in genital region
sum(data.si$contact1_sexual == 1 | data.si$contact2_sexual == 1 | data.si$contact3_sexual == 1, na.rm = T)
sum(data.si$n.lesion.genital > 0, na.rm = T)
data.si$genital_lesion <- ifelse(data.si$n.lesion.genital > 0, 1, 0)
table(data.si$genital_lesion)

##-------------------------------------------------------------------------------
## Create dataframe with infector-infectee and characteristics

## Most likely network

tab <- data.frame(table(out$networkIDs)/sum(table(out$networkIDs)))
net.id <- as.numeric(as.character(tab$Var1[tab$Freq == max(tab$Freq)]))
# load('code/trees_211024.RData')
# load('code/onsets_211024.RData')
# load('code/routes_211024.RData')
Network <- trees[net.id, c(-1)]
Time <- onsets[net.id, c(-1)]
Route <- transroutes[net.id, c(-1)]
table(Route)/sum(table(Route)) # 37.23% sexual transmission

case.ids <- data.si$ID
IsContributorToLikel <- case.ids[Network != 0]
IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
SerialInterval <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]]
hist(SerialInterval)
summary(SerialInterval)

## Transmission pair characteristics
Infector <- Network
Infectee <- 1:length(Network)
Type <- Route
dat <- data.frame(Infectee, Infector, Type)
tree <- cbind(Infector, Infectee, Type)
tree <- tree[which(Infector != 0), ]

tree <- as.data.frame(tree)
# Serial interval most likely network
tree$SerialInterval <- NA
for(i in 1:dim(tree)[1]){
  infector <- tree$Infector[i]
  infectee <- tree$Infectee[i]
  tree$SerialInterval[i] <- Time[infectee] - Time[infector]
}
hist(tree$SerialInterval); summary(tree$SerialInterval); summary(SerialInterval)

# Demographics
tree$age.infector <- NA; tree$age.infectee <- NA
tree$gender.infector <- NA; tree$gender.infectee <- NA
tree$occupation.infector <- NA; tree$occupation.infectee <- NA

for(i in 1:dim(tree)[1]){
  infector <- tree$Infector[i]
  infectee <- tree$Infectee[i]
  
  tree$age.infector[i] <- data.si$agenum[data.si$ID == infector]
  tree$age.infectee[i] <- data.si$agenum[data.si$ID == infectee]
  tree$gender.infector[i] <- data.si$gender[data.si$ID == infector]
  tree$gender.infectee[i] <- data.si$gender[data.si$ID == infectee]
  tree$occupation.infector[i] <- data.si$occupation_prim[data.si$ID == infector]
  tree$occupation.infectee[i] <- data.si$occupation_prim[data.si$ID == infectee]
}
tree$occupation.infectee <- ifelse(tree$occupation.infectee == 3, 1, ifelse(tree$occupation.infectee == 4, 2, 3))# 1 = sex worker, 2 = mine worker, 3 = other
tree$occupation.infector <- ifelse(tree$occupation.infector == 3, 1, ifelse(tree$occupation.infector == 4, 2, 3))# 1 = sex worker, 2 = mine worker, 3 = other

# Household status
table(data.si$contact1_rel) # 1 = HH member

NCases <- length(unique(data.si$ID))
Case <- data.si$ID
infectors <- list()
routes <- list()
PossibleInfector <- matrix(nrow = NCases, ncol = 4)
TransmissionRoutes <- matrix(NA, nrow = NCases, ncol = 4)
for(i in 1:NCases){
  inf <- as.numeric(unlist(strsplit(data.si$contacts[i], ",")))
  infectors[[i]] <- data.si$ID[data.si$ID_orig %in% inf]
  if(length(infectors[[i]]) > 0){
    PossibleInfector[i, 1:length(infectors[[i]])] <- infectors[[i]]
  }else{
    PossibleInfector[i, 1] <- 0 # no contacts among included cases
  }
  HHtrans <- numeric()
  # transmission routes: 1 = HH member, 2 = other (also if not specifically reported) (data columns: 24,25,26,27)
  if(length(inf) > 0){
    for(c in 1:length(inf)){
      HH <- ifelse(data.si[i, 23 + c] == 1 & !is.na(data.si[i, 23 + c]), 1, 2)
      HHtrans <- c(HHtrans, HH)
    }
  }
  routes[[i]] <- HHtrans
  if(length(routes[[i]]) > 0){
    TransmissionRoutes[i, 1:length(routes[[i]])] <- routes[[i]]
  }
}
inf.mat <- matrix(nrow = NCases, ncol = NCases)
trans.mat <- matrix(nrow = NCases, ncol = NCases)
for(i in 1:NCases){
  for(j in 1:NCases){
    if(i == j){
      inf.mat[i, j] <- 0
    }else{
      inf.mat[i, j] <- ifelse(j %in% PossibleInfector[i, ], 1, 0)
    }
    if(inf.mat[i, j] == 1){
      id <- j
      p <- which(PossibleInfector[i, ] == j)
      trans.mat[i, j] <- TransmissionRoutes[i, p]
    }
  }
}  
for(i in 1:NCases){
  for(j in 1:NCases){
    if(i != j){
      if(inf.mat[i, j] == 1){
        inf.mat[j, i] <- 1
        trans.mat[j, i] <- trans.mat[i, j]
      }
    }
  }
}

tree$HH.member <- NA
for(i in 1:dim(tree)[1]){
  infector <- tree$Infector[i]
  infectee <- tree$Infectee[i]
  if(!is.na(trans.mat[infector, infectee])){
    tree$HH.member[i] <- trans.mat[infector, infectee]
  }else{
    tree$HH.member[i] <- trans.mat[infectee, infector]
  }
}
table(tree$HH.member)

edge.mat <- tree
## Presence of genital lesions
edge.mat$genital_lesion <- NA
for(i in 1:dim(edge.mat)[1]){
  edge.mat$genital_lesion[i] <- ifelse(data.si$genital_lesion[data.si$ID == edge.mat$Infectee[i]] == 1,
                                       1, 0)
}
table(edge.mat$genital_lesion)
table(edge.mat$Type)

##-------------------------------------------------------------------------
## Regression model

edge.mat[,c(3,7:11)] <- lapply(edge.mat[,c(3,7:11)], factor)
edge.mat$SI.transformed <- edge.mat$SerialInterval + 6

mod <- glm(SI.transformed ~ Type*occupation.infector + age.infector + HH.member, family = Gamma(link='log'), data = edge.mat)
summary(mod)

summary(edge.mat$SerialInterval[edge.mat$Type == 1])
summary(edge.mat$SerialInterval[edge.mat$Type == 2])

summary(edge.mat$SerialInterval[edge.mat$genital_lesion == 1])
summary(edge.mat$SerialInterval[edge.mat$genital_lesion == 0])

##--------------------------------------------------------------------------
## Bayesian regression

library(mlbench)
library(rstanarm)
library(bayestestR)
library(bayesplot)
library(insight)
library(broom)

modB <- rstanarm::stan_glm(SerialInterval ~ Type + occupation.infector + age.infector + HH.member,
                           data = edge.mat,
                           family = gaussian(),
                           prior = NULL, #(for regression coefs, Default: normal); set to NULL for flat uniform prior
                           # prior_intercept = NULL,
                           # prior_aux = NULL, # auxiliary variables like std. error
                           algorithm = 'sampling', # MCMC
                           iter = 2000,
                           chains = 4, 
                           warmup = floor(2000/2), # default: burnin is half the number of iterations
                           seed = 123)
print(modB, digits = 3) # MAD_SD is the median absolute deviaiton
# Plot for each predictor
bayesplot::mcmc_dens(modB) #, pars = c("Type2"))
bayestestR::describe_posterior(modB)
# CI = highest density interval (HDI)
# pd = probability of direction: equivalent to p-value?
# ROPE = region of practical equivalence (i.e. no effect)
hdi(modB) # if 0 not in 95% HDI, the coefficient is significant
eti(modB)
# test significance using ROPE interval: if low % in ROPE, highly significant (i.e. probability of coefficient to be zero)

