
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
hist(data.si$serial.interval[data.si$contact1_sexual == 1]); summary(data.si$serial.interval[data.si$contact1_sexual == 1])
hist(data.si$serial.interval[data.si$contact1_sexual == 2]); summary(data.si$serial.interval[data.si$contact1_sexual == 2])

## Nonparametric estimation using EpiLPS



##--------------------------------------------------------------
## Reporting multiple contacts

data.si <- data[data$contacts != '', ]
table(data.si$contact1_sexual)
table(data.si$transm_sexual) # not specified which contact was sexual

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

## Initial infector-infectee matrix + transmission routes
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
  
  sextrans <- numeric()
  # transmission routes: 1 = sexual, 2 = other (also if not specifically sexual reported) (data columns: 36,37,38,39)
  if(length(inf) > 0){
    for(c in 1:length(inf)){
      sex <- ifelse(data.si[i, 35 + c] == 1 & !is.na(data.si[i, 35 + c]), 1, 2)
      sextrans <- c(sextrans, sex)
    }
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

rm(Case); rm(NCases); rm(i)
library(igraph)

# ## Sample networks
# source('code/functions/fun_network.R')
# 
# num.nets <- 1000
# trees <- matrix(NA, nrow = num.nets, ncol = length(unique(data.si$ID)) + 1)
# onsets <- matrix(NA, nrow = num.nets, ncol = length(unique(data.si$ID)) + 1)
# transroutes <- matrix(NA, nrow = num.nets, ncol = length(unique(data.si$ID)) + 1)
# 
# # ptm <- proc.time()
# for(i in 1:num.nets){
#   net <- setup_network(case.ids = data.si$ID,
#                        cluster = NA,
#                        contact.list = data.si$contact_list,
#                        transm.list = data.si$transmission_list,
#                        infector.mat = PossibleInfector,
#                        transm.route = TransmissionRoutes,
#                        symptom.onset = as.Date(data.si$symptom.onset),
#                        helper.date = as.Date(data.si$date), # date of questionaire
#                        min.si = min.si, # absolute value of the max. allowed negative serial interval
#                        max.si = 30
#   )
#   
#   trees[i, ] <- c(i, net$network)
#   onsets[i, ] <- c(i, net$onset.times)
#   transroutes[i, ] <- c(i, net$transm.route)
#   
#   if(i%%100 == 0){
#     print(i)
#   }
#   
# }
# # proc.time() - ptm
# save(trees, file = 'code/trees_211024.RData')
# save(onsets, file = 'code/onsets_211024.RData')
# save(transroutes, file = 'code/routes_211024.RData')

## MCMC to estimate serial interval using normal distribution
load('code/trees_211024.RData')
load('code/onsets_211024.RData')
load('code/routes_211024.RData')

source('code/functions/fun_mcmc_si_strat.R')

nrun <- 2000000
burnin <- 0.4
thin <- 100
updatefr <- 2 # update network every other run
out <- estimate_si(case.ids = data.si$ID,
                   networks = trees,
                   onsets = onsets,
                   routes = transroutes,
                   max.si = 30,
                   start.parms = c(1, 1, 1, 1), # mean1 = theta[1], sd1 = theta[2], mean2 = theta[3], sd2 = theta[4]
                   tuning.parms = c(0.2, 0.2, 0.2, 0.2),
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
colnames(chain1) = c("mean SI sexual","sd SI sexual","mean SI non-sexual","sd SI non-sexual")
autocorr.plot(chain1[,c(1,2,3,4)]) # see if thinning has to be increased
plot(chain1[,c(1,2,3,4)])
summary(chain1); sumstats <- summary(chain1)

par(mfrow=c(2,1))
plot(out$log_post[seq(1, (nrun-(burnin*nrun))/thin)], type="l", ylab="Posterior")
plot(out$log_lik[seq(1, (nrun-(burnin*nrun))/thin)], type="l", ylab="Likelihood")

## Network posterior probabilities
NCases <- length(unique(data.si$ID))
infector.mat <- matrix(0, NCases + 1, NCases)
for(i in 1:NCases){
  tab <- data.frame(table(out$network[i, seq(1, (nrun-(burnin*nrun))/thin)]) / sum(table(out$network[i, seq(1, (nrun-(burnin*nrun))/thin)])))
  for(k in 1:length(tab$Var1)){
    infector.mat[as.numeric(as.character(tab$Var1[k])) + 1, i] <- tab$Freq[k]
  }
  
}
library(plot.matrix)
library(RColorBrewer)
rownames(infector.mat) <- c(0:NCases); colnames(infector.mat) <- c(1:NCases)
par(mfrow = c(1,1))
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(infector.mat, xlab="Case", ylab="Infector", main="", fmt.key="%.1f", cex.axis=0.5,
     col=colorRampPalette(brewer.pal(5, "Oranges")))

##-----------------------------------------------------
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

## Plot estimated serial interval to observed 

# Sexual (1) vs non-sexual (2) transmission
case.ids <- data.si$ID
IsContributorToLikel <- case.ids[Network != 0]
IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
SerialInterval <- list()
for(r in c(1,2)){
  SerialInterval[[r]] <- Time[IsContributorToLikel[Route[IsContributorToLikel] == r]] - Time[Network[IsContributorToLikel[Route[IsContributorToLikel] == r]]]
}

jpeg('results/plotFitMostLikely.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(2, 1))
hist(SerialInterval[[1]], prob = T, xlim = c(-20, 40), breaks = 20, 
     main = paste0('Sexual transmission (n = ', length(SerialInterval[[1]]), ')'), xlab = 'Serial interval') 
curve(dnorm(x, mean = sumstats$quantiles[1, 3], sd = sumstats$quantiles[2, 3]), from = -20, to = 40,
      add = T, col = 2, lwd = 2)
hist(SerialInterval[[2]], prob = T, xlim = c(-20, 40), breaks = 20, 
     main = paste0('Non-sexual transmission (n = ', length(SerialInterval[[2]]), ')'), xlab = 'Serial interval') 
curve(dnorm(x, mean = sumstats$quantiles[3, 3], sd = sumstats$quantiles[4, 3]), from = -20, to = 40,
      add = T, col = 2, lwd = 2)
dev.off()

## Transmission pair characteristics
Infector <- Network
Infectee <- 1:length(Network)
Type <- Route
dat <- data.frame(Infectee, Infector, Type)
tree <- cbind(Infector, Infectee, Type)
tree <- tree[which(Infector != 0), ]
library(igraph)
g <- graph_from_edgelist(tree[,c(1,2)])
source('code/functions/fun_network.R')
length(FindCycles(g)) == 0

edge.mat <- tree

# case <- c(1:max(Infector))
case <- 1:length(Network)
age <- c() # 1 = >18y; 2 = 12-17y
gender <- c() # 1 = male
for(i in 1:length(case)){
  # age[i] <- floor(data.si$agenum[data.si$ID == i])
  age[i] <- data.si$agecat[data.si$ID == i]
  gender[i] <- data.si$gender[data.si$ID == i]
}

table(age)/sum(table(age))
table(gender)/sum(table(gender))

vertex.mat <- data.frame(case, age, gender)

net <- graph.data.frame(edge.mat, vertex.mat, directed = T)

## Plot network
library(extrafont)
library(RColorBrewer)
colr2 <- c("#B569DB", "#2A9832", "#F39110")
V(net)$color <- colr2[V(net)$age]
edge.col <- c("red","black")
E(net)$color <- edge.col[E(net)$Type]
v.shape <- c('circle', 'square')
V(net)$shape <- v.shape[V(net)$gender]

jpeg('results/plotNetwork.jpeg', width = 30, height = 30, units = 'cm', res = 300)
par(mfrow = c(1,1))
plot(net, vertex.size = 3, edge.color = E(net)$color, vertex.shape = V(net)$shape,
     vertex.label = '', layout=layout.fruchterman.reingold,
     vertex.label.cex = 1.2, edge.arrow.size = 0.5)
legend('topleft', c(">18y", "12-17y", "<12y", "male", "female", "sexual"),
       col = c("#B569DB", "#2A9832", "#F39110", 1, 1, "red"),
       lty = c(rep(NA, 5), 1),
       pch = c(rep(16, 3), 1, 0, NA),
       cex = 1.2)
dev.off()

##-------------------------------------------------------------------------------------------
##
## EpiLPS to obtain nonparametric estimate of non-sexual SI based on most likely network
## If bimodal -> further stratification needed?
##
##-------------------------------------------------------------------------------------------

library(EpiLPS)
source('code/functions/estimSI_boot.R')

## Non-sexual transmission
xNonSexual <- data.frame(sL = SerialInterval[[2]] - 0.5, sR = SerialInterval[[2]] + 0.5)
set.seed(2022)
fitNonSexual <- estimSI_boot(x = xNonSexual)
round(fitNonSexual$estim, 2)

# Plot cdf
bsLO <- min(xNonSexual$sL) - 0.5
bsRO <- max(xNonSexual$sR) + 0.5
sfineO <- seq(bsLO, bsRO, length = 100)
dsfineO <- sfineO[2] - sfineO[1]
BO <- nrow(fitNonSexual$bootsamples)
sboot_cdfO <- matrix(0, nrow = BO, ncol = 100)
for(b in 1:BO){
  fsbootO <- histosmooth(fitNonSexual$bootsamples[b,], xl = bsLO, xr = bsRO, K = 12)
  sboot_densO <- sapply(sfineO, fsbootO$fdens)
  sboot_densO <- sboot_densO/sum(sboot_densO * dsfineO)
  sboot_cdfO[b,] <- cumsum(sboot_densO * dsfineO)
  print(b)
}

jpeg('results/plotCDF.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Non-parametric CDF
plot(sfineO, sapply(sfineO, fitNonSexual$Fhat), type = 'l', col = 'darkblue',
     xlab = "Serial interval (days)",
     ylab = "Cumulative distribution function",
     xaxt = 'n', xlim = c(min(xNonSexual$sL)-1, max(xNonSexual$sR)+1),
     cex.lab = 2, cex.axis = 1.5, lwd = 3
)
title(main = 'Non-sexual transmission', cex.main = 2)
axis(1, at=seq(-3, 29, by = 2), cex.axis = 1.5)
# Normal CDF
lines(sfineO, sapply(sfineO, pnorm, mean = sumstats$quantiles[3, 3], sd = sumstats$quantiles[4, 3]),
      col = '#54C8F0', lwd = 3, lty = 2)
legend('topleft', c('Normal CDF','Non-parametric CDF'), col = c('#54C8F0','darkblue'), lwd = c(3,3), cex = 1.5)
# 95%CI for selected percentiles
low_95 <- fitNonSexual$estim$CI95p_l[3:7] 
up_95 <- fitNonSexual$estim$CI95p_r[3:7]
perc <- c(0.05,0.25,0.50,0.75,0.95)
perctxt <- paste0(perc * 100, "th percentile")
for(j in 1:length(up_95)){
  lines(x = c(low_95[j],up_95[j]), y = c(perc[j],perc[j]), type = "l",
        col = "darkblue", lwd = 3, lty = 1)
  lines(x = fitNonSexual$estim$Estim[3:7], perc, type = "p", pch = 16,
        col = "darkblue", cex = 2)
  text(x = up_95[j]+1.4, y = perc[j]-0.05, perctxt[j], offset = 1, cex = 1.5)
}
lines(x = c(sumstats$quantiles[3, 1], sumstats$quantiles[3,5]), y = c(0.5, 0.5),
      type = "l", col = "#54C8F0", lwd = 3, lty = 2)
lines(x = sumstats$quantiles[3, 3], 0.5, type = "p", pch = 16, col = "#54C8F0", cex = 2)
dev.off()
