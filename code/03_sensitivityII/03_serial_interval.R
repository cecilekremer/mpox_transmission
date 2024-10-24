
load('data/contact_data.RData')
data <- data.contact.clean

##--------------------------------------------------------------
## Data preparation

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
data.si <- data.si[,c(1:3,5,9:11,13:15,17,18,19,22,23,30:33,38:41,46:49,62:89,100,111:112)]

data.si$ID_orig <- data.si$ID # save original IDs
data.si$ID <- 1:dim(data.si)[1] # IDs need to be 1 to 171

table(data.si$contact1_rel)
# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)

##--------------------------------------------------------------
## Sexual transmission hypothesis

table(data.si$transm_sexual); table(data.si$transm_repeatedcontact)
table(data.si$transm_sexual, data.si$transm_repeatedcontact)
# table(data.si$contact1_sexual)

table(data.si$transm_sexual[data.si$contact2_id == '' & data.si$contact1_id != ''],
      data.si$contact1_sexual[data.si$contact2_id == '' & data.si$contact1_id != ''])
unique(data.si$transm_hyp[data.si$contact1_sexual == 2 & data.si$contact2_id == '' & data.si$contact1_id != '' & data.si$transm_sexual == 1])
# often hypothesis is sexual, but also possible from another repeated contact --> probably that one is then included in the contacts

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

##--------------------------------------------------------------
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

# ## Sample networks --> VSC
# source('code/03_sensitivityII/functions/fun_network.R')
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
#                        max.si = max.si
#   )
# 
#   trees[i, ] <- c(i, net$network)
#   onsets[i, ] <- c(i, net$onset.times)
#   transroutes[i, ] <- c(i, net$transm.route)
# 
#   # if(i%%100 == 0){
#     print(i)
#   # }
# 
# }
# # proc.time() - ptm
# save(trees, file = 'code/trees_211024.RData')
# save(onsets, file = 'code/onsets_211024.RData')
# save(transroutes, file = 'code/routes_211024.RData')

##-----------------------------------------------------------
## MCMC to estimate serial interval using normal distribution

# Change network when changing min and max si !!!
load('code/03_sensitivityII/trees_221024_sens2.RData')
load('code/03_sensitivityII/onsets_221024_sens2.RData')
load('code/03_sensitivityII/routes_221024_sens2.RData')

# # check if all networks are unique
# library(mgcv)
# dim(uniquecombs(trees[1:1000,c(-1)]))

trees <- trees[1:1000, ]
onsets <- onsets[1:1000, ]
transroutes <- transroutes[1:1000, ]

trees[,1] <- 1:dim(trees)[1]
onsets[,1] <- 1:dim(onsets)[1]
transroutes[,1] <- 1:dim(transroutes)[1]

source('code/03_sensitivityII/functions/fun_mcmc_si_strat.R')

nrun <- 5000000
burnin <- 0.4
thin <- 200
updatefr <- 2 # update network every other run
out <- estimate_si(case.ids = data.si$ID,
                   networks = trees,
                   onsets = onsets,
                   routes = transroutes,
                   max.si = max.si, ## Change network when changing min and max si !!!
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
jpeg('results/sensitivity2/plotConvergence.jpeg', width = 30, height = 30, units = 'cm', res = 300)
plot(chain1[,c(1,2,3,4)])
dev.off()
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

jpeg('results/sensitivity2/plotFitMostLikely.jpeg', width = 30, height = 20, units = 'cm', res = 300)
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
source('code/03_sensitivityII/functions/fun_network.R')
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
occupation <- c() # occupation_prim = 3 = sex worker' ; occupation_prim = 4 = mine worker
for(i in 1:length(case)){
  occupation[i] <- ifelse(data.si$occupation_prim[data.si$ID == i] == 3, 1, ifelse(data.si$occupation_prim[data.si$ID == i] == 4, 2, 3))  # 1 = sex worker, 2 = mine worker, 3 = other
}
occupation <- ifelse(is.na(occupation), 3, occupation)

table(age)/sum(table(age))
table(gender)/sum(table(gender))
table(occupation)/sum(table(occupation))

vertex.mat <- data.frame(case, age, gender, occupation)
vertex.mat$gender.occupation <- NA
for(i in 1:dim(vertex.mat)[1]){
  if(vertex.mat$gender[i] == 2){
    if(vertex.mat$occupation[i] == 1){
      vertex.mat$gender.occupation[i] <- 1 # female sex worker
    }else{
      vertex.mat$gender.occupation[i] <- 3 # female other
    }
  }else if(vertex.mat$gender[i] == 1){
    if(vertex.mat$occupation[i] == 2){
      vertex.mat$gender.occupation[i] <- 2 # male miner
    }else{
      vertex.mat$gender.occupation[i] <- 4 # male other
    }
  }
}

## Transmission matrix by age
head(vertex.mat)
head(edge.mat)

vertex.matrix <- vertex.mat
vertex.matrix$agenum <- NA
for(i in 1:dim(vertex.matrix)[1]){
  vertex.matrix$agenum[i] <- data.si$agenum[data.si$ID == i]
}
net.matrix <- as.data.frame(edge.mat)

net.matrix$source_age <- NA
net.matrix$case_age <- NA
for(i in 1:dim(net.matrix)[1]){
  net.matrix$source_age[i] <- vertex.matrix$agenum[vertex.matrix$case == net.matrix$Infector[i]]
  net.matrix$case_age[i] <- vertex.matrix$agenum[vertex.matrix$case == net.matrix$Infectee[i]]
}

net.matrix["agegroup_source"] <- cut(net.matrix$source_age, c(0, 2, 5, 10, 15, 20, 25, 30, 40, 50, 65),
                                     c("0-2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-40", "41-50", "51-65"),
                                     include.lowest = TRUE)
table(net.matrix$agegroup_source)
net.matrix["agegroup_case"] <- cut(net.matrix$case_age, c(0, 2, 5, 10, 15, 20, 25, 30, 40, 50, 65),
                                   c("0-2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-40", "41-50", "51-65"),
                                   include.lowest = TRUE)
table(net.matrix$agegroup_case)

age_pairs <- data.frame(xtabs(~agegroup_source + agegroup_case, net.matrix))
trans.mat = matrix(age_pairs$Freq, ncol=10, nrow=10) # row = index, col = contact
rownames(trans.mat) = c("0-2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-40", "41-50", "51-65")
colnames(trans.mat) = c("0-2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-40", "41-50", "51-65")
trans.mat = apply(t(trans.mat),2,rev) # row = contact, col = index
# index age groups
N.agegroup = apply(trans.mat,2,sum) # total number of index cases in age group i
N.mat = matrix(rep(N.agegroup,10), nrow=10, ncol=10, byrow=T)
mean.trans.mat = trans.mat / N.mat
apply(mean.trans.mat,2,sum) # should sum to 1

library(plot.matrix)
library(RColorBrewer)
jpeg("results/sensitivity2/agemat.jpeg", width=10, height=10, units="cm", res=300)
# par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(mean.trans.mat, xlab="Source case", ylab="Case", main="", col=colorRampPalette(brewer.pal(5, "GnBu")),
     fmt.key="%.0f", border=NA, asp=T, cex.axis=0.7, 
     axis.col=list(side=1, las=2), 
     axis.row=list(side=2, las=2),
     fmt.cell="%.2f", text.cell=list(cex=0.5), key=NULL,#  key=list(tick=F, at=c(0,2,4,6), labels=c(0,2,4,6)), 
     breaks=seq(0,1,0.1)
)
dev.off()

##---------------------------------------------------------------------
## Plot network

net <- graph.data.frame(edge.mat, vertex.mat, directed = T)

library(extrafont)
library(RColorBrewer)
colr2 <- c("#B569DB", "#2A9832", "#F39110")
V(net)$color <- colr2[V(net)$age]
edge.col <- c("red","black")
E(net)$color <- edge.col[E(net)$Type]
# v.shape <- c('csquare', 'crectangle', 'square', 'rectangle')
v.shape <- c('circle', 'square')
# v.shape.occ <- c('filled square', 'filled circle', NA)
# V(net)$shape <- v.shape[V(net)$gender.occupation]
V(net)$shape <- v.shape[V(net)$gender]
v.size <- c(5,5,3,3)
V(net)$size <- v.size[V(net)$gender.occupation]

jpeg('results/sensitivity2/plotNetwork.jpeg', width = 50, height = 50, units = 'cm', res = 300)
par(mfrow = c(1,1))
plot(net, vertex.size = V(net)$size, edge.color = E(net)$color, vertex.shape = V(net)$shape,
     vertex.label = '', layout=layout.fruchterman.reingold,
     vertex.label.cex = 1.2, edge.arrow.size = 0.5)
legend('topleft', c(">18y", "12-17y", "<12y", "male miner", "male other", "female sexworker", "female other", "sexual"),
       col = c("#B569DB", "#2A9832", "#F39110", 1, 1, 1, 1, "red"),
       lty = c(rep(NA, 7), 1),
       pch = c(rep(16, 3), 1, 1, 0, 0, NA),
       pt.cex = c(3,3,3,5,3,5,3,NA),
       cex = 1.2,
       y.intersp = 2)
dev.off()

##-------------------------------------------------------------------------------------------
##
## EpiLPS to obtain nonparametric estimate of non-sexual SI based on most likely network
## If bimodal -> further stratification needed?
##
##-------------------------------------------------------------------------------------------

library(EpiLPS)
source('code/03_sensitivityII/functions/estimSI_boot.R')

xS <- data.frame(sL = unlist(SerialInterval) - 0.5, sR = unlist(SerialInterval) + 0.5)
set.seed(2022)
fitS <- estimSI_boot(x = xS)
round(fitS$estim, 2)

# Plot cdf
jpeg('results/sensitivity2/plotCDF.jpeg', width = 30, height = 20, units = 'cm', res = 300)
bsLO <- min(xS$sL) - 0.5
bsRO <- max(xS$sR) + 0.5
sfineO <- seq(bsLO, bsRO, length = 100)
dsfineO <- sfineO[2] - sfineO[1]
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(unlist(SerialInterval))
plot(ecdf_obs, main = '', xlab = '', ylab = '', xaxt = 'n', col = 'grey')
# Non-parametric CDF
lines(sfineO, sapply(sfineO, fitS$Fhat), type = 'l', col = 'darkblue',
      xlab = "Serial interval (days)",
      ylab = "Cumulative distribution function",
      xaxt = 'n', xlim = c(min(xS$sL)-1, max(xS$sR)+1),
      cex.lab = 2, cex.axis = 1.5, lwd = 3
)
title(main = 'All transmission', cex.main = 2)
axis(1, at=seq(bsLO, bsRO, by = 2), cex.axis = 1.5)
# # Normal CDF
# lines(sfineO, sapply(sfineO, pnorm, mean = sumstats$quantiles[3, 3], sd = sumstats$quantiles[4, 3]),
#       col = '#54C8F0', lwd = 3, lty = 2)
# legend('topleft', c('Normal CDF','Non-parametric CDF'), col = c('#54C8F0','darkblue'), lwd = c(3,3), cex = 1.5)
# 95%CI for selected percentiles
low_95 <- fitS$estim$CI95p_l[3:7] 
up_95 <- fitS$estim$CI95p_r[3:7]
perc <- c(0.05,0.25,0.50,0.75,0.95)
perctxt <- paste0(perc * 100, "th percentile")
for(j in 1:length(up_95)){
  lines(x = c(low_95[j],up_95[j]), y = c(perc[j],perc[j]), type = "l",
        col = "darkblue", lwd = 3, lty = 1)
  lines(x = fitS$estim$Estim[3:7], perc, type = "p", pch = 16,
        col = "darkblue", cex = 2)
  text(x = up_95[j]+1.4, y = perc[j]-0.05, perctxt[j], offset = 1, cex = 1.5)
}
# lines(x = c(sumstats$quantiles[3, 1], sumstats$quantiles[3,5]), y = c(0.5, 0.5),
#       type = "l", col = "#54C8F0", lwd = 3, lty = 2)
# lines(x = sumstats$quantiles[3, 3], 0.5, type = "p", pch = 16, col = "#54C8F0", cex = 2)
dev.off()

##-----------------------------------------------------------------------------------
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

jpeg('results/sensitivity2/plotCDFnonSexual.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(SerialInterval[[2]])
plot(ecdf_obs, main = '', xlab = '', ylab = '', xaxt = 'n', col = 'grey')
# Non-parametric CDF
lines(sfineO, sapply(sfineO, fitNonSexual$Fhat), type = 'l', col = 'darkblue',
      xlab = "Serial interval (days)",
      ylab = "Cumulative distribution function",
      xaxt = 'n', xlim = c(min(xNonSexual$sL)-1, max(xNonSexual$sR)+1),
      cex.lab = 2, cex.axis = 1.5, lwd = 3
)
title(main = 'Non-sexual transmission', cex.main = 2)
axis(1, at=seq(bsLO, bsRO, by = 2), cex.axis = 1.5)
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

##-----------------------------------------------------------------------------------
## Sexual transmission
xSexual <- data.frame(sL = SerialInterval[[1]] - 0.5, sR = SerialInterval[[1]] + 0.5)
set.seed(2022)
fitSexual <- estimSI_boot(x = xSexual)
round(fitSexual$estim, 2)

# Plot cdf
bsLO <- min(xSexual$sL) - 0.5
bsRO <- max(xSexual$sR) + 0.5
sfineO <- seq(bsLO, bsRO, length = 100)
dsfineO <- sfineO[2] - sfineO[1]
BO <- nrow(fitSexual$bootsamples)
sboot_cdfO <- matrix(0, nrow = BO, ncol = 100)
for(b in 1:BO){
  fsbootO <- histosmooth(fitSexual$bootsamples[b,], xl = bsLO, xr = bsRO, K = 12)
  # fsbootO <- histosmooth(fitSexual$bootsamples[b,], K = 12)
  sboot_densO <- sapply(sfineO, fsbootO$fdens)
  sboot_densO <- sboot_densO/sum(sboot_densO * dsfineO)
  sboot_cdfO[b,] <- cumsum(sboot_densO * dsfineO)
  print(b)
}

jpeg('results/sensitivity2/plotCDFSexual.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(SerialInterval[[1]])
plot(ecdf_obs, main = '', xlab = '', ylab = '', xaxt = 'n', col = 'grey')
# Non-parametric CDF
lines(sfineO, sapply(sfineO, fitSexual$Fhat), type = 'l', col = 'darkblue',
      xlab = "Serial interval (days)",
      ylab = "Cumulative distribution function",
      xaxt = 'n', xlim = c(min(xSexual$sL)-1, max(xSexual$sR)+1),
      cex.lab = 2, cex.axis = 1.5, lwd = 3
)
title(main = 'Sexual transmission', cex.main = 2)
axis(1, at=seq(bsLO, bsRO, by = 2), cex.axis = 1.5)
# Normal CDF
lines(sfineO, sapply(sfineO, pnorm, mean = sumstats$quantiles[1, 3], sd = sumstats$quantiles[1, 3]),
      col = '#54C8F0', lwd = 3, lty = 2)
legend('topleft', c('Normal CDF','Non-parametric CDF'), col = c('#54C8F0','darkblue'), lwd = c(3,3), cex = 1.5)
# 95%CI for selected percentiles
low_95 <- fitSexual$estim$CI95p_l[3:7] 
up_95 <- fitSexual$estim$CI95p_r[3:7]
perc <- c(0.05,0.25,0.50,0.75,0.95)
perctxt <- paste0(perc * 100, "th percentile")
for(j in 1:length(up_95)){
  lines(x = c(low_95[j],up_95[j]), y = c(perc[j],perc[j]), type = "l",
        col = "darkblue", lwd = 3, lty = 1)
  lines(x = fitSexual$estim$Estim[3:7], perc, type = "p", pch = 16,
        col = "darkblue", cex = 2)
  text(x = up_95[j]+1.4, y = perc[j]-0.05, perctxt[j], offset = 1, cex = 1.5)
}
lines(x = c(sumstats$quantiles[1, 1], sumstats$quantiles[1,5]), y = c(0.5, 0.5),
      type = "l", col = "#54C8F0", lwd = 3, lty = 2)
lines(x = sumstats$quantiles[1, 3], 0.5, type = "p", pch = 16, col = "#54C8F0", cex = 2)
dev.off()
