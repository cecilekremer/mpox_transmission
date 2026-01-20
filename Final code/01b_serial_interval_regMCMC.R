
# https://khayatrayen.github.io/MCMC.html

# https://medium.com/@tinonucera/bayesian-linear-regression-from-scratch-a-metropolis-hastings-implementation-63526857f191

# setwd('./Final code')

load('./Final code/data/contact_data_all.RData')
data <- data.contact.clean

data.si <- data[data$contacts != '', ]
table(data.si$contact1_sexual)
table(data.si$transm_sexual) # not specified which contact was sexual

# Make sure all contacts are included as cases
all.contact.ids <- unique(c(as.numeric(gsub("SCREEN", "", data.si$contact1_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact2_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact3_id)),
                            as.numeric(gsub("SCREEN", "", data.si$contact4_id))
))
table(all.contact.ids %in% data.si$ID)
all.case.ids <- unique(c(all.contact.ids, unique(data$ID[data$contacts != ''])))

data.si <- data[data$ID %in% all.case.ids, ]
# data.si <- data.si[,c(1:3,5,9:11,13:15,17,18,19,22,23,30:33,38:41,46:49,62:89,100)]
dim(data.si)

data.si$ID_orig <- data.si$ID # save original IDs
data.si$ID <- 1:dim(data.si)[1] # IDs need to be 1 to N

table(data.si$contact1_rel)
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
  # transmission routes: 1 = sexual, 2 = other (also if not specifically sexual reported) (data columns: 40,41,42,43)
  if(length(inf) > 0){
    for(c in 1:length(inf)){
      # sex <- ifelse(data.si[i, 82 + c] == 1 & !is.na(data.si[i, 82 + c]), 1, 2) # Scenario 2: NA = non-sexual
      # Scenario 1: NA = sexual if not family
      sex <- ifelse((is.na(data.si[i, 39 + c]) & data.si[i, 31 + c] %in% c(1,2)) | data.si[i, 39 + c] == 2, 2, 1)
      
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
      HH <- ifelse(data.si[i, 31 + c] == 1 & !is.na(data.si[i, 31 + c]), 1, 2)
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
load('./Final code/results/trees_100325_base.RData')
load('./Final code/results/onsets_100325_base.RData')
# load('code/serialinterval_final/routes_100325_base.RData')

# # check if all networks are unique
# library(mgcv)
# dim(uniquecombs(trees[,c(-1)]))

trees <- trees[!is.na(trees[,1]),]
onsets <- onsets[trees[,1],]
# transroutes <- transroutes[trees[,1], ]

# trees <- trees[1:1000, ]
# onsets <- onsets[1:1000, ]
# transroutes <- transroutes[1:1000, ]

trees[,1] <- 1:dim(trees)[1]
onsets[,1] <- 1:dim(onsets)[1]
# transroutes[,1] <- 1:dim(transroutes)[1]

source('./Final code/fun_reg_si.R')

nrun <- 5000000
burnin <- 0.4
thin <- 500
updatefr <- 2 # update network every other run
out <- estimate_si_reg(case.ids = data.si$ID,
                   networks = trees,
                   onsets = onsets,
                   # routes = transroutes,
                   data = data.si,
                   max.si = max.si, ## Change network when changing min and max si !!!
                   start.parms = rep(1, 5),
                   tuning.parms = rep(0.1, 5),
                   mcmc.runs = nrun,
                   burnin = burnin,
                   thin = thin,
                   update.freq = updatefr
)

# load('./results/SerialIntervalReg_Sex_HH_020425.RData')
# load('./results/SerialIntervalReg_age_HH_310325.RData') # with HH x sexual interaction

## Final result in manuscript (model incl. sexual, age group, and HH status)
load('./Final code/results/SerialIntervalReg_age_HH_280325.RData')
# ## Alternative model with age cutoff at 15y
# load('./results/SerialIntervalReg_Sex_HH_age_280725.RData')

length(unique(out$networkIDs))

## Check MCMC convergence
library(coda)
# thin <- 1
chain1 <- coda::mcmc(out$parms[seq(1, (nrun-(burnin*nrun))/thin), ])
# colnames(chain1) = c("Intercept", "b.sexual", "b.household", "b.interaction", "SD serial interval")
colnames(chain1) = c("Intercept", "b.sexual", "b.age", "b.household", "sigma")#, "b.hh.sexual")
# colnames(chain1) = c("Intercept", "b.sexual", "b.hh.sexual",  "b.household", "sigma")#, "b.hh.sexual")
# colnames(chain1) = c("Intercept", "b.sexual",  "b.household", "sigma")#, "b.hh.sexual")
autocorr.plot(chain1[,c(1:4)]) # see if thinning has to be increased
# jpeg('results/baseline/plotConvergenceReg.jpeg', width = 30, height = 30, units = 'cm', res = 300)
# par(mfrow = c(5, 2))
# plot(chain1[,c(1:4)])
# dev.off()
summary(chain1); sumstats <- summary(chain1)

# jpeg('./Final code/results/FigureS3.jpeg', width = 40, height = 30, units = 'cm', res = 300)
# par(mfrow = c(3, 2))
# plot(out$parms[,1], type = 'l', ylab = 'Intercept')
# plot(out$parms[,2], type = 'l', ylab = 'Sexual')
# plot(out$parms[,3], type = 'l', ylab = 'Age')
# plot(out$parms[,4], type = 'l', ylab = 'Household')
# # plot(out$parms[,3], type = 'l', ylab = 'HH x Sexual')
# plot(out$parms[,5], type = 'l', ylab = 'Sigma')
# dev.off()
library(ggplot2)
library(dplyr)
library(tidyr)
param_df <- as.data.frame(out$parms)
colnames(param_df) <- c("Intercept","Sexual","Age","Household","sigma")
param_df_long <- param_df %>%
  mutate(Iteration = row_number()) %>%
  pivot_longer(-Iteration, names_to = "Parameter", values_to = "Value")
ggplot(param_df_long, aes(x = Iteration, y = Value)) +
  geom_line() +
  facet_wrap(~Parameter, scales = 'free_y', ncol = 2) +
  theme_minimal(base_size = 18) +
  ylab("Parameter value") + xlab("Iteration")
# ggsave('./results/FigureS3.jpeg', width = 30, height = 20, units = 'cm', dpi = 300)

par(mfrow=c(2,1))
plot(out$log_post[seq(1, (nrun-(burnin*nrun))/thin)], type="l", ylab="Posterior")
plot(out$log_lik[seq(1, (nrun-(burnin*nrun))/thin)], type="l", ylab="Likelihood")

## Test significance of predictors: 0 in 95% CrI?
sumstats$quantiles[,c(1,3,5)]

## Posteriors
post.overall <- out$parms[seq(1, (nrun-(burnin*nrun))/thin), 1]
quantile(post.overall, c(0.025, 0.5, 0.975))
# Sexual transmission
post.sexual <- out$parms[seq(1, (nrun-(burnin*nrun))/thin), 1] + out$parms[seq(1, (nrun-(burnin*nrun))/thin), 2]
quantile(post.sexual, c(0.025, 0.5, 0.975))
# Household transmission
post.household <- out$parms[seq(1, (nrun-(burnin*nrun))/thin), 1] + out$parms[seq(1, (nrun-(burnin*nrun))/thin), 4]
quantile(post.household, c(0.025, 0.5, 0.975))
# Adult secondary case
post.household.sexual <- out$parms[seq(1, (nrun-(burnin*nrun))/thin), 1] + out$parms[seq(1, (nrun-(burnin*nrun))/thin), 3]
quantile(post.household.sexual, c(0.025, 0.5, 0.975))

## Network posterior probabilities
NCases <- length(unique(data.si$ID))
infector.mat <- matrix(0, NCases + 1, NCases)
for(i in 1:NCases){
  tab <- data.frame(table(out$network[i, seq(1, (nrun-(burnin*nrun))/thin)]) / sum(table(out$network[i, seq(1, (nrun-(burnin*nrun))/thin)])))
  for(k in 1:length(tab$Var1)){
    infector.mat[as.numeric(as.character(tab$Var1[k])) + 1, i] <- tab$Freq[k]
  }
  
}
# library(plot.matrix)
# library(RColorBrewer)
# rownames(infector.mat) <- c(0:NCases); colnames(infector.mat) <- c(1:NCases)
# par(mfrow = c(1,1))
# par(mar=c(5.1, 4.1, 4.1, 4.1))
# plot(infector.mat, xlab="Case", ylab="Infector", main="", fmt.key="%.1f", cex.axis=0.5,
#      col=colorRampPalette(brewer.pal(5, "Oranges")))

### Most likely network
tab <- data.frame(table(out$networkIDs)/sum(table(out$networkIDs)))
net.id <- as.numeric(as.character(tab$Var1[tab$Freq == max(tab$Freq)]))
# load('code/trees_211024.RData')
# load('code/onsets_211024.RData')
# load('code/routes_211024.RData')
Network <- trees[net.id, c(-1)]
Time <- onsets[net.id, c(-1)]

## Add transmission route
## Household / Sexual transmission yes/no
hh.trans <- c()
sex.trans <- c()
age <- c()
for(i in 1:NCases){
  age[i] <- ifelse(data.si$agecat[i] == 3, 0, 1) # 0 = child, 1 = adult
  inf.id <- Network[i]
  if(inf.id == 0){
    hh.trans[i] <- NA
    sex.trans[i] <- NA
  }else{
    inf.id.orig <- data.si$ID_orig[data.si$ID == inf.id]
    c <- which(as.numeric(unlist(strsplit(data.si$contacts[i], ","))) == inf.id.orig)
    if(length(c) == 0){
      id.orig <- data.si$ID_orig[data.si$ID == i]
      c <- which(as.numeric(unlist(strsplit(data.si$contacts[inf.id], ","))) == id.orig)
      hh.trans[i] <- ifelse(data.si[inf.id, 31 + c] == 1 & !is.na(data.si[inf.id, 31 + c]), 1, 0)
      # sex.trans[i] <- ifelse(data.si[inf.id, 82 + c] == 1 & !is.na(data.si[inf.id, 82 + c]), 1, 0) # Scenario 2
      sex.trans[i] <- ifelse((is.na(data.si[inf.id, 39 + c]) & data.si[inf.id, 31 + c] %in% c(1,2)) | (!is.na(data.si[inf.id, 39 + c]) & data.si[inf.id, 39 + c] == 2), 0, 1) # Scenario 1
    }else{
      hh.trans[i] <- ifelse(data.si[i, 31 + c] == 1 & !is.na(data.si[i, 31 + c]), 1, 0)
      # sex.trans[i] <- ifelse(data.si[i, 82 + c] == 1 & !is.na(data.si[i, 82 + c]), 1, 0) # Scenario 2
      sex.trans[i] <- ifelse((is.na(data.si[i, 39 + c]) & data.si[i, 31 + c] %in% c(1,2)) | (!is.na(data.si[i, 39 + c]) & data.si[i, 39 + c] == 2), 0, 1) # Scenario 1
    }
  }
}

table(hh.trans[Network!=0])
table(sex.trans[Network!=0])
table(age) # 0 = child, 1 = adult

##--------------------------------------------------
## Transmission pair characteristics
Infector <- Network
Infectee <- 1:length(Network)
# Type <- Route
HH <- hh.trans
Sexual <- sex.trans

table(HH[Network!=0])
table(Sexual[Network!=0])

table(HH)
table(Sexual)
table(HH, Sexual)

dat <- data.frame(Infectee, Infector, HH, Sexual)
tree <- cbind(Infector, Infectee, HH, Sexual)
tree <- tree[which(Infector != 0), ]
library(igraph)
g <- graph_from_edgelist(tree[,c(1,2)])
source('./Final code/fun_network.R')
length(FindCycles(g)) == 0

edge.mat <- tree

## Infectees
case <- 1:length(Network)
gender <- c()
occupation <- c()
for(i in 1:length(case)){
  gender[i] <- ifelse(data.si$gender[i] == 1, 1, 0) # 1 = male, 0 = female
  occupation[i] <- ifelse(data.si$occupation_prim[i] == 3, 1,
                          ifelse(data.si$occupation_prim[i] == 4, 2, 3)) # 1 = SW, 2 = mine worker, 3 = other
}
table(gender)
table(occupation)
occupation <- ifelse(is.na(occupation), 3, occupation)

# Infectee characteristics
table(age)/sum(table(age))
table(gender)/sum(table(gender))
table(occupation)/sum(table(occupation))

vertex.mat <- data.frame(case, age, gender, occupation)
vertex.mat$gender.occupation <- NA
for(i in 1:dim(vertex.mat)[1]){
  if(vertex.mat$gender[i] == 0){
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
table(vertex.mat$gender.occupation)

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
# jpeg("results/agemat.jpeg", width=10, height=10, units="cm", res=300)
# # par(mar=c(5.1, 4.1, 4.1, 4.1))
# plot(mean.trans.mat, xlab="Source case", ylab="Case", main="", col=colorRampPalette(brewer.pal(5, "GnBu")),
#      fmt.key="%.0f", border=NA, asp=T, cex.axis=0.7, 
#      axis.col=list(side=1, las=2), 
#      axis.row=list(side=2, las=2),
#      fmt.cell="%.2f", text.cell=list(cex=0.5), key=NULL,#  key=list(tick=F, at=c(0,2,4,6), labels=c(0,2,4,6)), 
#      breaks=seq(0,1,0.1)
# )
# dev.off()

## Transmission pair characteristics
edgedat <- as.data.frame(edge.mat)
vertexdat <- as.data.frame(vertex.mat)
inf.ids <- unique(edgedat$Infector); length(inf.ids)
table(vertexdat$age[vertexdat$case %in% inf.ids])
table(vertexdat$occupation[vertexdat$case %in% inf.ids])
table(vertexdat$occupation[vertexdat$case %in% inf.ids], vertexdat$age[vertexdat$case %in% inf.ids])
table(vertexdat$gender[vertexdat$case %in% inf.ids])
table(edgedat$Sexual)
table(edgedat$HH)
table(edgedat$HH, edgedat$Sexual)

##------------------------------------------------------------
## Plot network

# vertex.mat$age <- ifelse(data.si$agecat == 1, 1, 
#                          ifelse(data.si$agecat == 2, 2, 3)) # 1 = adult, 2 = 12-17, 3 = <12
# edge.mat[,3] <- ifelse(edge.mat[,3] == 1, 1, 2) # 1 = HH
# edge.mat[,4] <- ifelse(edge.mat[,4] == 1, 1, 2) # 1 = Sexual
# vertex.mat$gender <- ifelse(data.si$gender == 1, 1, 2) # 1 = male
# 
# edge.mat <- cbind(edge.mat, rep(NA, dim(edge.mat)[1]))
# edge.mat[,5] <- ifelse(edge.mat[,3] == 1 & edge.mat[,4] == 2, '1', # HH non-sexual
#                              ifelse(edge.mat[,3] == 1 & edge.mat[,4] == 1, 2, # HH sexual
#                                     ifelse(edge.mat[,3] == 2 & edge.mat[,4] == 1, 3, # non-HH sexual
#                                            ifelse(edge.mat[,3] == 2 & edge.mat[,4] == 2, 4, NA)))) # non-HH non-sexual
# colnames(edge.mat)[5] <- "HH_sexual"
# edge.mat[,5] <- as.numeric(edge.mat[,5])
# 
# net <- graph.data.frame(edge.mat, vertex.mat, directed = T)
# 
# library(extrafont)
# library(RColorBrewer)
# colr2 <- c("#B569DB", "#2A9832", "#F39110")
# V(net)$color <- colr2[V(net)$age]
# edge.col <- c("red","black","green","blue")
# E(net)$color <- edge.col[as.numeric(E(net)$HH_sexual)]
# # v.shape <- c('csquare', 'crectangle', 'square', 'rectangle')
# v.shape <- c('circle', 'square')
# # v.shape.occ <- c('filled square', 'filled circle', NA)
# # V(net)$shape <- v.shape[V(net)$gender.occupation]
# V(net)$shape <- v.shape[V(net)$gender]
# v.size <- c(5,5,3,3)
# V(net)$size <- v.size[V(net)$gender.occupation]
# # edge.lty <- c(2,1)
# # E(net)$lty <- edge.lty[E(net)$HH]
# 
# # jpeg('results/Figure1_final.jpeg', width = 180, height = 180, units = 'cm', res = 300)
# pdf('results/Figure1_final.pdf')#, res = 300)
# # par(mar = c(5,5,5,15))
# layout(matrix(c(2,1), nrow = 2), heights = c(5, 20))
# # net.layout <- layout_with_kk(net)
# plot(net, 
#      # layout = net.layout, 
#      vertex.size = V(net)$size, edge.color = E(net)$color, vertex.shape = V(net)$shape, #edge.lty = E(net)$lty,
#      vertex.label = '', layout=layout.fruchterman.reingold,
#      vertex.label.cex = 2, edge.arrow.size = 1.5, edge.width = 3)
# plot.new()
# legend('top', c(">18y", "12-17y", "<12y", "male miner", "male other", "female sexworker", "female other", "HH non-sexual", "HH sexual", "non-HH sexual","non-HH non-sexual"),
#        col = c("#B569DB", "#2A9832", "#F39110", rep('black', 4), "red","black","green","blue"),
#        lty = c(rep(NA, 7), 1,1,1,1),
#        pch = c(rep(19, 3), 19,19, 15,15, NA,NA,NA,NA),
#        pt.cex = c(10,10,10,15,10,15,10,NA,NA,NA,NA),
#        lwd = c(rep(NA,7), 10, 10,10,10),
#        cex = 9, ncol = 3,
#        # inset = c(-0.4, 0), xpd = NA,
#        y.intersp = 1, x.intersp = 3)
# dev.off()
# 
# jpeg('results/plotNetworkNoLegend.jpeg', width = 180, height = 180, units = 'cm', res = 300)
# par(mfrow = c(1,1))
# plot(net, vertex.size = V(net)$size, edge.color = E(net)$color, vertex.shape = V(net)$shape, #edge.lty = E(net)$lty,
#      vertex.label = '', layout=layout.fruchterman.reingold,
#      vertex.label.cex = 1.2, edge.arrow.size = 1)
# dev.off()   

# ---------- FINAL FIGURE----------
# convert cm to inches:
cm2in <- function(x) x / 2.54

# Desired figure size in cm:
fig_w_cm <- 18
fig_h_cm <- 18

# Use cairo_pdf to better embed fonts
cairo_pdf("./Final code/results/Figure1_final.pdf",
          width  = cm2in(fig_w_cm),
          height = cm2in(fig_h_cm),
          family = "sans")   # change to "Arial" if required

# Use layout: 2 rows (network big, legend small)
# heights are relative: give network much more vertical space than the legend
layout(matrix(c(1,2), nrow = 2, byrow = TRUE), heights = c(20, 5))

## Panel 1: network plot
# minimal margins for the network
par(mar = c(0, 0, 0, 0))   # bottom, left, top, right (in lines)
# fix the layout reproducibly if needed:
set.seed(42)
# If you have a precomputed layout, use it; otherwise:
# net.layout <- layout_with_fr(net)  # or layout_fruchterman_reingold, layout_with_kk, etc.
# plot:
plot(net,
     layout = layout.fruchterman.reingold,   # or net.layout
     vertex.size = V(net)$size,
     vertex.color = V(net)$color,
     vertex.shape = V(net)$shape,
     vertex.label = NA,        # no labels on main panel for clarity
     edge.color = E(net)$color,
     edge.width = ifelse(is.null(E(net)$width), 0.5, E(net)$width),
     edge.arrow.size = 0.2)    # scale down arrow size for publication

## Panel 2: legend panel
# small margins but allow drawing anywhere inside panel
par(mar = c(1, 0.5, 0, 0.5), xpd = NA)  # leave a little space
plot.new()

legend(x = "center", 
       legend = c(">18y", "12-17y", "<12y", 
                  "male miner", "male other", "female sexworker", "female other",
                  "HH non-sexual", "HH sexual", "non-HH sexual", "non-HH non-sexual"),
       col = c("#B569DB", "#2A9832", "#F39110", rep('black', 4), "red","black","green","blue"),
       pch = c(19,19,19, 19,19, 15,15, NA, NA, NA, NA),
       pt.cex = c(1.2, 1.2, 1.2, 1.5, 1.2, 1.5, 1.2, NA, NA, NA, NA),
       lty = c(rep(NA,7), 1,1,1,1),
       lwd = c(rep(NA,7), 2, 2, 2, 2),
       cex = 1.0,
       ncol = 3,
       bty = "n",
       y.intersp = 1.2, x.intersp = 1.2)

dev.off()
# ----------------------------------------

