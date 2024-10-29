
load('data/contact_data.RData')
data <- data.contact.clean

##---------------------------------------------------------------
## Individuals reporting only one contact (N = 86)

data.si <- (data[data$contact1_id != '' & data$contact2_id == '' & data$contact3_id == '', ])
data.si <- data.si[data.si$contacts != '', ]
all(nchar(data.si$contacts) <= 3 )
data.si$contacts <- as.numeric(data.si$contacts)

summary(data.si$symptom.onset)
data.si$symptom.onset.index <- NA
for(i in 1:dim(data.si)[1]){
  idcon <- data.si$contacts[i] # probable index case
  dt <- data$symptom.onset[data$ID == idcon] # take from full dataset
  data.si$symptom.onset.index[i] <- ifelse(is.na(dt), NA, as.Date(dt))
}
data.si$symptom.onset.index <- as.Date(data.si$symptom.onset.index)

data.si$serial.interval <- as.numeric(data.si$symptom.onset - data.si$symptom.onset.index)

summary(data.si$serial.interval)
hist(data.si$serial.interval)
hist(data.si$serial.interval[data.si$contact1_sexual == 1]); summary(data.si$serial.interval[data.si$contact1_sexual == 1])
hist(data.si$serial.interval[data.si$contact1_sexual == 2]); summary(data.si$serial.interval[data.si$contact1_sexual == 2])

# data.si$ID[data.si$serial.interval == 70]
# i <- 477
# idcon <- data.si$contacts[data.si$ID == 477]
# data$symptom.onset[data$ID == i]; data$symptom.onset[data$ID == idcon]
# data$days_since_last_contact1[data$ID == i]; data$date[data$ID == i]
# data$date[data$ID == i] - data$days_since_last_contact1[data$ID == i]

## Exclude outlier (last contact was 45 days after index symptom onset)
data.si <- data.si[data.si$serial.interval < 70, ]
data.si <- data.si[!is.na(data.si$serial.interval), ]
dim(data.si)
table(data.si$contact1_sexual)

hist(data.si$days_since_last_contact1)

# Check whether contact != ID
sum(data.si$contacts == data.si$ID)
data.si <- data.si[data.si$contacts != data.si$ID, ]

# Contacts reporting each other?
ddd <- data.si[data.si$contacts %in% data.si$ID, ]
ids <- c(ddd$contacts, ddd$ID)
View(data.si[data.si$ID %in% ids,][,c("ID","contacts")])

##-------------------------------------------------------------------------------------------
## Characteristics of the pairs with sexual transmission

table(data.si$contact1_sexual, data.si$n.lesion.genital)
table(data.si$n.lesion.genital[is.na(data.si$contact1_sexual)])
table(data.si$n.lesion.anal[is.na(data.si$contact1_sexual)])

# HH members
table(data.si$contact1_rel[data.si$contact1_sexual == 1])
# 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)

## plot matrix?
table(data.si$occupation_prim[data.si$contact1_sexual == 1], data.si$contact1_rel[data.si$contact1_sexual == 1])
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
table(data.si$contact1_sexual[data.si$serial.interval < 0])
sum(data.si$serial.interval[is.na(data.si$contact1_sexual)] < 0)
sum(data.si$serial.interval[!is.na(data.si$contact1_sexual)] < 0)


##-------------------------------------------------------------------------------------------
##
## EpiLPS to obtain nonparametric estimate of the serial interval
##
##-------------------------------------------------------------------------------------------

library(EpiLPS)
source('code/functions/estimSI_boot.R')

hist(data.si$serial.interval, breaks = 10)

xS <- data.frame(sL = data.si$serial.interval - 0.5, sR = data.si$serial.interval + 0.5)
set.seed(2022)
fitS <- estimSI_boot(x = xS)
round(fitS$estim, 2)

# Plot cdf
jpeg('results/SI_observed/plotCDF.jpeg', width = 30, height = 20, units = 'cm', res = 300)
bsLO <- min(xS$sL) - 0.5
bsRO <- max(xS$sR) + 0.5
sfineO <- seq(bsLO, bsRO, length = 100)
dsfineO <- sfineO[2] - sfineO[1]
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.si$serial.interval)
plot(ecdf_obs, main = '', xlab = '', ylab = '', xaxt = 'n', col = 'grey')
# Non-parametric CDF
lines(sfineO, sapply(sfineO, fitS$Fhat), type = 'l', col = 'darkblue',
      xlab = "Serial interval (days)",
      ylab = "Cumulative distribution function",
      xaxt = 'n', xlim = c(min(xS$sL)-1, max(xS$sR)+1),
      cex.lab = 2, cex.axis = 1.5, lwd = 3
)
title(main = 'All observed transmission pairs', cex.main = 2)
axis(1, at=seq(bsLO, bsRO, by = 2), cex.axis = 1.5)
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
dev.off()

##---------------------------------------------------------------------------------------
## Sexual vs non-sexual transmission: assuming NAs are sexual

data.sens1 <- data.si
data.sens1$contact1_sexual <- ifelse(is.na(data.sens1$contact1_sexual) | data.sens1$contact1_sexual == 1, 1, 2)
table(data.sens1$contact1_sexual)

## Non-sexual transmission
xNonSexual <- data.frame(sL = data.sens1$serial.interval[data.sens1$contact1_sexual == 2] - 0.5, 
                         sR = data.sens1$serial.interval[data.sens1$contact1_sexual == 2] + 0.5)
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

jpeg('results/SI_observed/plotCDFnonSexual.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.sens1$serial.interval[data.sens1$contact1_sexual == 2])
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
dev.off()

## Sexual transmission
xSexual <- data.frame(sL = data.sens1$serial.interval[data.sens1$contact1_sexual == 1] - 0.5, 
                      sR = data.sens1$serial.interval[data.sens1$contact1_sexual == 1] + 0.5)
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

jpeg('results/SI_observed/plotCDFSexual.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.sens1$serial.interval[data.sens1$contact1_sexual == 1])
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
dev.off()

##---------------------------------------------------------------------------------------
## Sexual vs non-sexual transmission: assuming NAs are non-sexual

data.sens2 <- data.si
data.sens2$contact1_sexual <- ifelse(is.na(data.sens2$contact1_sexual) | data.sens2$contact1_sexual == 2, 2, 1)
table(data.sens2$contact1_sexual)

## Non-sexual transmission
xNonSexual <- data.frame(sL = data.sens2$serial.interval[data.sens2$contact1_sexual == 2] - 0.5, 
                         sR = data.sens2$serial.interval[data.sens2$contact1_sexual == 2] + 0.5)
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

jpeg('results/SI_observed/plotCDFnonSexual_scen2.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.sens2$serial.interval[data.sens2$contact1_sexual == 2])
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
dev.off()

## Sexual transmission
xSexual <- data.frame(sL = data.sens2$serial.interval[data.sens2$contact1_sexual == 1] - 0.5, 
                      sR = data.sens2$serial.interval[data.sens2$contact1_sexual == 1] + 0.5)
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

jpeg('results/SI_observed/plotCDFSexual_scen2.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.sens2$serial.interval[data.sens2$contact1_sexual == 1])
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
dev.off()


##-------------------------------------------------------------------------------------------------------
## Sexual vs non-sexual transmission: assuming presence of genital/anal lesions indicate sexual contact

data.sens3 <- data.si
for(i in 1:dim(data.sens3)[1]){
  if(is.na(data.sens3$contact1_sexual[i])){
    if(is.na(data.sens3$n.lesion.anal[i]) || is.na(data.sens3$n.lesion.genital[i])){
      data.sens3$contact1_sexual[i] <- 2
    }else{ 
      if(data.sens3$n.lesion.anal[i] > 0 | data.sens3$n.lesion.genital[i] > 0){
        data.sens3$contact1_sexual[i] <- 1
      }else{
        data.sens3$contact1_sexual[i] <- 2
      }
    }
  }
}
table(data.sens3$contact1_sexual)

## Non-sexual transmission
xNonSexual <- data.frame(sL = data.sens3$serial.interval[data.sens3$contact1_sexual == 2] - 0.5, 
                         sR = data.sens3$serial.interval[data.sens3$contact1_sexual == 2] + 0.5)
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

jpeg('results/SI_observed/plotCDFnonSexual_scen3.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.sens3$serial.interval[data.sens3$contact1_sexual == 2])
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
dev.off()

## Sexual transmission
xSexual <- data.frame(sL = data.sens3$serial.interval[data.sens3$contact1_sexual == 1] - 0.5, 
                      sR = data.sens3$serial.interval[data.sens3$contact1_sexual == 1] + 0.5)
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

jpeg('results/SI_observed/plotCDFSexual_scen3.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,1))
# Empirical CDF
ecdf_obs <- ecdf(data.sens3$serial.interval[data.sens3$contact1_sexual == 1])
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
dev.off()






