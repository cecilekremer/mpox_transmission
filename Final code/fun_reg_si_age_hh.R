
estimate_si_reg <- function(case.ids,
                            networks,
                            onsets,
                            # routes,
                            data,
                            max.si,
                            start.parms = rep(1, 6), 
                            tuning.parms = rep(0.1, 6),
                            mcmc.runs = 1000000,
                            burnin = 0.3,
                            thin = 100,
                            update.freq = 1
){
  
  ncases <- length(case.ids)
  
  ## Sample initial network
  t <- sample(1:dim(networks)[1], 1)
  Network <- networks[t, c(-1)]
  Network.ID <- networks[t, 1]
  # net <- setup_network(case.ids = case.ids,
  #                      cluster = NA,
  #                      contact.list = contact.list,
  #                      infector.mat = infector.mat,
  #                      symptom.onset = symptom.onset,
  #                      helper.date = helper.date, # date of questionaire
  #                      min.si = min.si, # absolute value of the max. allowed negative serial interval
  #                      max.si = max.si)
  # Network <- net$network
  AcceptedNetwork <- Network
  
  Time <- onsets[t, c(-1)]
  # Time <- net$onset.times
  AcceptedTime <- Time
  
  # Route <- routes[t, c(-1)]
  # AcceptedRoute <- Route
  
  ## Household / Sexual transmission yes/no
  hh.trans <- c();
  sex.trans <- c()
  age <- c()
  for(i in 1:ncases){
    # age[i] <- ifelse(data$agecat[i] == 3, 0, 1) # 0 = child, 1 = adult
    age[i] <- ifelse(data$agenum[i] < 15, 0, 1)
    inf.id <- Network[i]
    if(inf.id == 0){
      hh.trans[i] <- NA
      sex.trans[i] <- NA
    }else{
      inf.id.orig <- data$ID_orig[data$ID == inf.id]
      c <- which(as.numeric(unlist(strsplit(data$contacts[i], ","))) == inf.id.orig)
      if(length(c) == 0){
        id.orig <- data$ID_orig[data$ID == i]
        c <- which(as.numeric(unlist(strsplit(data$contacts[inf.id], ","))) == id.orig)
        hh.trans[i] <- ifelse(data[inf.id, 58 + c] == 1 & !is.na(data[inf.id, 58 + c]), 1, 0)
        # sex.trans[i] <- ifelse(data[inf.id, 82 + c] == 1 & !is.na(data[inf.id, 82 + c]), 1, 0) # Scenario 2
        sex.trans[i] <- ifelse((is.na(data[inf.id, 82 + c]) & data[inf.id, 58 + c] %in% c(1,2)) | (!is.na(data[inf.id, 82 + c]) & data[inf.id, 82 + c] == 2), 0, 1) # Scenario 1
      }else{
        hh.trans[i] <- ifelse(data[i, 58 + c] == 1 & !is.na(data[i, 58 + c]), 1, 0)
        # sex.trans[i] <- ifelse(data[i, 82 + c] == 1 & !is.na(data[i, 82 + c]), 1, 0) # Scenario 2
        sex.trans[i] <- ifelse((is.na(data[i, 82 + c]) & data[i, 58 + c] %in% c(1,2)) | (!is.na(data[i, 82 + c]) & data[i, 82 + c] == 2), 0, 1) # Scenario 1
      }
    }
  }
  
  AcceptedHH <- hh.trans
  AcceptedSex <- sex.trans
  
  IsContributorToLikel <- case.ids[Network != 0]
  IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
  
  SerialInterval <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]]
  
  ## Likelihood
  likelihood <- function(){
    
    int <- theta[1] # intercept
    s <- theta[2] # main effect sexual
    h <- theta[3] # main effect age
    # sh <- theta[3] # interaction sexual-hh
    hh <- theta[4] # main effect household
    
    pred <- int + s*sex.trans[IsContributorToLikel] + hh*hh.trans[IsContributorToLikel] + h*age[IsContributorToLikel] #+ 
      # sh*sex.trans[IsContributorToLikel]*hh.trans[IsContributorToLikel]
    
    L <- c()
    for(i in 1:length(IsContributorToLikel)){
      # L[i] <- dnorm(SerialInterval[i], mean = pred[i], sd = theta[5], log = FALSE)
      L[i] <- dnorm(SerialInterval[i], mean = pred[i], sd = theta[5], log = FALSE)
    
    }
    
    return(sum(log(1e-50 + L)))
    
  }
  
  ## Prior
  prior <- function(){
    
    int.prior <- dnorm(theta[1], mean = 0, sd = 10)
    s.prior <- dnorm(theta[2], mean = 0, sd = 10)
    h.prior <- dnorm(theta[3], mean = 0, sd = 10)
    # sh.prior <- dnorm(theta[3], mean = 0, sd = 10)
    hh.prior <- dnorm(theta[4], mean = 0, sd = 10)
    # mu.prior <- dunif(theta[1], 0, max.si)
    # sd.prior <- dunif(theta[5], 0, 50) # Use normal prior on log(1/sd) ?? normal prior is the conjugate for Bayesian regression assuming a Normal distribution
    sd.prior <- dnorm(log(1/theta[5]), 0, 1)
    
    # return(log(1e-50+mu.prior) + log(1e-50+sd.prior))
    return(log(1e-50+int.prior) + log(1e-50+s.prior) + log(1e-50+hh.prior) + log(1e-50+h.prior) + 
             # log(1e-50+sh.prior) +
             log(1e-50+sd.prior))
    
  }
  
  ## Posterior
  posterior <- function(){
    
    return(likelihood() + prior())
    
  }
  
  ##---------------------------------------------------------------------------##
  ##---------------------------- MCMC algorithm -------------------------------##
  
  ## Initial values
  AcceptedTheta <- theta <- start.parms
  
  P <- posterior()
  AcceptedP <- P
  
  NRuns <- mcmc.runs
  Burnin <- round(burnin*NRuns)
  Thinning <- thin
  
  SaveP <- numeric()
  SaveL <- numeric()
  SaveNetwork <- matrix(nrow = ncases, ncol = (NRuns - Burnin)/Thinning)
  SaveNetworkID <- numeric()
  Net.ID <- numeric()
  SaveTimes <- matrix(nrow = ncases, ncol = (NRuns - Burnin)/Thinning)
  # SaveRoutes <- matrix(nrow = ncases, ncol = (NRuns - Burnin)/Thinning)
  Savetheta <- matrix(nrow = (NRuns - Burnin)/Thinning, ncol = length(theta))
  
  anetwork <- asd <- 0
  tuning <- tuning.parms
  a <- 0
  sd <- 0.5
  
  progressbar <- txtProgressBar(min = 0, max = NRuns, style = 3)
  
  for(b in 1:NRuns){
    
    # Update network
    if(b%%update.freq == 0){
      
      theta <- AcceptedTheta
      
      t <- sample(1:dim(networks)[1], 1)
      Network <- networks[t, c(-1)]
      Network.ID <- networks[t, 1]
      # net <- setup_network(case.ids = case.ids,
      #                      cluster = NA,
      #                      contact.list = contact.list,
      #                      infector.mat = infector.mat,
      #                      symptom.onset = symptom.onset,
      #                      helper.date = helper.date, # date of questionaire
      #                      min.si = min.si, # absolute value of the max. allowed negative serial interval
      #                      max.si = max.si)
      # Network <- net$network
      Time <- onsets[t, c(-1)]
      # Time <- net$onset.times
      # Route <- routes[t, c(-1)]
      ## Household / Sexual transmission yes/no
      hh.trans <- c();
      sex.trans <- c(); age <- c()
      for(i in 1:ncases){
        # age[i] <- ifelse(data$agecat[i] == 3, 0, 1) # 0 = child, 1 = adult
        age[i] <- ifelse(data$agenum[i] < 15, 0, 1)
        inf.id <- Network[i]
        if(inf.id == 0){
          hh.trans[i] <- NA
          sex.trans[i] <- NA
        }else{
          inf.id.orig <- data$ID_orig[data$ID == inf.id]
          c <- which(as.numeric(unlist(strsplit(data$contacts[i], ","))) == inf.id.orig)
          if(length(c) == 0){
            id.orig <- data$ID_orig[data$ID == i]
            c <- which(as.numeric(unlist(strsplit(data$contacts[inf.id], ","))) == id.orig)
            hh.trans[i] <- ifelse(data[inf.id, 58 + c] == 1 & !is.na(data[inf.id, 58 + c]), 1, 0)
            # sex.trans[i] <- ifelse(data[inf.id, 82 + c] == 1 & !is.na(data[inf.id, 82 + c]), 1, 0) # Scenario 2
            sex.trans[i] <- ifelse((is.na(data[inf.id, 82 + c]) & data[inf.id, 58 + c] %in% c(1,2)) | (!is.na(data[inf.id, 82 + c]) & data[inf.id, 82 + c] == 2), 0, 1) # Scenario 1
          }else{
            hh.trans[i] <- ifelse(data[i, 58 + c] == 1 & !is.na(data[i, 58 + c]), 1, 0)
            # sex.trans[i] <- ifelse(data[i, 82 + c] == 1 & !is.na(data[i, 82 + c]), 1, 0) # Scenario 2
            sex.trans[i] <- ifelse((is.na(data[i, 82 + c]) & data[i, 58 + c] %in% c(1,2)) | (!is.na(data[i, 82 + c]) & data[i, 82 + c] == 2), 0, 1) # Scenario 1
          }
        }
      }
      
      
      IsContributorToLikel <- case.ids[Network != 0]
      IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
      
      SerialInterval <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]]
      
    }
    
    if(b%%update.freq != 0){
      
      Network <- AcceptedNetwork
      Time <- AcceptedTime
      # Route <- AcceptedRoute
      hh.trans <- AcceptedHH
      sex.trans <- AcceptedSex
      
      IsContributorToLikel <- case.ids[Network != 0]
      IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
      
      SerialInterval <- Time[IsContributorToLikel] - Time[Network[IsContributorToLikel]]
      
      # theta[1] <- runif(1, AcceptedTheta[1] - tuning[1], AcceptedTheta[1] + tuning[1])
      # # if(theta[1] < 0) theta[1] <- AcceptedTheta[1]
      # theta[2] <- runif(1, AcceptedTheta[2] - tuning[2], AcceptedTheta[2] + tuning[2])
      # theta[3] <- runif(1, AcceptedTheta[3] - tuning[3], AcceptedTheta[3] + tuning[3])
      # theta[4] <- runif(1, AcceptedTheta[4] - tuning[4], AcceptedTheta[4] + tuning[4])
      # # theta[6] <- runif(1, AcceptedTheta[6] - tuning[6], AcceptedTheta[6] + tuning[6])
      
      theta[1] <- rnorm(1, AcceptedTheta[1], tuning[1])
      theta[2] <- rnorm(1, AcceptedTheta[2], tuning[2])
      theta[3] <- rnorm(1, AcceptedTheta[3], tuning[3])
      theta[4] <- rnorm(1, AcceptedTheta[4], tuning[4])
      # theta[6] <- rnorm(1, AcceptedTheta[6], tuning[6])
      
      # Variance
      # theta[4] <- runif(1, AcceptedTheta[4] - tuning[4], AcceptedTheta[4] + tuning[4])
      # if(theta[4] <= 0) theta[4] <- AcceptedTheta[4]
      theta[5] <- runif(1, AcceptedTheta[5] - tuning[5], AcceptedTheta[5] + tuning[5])
      if(theta[5] <= 0) theta[5] <- AcceptedTheta[5]
      
      
    }
    
    P <- posterior()
    Lik <- likelihood()
    
    ## Accept/reject new parameter values / new network
    AcceptYN <- runif(1, min = 0, max = 1) <= exp(P - AcceptedP)
    
    if(AcceptYN){
      
      if(b%%update.freq == 0){
        anetwork <- anetwork + 1
        AcceptedNetwork <- Network
        AcceptedNetworkID <- Network.ID
        AcceptedTime <- Time
        # AcceptedRoute <- Route
        AcceptedHH <- hh.trans
        AcceptedSex <- sex.trans
      }
      if(b%%update.freq != 0){
        asd <- asd + 1
        AcceptedTheta <- theta
      }
      
      AcceptedP <- P
      AcceptedL <- Lik
      
    }
    
    ## Save output
    if(b%%thin == 0 && b > Burnin){
      
      a <- a + 1
      Savetheta[a, ] <- AcceptedTheta
      SaveNetwork[, a] <- AcceptedNetwork
      SaveNetworkID[a] <- AcceptedNetworkID
      SaveTimes[, a] <- AcceptedTime
      # SaveRoutes[, a] <- AcceptedRoute
      SaveP[a] <- AcceptedP
      SaveL[a] <- AcceptedL
      
      # par(mfrow = c(3, 2))
      # plot(Savetheta[,1], type = 'l', ylab = 'Intercept')
      # plot(Savetheta[,2], type = 'l', ylab = 'Sexual')
      # plot(Savetheta[,3], type = 'l', ylab = 'Age')
      # plot(Savetheta[,4], type = 'l', ylab = 'Interaction')
      # 
      # plot(Savetheta[,5], type = 'l', ylab = 'SD serial interval')
      # 
      # plot(SaveP, type = 'l', ylab = 'Posterior')
      
    }
    # 
    # print(asd / b)
    # print(anetwork / b)
    
    setTxtProgressBar(progressbar, b)
    
  }
  
  close(progressbar)
  
  return(list(parms = Savetheta,
              network = SaveNetwork,
              networkIDs = SaveNetworkID,
              times = SaveTimes,
              # transroutes = SaveRoutes,
              log_post = SaveP,
              log_lik = SaveL,
              accept.prob.theta = asd/NRuns,
              accept.prob.network = anetwork/NRuns
  ))
  
}