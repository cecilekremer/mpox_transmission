
estimate_si <- function(case.ids,
                        networks,
                        onsets,
                        routes,
                        max.si,
                        start.parms = c(1, 1, 1, 1, 1, 1), # mean1 = theta[1], sd1 = theta[2], mean2 = theta[3], sd2 = theta[4]
                        tuning.parms = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
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
  
  Route <- routes[t, c(-1)]
  AcceptedRoute <- Route
  
  IsContributorToLikel <- case.ids[Network != 0]
  IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
  
  ## Transmission by age of source case: 1 = 18+, 2 = 12-17, 3 = <12
  SerialInterval <- list()
  for(r in c(1,2,3)){
    SerialInterval[[r]] <- Time[IsContributorToLikel[Route[IsContributorToLikel] == r]] - Time[Network[IsContributorToLikel[Route[IsContributorToLikel] == r]]]
  }
  
  ## Likelihood
  likelihood <- function(){
    
    L1 <- c(); L2 <- c(); L3 <- c()

    # Age 18+
    for(i in 1:length(IsContributorToLikel[Route[IsContributorToLikel] == 1])){
      L1[i] <- dnorm(SerialInterval[[1]][i], mean = theta[1], sd = theta[2], log = FALSE)
    }
    
    # Age 12-17
    for(i in 1:length(IsContributorToLikel[Route[IsContributorToLikel] == 2])){
      L2[i] <- dnorm(SerialInterval[[2]][i], mean = theta[3], sd = theta[4], log = FALSE)
    }
    
    # Age <12
    for(i in 1:length(IsContributorToLikel[Route[IsContributorToLikel] == 3])){
      L3[i] <- dnorm(SerialInterval[[3]][i], mean = theta[5], sd = theta[6], log = FALSE)
    }
    
    return(sum(log(1e-50 + L1)) + sum(log(1e-50 + L2)) + sum(log(1e-50 + L3)))
    
  }
  
  ## Prior
  prior <- function(){
    
    mu1.prior <- dunif(theta[1], 0, max.si)
    sd1.prior <- dunif(theta[2], 0, 50)
    
    mu2.prior <- dunif(theta[3], 0, max.si)
    sd2.prior <- dunif(theta[4], 0, 50)
    
    mu3.prior <- dunif(theta[5], 0, max.si)
    sd3.prior <- dunif(theta[6], 0, 50)
    
    return(log(1e-50+mu1.prior) + log(1e-50+sd1.prior) + log(1e-50+mu2.prior) + log(1e-50+sd2.prior) + log(1e-50+mu3.prior) + log(1e-50+sd3.prior))
    
  }
  
  ## Posterior
  posterior <- function(){
    
    return(likelihood() + prior())
    
  }
  
  ##---------------------------------------------------------------------------##
  ##---------------------------- MCMC algorithm -------------------------------##
  
  ## Initial values
  AcceptedTheta <- theta <- c(start.parms[1], start.parms[2], start.parms[3], start.parms[4], start.parms[5], start.parms[6])
  
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
  SaveRoutes <- matrix(nrow = ncases, ncol = (NRuns - Burnin)/Thinning)
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
      Route <- routes[t, c(-1)]
      
      
      IsContributorToLikel <- case.ids[Network != 0]
      IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
      
      ## Sexual (1) vs non-sexual (2) transmission
      SerialInterval <- list()
      for(r in c(1,2,3)){
        SerialInterval[[r]] <- Time[IsContributorToLikel[Route[IsContributorToLikel] == r]] - Time[Network[IsContributorToLikel[Route[IsContributorToLikel] == r]]]
      }
            
    }
    
    if(b%%update.freq != 0){
      
      Network <- AcceptedNetwork
      Time <- AcceptedTime
      Route <- AcceptedRoute
      
      IsContributorToLikel <- case.ids[Network != 0]
      IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
      
      ## Sexual (1) vs non-sexual (2) transmission
      SerialInterval <- list()
      for(r in c(1,2,3)){
        SerialInterval[[r]] <- Time[IsContributorToLikel[Route[IsContributorToLikel] == r]] - Time[Network[IsContributorToLikel[Route[IsContributorToLikel] == r]]]
      }
      
      theta[1] <- runif(1, AcceptedTheta[1] - tuning[1], AcceptedTheta[1] + tuning[1])
      # if(theta[1] < 0) theta[1] <- AcceptedTheta[1]
      
      theta[2] <- runif(1, AcceptedTheta[2] - tuning[2], AcceptedTheta[2] + tuning[2])
      if(theta[2] < 0) theta[2] <- AcceptedTheta[2]
      
      theta[3] <- runif(1, AcceptedTheta[3] - tuning[3], AcceptedTheta[3] + tuning[3])
      # if(theta[1] < 0) theta[1] <- AcceptedTheta[1]
      
      theta[4] <- runif(1, AcceptedTheta[4] - tuning[4], AcceptedTheta[4] + tuning[4])
      if(theta[4] < 0) theta[4] <- AcceptedTheta[4]
      
      theta[5] <- runif(1, AcceptedTheta[5] - tuning[5], AcceptedTheta[5] + tuning[5])
      # if(theta[1] < 0) theta[1] <- AcceptedTheta[1]
      
      theta[6] <- runif(1, AcceptedTheta[6] - tuning[6], AcceptedTheta[6] + tuning[6])
      if(theta[6] < 0) theta[6] <- AcceptedTheta[6]
      
      
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
        AcceptedRoute <- Route
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
      SaveRoutes[, a] <- AcceptedRoute
      SaveP[a] <- AcceptedP
      SaveL[a] <- AcceptedL
      
      par(mfrow = c(3, 2))
      plot(Savetheta[,1], type = 'l', ylab = 'Mean SI 18+')
      plot(Savetheta[,2], type = 'l', ylab = 'SD SI 18+')
      plot(Savetheta[,3], type = 'l', ylab = 'Mean SI 12-17')
      plot(Savetheta[,4], type = 'l', ylab = 'SD SI 12-17')
      plot(Savetheta[,5], type = 'l', ylab = 'Mean SI <12')
      plot(Savetheta[,6], type = 'l', ylab = 'SD SI <12')
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
              transroutes = SaveRoutes,
              log_post = SaveP,
              log_lik = SaveL,
              accept.prob.theta = asd/NRuns,
              accept.prob.network = anetwork/NRuns
  ))
  
}