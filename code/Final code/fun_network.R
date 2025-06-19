
###########################
### FIND CYCLES
# Function to find cycles in network
FindCycles = function(g) {
  Cycles = NULL
  for(v1 in V(g)) {
    if(degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}

##############################
### CREATE NETWORK
# Function to set up a network based on symptom onset dates and contact info

setup_network <- function(case.ids, # ordered case IDs
                          cluster, # cluster IDs ordered by case ID
                          contact.list, 
                          transm.list,
                          infector.mat, # initial infector-infectee matrix
                          transm.route, # transmission routes for infector-intectee matrix
                          symptom.onset, # symptom onset dates by ordered case ID
                          helper.date, # date of notification / last present at collectivity - used to impute symptom onset dates if missing
                          min.si = 5, # absolute value of the max. allowed negative serial interval
                          max.si = 200 # suggestion of the max. positive serial interval (to assign index cases to clusters)
){
  
  NCases <- length(unique(case.ids))
  
  # Sample missing symptom onset dates
  minDate <- min(c(symptom.onset, helper.date), na.rm = T)-1
  Time <- as.numeric(symptom.onset - minDate) # symptom onset
  Time2 <- as.numeric(helper.date - minDate) # questionnaire
  timeLapsed <- max(as.numeric(helper.date - symptom.onset), na.rm = T) # max. time between questionnaire and symptom onset
  
  for(i in case.ids){
    
    if(is.na(Time[i])){
      
      # sample symptom onset times if no contact given
      min.date <- helper.date[i] - timeLapsed
      maxDate <- helper.date[i]
      Time[i] <- as.numeric(as.Date(round(runif(1, min.date, maxDate)), origin="1970-01-01") - minDate)
      
    }
    
  }
  
  ## Constrain negative serial interval to be max -X days
  # infection matrix based on observed pairs
  inf.mat <- matrix(nrow = NCases, ncol = NCases)
  trans.mat <- matrix(nrow = NCases, ncol = NCases)
  for(i in 1:NCases){
    for(j in 1:NCases){
      if(i == j){
        inf.mat[i, j] <- 0
      }else{
        inf.mat[i, j] <- ifelse(j %in% infector.mat[i, ], 1, 0)
      }
      
      if(inf.mat[i, j] == 1){
        id <- j
        p <- which(PossibleInfector[i, ] == j)
        trans.mat[i, j] <- TransmissionRoutes[i, p]
      }
    }
  }
  # reverse directionality also
  # i = infectee; j = infector
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
  
  # adapt infector-infectee matrix based on possible serial intervals
  time.mat <- matrix(nrow = NCases, ncol = NCases)
  for(i in 1:NCases){
    for(j in 1:NCases){
      time.mat[i, j] <- Time[i] - Time[j]
    }
  }
  # SI.mat <- inf.mat * time.mat
  for(i in 1:NCases){
    for(j in 1:NCases){
      inf.mat[i, j] <- ifelse((time.mat[i ,j] < -min.si || time.mat[i, j] > max.si), 0, inf.mat[i, j])
    }
  }
  
  # Update infector and transmission route matrix
  vi.list <- rep(list(NULL), NCases)
  transm.list <- rep(list(NULL), NCases)
  for(i in case.ids){
    for(j in case.ids){
      if(inf.mat[i, j] == 1){
        vi.list[[i]] <- c(vi.list[[i]], case.ids[j])
        transm.list[[i]] <- c(transm.list[[i]], trans.mat[i, j])
      }
    }
    # all cases possible index
    # vi.list[[i]] <- c(0, vi.list[[i]])
  }
  
  PossibleInfector2 <- matrix(nrow = NCases, ncol = max(lengths(vi.list)))
  for(i in 1:NCases){
    if(!is.null(vi.list[[i]])){
      PossibleInfector2[i, 1:lengths(vi.list)[i]] <- vi.list[[i]]
    }else{
      PossibleInfector2[i, 1] <- 0 # index case if no contacts
    }
  }
  NPossibleInfector2 <- rowSums(!is.na(PossibleInfector2))
  
  ## When 2 cases only contacted each other (i.e. pairs), randomly select index (both are possible based on SI constraint)
  for(i in case.ids){
    for(j in case.ids){
      if((NPossibleInfector2[i] == 1 & NPossibleInfector2[j] == 1 &
          PossibleInfector2[i, 1] == j & PossibleInfector2[j, 1] == i)[1]){
        PossibleInfector2[sample(c(i,j), 1), 1] <- 0
      }
    }
  }
  NPossibleInfector2    <- rowSums(!is.na(PossibleInfector2))
  
  ###-------------------------------------------------------------------------------------
  ### Additional checks to avoid cycles
  
  # ## Add cluster ID
  # c <- 1; cluster <- rep(NA, length(case.ids))
  # for(id in case.ids){
  #   
  #   contacts <- which(inf.mat[id, ] == 1 | inf.mat[,id] == 1)
  #   contacts <- contacts[!is.na(contacts)]
  #   contacts <- contacts[!is.na(cluster[contacts])]
  #   
  #   if(length(contacts) > 0){
  #     clust.id <- unique(cluster[contacts[!is.na(cluster[contacts])]])
  #     if(length(clust.id) > 0 && !is.na(unique(clust.id[1]))){
  #       cluster[id] <- clust.id[1]
  #       cluster[contacts] <- clust.id[1]
  #     }else{
  #       cluster[id] <- c
  #       c <- c + 1
  #     }
  #   }else{
  #     cluster[id] <- c
  #     c <- c + 1
  #   }
  # }
  # ids <- 1
  # while(length(unique(ids)) > 0){
  #   ids <- c()
  #   for(i in case.ids){
  #     for(j in case.ids){
  #       if(inf.mat[i,j] == 1 || inf.mat[j,i] == 1){
  #         if(cluster[i] != cluster[j]){
  #           ids <- unique(c(ids, i, j))
  #           cluster[i] <- cluster[j]
  #         }
  #       }
  #     }
  #   }
  #   # print(ids)
  # }
  
  ##---------------------------------------------------------------------
  ## Add cluster ID
  c <- 1; cluster <- NA
  for(id in 1:length(case.ids)){
    contacts <- PossibleInfector2[id, ]
    contacts <- contacts[!is.na(contacts)]
    # if((length(contacts) == 1 && contacts == 0)  || (length(contacts) == 0)){
    if(all(is.na(cluster[contacts]))){
      cluster[case.ids[case.ids == id]] <- c
      c <- c + 1
    }else{
      clust.id <- cluster[contacts[!is.na(contacts)]]
      if(length(unique(clust.id == 1))){
        cluster[case.ids[case.ids == id]] <- clust.id[1]
      }else{stop()}
    }
  }
  # cluster[c(167,168,170)] <- c
  # cluster[c(60,61,65,66,77)] <- c + 1
  # cluster[c(72,74,75)] <- c + 2
  # cluster[c(57,68,98,58,78,49,54,81,82)] <- 27
  cluster[c(49,54,57,58,68,78,81,82,98)] <- 28
  cluster[c(244,245,246,249)] <- 124
  # while(NA %in% cluster){
  #   for(i in case.ids){
  #     for(j in case.ids){
  #       if(inf.mat[i, j] == 1){
  #         if(is.na(cluster[i])){
  #           cluster[i] <- cluster[j]
  #         }
  #         if(is.na(cluster[j])){
  #           cluster[j] <- cluster[i]
  #         }
  #       }
  #     }
  #   }
  # }
  ids <- 1
  while(length(unique(ids)) > 0){
    ids <- c()
    for(i in case.ids){
      for(j in case.ids){
        if(inf.mat[i,j] == 1 || inf.mat[j,i] == 1){
          if(cluster[i] != cluster[j]){
            ids <- unique(c(ids, i, j))
            cluster[i] <- cluster[j]
          }
        }
      }
    }
    print(ids)
  }
  
  ##-------------------------------------------------------------------------
  
  # Resample network until there are no cycles
  PossibleInfector3 <- PossibleInfector2
  n.cycles <- 1
  n.try <- 1 # set limit on number of times to resample
  while(n.cycles > 0 || (NA %in% PossibleInfector3[,1])){
    
    PossibleInfector3 <- PossibleInfector2
    
    ## For clusters with unknown index, randomly select one
    for(c in unique(cluster)){
      if(length(which(cluster == c)) == 1){ # if only one case in the cluster
        PossibleInfector3[case.ids[cluster == c], ] <- c(0, rep(NA, length(PossibleInfector3[case.ids[cluster == c], ]) - 1))
      }else{
        if(!(0 %in% PossibleInfector3[case.ids[cluster == c], ])){
          if(length(case.ids[cluster == c]) > 1){
            temp <- c()
            for(i in case.ids[cluster == c]){
              for(j in case.ids[cluster == c]){
                if(j %in% PossibleInfector3[i, ] & i %in% PossibleInfector3[j, ] & (Time[i]-Time[j]) %in% seq(-min.si, max.si, 1)){
                  temp <- c(temp, i, j)
                }
              }
            }
            if(!is.null(temp)){
              ind <- sample(unique(temp), 1)
              PossibleInfector3[ind, ] <- c(0, rep(NA, length(PossibleInfector3[ind, ]) - 1))
            }
          }else{ # if only one case in the cluster
            id <- case.ids[cluster == c]
            PossibleInfector3[id, ] <- c(0, rep(NA, length(PossibleInfector3[id, ]) - 1))
          }
        }
      }
    }
    NPossibleInfector3    <- rowSums(!is.na(PossibleInfector3))
    
    ## When case i has only one contact, j, and j is definitely the infector (based on SI constraint) --> i should be removed from j's infector list
    # Shouldn't be possible because SI constraint already accounted for?
    for(i in case.ids){
      for(j in case.ids){
        if(NPossibleInfector3[i] == 1 & PossibleInfector3[i, 1] == j & (i %in% PossibleInfector3[j, ])){
          if((Time[i] - Time[j]) >= -min.si){ # sample infector if serial interval larger than negative constraint (presymptomatic transmission possible)
            ind <- sample(c(i,j), 1)
            if(ind == i){
              PossibleInfector3[i, 1] <- 0
            }else{
              temp <- PossibleInfector3[j, ]
              PossibleInfector3[j, ] <- c(temp[temp != case.ids[i]], NA)
            }
          }else{ # if t_i more than X days before t_j, i is infector (i.e. outside of serial interval constraint)
            temp <- PossibleInfector3[i, ]
            PossibleInfector3[i, ] <- c(0, rep(NA, length(PossibleInfector3[i, ])-1))
          }
        }
      }
    }
    NPossibleInfector3   <- rowSums(!is.na(PossibleInfector3))
    
    ###-----------------------------------------------------------------------------------------------
    
    infector.mat.final <- PossibleInfector3
    transm.mat.final <- matrix(nrow = NCases, ncol = NCases)
    for(i in 1:NCases){
      for(j in 1:NCases){
        if(i == j){
          transm.mat.final[i, j] <- NA
        }else{
          transm.mat.final[i, j] <- ifelse(j %in% infector.mat.final[i, ], trans.mat[i, j], NA)
        }
      }
    }
    
    IsNotContributorToLikelorg <- c(which(infector.mat.final[, 1] == 0)) 
    IsContributorToLikelorg <- case.ids[!case.ids%in%IsNotContributorToLikelorg]
    
    ### Sample network
    Network <- numeric(NCases)+0
    Update <- IsContributorToLikelorg
    
    # # Resample network until there are no cycles
    # n.cycles <- 1
    # n.try <- 1 # set limit on number of times to resample
    # ptm = proc.time()
    # while(n.cycles > 0){
    Draw <- round(runif(length(Update), min = 0.6, max = NPossibleInfector3[Update]+0.4)) # sample infector
    for(i in 1:length(IsContributorToLikelorg)){
      Network[Update[i]] <- infector.mat.final[Update[i], Draw[i]]
    }
    # AcceptedNetwork <Network
    Infector <- Network[IsContributorToLikelorg]
    Infectee <- case.ids[IsContributorToLikelorg]
    
    tree <- cbind(Infector, Infectee); tree <- tree[which(Infector != 0), ]
    g <- graph_from_edgelist(tree)
    n.cycles <- length(FindCycles(g))
    print(n.cycles)
    n.try <- n.try + 1
    if(n.try > 100000) stop()
  }
  # proc.time()-ptm
  
  IsContributorToLikel <- case.ids[Network != 0]
  IsNotContributorToLikel <- case.ids[!(case.ids %in% IsContributorToLikel)]
  
  TransRoute <- c()
  for(i in 1:NCases){
    TransRoute[i] <- ifelse(Network[i] != 0, transm.mat.final[i, Network[i]], NA)
  }
  
  # if(NA %in% TransRoute[Network != 0]){
  #   ret <- list(net = Network, trans = TransRoute, transmat = transm.mat.final, infect.mat = infector.mat.final, ninf = NPossibleInfector2)
  #   save(ret, file = 'error.RData')
  #   stop('no transmission route')
  # }
  
  return(list(network = Network,
              onset.times = Time,
              transm.route = TransRoute))
  
}
