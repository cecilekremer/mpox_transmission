
data <- read.csv('data/rawData.csv', header = T, sep = ',')

## Occupation
data$occupation1 <- apply(data[, paste0("occupation___", 1:19)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[1])
  }
})
data$occupation1 <- factor(data$occupation1, levels = 1:19)
table(data$occupation1)
# some people have multiple occupations
data$occupation2 <- apply(data[, paste0("occupation___", 1:19)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[2])
  }
})
data$occupation2 <- factor(data$occupation2, levels = 1:19)
table(data$occupation2)
# other occupation
table(data$otr_occp)

## Hunting
table(data$part_chas) # participated in hunting session; 1 = yes, 2 = no, 3 = no response
data$hunting_species1 <- apply(data[, paste0("espac_attr___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[1])
  }
})
table(data$hunting_species1) # 1 = monkey, 3 = rat, 9 = other
data$hunting_species2 <- apply(data[, paste0("espac_attr___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[2])
  }
})
table(data$hunting_species2) # 2 = squirrel
table(data$otr_espatrap)

table(data$consronger) # rodents caught or consumed; 1 = yes, 2 = no, 3 = no response
data$consronger_spec1 <- apply(data[, paste0("kelespec___", 1:5)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[1])
  }
})
table(data$consronger_spec1) # 1 = rat; 2 = Gambian rat; 3 = squirrel; 4 = other
data$consronger_spec2 <- apply(data[, paste0("kelespec___", 1:5)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[2])
  }
})
table(data$consronger_spec2) # 2 = Gambian rat; 3 = squirrel; 4 = other
table(data$otr_espattra)

## Bushmeat consumption
table(data$viandbrous) # 1 = never, 2 = <2x per week, 3 = >2x per week, 4 = no response
data$viandbrous_spec1 <- apply(data[, paste0("espconsm___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[1])
  }
})
table(data$viandbrous_spec1) # 1 = monkey, 2 = squirrel, 4 = bat, 5 = pangolin, 8 = antilope
data$viandbrous_spec2 <- apply(data[, paste0("espconsm___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[2])
  }
})
table(data$viandbrous_spec2)# 3 = rat, 9 = other
data$viandbrous_spec3 <- apply(data[, paste0("espconsm___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[3])
  }
})
table(data$viandbrous_spec3) # 7 = gazelle
data$viandbrous_spec4 <- apply(data[, paste0("espconsm___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[4])
  }
})
table(data$viandbrous_spec4)
table(data$otr_avicos)

## Touching or eating of animals that died naturally?
table(data$manipulmort) # 1 = yes, 2 = no, 3 = no response
data$manipulmort_spec1 <- apply(data[, paste0("espcmort___", 1:10)], 1, function(row){
  idx <- which(row == 1)
  # return first index if 1 exists, otherwise return NA
  if(length(idx) == 0){
    return(NA)
  }else{
    return(idx[1])
  }
})
table(data$manipulmort_spec1) # 9 = other
table(data$otr_vimr) # pork

##-----------------------------------------------------------------------------
## PCR confirmed cases
table(data$rsltlc) # 1 = positive, 2 = negative, 3 = not analysed yet

table(data$rsltlc[data$sexe==2 & (data$cathgr==1|data$cathgr==2)])
## Females +12y (N = 165?)
table(data$codeadultefemmeconfirmeecorrecte) # correct adult code, 1 = yes
table(data$codeenfantconfirmeecorrecte_2)
sum(!is.na(as.Date(data$dat_invesc)))

# table(data$rsltlc[data$sexe==1 & (data$cathgr==1|data$cathgr==2)])
# ## Males +12y (N = 171?)

##----------------------------------------------------------------------------------------------
## Longitudinal data (hospital follow-up)

sum(data$codeadulte != ""); sum(data$codeenfant != "") # follow-up IDs
length(unique(data$codeadulte)); length(unique(data$codeenfant)) 
length(unique(data$codeadulte[data$suivi_fait==1])); length(unique(data$codeenfant[data$suivi_fait==1])) 
# 417 patients with follow-up data, but almost none 'suivi_fait'?

table(data$redcap_repeat_instrument)

data$visit <- ifelse(data$redcap_repeat_instrument == '', 'baseline', data$redcap_repeat_instrument)
table(data$visit, data$redcap_repeat_instance)

length(unique(data$isac[data$suivi_fait==1])) # ID's of patients that follow-up during hospitalisation

table(data$isd_part) # hospital outcome: 1 = released, 2 = withdrawal of consent, 3 = deceased, 4 = LTFU

##------------------------------------------------------------------------------------------------
## Longitudinal data (regular follow-up)

table(data$jr_delavisite) # 1 = day 29; 2 = day 59


##------------------------------------------------------------------------------------------------
## Longitudinal data (pregnancy follow-up)




















