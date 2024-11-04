
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

## Cases that meet inclusion criteria
ids.incl <- data$isac[data$critcompl_incl == 1]
data <- data[data$isac %in% ids.incl, ]

##-----------------------------------------------------------------------------
## PCR confirmed cases
table(data$rsltlc) # 1 = positive, 2 = negative, 3 = not analysed yet; only positive cases are retained in longitudinal part of the study

table(data$rsltlc[data$sexe==2 & (data$cathgr==1|data$cathgr==2)])
## Females +12y (N = 165?)
table(data$codeadultefemmeconfirmeecorrecte) # correct adult code, 1 = yes
table(data$codeenfantconfirmeecorrecte_2)
sum(!is.na(as.Date(data$dat_invesc)))

# table(data$rsltlc[data$sexe==1 & (data$cathgr==1|data$cathgr==2)])
# ## Males +12y (N = 171?)

## IDs of confirmed cases at baseline visit
ids.confirmed <- data$isac[data$rsltlc == 1 & data$redcap_repeat_instrument == '']
data <- data[data$isac %in% ids.confirmed, ]
dim(data)

length(unique(data$isac))
table(data$sexe) # 1 = male, 2 = female
table(data$cathgr) # age category: 1 = adult (>18y), 2 = minor (12-17y), 3 = minor (<12y)

##----------------------------------------------------------------------------------------------
## Longitudinal data (hospital follow-up)

sum(data$codeadulte != ""); sum(data$codeenfant != "") # follow-up IDs
length(unique(data$codeadulte)); length(unique(data$codeenfant)) 
length(unique(data$codeadulte[data$suivi_fait==1])); length(unique(data$codeenfant[data$suivi_fait==1])) 
# 417 patients with follow-up data, but almost none 'suivi_fait'?

table(data$redcap_repeat_instrument)

data$visit <- ifelse(data$redcap_repeat_instrument == '', 'baseline', data$redcap_repeat_instrument)
table(data$visit, data$redcap_repeat_instance)

length(unique(data$isac[data$suivi_fait==1])) # ID's of patients that had follow-up during hospitalisation

table(data$isd_part) # hospital outcome: 1 = released, 2 = withdrawal of consent, 3 = deceased, 4 = LTFU

##------------------------------------------------------------------------------------------------
## Longitudinal data (regular follow-up)

table(data$jr_delavisite) # 1 = day 29; 2 = day 59


##------------------------------------------------------------------------------------------------
## Longitudinal data (pregnancy follow-up)







##------------------------------------------------------------------------------------------------
## Rename variables

data$agecat <- data$cathgr # 1 = >18y; 2 = 12-17y
data$ageyear <- data$ag_an # NA if agecat = 3 & under 2y old
data$agemonth <- data$ag_mois
data$ageweek <- data$ag_seme
data$gender <- data$sexe

data$agenum[data$redcap_repeat_instrument==''] <- ifelse(!is.na(data$ageyear[data$redcap_repeat_instrument=='']),
                                                         data$ageyear[data$redcap_repeat_instrument==''],
                                                         ifelse(!is.na(data$agemonth[data$redcap_repeat_instrument=='']),
                                                                data$agemonth[data$redcap_repeat_instrument=='']/12,
                                                                data$ageweek[data$redcap_repeat_instrument=='']/52))
hist(data$agenum[data$redcap_repeat_instrument==''])

data$HZ <- data$zs # health zone
data$air_HZ <- data$air_st # aire de sante

data$occupation_prim <- data$occupation1
data$occupation_sec <- data$occupation2

data$HH_size <- data$tail_men
data$HH_children <- data$nbr_inf18
data$HH_mpox <- data$nbr_menage # 483 patients reported that at least one other HH member has mpox

data$travel <- data$partvoyage # 1 = yes, 2 = no

## Information on sexual contacts
data$intercourse <- data$rapsex # intercourse in the last 3 weeks; 1 = yes
data$num_partners <- data$cbpart # number of different sexual partners in the last 3 weeks
data$intercourse_SW <- data$rapprofsex # intercourse with professional sex worker in the last 3 weeks; 1 = yes
data$condomuse <- data$prsvraport # condom use: 1 = (almost) never, 2 = sometimes, 3 = (almost) always

##-------------------------------------------------------------------------------------------------
## Contact investigation
# unique(data$hyp_com) # transmission hypothesis
# sum(grepl('sexuel', data$hyp_com, ignore.case = T) | grepl('sexual', data$hyp_com, ignore.case = T))
# sum(grepl('non sexuel', data$hyp_com, ignore.case = T))
# sum(grepl('contact direct répété', data$hyp_com, ignore.case = T))
data$transm_hyp <- data$hyp_com

data$contact_mpox <- data$pesattient # contact w/ (suspected) mpox case: 1 = yes, 2 = no, 3 = don't know
data$num_contact_mpox <- data$cbinpers # if contact, w/ how many different persons
## Information recorded for up to 4 contacts:
data$contact1_included <- data$ptatmpox # 1 = included in this study
data$contact2_included <- data$ptatmpox_2
data$contact3_included <- data$ptatmpox_3
data$contact4_included <- data$ptatmpox_4
# If contact included in study:
data$contact1_id <- data$screnprs
data$contact2_id <- data$screnprs_2
data$contact3_id <- data$screnprs_3
data$contact4_id <- data$screnprs_4
data$contact1_FUid <- ifelse(data$numetdlg == 1 & data$numetdlg != 3, data$numadult, data$numenft) # follow-up ID of contact
data$contact2_FUid <- ifelse(data$numetdlg_2 == 1 & data$numetdlg_2 != 3, data$numadult_2, data$numenft_2)
data$contact3_FUid <- ifelse(data$numetdlg_3 == 1 & data$numetdlg_3 != 3, data$numadult_3, data$numenft_3)
data$contact4_FUid <- ifelse(data$numetdlg_4 == 1 & data$numetdlg_4 != 3, data$numadult_4, data$numenft_4)
data$contact1_type <- data$typconct # 1 = single contact, 2 = multiple occasions / ongoing
data$contact2_type <- data$typconct_2
data$contact3_type <- data$typconct_3
data$contact4_type <- data$typconct_4
data$days_since_last_contact1 <- data$dernct
data$days_since_last_contact2 <- data$dernct_2
data$days_since_last_contact3 <- data$dernct_3
data$days_since_last_contact4 <- data$dernct_4

data$contact1_rel <- ifelse(grepl('voisin', data$otr_typrel, ignore.case = T), 9, data$typrltpatien) # 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)  
data$contact1_rel <- ifelse(grepl('epoux', data$otr_typrel, ignore.case = T), 1, data$contact1_rel)
data$contact1_rel <- ifelse(grepl('glise', data$otr_typrel, ignore.case = T), 10, data$contact1_rel)
data$contact1_rel <- ifelse(grepl('patient', data$otr_typrel, ignore.case = T), 11, data$contact1_rel)
data$contact1_rel <- ifelse(grepl('mari', data$otr_typrel, ignore.case = T), 1, data$contact1_rel)
data$contact1_rel <- ifelse(grepl('conjointe', data$otr_typrel, ignore.case = T), 1, data$contact1_rel)
data$contact1_rel <- ifelse(grepl('sexuel', data$otr_typrel, ignore.case = T), 12, data$contact1_rel)

data$contact2_rel <- ifelse(grepl('voisin', data$otr_typrel_2, ignore.case = T), 9, data$typrltpatien_2)
data$contact2_rel <- ifelse(grepl('patient', data$otr_typrel_2, ignore.case = T), 11, data$contact2_rel)
data$contact2_rel <- ifelse(grepl('client', data$otr_typrel_2, ignore.case = T), 13, data$contact2_rel)

data$contact3_rel <- ifelse(grepl('voisin', data$otr_typrel_3, ignore.case = T), 9, data$typrltpatien_3)
data$contact3_rel <- ifelse(grepl('patient', data$otr_typrel_3, ignore.case = T), 11, data$contact3_rel)

data$contact4_rel <- ifelse(grepl('patient', data$otr_typrel_4, ignore.case = T), 11, data$typrltpatien_4)

data$contact1_famtype <- data$ty_mbre # 1 = parent, 2 = child, 3 = sibling, 4 = other ('data$otr_mbfam1')
data$contact2_famtype <- data$ty_mbre_2
data$contact3_famtype <- data$ty_mbre_3
data$contact4_famtype <- data$ty_mbre_4

# zone de sante where contact took place (if other: 'data$ot_aisnt')
data$contact1_zs <- data$airsntpatien 
data$contact2_zs <- data$airsntpatien_2
data$contact3_zs <- data$airsntpatien_3
data$contact4_zs <- data$airsntpatien_4
data$contact1_bar <- data$nmbar # name of bar where contact took place
data$contact2_bar <- data$nmbar_2
data$contact3_bar <- data$nmbar_3
data$contact4_bar <- data$nmbar_4

# Details on contact duration
data$contact1_prolonged <- data$causprolong # >15min: 1 = yes, 2 = no, 3 = don't remember
data$contact2_prolonged <- data$causprolong_2
data$contact3_prolonged <- data$causprolong_3
data$contact4_prolonged <- data$causprolong_4
data$contact1_direct <- data$tochdirect # direct touch
data$contact2_direct <- data$tochdirect_2
data$contact3_direct <- data$tochdirect_3
data$contact4_direct <- data$tochdirect_4
data$contact1_sexual <- data$contsex # sexual contact
data$contact2_sexual <- data$contsex_2
data$contact3_sexual <- data$contsex_3
data$contact4_sexual <- data$contsex_4
data$contact1_care <- data$prisoin # regularly took care of mpox patient
data$contact2_care <- data$prisoin_2
data$contact3_care <- data$prisoin_3
data$contact4_care <- data$prisoin_4
data$contact1_food <- data$mngemasit # eaten from the same plate
data$contact2_food <- data$mngemasit_2
data$contact3_food <- data$mngemasit_3
data$contact4_food <- data$mngemasit_4
data$contact1_bed <- data$mmlit # shared the same bed
data$contact2_bed <- data$mmlit_2
data$contact3_bed <- data$mmlit_3
data$contact4_bed <- data$mmlit_4
data$contact1_clothes <- data$mmvetm # shared the same clothes
data$contact2_clothes <- data$mmvetm_2
data$contact3_clothes <- data$mmvetm_3
data$contact4_clothes <- data$mmvetm_4

## Vaccination
data$smallpox_vacc <- data$partvaccin # vaccinated against smallpox before 1980; 1 = yes, 2 = no, 3 = NA
data$mpox_vacc_recent <- data$rsenvacci # recent mpox vaccination (since 2024) --> none

## Previous infection
data$previous_mpox <- data$ftmpox # already had mpox --> N = 9
library(lubridate)
data$date_previous_mpox <- ifelse(data$temps == 1,
                                      as.Date(data$date) - weeks(data$tpjour), # X weeks ag0
                                      as.Date(data$date) - data$tpmois*30.437 # X months ago
)
data$already_included <- data$patinclu # 3 patients already included
data$already_included_id <- data$nmpatincl
data$already_included_id[data$already_included_id!='']
data$already_included_id_long <- ifelse(data$nmetdlgnal == 1, data$nmadltetlnal, data$nmanftetlnal)

## Other infections (also included: TBC, diabetis, kidney, liver, other)
data$HIV <- data$antmedsvt___1
data$ART <- ifelse(data$HIV == 1 & data$siserop == 1, 1, 0)
data$HIV_test_during_consult <- data$tstvih # none

## Symptoms
data$symptom.onset <- as.Date(data$dtenqt) - data$dpdsympt
summary(data$symptom.onset)
# info on which symptoms

## Clinical presentation: lesions on mucous membranes
data$oral_lesion <- data$orale # 1 = yes, 2 = no
data$tonsil_lesion <- data$tonsil
data$penis_lesion <- data$glpenis
data$vagina_lesion <- data$vulv_gen

## Diagnostic sample results
data$diagnostics_complete <- data$resultats_des_echantillons_diagnostique_complete # 0 = incomplete, 2 = complete
data$mpox_pcr <- data$rsltlc # 1 = positive, 2 = negative, 3 = not analysed yet
data$sample_lesion <- data$base_quels_echant___1
data$sample_oropharyngal <- data$base_quels_echant___2
# longitudinal ID assigned if positive PCR result
data$id_adult <- data$codeadulte
data$id_child <- data$codeenfant

## Hospital follow-up 
summary(data$dtvisit); table(data$visit); table(data$suivi_fait) 

## Extract probable transmission route
# unique(data$transm_hyp)

data$transm_sexual <- ifelse(
  (grepl('sexuel', data$transm_hyp, ignore.case = T) | grepl('sexual', data$transm_hyp, ignore.case = T)) 
  & !grepl('non sexuel', data$transm_hyp, ignore.case = T),
  1, 0)
table(data$transm_sexual)
data$transm_repeatedcontact <- ifelse(
  grepl('contact direct répété', data$transm_hyp, ignore.case = T),
  1, 0
)
table(data$transm_repeatedcontact)
table(data$transm_repeatedcontact, data$transm_sexual)

table(data$contact1_rel) # 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)  
table(data$contact1_rel, data$transm_sexual)
table(data$contact1_rel, data$transm_repeatedcontact)

##------------------------------------------------------------------------------------------------
## Save clean dataset

save(data, file = 'data/clean_data.RData')

##-----------------------------------------------------------------------------
## Longitudinal data

load('data/clean_data.RData')
names(data); dim(data)

## PCR / Ct values


## Follow-up day 29/59
table(data$patien_lesicuta, data$visit) # lesions still present? 1 = yes


## Pregnancy (N = 20, 10% of women) --> vertical transmission?
table(data$pti_encient3m) # pregnant: 1 = yes, 2 = no, 3 = don't know
table(data$rslr_test) # if not pregnant or unkown, pregnancy test is done (1 = pos, 2 = neg, 3 = not done)

table(data$encetactuel) # is the woman really pregnant -> N = 20
table(data$pt_etaienct) # if not really pregnant, have they been pregnant in the 3 weeks before mpox?

summary(as.Date(data$datenaissprevue))
sum(as.Date(data$datenaissprevue) < '2024-10-01', na.rm = T)
summary(as.Date(data$dtissu))
sum(as.Date(data$dtissu) < '2024-10-01', na.rm = T)

table(data$issu_gross) # pregnancy outcome: 1 = >=37w, 2 = prem (<37w), 3 = stillbirth (>20w), 4 = miscarriage (<20w), 5 = abortion
table(data$issu_gross[as.Date(data$dtissu) < '2024-10-01'])
table(data$vnbmenfconu) # number of fetuses
table(data$signnass) # signs of mpox infection at birth: 1 = yes, 2 = no, 3 = unkown
table(data$descrextfoet)

