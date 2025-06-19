
load('data/clean_data_100325.RData')

names(data); dim(data)

## Hospital outcome
data$hospital_outcome <- data$isd_part # 1 = exit, 2 = consent withdrawn, 3 = deceased, 4 = LTFU
data$hospital_discharge <- as.Date(data$dt_dissue)
data$lesions_resolved <- data$si_pde___1
data$symptoms_resolved <- data$si_pde___2
data$home_care <- data$si_pde___3
data$other_hosp <- data$si_pde___4

##--------------------------------------------------
## Baseline visit: patient characteristics

data_base <- data[data$visit == 'baseline', ]
dim(data_base)
data_base$date <- as.Date(data_base$dtenqt) # inclusion date
summary(data_base$date)

table(data_base$critcompl_incl)
table(data_base$gender)/sum(table(data_base$gender)) # 1 = male, 2 = female
table(data_base$agecat)/sum(table(data_base$agecat)) # 1 = >18y; 2 = 12-17y; 3 = <12y
summary(data_base$agenum)

round(table(data_base$occupation_prim[data_base$agecat == 1])/sum(table(data_base$occupation_prim[data_base$agecat == 1])), 4)*100
table(data_base$occupation_sec)

table(data_base$HH_size); median(data_base$HH_size, na.rm = T); median(data_base$HH_size[data_base$HH_size!=132], na.rm = T)
summary(data_base$HH_size[data_base$HH_size!=132])
table(data_base$HH_mpox); table(data_base$HH_mpox[data_base$HH_size > 1])
table(data_base$HH_children); summary(data_base$HH_children)
sum(data_base$HH_children[data_base$HH_size > 0] == 0, na.rm = T)

table(data_base$HH_mpox)

table(data_base$travel) # 1 = yes, 2 = no
table(data_base$provnac)

## Behavior related to animals
table(data_base$part_chas) # participated in hunting session; 1 = yes, 2 = no, 3 = no response
table(data_base$consronger) # rodents caught or consumed; 1 = yes, 2 = no, 3 = no response
sum(data_base$contact1_id[data_base$consronger == 1] != "" , na.rm = T)
table(data_base$consronger_spec1)
table(data_base$viandbrous) # bushmeat consumption: 1 = never, 2 = <2x per week, 3 = >2x per week, 4 = no response
table(data_base$viandbrous_spec1) # 1 = monkey, 2 = squirrel, 3 = rat, 4 = bat, 5 = pangolin, 8 = antilope, 9 = other
table(data_base$viandbrous_spec2)
table(data_base$manipulmort) # touching or eating of animals that died naturally? 1 = yes, 2 = no, 3 = no response

## Previous vaccination
table(data_base$mpox_vacc_recent) # 1 = yes, 2 = no, 3 = NA
table(data_base$smallpox_vacc) # 1 = yes, 2 = no, 3 = NA (vaccinated before 1980, verified by skin mark)
data_base$contact1_id[data_base$smallpox_vacc==1] # only some of those that are vaccinated have identified contacts in the study

## Previous infection
table(data_base$previous_mpox)
summary(as.Date(data_base$date_previous_mpox))
table(data_base$already_included)
data_base$already_included_id[data_base$already_included==1]

## Other infections
table(data_base$HIV)
table(data_base$ART[data_base$HIV == 1])

## Information on sexual contacts
table(data_base$intercourse) # intercourse in the last 3 weeks; 1 = yes, 2 = no, 3 = no response
summary(data_base$num_partners) # number of different sexual partners in the last 3 weeks
sum(!is.na(data_base$num_partners[data_base$agecat != 3]))
hist(data_base$num_partners)
table(data_base$intercourse_SW) # intercourse with professional sex worker in the last 3 weeks; 1 = yes, 2 = no, 3 = no response, 4 = participant is SW
table(data_base$condomuse) # condom use: 1 = (almost) never, 2 = sometimes, 3 = (almost) always


##-------------------------------------------------------------------------------------------------
## Contact investigation

table(data_base$contact_mpox)
table(data_base$num_contact_mpox)
table(data_base$contact1_included) # mpox contact included in present study?
table(data_base$contact1_type) # 1 = single contact, 2 = multiple occasions / ongoing
table(data_base$contact1_rel) # 1 = family member in HH, 
# 2 = other family member, 3 = colleague, 4 = client (SW),
# 5 = SW visited, 6 = patient (participant is a HCW), 7 = deceased person, 8 =friend,
# 9 = neighbor, 10 = church, 11 = co-patient , 12 = sexual partner, 13 = client (shop)  
table(data_base$contact1_famtype) # 1 = parent, 2 = child, 3 = sibling, 4 = other ('data$otr_mbfam1')
table(data_base$contact1_zs)
table(data_base$contact1_bar)
table(data_base$contact1_prolonged) # >15min: 1 = yes, 2 = no, 3 = don't remember
table(data_base$contact1_direct) # direct touch
table(data_base$contact1_sexual)
table(data_base$contact1_care) # regularly took care of mpox patient
table(data_base$contact1_food)
table(data_base$contact1_bed)
table(data_base$contact1_clothes)

##-------------------------------------------------------------------------------------------------
## Clinical presentation: lesions on mucous membranes
table(data_base$oral_lesion)
table(data_base$tonsil_lesion)
table(data_base$penis_lesion)
table(data_base$vagina_lesion)

table(data_base$transm_sexual, data_base$penis_lesion)
table(data_base$transm_sexual, data_base$vagina_lesion)

## Symptoms: 1 = yes, 2 = no, 3 = don't know/no response
data_base$fatigue <- data_base$faigue # fatigue
# table(data_base$malaise) # malaise
# table(data_base$myalgie)
data_base$malaise <- data_base$malaise
data_base$myalgie <- data_base$myalgie
data_base$lesion_pain <- data_base$dassasio # pain associated to lesion
data_base$itching <- data_base$pruritcorps # itching on the body
data_base$cough <- data_base$toux 
data_base$dyspnea <- data_base$dyspn_e
data_base$anorexia <- data_base$anorexi
data_base$sore_throat <- data_base$malgorge
data_base$abdominal_pain <- data_base$dabdomie
data_base$vomiting <- data_base$vomissement
data_base$diarrhea <- data_base$diarrh_e
data_base$headache <- data_base$matete
data_base$confusion <- data_base$trblla_cns
data_base$seizures <- data_base$crise_clsives
data_base$eyepain <- data_base$dlr_oculaire
data_base$vision <- data_base$prt_vision
data_base$hemorrhage <- data_base$sagne
data_base$rectal_pain <- data_base$dlr_rectale
data_base$pain_urinating <- data_base$doulurin
data_base$genital_pain <- data_base$doulvagpeni
data_base$other_symptoms <- data_base$autre

## Follow-up ID
sum(data_base$codeadulte != ""); sum(data_base$codeenfant != "") # follow-up IDs
data_base$followupID <- ifelse(data_base$agecat == 1, data_base$codeadulte, data_base$codeenfant)
sum(data_base$followupID != '')

## Diagnostic sample results
table(data_base$diagnostics_complete) # 0 = incomplete, 2 = complete
table(data_base$mpox_pcr) # 1 = positive, 2 = negative, 3 = not analysed yet
table(data_base$sample_lesion)
table(data_base$sample_oropharyngal)
# longitudinal ID assigned if positive PCR result
data_base$id_adult 
data_base$id_child 
data_base$followupID

##----------------------------------------------------------------------------------------------
## Longitudinal data (hospital follow-up)
summary(as.Date(data$dtvisit))
table(data$visit[data$dtvisit == ''])
data$dtvisit <- as.Date(data$dtvisit)
summary(data$dtvisit)

data_hospital <- data[!is.na(data$dtvisit), ]
data_hospital <- data_hospital[data_hospital$suivi_fait == 1,]
data_hospital <- data_hospital[!is.na(data_hospital$isac), ]
dim(data_hospital); length(unique(data_hospital$isac))

# Add visit number
library(dplyr)
data_hospital <- data_hospital %>%
  dplyr::distinct(isac, dtvisit) %>%
  group_by(isac) %>%
  arrange(dtvisit, .by_group = TRUE) %>%
  dplyr::mutate(n_visit = row_number()) %>%
  right_join(data_hospital, by = c('isac', 'dtvisit'))
dim(data_hospital)
table(data_hospital$n_visit)

data_hospital$fatigue <- data_hospital$fatigue
data_hospital$malaise <- data_hospital$malaise_2b170b
data_hospital$myalgie <- data_hospital$myalgie_art
data_hospital$lesion_pain <- data_hospital$douleur_as # pain associated to lesion
data_hospital$itching <- data_hospital$prurit_a # itching on the body
data_hospital$cough <- data_hospital$toux_5c4877 
data_hospital$dyspnea <- data_hospital$dyspn_e_2dfa93
data_hospital$anorexia <- data_hospital$anorexie
data_hospital$sore_throat <- data_hospital$mal_d
data_hospital$abdominal_pain <- data_hospital$douleur_a
data_hospital$vomiting <- data_hospital$vomisse
data_hospital$diarrhea <- data_hospital$diarrh_e_19dff4
data_hospital$headache <- data_hospital$maux_dt
data_hospital$confusion <- data_hospital$troubles_d
data_hospital$seizures <- data_hospital$crises_con
data_hospital$eyepain <- data_hospital$douleur_oc
data_hospital$vision <- data_hospital$perte_d
data_hospital$hemorrhage <- data_hospital$saignement_h
data_hospital$rectal_pain <- data_hospital$douleur_rc
data_hospital$pain_urinating <- data_hospital$douleur_uri
data_hospital$genital_pain <- data_hospital$douleur_gnit
data_hospital$other_symptoms <- data_hospital$autre_fc4d0d


##---------------------------------------------------------------  
## Follow-up visit  

data_fu <- data[data$visit == 'visite_de_suivi', ]
dim(data_fu)
data_fu$dt_visite <- as.Date(data_fu$dt_visite)
data_fu <- data_fu[!is.na(data_fu$dt_visite), ]
table(data_fu$jr_delavisite) # 1 = day 29, 2 = day 59

table(data_fu$nbr_prsatteintmpox) # number of HH members that have acquired mpox

data_fu$fatigue <- data_fu$fatigue1
data_fu$malaise <- data_fu$malaise1
data_fu$myalgie <- data_fu$myalgie_arthralgie
data_fu$lesion_pain <- data_fu$dlr_associ1 # pain associated to lesion
data_fu$itching <- data_fu$prurites1 # itching on the body
data_fu$cough <- data_fu$toux1 
data_fu$dyspnea <- data_fu$dyspn_e1
data_fu$anorexia <- data_fu$anorexie1
data_fu$sore_throat <- data_fu$mal_dgorge1
data_fu$abdominal_pain <- data_fu$douleur_abdom1
data_fu$vomiting <- data_fu$vomissements1
data_fu$diarrhea <- data_fu$diarrh_e1
data_fu$headache <- data_fu$maux_tte1
data_fu$confusion <- data_fu$trbls_consc

data_fu$depression <- data_fu$d_pression
data_fu$anguish <- data_fu$anguish

data_fu$seizures <- data_fu$crises_convulsives
data_fu$eyepain <- data_fu$douleur_oculaire
data_fu$vision <- data_fu$perte_vision1

data_fu$decreased_vision <- data_fu$dimi_acuitvisuell
data_fu$itchy_eye <- data_fu$d_mangeaison

data_fu$photophobia <- data_fu$photophobie1

data_fu$hemorrhage <- data_fu$saignement1
data_fu$rectal_pain <- data_fu$douleur_rectale1
data_fu$pain_urinating <- data_fu$douleur_urinant
data_fu$genital_pain <- data_fu$douleur_gnitale
data_fu$painful_sex <- data_fu$douleur_contsex
data_fu$vaginal_discharge <- data_fu$pertes_vaginales
data_fu$vaginal_dryness <- data_fu$s_cheresse_vaginale

data_fu$other_symptoms <- data_fu$autre1

data_fu$lesions_present <- data_fu$patien_lesicuta


data_fu$HH_mpox <- data_fu$nbr_prsatteintmpox
data_fu$travel <- data_fu$partic_voyage
data_fu$intercourse <- data_fu$raport_sexdps # had intercourse since last visit (hospital discharge or day 29) 1 = yes, 2 = no, 3 = no response
data_fu$num_partners <- data_fu$avc_cbipartic
data_fu$intercourse_SW <- data_fu$rappsex_ps
data_fu$condomuse <- data_fu$utilis_presevatif
data_fu$oral_lesion <- data_fu$orale3
data_fu$tonsil_lesion <- data_fu$tonsilles
data_fu$penis_lesion <- data_fu$gland_pnis
data_fu$vagina_lesion <- data_fu$vulvo_vaginale


##---------------------------------------------------------------
## Combine symptom data

names(data_base)
data_base <- data_base[,c(1,7,21:26,234,235,449,881:1027)]
data_base <- data_base[,c(1,2,23,9,10,11,24:158)]
# data_base <- data_base[,c(1,2,6, 7:120, 122, 3:5, 123:141)]
dim(data_base)
names(data_base)[2:3] <- c('date','visit')
data_base$date <- as.Date(data_base$date)
data_base <- data_base[,c(1:120, 122:141)]

names(data_hospital)
data_hospital <- data_hospital[,c(1,2,3,9,236,237,450,894:1026)]
data_hospital$visit <- paste0('hosp_', data_hospital$n_visit)
data_hospital$followupID <- NA
data_hospital <- data_hospital[,c(1,2,141,5,6,7,8:140,142)]
dim(data_hospital)
names(data_hospital)[2:3] <- c('date','visit')

names(data_fu)
data_fu <- data_fu[,c(1,728,729,234,235,449,881:1033)]
data_fu$visit <- data_fu$jr_delavisite
data_fu <- data_fu[,c(1:6, 18:143, 145,146,147,151,152,153,154,158)]
data_fu$followupID <- NA
dim(data_fu)
data_fu <- data_fu[,c(1:6,8:141)]
names(data_fu)[2:3] <- c('date','visit')
data_fu$visit <- ifelse(data_fu$visit == 1, 'day29', 'day59')

names(data_fu) == names(data_base)
names(data_base) == names(data_hospital)

save(data_hospital, file = 'data/hospital_data_100325.RData')
save(data_base, file = 'data/baseline_data_100325.RData')
save(data_fu, file = 'data/followup_data_100325.RData')

data_symptoms <- rbind(data_hospital, data_base, data_fu)
dim(data_symptoms); length(unique(data_symptoms$isac))
table(data_symptoms$visit)

save(data_symptoms, file = 'data/symptom_data_100325.RData')

##----------------------------------------------------------------------------------------------
## Duration of hospitalization

load('data/symptom_data_100325.RData')

names(data_symptoms[grepl('hosp', data_symptoms$visit), ])
data_symptoms$hospital_duration <- as.numeric(data_symptoms$hospital_discharge - data_symptoms$date)
# summary(as.numeric(data_base$hospital_discharge - data_base$date))
summary(data_symptoms$hospital_duration); sum(!is.na(data_symptoms$hospital_duration))

length(unique(data_symptoms$isac[data_symptoms$hospital_duration < 0]))

##-----------------------------------------------------------------------------------------------
## Duration of symptoms

symptom.dat <- data_symptoms[,c(1,2,3,4,5,6,102:106,122:139)]
symptom.dat <- symptom.dat[,c(1:3,7,4:6, 8:29)]
symptom.dat <- data.frame(symptom.dat)
prop.included <- as.data.frame(table(symptom.dat$visit))
prop.included$Var1 <- factor(as.character(prop.included$Var1), levels = c('baseline', paste0('hosp_',c(1:28)), 'day29', 'day59'))
prop.included <- prop.included[order(prop.included$Var1), ]
prop.included$Prop <- prop.included$Freq / 590
symp <- c()
cols <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", 
  "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", 
  "khaki2",
  "maroon", 
  "orchid1", 
  "deeppink1", 
  "blue1", "steelblue4",
  "darkturquoise",
  "green1",
  "yellow4", "yellow3",
  "darkorange4",
  "brown"
)
c <- 1
symp.evol <- matrix(NA, nrow = length(unique(symptom.dat$visit)), ncol = length(names(symptom.dat)[5:29]))
jpeg('results/symptoms.jpeg', width = 30, height = 20, units = 'cm', res = 300)
for(s in 5:29){
  dat <- as.data.frame(table(symptom.dat[,3], symptom.dat[,s]))
  dat <- dat[dat$Var2 == 1, ]
  dat$Var1 <- factor(as.character(dat$Var1), levels = c('baseline', paste0('hosp_',c(1:28)), 'day29', 'day59'))
  dat <- dat[order(dat$Var1), ]
  dat$Prop <- dat$Freq / prop.included$Freq
  dat$Freq <- dat$Freq / 590
  dat$visit <- 1:dim(dat)[1]
  symp.evol[,c] <- dat$Prop
  symp <- c(symp, names(symptom.dat)[s])
  if(s == 5){
    plot(dat$visit, dat$Freq, type = 'l', xlab = '', ylab = 'Proportion of included patients', xaxt = 'n', ylim = c(0, 1))
  }else{
    lines(dat$visit, dat$Freq, type = 'l', xlab = '', ylab = 'Proportion of included patients', xaxt = 'n', ylim = c(0, 1), col = cols[c])
  }
  c <- c + 1
}
axis(side = 1, las = 2, labels = dat$Var1, at = seq(1,dim(dat)[1]))
lines(c(1:dim(prop.included)[1]), prop.included$Prop, lty = 2, col = 'darkgrey', lwd = 2)
legend('topright', c(symp, 'prop of cases included'), col = c(cols, 'darkgrey'), lty = c(rep(1, 20), 2), lwd = c(rep(1, 20), 2))
dev.off()

# Proportion of patients in hospital on that visit that report a certain symptom
colnames(symp.evol) <- names(symptom.dat)[5:29]
rownames(symp.evol) <- dat$Var1
brk <- 50
library(plot.matrix)

jpeg('results/symptom_evolution.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mar=c(7.1, 5.1, 4.1, 4.1))
plot(symp.evol, xlab = '', ylab = '', main = '',
     axis.col = list(side = 1, las = 2),
     axis.row = list(side = 2, las = 1),
     breaks = brk,
     col = rev(heat.colors(brk)),
     key = list(side = 3, cex.axis = 0.75), fmt.key = '%.2f',
     polygon.key = NULL, axis.key = NULL, #spacing.key = c(3,2,2),
     border = NA
     )
dev.off()

##-----------------------------------------------------------------------------------------------
## Duration of lesions

table(data$oral_lesion)
table(data$tonsil_lesion)
table(data$penis_lesion)
table(data$vagina_lesion)

symptom.dat <- data_symptoms[,c(1,2,3,4,5,6,102:106,122:139)]
symptom.dat <- symptom.dat[,c(1:3,7,4:6, 8:29)]
symptom.dat <- data.frame(symptom.dat)
prop.included <- as.data.frame(table(symptom.dat$visit))
prop.included$Var1 <- factor(as.character(prop.included$Var1), levels = c('baseline', paste0('hosp_',c(1:28)), 'day29', 'day59'))
prop.included <- prop.included[order(prop.included$Var1), ]
prop.included$Prop <- prop.included$Freq / 590
symp <- c()
cols <- c(
  "dodgerblue2", "#E31A1C", 
  "green4",
  "#FF7F00" 
)
c <- 1
symp.evol <- matrix(NA, nrow = length(unique(symptom.dat$visit)), ncol = length(names(symptom.dat)[8:11]))
jpeg('results/lesions.jpeg', width = 30, height = 20, units = 'cm', res = 300)
for(s in 8:11){
  dat <- as.data.frame(table(symptom.dat[,3], symptom.dat[,s]))
  dat <- dat[dat$Var2 == 1, ]
  dat$Var1 <- factor(as.character(dat$Var1), levels = c('baseline', paste0('hosp_',c(1:28)), 'day29', 'day59'))
  dat <- dat[order(dat$Var1), ]
  dat$Prop <- dat$Freq / prop.included$Freq
  dat$Freq <- dat$Freq / 590
  dat$visit <- 1:dim(dat)[1]
  symp.evol[,c] <- dat$Prop
  symp <- c(symp, names(symptom.dat)[s])
  if(s == 8){
    plot(dat$visit, dat$Freq, type = 'l', xlab = '', ylab = 'Proportion of included patients', xaxt = 'n', ylim = c(0, 1))
  }else{
    lines(dat$visit, dat$Freq, type = 'l', xlab = '', ylab = 'Proportion of included patients', xaxt = 'n', ylim = c(0, 1), col = cols[c])
  }
  c <- c + 1
}
axis(side = 1, las = 2, labels = dat$Var1, at = seq(1,dim(dat)[1]))
lines(c(1:dim(prop.included)[1]), prop.included$Prop, lty = 2, col = 'darkgrey', lwd = 2)
legend('topright', c(symp, 'prop of cases included'), col = c(cols, 'darkgrey'), lty = c(rep(1, 4), 2), lwd = c(rep(1, 4), 2))
dev.off()

# Proportion of patients in hospital on that visit that report a certain symptom
colnames(symp.evol) <- names(symptom.dat)[8:11]
rownames(symp.evol) <- dat$Var1
brk <- 50
library(plot.matrix)

jpeg('results/lesion_evolution.jpeg', width = 30, height = 20, units = 'cm', res = 300)
par(mar=c(7.1, 5.1, 4.1, 4.1))
plot(symp.evol, xlab = '', ylab = '', main = '',
     axis.col = list(side = 1, las = 2),
     axis.row = list(side = 2, las = 1),
     breaks = brk,
     col = rev(heat.colors(brk)),
     key = list(side = 3, cex.axis = 0.75), fmt.key = '%.2f',
     polygon.key = NULL, axis.key = NULL, #spacing.key = c(3,2,2),
     border = NA
)
dev.off()
