
load('data/clean_data.RData')

dim(data); names(data)

##--------------------------------------------------
## Baseline visit: patient characteristics

data_base <- data[data$visit == 'baseline', ]
data_base$date <- as.Date(data_base$dtenqt) # inclusion date
summary(data_base$date)

table(data_base$critcompl_incl)
table(data_base$gender)/sum(table(data_base$gender)) # 1 = male, 2 = female
table(data_base$agecat)/sum(table(data_base$agecat)) # 1 = >18y; 2 = 12-17y; 3 = <12y
summary(data_base$agenum)

round(table(data_base$occupation_prim[data_base$agecat == 1])/sum(table(data_base$occupation_prim[data_base$agecat == 1])), 4)*100
table(data_base$occupation_sec)

table(data_base$HH_size); median(data_base$HH_size, na.rm = T); median(data_base$HH_size[data_base$HH_size!=132], na.rm = T)
table(data_base$HH_mpox); table(data_base$HH_mpox[data_base$HH_size > 1])
table(data_base$HH_children); summary(data_base$HH_children)
sum(data_base$HH_children[data_base$HH_size > 1] == 0, na.rm = T)

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
data_base$contact1_id[data_base$smallpox_vacc==1] # none of those that are vaccinated have identified contacts in the study

## Previous infection
table(data_base$previous_mpox)
summary(data_base$date_previous_mpox)
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


##----------------------------------------------------------------------------------------------
## Longitudinal data (hospital follow-up)

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


## All data
load('data/clean_data.RData')
dim(data); names(data)
data$followupID <- ifelse(data$agecat == 1, data$codeadulte, data$codeenfant)
sum(data$followupID != '', na.rm = T)

table(data$suivi_fait)
length(unique(data$isac[data$suivi_fait==1])) # ID's of patients that had follow-up during hospitalisation

## Hospital follow-up 
summary(data$dtvisit); table(data$visit); table(data$suivi_fait) 

