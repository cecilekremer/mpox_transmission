
load('data/clean_data.RData')

##------------------------------------------------------
## Only baseline visit data is needed (N = 415)

data.contact <- data[data$visit == 'baseline', ]
data.contact$date <- as.Date(data.contact$dtenqt) # inclusion date

## Number of mpox cases contacted
table(data.contact$contact_mpox)
sum(data.contact$num_contact_mpox, na.rm = T)

table(data.contact$contact1_included)
sum(data.contact$contact1_included!=1, na.rm = T)/sum(table(data.contact$contact1_included))
sum(data.contact$contact2_included!=1, na.rm = T)/sum(table(data.contact$contact2_included))
sum(data.contact$contact3_included!=1, na.rm = T)/sum(table(data.contact$contact3_included))

num.contacts <- sum(data.contact$num_contact_mpox, na.rm = T)
num.contacts.included <- sum(data.contact$contact1_included == 1, na.rm = T) + sum(data.contact$contact2_included == 1, na.rm = T) + sum(data.contact$contact3_included == 1, na.rm = T)
(num.contacts - num.contacts.included)/num.contacts

## Info available for contacts not included?
table(data.contact$contact1_rel[data.contact$contact1_included != 1])

##------------------------------------------------------
## Contacts that are included in the study

table(data.contact$contact4_included) # 1 = yes
sum(data.contact$contact4_id != '')

data.contact$contact1_id <- ifelse(grepl('SCREEN', data.contact$contact1_id, ignore.case = T), data.contact$contact1_id,
                                   ifelse(data.contact$contact1_id != '',
                                          ifelse(nchar(data.contact$contact1_id) == 3, paste0('SCREEN', data.contact$contact1_id), paste0('SCREEN0', data.contact$contact1_id)),
                                          data.contact$contact1_id)
)
data.contact$contact1_id <- ifelse(grepl('SCREENE', data.contact$contact1_id, ignore.case = T), gsub('SCREENE', 'SCREEN', data.contact$contact1_id),
                                   data.contact$contact1_id)
data.contact$contact1_id <- gsub(' ', '', data.contact$contact1_id)
data.contact$contact1_id <- gsub('Screen', 'SCREEN', data.contact$contact1_id)
data.contact$contact1_id <- gsub('SCREEN88', 'SCREEN088', data.contact$contact1_id)
all(nchar(data.contact$contact1_id[data.contact$contact1_id!='']) == 9) 

data.contact$contact2_id <- ifelse(grepl('SCREEN', data.contact$contact2_id, ignore.case = T), data.contact$contact2_id,
                                   ifelse(data.contact$contact2_id != '',
                                          ifelse(nchar(data.contact$contact2_id) == 3, paste0('SCREEN', data.contact$contact2_id), paste0('SCREEN0', data.contact$contact2_id)),
                                          data.contact$contact2_id)
)
all(nchar(data.contact$contact2_id[data.contact$contact2_id!='']) == 9) 

data.contact$contact3_id <- ifelse(grepl('SCREEN', data.contact$contact3_id, ignore.case = T), data.contact$contact3_id,
                                   ifelse(data.contact$contact3_id != '',
                                          ifelse(nchar(data.contact$contact3_id) == 3, paste0('SCREEN', data.contact$contact3_id), paste0('SCREEN0', data.contact$contact3_id)),
                                          data.contact$contact3_id)
)
all(nchar(data.contact$contact3_id[data.contact$contact3_id!='']) == 9) 

data.contact$ID <- as.numeric(unlist(regmatches(data.contact$isac, gregexpr("[0-9]+", data.contact$isac))))

data.contact$contacts <- NA
for(i in 1:dim(data.contact)[1]){
  
  # IDs of contacts that are included in the study
  contacts <- c(as.numeric(unlist(regmatches(data.contact$contact1_id[i], gregexpr("[0-9]+", data.contact$contact1_id[i])))),
                as.numeric(unlist(regmatches(data.contact$contact2_id[i], gregexpr("[0-9]+", data.contact$contact2_id[i])))),
                as.numeric(unlist(regmatches(data.contact$contact3_id[i], gregexpr("[0-9]+", data.contact$contact3_id[i]))))
  )
  contacts <- contacts[contacts %in% unique(data.contact$ID)]
  
  data.contact$contacts[i] <- paste(contacts, collapse = ',')
  
}
sum(data.contact$contacts != '')

##-------------------------------------------------------------------------
## Rash onset

table(data.contact$lesion_cuta)
table(data.contact$dptemps) # since how long is there a rash (in days)

data.contact$rash_onset <- as.Date(data.contact$date) - data.contact$dptemps
summary(data.contact$rash_onset)

table(data.contact$lsprof)
table(data.contact$reg_genital) # how many lesions on genital region
data.contact$n.lesion.genital <- data.contact$reg_genital
table(data.contact$reg_anl) # how many lesions on anal region
data.contact$n.lesion.anal <- data.contact$reg_anl

##-------------------------------------------------------------------------
## Save contact dataset

# Variables to keep:

data.contact.clean <- data.contact[, c(987:989,
                                       878:986,
                                       867:877,
                                       990:992
)]
dim(data.contact.clean)
save(data.contact.clean, file = "data/contact_data.RData")

