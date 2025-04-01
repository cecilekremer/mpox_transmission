
##--------------------------------
## Goma
load('./data/clean_data_Goma_100325.RData')
ct_goma <- data[,c(1,370:403,523:526,365)]
ct_goma <- ct_goma[,c(1,5:9,13:17,21:25,29:40)]
# HEX or OPX ct values
ct_goma <- ct_goma[,c(1,2,5,6,7,10,11,12,15,16,17,20,21,22:28)]
names(ct_goma) <- c('isac', 'result_lesion', 'HEX_ct_lesion', 'OPX_ct_lesion', 'result_oro', 'HEX_ct_oro', 'OPX_ct_oro',
                    'result_blood', 'HEX_ct_blood', 'OPX_ct_blood',
                    'result_saliva', 'HEX_ct_saliva', 'OPX_ct_saliva', 'concl_result', 'diagnostics_complete',
                    names(ct_goma)[16:19],'sample_date')

## Take lesion Ct values (to be same as Kamituga)
ct_goma$HEX_ct_lesion <- as.numeric(ct_goma$HEX_ct_lesion)
summary(ct_goma$HEX_ct_lesion)
ct_goma$OPX_ct_lesion <- as.numeric(ct_goma$OPX_ct_lesion)
summary(ct_goma$OPX_ct_lesion)
sum(is.na(ct_goma$HEX_ct_lesion) & is.na(ct_goma$OPX_ct_lesion))

ct_goma <- ct_goma[!(is.na(ct_goma$HEX_ct_lesion) & is.na(ct_goma$OPX_ct_lesion)), ]
length(unique(ct_goma$isac))

## Include other data: transmission route, age, gender
demo_goma <- data[,c(1,422:427,430,431:442,507,508,509,513,514,517:528,516)]
names(demo_goma)[1] <- 'isac'

## To combine with Kamituga data
ct_goma <- ct_goma[,c(1,20,16,2,3,4,14,15)]
ct_goma <- merge(ct_goma, demo_goma, by.x = 'isac', by.y = 'isac', all.x = TRUE, all.y = FALSE)
names(ct_goma)[c(3,8)] <- c('sample_lesion','diagnostics_complete')
ct_goma <- ct_goma[,-c(37,39)]
length(unique(ct_goma$isac))
names(ct_goma)

##--------------------------------
## Kamituga
ct_kamituga <- read.csv('./data/kamituga_ct.csv', header = TRUE, sep = ';')
dim(ct_kamituga); length(unique(ct_kamituga$isac)); length(unique(ct_kamituga$id.conf))

sum(ct_kamituga$clade.Ib.ct != '')
ct_kamituga <- ct_kamituga[ct_kamituga$clade.Ib.ct != '', ]
dim(ct_kamituga); length(unique(ct_kamituga$isac)); length(unique(ct_kamituga$id.conf))
# 219 Ct values for 147 individuals?

table(ct_kamituga$sample.type)
table(ct_kamituga$sample.n)

table(ct_kamituga$xpert.result)
table(ct_kamituga$xpert.ct[ct_kamituga$xpert.result == 'POSITIF'])

table(ct_kamituga$clade.Ib.result)
table(ct_kamituga$clade.Ib.ct[ct_kamituga$clade.Ib.result == 'positif'])

ct_kamituga$CT <- as.numeric(gsub(',','.',ct_kamituga$clade.Ib.ct))
summary(ct_kamituga$CT) # cutoff = 38? (above = negative)
summary(ct_kamituga$CT[ct_kamituga$clade.Ib.result == 'positif'])
summary(ct_kamituga$CT[ct_kamituga$clade.Ib.result == 'negatif'])

### Individuals with more than one Ct value?
ids.mult.ct <- c()
for(i in ct_kamituga$isac[duplicated(ct_kamituga$isac)]){
  cts <- ct_kamituga$CT[ct_kamituga$isac == i]
  
  if(max(cts) != min(cts)){
    ids.mult.ct <- c(ids.mult.ct, i)
  }
}
ids.mult.ct # all equal

ct_kamituga <- ct_kamituga[,c(12,15,5,19,27,23,21)]
names(ct_kamituga)
length(unique(ct_kamituga$isac))

library(dplyr)
ct_kamituga <- ct_kamituga %>% 
  distinct(isac, .keep_all = TRUE)

## Add demographics
load('./data/clean_data_100325.RData')
data <- data[data$visit=='baseline',]
demo_kamituga <- data[,c(1,993,893:913,978:980,985,986,989:992,994:996,999,1000,988)]

ct_kamituga <- merge(ct_kamituga, demo_kamituga, by.x = 'isac', by.y = 'isac', all.x = TRUE, all.y = FALSE)
ct_kamituga <- ct_kamituga[,c(1:14,17:39,41,40,42:44)]
ct_kamituga$sample_blood <- NA
ct_kamituga$sample_saliva<-NA

ct_kamituga <- ct_kamituga[,c(1,2,39,4,5,6,7,8,9:38,42,43,40,41,42)]

### Combine datasets
names(ct_kamituga) <- names(ct_goma)
ct_all <- rbind(ct_kamituga, ct_goma)
dim(ct_all); length(unique(ct_all$isac))

save(ct_all, file = 'CT_data_170325.RData')

