
# load('./data/CT_data_230525.RData')
load('./Final code/data/ct_data.RData')
names(ct_all)
length(unique(ct_all$isac)); sum(grepl('SCREEN', ct_all$isac)) # Kamituga cases
table(ct_all$site)
summary(ct_all$Ct.value)

ct_all <- ct_all[!is.na(ct_all$Ct.value), ]
table(ct_all$site)

ct_pos <- ct_all
# ct_pos <- ct_pos[!is.na(ct_pos$HEX_ct_lesion), ]
dim(ct_pos); length(unique(ct_pos$isac))
# sum(grepl('SCREEN', ct_pos$isac))
# sum(!is.na(ct_pos$OPX_ct_lesion[!grepl('SCREEN',ct_pos$isac)]) & ct_pos$HEX_ct_lesion[!grepl('SCREEN',ct_pos$isac)] != 0, na.rm = T)

summary(ct_pos$Ct.value)
hist(ct_pos$Ct.value)

par(mfrow=c(2,1))
hist(ct_pos$Ct.value[ct_pos$site == 'Kamituga'], breaks = 20, main = '', xlab = 'OPX Ct (Kamituga)')
hist(ct_pos$Ct.value[ct_pos$site == 'Goma'], breaks = 20, main = '', xlab = 'OPX Ct (Goma)')

###-----------------------------------------------------------------
### Based on rash onset
load('./Final code/data/data_GM.RData')
data_Goma <- data_Goma
data_Goma$rash_onset <- as.Date(data_Goma$dtenqt) - data_Goma$dptemps
data_Goma$isac <- data_Goma$record_id
summary(data_Goma$rash_onset)
load('./Final code/data/data_KT.RData')
# data_Kamituga <- data_Kamituga[data_Kamituga$visit == 'baseline', ]
data_Kamituga$rash_onset <- as.Date(data_Kamituga$dtenqt) - data_Kamituga$dptemps
summary(data_Kamituga$rash_onset)

ids <- c(data_Goma$record_id, data_Kamituga$isac)
sum(ct_pos$isac %in% ids)

data_all <- rbind(data_Goma[,c('isac', 'rash_onset')],
                  data_Kamituga[,c('isac', 'rash_onset')])
summary(data_all$rash_onset)

ct_pos <- merge(ct_pos, data_all, by.x = 'isac', by.y = 'isac', all.x = TRUE)
summary(ct_pos$rash_onset)

ct_pos$rashOnset.to.sample <- as.numeric(as.Date(ct_pos$sample_date) - as.Date(ct_pos$rash_onset))
summary(ct_pos$rashOnset.to.sample)
par(mfrow = c(1,1))
hist(ct_pos$rashOnset.to.sample)
sum(ct_pos$rashOnset.to.sample <= 3, na.rm = T)

aggregate(ct_pos$Ct.value[ct_pos$agecat == 1 & ct_pos$rashOnset.to.sample <= 3], 
          by = list(ct_pos$transm_sexual[ct_pos$agecat == 1 & ct_pos$rashOnset.to.sample <= 3]), mean)

library(ggplot2)
ggplot(ct_pos[ct_pos$rashOnset.to.sample<=3, ], aes(x = factor(agecat), y = Ct.value, fill = factor(transm_sexual))) +
  geom_boxplot() +
  xlab('Age group (1 = adults, 2 = 12-17y, 3 = <12y)') + ylab('Ct value') + labs(fill = 'Sexual')

# mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual*agenum, data = ct_pos[grepl('SCREEN', ct_pos$isac), ])
mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual + agenum, data = ct_pos[ct_pos$site == 'Kamituga', ])
summary(mod)
modG <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual * agenum, data = ct_pos[ct_pos$site == 'Goma', ])
summary(modG)

mod <- lm(Ct.value ~ rashOnset.to.sample + factor(transm_sexual)*agenum + factor(transm_sexual)*factor(site) , data = ct_pos)
summary(mod)
mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual + agenum + site, data = ct_pos)
summary(mod)

# summary(lm(Ct.value ~ transm_sexual + agenum, data = ct_pos))
# 
# summary(lm(Ct.value ~ rashOnset.to.sample + transm_sexual + agecat, data = ct_pos[ct_pos$agecat != 3, ]))

ct_pos$date.breaks <- ifelse(ct_pos$rashOnset.to.sample < 3, 1,
                             ifelse(ct_pos$rashOnset.to.sample < 7, 2, 
                                    ifelse(ct_pos$rashOnset.to.sample <= 10, 3, 4)))
table(ct_pos$date.breaks)

ct_pos$date.breaks <- factor(ct_pos$date.breaks,
                             labels = c('Sample < 3 days after rash', 'Sample 3-6 days after rash', 
                                        'Sample 7-10 days after rash', 'Sample > 10 days after rash'))
# names(facet.labs) <= c('1','2','3')
ggplot(ct_pos[!is.na(ct_pos$date.breaks), ], aes(x = factor(agecat), y = Ct.value, fill = factor(transm_sexual))) +
  facet_wrap(site~date.breaks, nrow = 2) +
  geom_boxplot() +
  # xlab('Age group (1 = adults, 2 = 12-17y, 3 = <12y)') + 
  xlab('Age group') +
  ylab('OPX Ct value') + 
  # labs(fill = 'Sexual') +
  labs(fill = 'Transmission mode') +
  scale_x_discrete(labels = c('1'='Adults', '2'='12-17y', '3'='<12y')) +
  scale_fill_discrete(labels = c('0' = 'Non-sexual', '1'='Sexual')) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text.x = element_text(size = 18)
        )
ggsave('./Final code/results/FigureS2.jpeg', width = 35, height = 20, units = 'cm')

###---------------------------------------------------------------------------------
### Cases included in serial interval estimation

ct_ids <- ct_pos[grepl('SCREEN', ct_pos$isac), ]
ct_ids$IDnum <- as.numeric(unlist(regmatches(ct_ids$isac, gregexpr("[0-9]+", ct_ids$isac))))
load('./Final code/data/data_obs_SI.RData')
# remove duplicate
data.si <- data.si.id
data.si <- data.si[-which(data.si$ID==839)[1], ]
dim(data.si)

sum(ct_ids$IDnum %in% unique(data.si$ID))
sum(!(ct_ids$IDnum %in% unique(data.si$ID)))

sum(data.si$ID %in% ct_ids$IDnum)

ct_ids <- ct_ids[ct_ids$IDnum %in% data.si$ID, ]
mod <- lm(Ct.value ~ rashOnset.to.sample + factor(transm_sexual) + agenum, data = ct_ids)
# mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual + factor(agecat), data = ct_ids)
# mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual, data = ct_ids) # not significant when excluding age...
summary(mod)



## Ct values available for infector as well as infectee?
sum(ct_ids$IDnum %in% data.si$contacts)

###------------------------------------------------------------------------------------
### Cases included in incubation period estimation?

load("./Final code/data/ids_incubation.RData")
sum(ids.incubation %in% ct_pos$isac)
sum(ct_pos$isac %in% ids.incubation)

ct_ip <- ct_pos[ct_pos$isac %in% ids.incubation, ]
dim(ct_ip)
table(ct_ip$site)

hist(ct_ip$Ct.value[ct_ip$site=='Kamituga'])
hist(ct_ip$Ct.value)

# mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual + agenum, data = ct_ip[ct_ip$site == 'Kamituga', ])
# mod <- lm(Ct.value ~ rashOnset.to.sample + transm_sexual, data = ct_ip[ct_ip$location == 'Kamituga', ]) # not significant when excluding age
mod <- lm(Ct.value ~ rashOnset.to.sample + factor(transm_sexual) + agenum + factor(site), data = ct_ip)
summary(mod)


