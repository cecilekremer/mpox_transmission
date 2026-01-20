
load('./Final code/data/contact_data.RData')

summary(data.all$date)

##------------------------------------
## Demographics
table(data.all$gender) # 2 = female
table(data.all$agecat) # 1 = >=18y, 2 = 12-17y, 3 = <12y
summary(data.all$agenum)
# h <- hist(data.all$agenum, breaks = 80, plot = F)
# hist(data.all$agenum, xlab = 'Age (years)', ylab = 'Number of participants', main = '', breaks = 80)
library(ggplot2)
ggplot(data = data.all, aes(x = agenum)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = 'steelblue', color = 'white', alpha = 0.7) +
  geom_density(size=2, color='red') +
  # geom_density(color = 'darkred', size = 1) + # only if usiing geom_histogram(aes(y = ...density...))
  scale_x_continuous(breaks = c(0,5,12,18,30,40,50,60,70,80)) +
  geom_vline(xintercept = c(12, 18), color = 'red', linetype = "dashed", size = 1) +
  xlab('Age (years)') + ylab('Relative frequency') +
  theme_minimal(base_size = 14)
# ggsave('./Final code/results/FigureS1.jpeg', width = 30, height = 20, units = 'cm', dpi = 300)

##-----------------------------------
## Sexual behavior
table(data.all$n.lesion.anal)
sum(data.all$n.lesion.genital[data.all$gender==2] > 0, na.rm = T) # genital lesions women
sum(data.all$n.lesion.genital[data.all$gender==1] > 0, na.rm = T) # genital lesions men

table(data.all$intercourse[data.all$agecat != 3]) # 1 = had intercourse in the last 3 weeks
table(data.all$intercourse[data.all$agecat != 3], data.all$agecat[data.all$agecat != 3])
summary(data.all$agenum[data.all$intercourse==1 & data.all$agecat!=3])
sum(data.all$agenum[data.all$intercourse==1 & data.all$agecat!=3] < 15, na.rm = T)

table(data.all$intercourse_SW[data.all$agecat != 3 & data.all$intercourse == 1]) # 1 = had intercourse with SW in the last 3 weeks; 4 = person is SW themselves
sum(data.all$agenum[data.all$intercourse_SW==1] < 15, na.rm = T)

##-----------------------------------
## Occupation
# load('./data/clean_data_100325.RData') # Kamituga
load('./Final code/data/data_KT.RData')
data_KT <- data_Kamituga
# data_KT <- data_KT[data_KT$visit == 'baseline', ] # only including baseline visit
# load('./data/clean_data_Goma_100325.RData') # Goma
load('./Final code/data/data_GM.RData')
data_GM <- data_Goma

table(data_KT$occupation_prim)
# 1=HCW, 2=traditional patrician, 3=SW, 4=mine worker, 5=barman, 6=hunter, 7=fisherman, 8=farmer, 9=trader, 10=biker
# 11=carpenter, 12=trucker, 13=law enforcement, 14=teacher, 15=pastor, 16=student, 17=child(not school-going), 18=other, 19=no response
table(data_GM$occupation_prim)
## Among adults
table(c(data_GM$occupation_prim[data_GM$agecat==1], data_KT$occupation_prim[data_KT$agecat==1]))


## HH mpox
table(c(data_KT$HH_mpox, data_GM$HH_mpox))
sum(c(data_KT$HH_mpox, data_GM$HH_mpox) > 0, na.rm = T)
