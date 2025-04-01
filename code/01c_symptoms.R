
load('data/symptom_data_100325.RData')

table(data_symptoms$transm_sexual[data_symptoms$visit!='baseline'])
data_symptoms$transm_sexual[data_symptoms$visit != 'baseline'] <- NA

table(data_symptoms$transm_sexual)

## Add transmission route to follow-up data
library(dplyr)
library(tidyr)
data_symptoms <- data_symptoms %>%
  group_by(isac) %>%
  fill(agecat, gender, transm_sexual, .direction = 'updown')
table(data_symptoms$transm_sexual)

table(data_symptoms$transm_sexual, data_symptoms$penis_lesion)
table(data_symptoms$transm_sexual, data_symptoms$vagina_lesion)

##-----------------------------------------------------------------------------------------------
## Duration of symptoms

data_plot <- data_symptoms[data_symptoms$visit %in% c('baseline', paste0('hosp_',1:6), 'day29', 'day59'), ]
table(data_plot$visit)
data_plot$visit <- factor(as.character(data_plot$visit), levels = c('baseline', paste0('hosp_',1:6), 'day29', 'day59'))
library(tidyverse)

data_plot <- data_plot %>%
  mutate(across(c(4:6, 103:106, 121:138), ~ ifelse(as.numeric(.) == 1 & !is.na(as.numeric(.)), 1, 0)))
# data_plot <- data_plot %>%
#   mutate(across(c(4,103,121), ~ ifelse(. == 1 & !is.na(.), 1, 0)))
data_plot$lesion_pain <- ifelse(data_plot$lesion_pain == 1 & !is.na(data_plot$lesion_pain), 1, 0)
data_plot$malaise <- ifelse(data_plot$malaise == 1 & !is.na(data_plot$malaise), 1, 0)
data_plot$oral_lesion <- ifelse(data_plot$oral_lesion == 1 & !is.na(data_plot$oral_lesion), 1, 0)
df_plot <- data_plot %>%
  pivot_longer(cols = c(4:6, 103:106, 121:138),
               names_to = 'symptom',
               values_to = 'frequency') %>%
  mutate(frequency = replace_na(frequency, 0))
group_sizes <- data_plot[data_plot$visit == 'baseline', ] %>%
  group_by(transm_sexual) %>%
  summarise(n_individuals = n(), .groups = 'drop')
df_plot <- df_plot %>%
  group_by(visit, transm_sexual, symptom) %>%
  summarise(frequency = sum(frequency, na.rm = T)) %>%
  left_join(group_sizes, by = 'transm_sexual') %>%
  mutate(proportion = frequency / n_individuals)
head(df_plot)
summary(df_plot$proportion); table(df_plot$symptom[df_plot$proportion>1])

library(ggplot2)
library(grid)
ggplot(df_plot, aes(x = symptom, y = proportion, fill = factor(transm_sexual))) +
  facet_wrap(~visit) +
  geom_bar(stat = 'identity', position = 'dodge') +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1),
  #       legend.position = 'top',
  #       axis.text.x.top = element_text(angle = 90, hjust = 1),
  #       axis.title.x.top = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        legend.position = 'top',
        axis.text.y = element_text(size = 14),  # Increase size of y-axis text
        axis.title.x = element_text(size = 16),  # Increase size of x-axis title
        axis.title.y = element_text(size = 16),  # Increase size of y-axis title
        strip.text = element_text(size = 16),  # Increase size of facet labels (visit labels)
        legend.text = element_text(size = 14),  # Increase size of legend text
        legend.title = element_text(size = 16)  # Increase size of legend title) 
  ) +
  # scale_x_discrete(sec.axis = dup_axis(name = NULL)) +
  labs(x = '', y = 'Proportion of cases', fill = 'Sexual transmission') 
# +
# geom_text(aes(x = symptom, y = max(proportion) + 0.05, label = symptom), 
# angle = 90, hjust = 1, vjust = 0, size = 3, inherit.aes = FALSE)
ggsave(filename = './results/symptoms_by_transm.jpeg', width = 40, height = 30, units = 'cm')

baseline_freq <- df_plot %>%
  filter(visit == 'baseline') %>%
  group_by(symptom, transm_sexual) %>%
  summarize(frequency = sum(frequency)) %>%
  arrange(desc(frequency))
baseline_freq[baseline_freq$transm_sexual == 1, ]

df_plot[df_plot$symptom %in% c('oral_lesion', 'tonsil_lesion', 'penis_lesion', 'vagina_lesion') & df_plot$visit == 'baseline', ]


##-----------------------------------------------------------------------------------------------

# Sexual transmission

symptom.dat <- data_symptoms[data_symptoms$transm_sexual == 0, ]
symptom.dat <- symptom.dat[,c(1,2,3,4,5,6,102:106,122:139)]
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
  "gray70", "khaki2", "maroon", "orchid1",  "deeppink1", "blue1", "steelblue4","darkturquoise", "green1",  "yellow4", "yellow3",  "darkorange4",  "brown")

c <- 1
symp.evol <- matrix(NA, nrow = length(unique(symptom.dat$visit)), ncol = length(names(symptom.dat)[c(8:11,26,27,28)]))
jpeg('results/symptoms_nonsexual.jpeg', width = 30, height = 20, units = 'cm', res = 300)
for(s in c(8:11,26,27,28)){
  dat <- as.data.frame(table(symptom.dat[,3], symptom.dat[,s]))
  dat <- dat[dat$Var2 == 1, ]
  dat$Var1 <- factor(as.character(dat$Var1), levels = c('baseline', paste0('hosp_',c(1:28)), 'day29', 'day59'))
  dat <- dat[order(dat$Var1), ]
  dat$Prop <- dat$Freq / prop.included$Freq
  dat$Freq <- dat$Freq / length(unique(symptom.dat$isac))
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
# lines(c(1:dim(prop.included)[1]), prop.included$Prop, lty = 2, col = 'darkgrey', lwd = 2)
# legend('topright', c(symp, 'prop of cases included'), col = c(cols, 'darkgrey'), lty = c(rep(1, 20), 2), lwd = c(rep(1, 20), 2))
legend('topright', symp, col = cols, lty = c(rep(1, 20)), lwd = c(rep(1, 20)))
dev.off()

# Proportion of patients in hospital on that visit that report a certain symptom
colnames(symp.evol) <- names(symptom.dat)[c(8:11,26,27,28)]
rownames(symp.evol) <- dat$Var1
brk <- 50
library(plot.matrix)

jpeg('results/symptom_evolution_nonsexual.jpeg', width = 30, height = 20, units = 'cm', res = 300)
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
