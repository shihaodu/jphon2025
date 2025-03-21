#This R script generates the figures and stats tables used in the paper

library(ggpubr)
library(patchwork)
library(ggplot2)
library(RGraphics)
library(R.matlab)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(modelsummary)
library(emmeans)

##################################Load data##################################
data = read.table(file.choose(), header=TRUE)
data$cluster <- factor(data$cluster, levels = c('bl','pl','gl','kl'))

DE = subset(data, language=="German")
SP = subset(data, language=="Spanish")
SG = subset(data, language=="Swiss German")

##################################Figure 2#################################
token = 'klage_w_05'
subj = 'DE5'
language = 'German'

matfile = readMat(file.choose()) #load the .mat file named 'klage_w_05.mat' in the repository

signal = matfile[[1]]

for (j in 1:length(signal)){
  if (is.character(signal[[j]]) & length(signal[[j]])==1){
    if (signal[[j]] == 'TT'){TT_ind = j + 2}
    if (signal[[j]] == 'TB'){TB_ind = j + 2}
  }
}

TT = as.data.frame(signal[[TT_ind]])
TT = mutate(TT, time = row_number())
TB = as.data.frame(signal[[TB_ind]])
TB = mutate(TB, time = row_number())

C1_start = data[data$SOURCE==token&data$subj==subj,]$C1_GONS..ms. / 4
C2_end = data[data$SOURCE==token&data$subj==subj,]$C2_GOFFS..ms. / 4

C1 = TB[(C1_start-26):(C2_end+25),]
C2 = TT[(C1_start-26):(C2_end+25),]

C1$time = C1$time * 4
C2$time = C2$time * 4

C1$dist = rep(NA,nrow(C1))
C2$dist = rep(NA,nrow(C2))
C1$vel = rep(NA,nrow(C1))
C2$vel = rep(NA,nrow(C2))

for (i in 2:nrow(C1)){
  C1[i,]$dist = sqrt((C1[i,]$V1-C1[i-1,]$V1)^2+
                       (C1[i,]$V2-C1[i-1,]$V2)^2+
                       (C1[i,]$V3-C1[i-1,]$V3)^2)
  C1[i,]$vel = C1[i,]$dist*100/4
}

for (i in 2:nrow(C2)){
  C2[i,]$dist = sqrt((C2[i,]$V1-C2[i-1,]$V1)^2+
                       (C2[i,]$V2-C2[i-1,]$V2)^2+
                       (C2[i,]$V3-C2[i-1,]$V3)^2)
  C2[i,]$vel = C2[i,]$dist*100/4
}

C1_pos = gather(C1[-1,],dim,value,c(V1,V3,vel))
C1_pos$dim = c(rep('X', 130), rep('Y', 130), rep('Velocity', 130))
C1_pos$type = c(rep('Position (mm)', 260), rep('Velocity (cm/s)',130))

C2_pos = gather(C2[-1,],dim,value,c(V1,V3,vel))
C2_pos$dim = c(rep('X', 130), rep('Y', 130), rep('Velocity', 130))
C2_pos$type = c(rep('Position (mm)', 260), rep('Velocity (cm/s)',130))

p1 <- ggplot(C1_pos, aes(x = time, y = value)) +
  geom_line(linewidth = 1, aes(col = dim)) +
  facet_grid(rows = vars(type), scales = 'free_y')+
  labs(
    x = 'Time (ms)',
    y = 'TB',
    title = '',
    col = '') +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12))

ann_text <- data.frame(time = c(1804,1840,1888,1924,1964,2012,2100),
                       value = c(0,25,0,3,0,18,0), lab = c("onset","pvel_to","target","mvel","release","pvel_fro","offset"),
                       type = c(factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Velocity (cm/s)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Velocity (cm/s)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Velocity (cm/s)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)'))))

p_TB = p1 + geom_vline(xintercept=c(1804,1840,1888,1924,1964,2012,2100), color = 'black', linewidth = 0.5, linetype = 'dashed') + 
  geom_text(data = ann_text, size = 5,
            label = c("onset","pvel_to","target","mvel","release","pvel_fro","offset"))

p2 <- ggplot(C2_pos, aes(x = time, y = value)) +
  geom_line(linewidth = 1, aes(col = dim)) +
  facet_grid(rows = vars(type), scales = 'free_y')+
  labs(
    x = 'Time (ms)',
    y = 'TT',
    title = '',
    col = '') +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12))

ann_text <- data.frame(time = c(1940,1976,1992,2008,2028,2084,2124),
                       value = c(-5,20,-5,10,-5,25,-5), lab = c("onset","pvel_to","target","mvel","release","pvel_fro","offset"),
                       type = c(factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Velocity (cm/s)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Velocity (cm/s)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Velocity (cm/s)',levels = c('Position (mm)','Velocity (cm/s)')),
                                factor('Position (mm)',levels = c('Position (mm)','Velocity (cm/s)'))))

p_TT = p2 + geom_vline(xintercept=c(1940,1976,2000,2008,2020,2080,2124), color = 'black', linewidth = 0.5, linetype = 'dashed') + 
  geom_text(data = ann_text, size = 5,
            label = c("onset","pvel_to","target","mvel","release","pvel_fro","offset"))


figure <- ggarrange(p_TB, p_TT, nrow = 2, ncol = 1)

annotate_figure(figure,
                bottom = text_grob("Time (msec)", 
                                   color = "black", face = "bold", size = 16))


##################################Figure 5##################################
#Top panel: C1 opening sitffness
scatterplot <- function(data, title){
  plot1 <- ggplot(data, aes(x = recC1stiffness2, y = C2onset_C1offset)) +
    geom_point(aes(color = cluster), shape = 20, size = 2, alpha = 0.7) +
    geom_smooth(aes(color = cluster),
                method = "gam", formula = y ~ 0 + x,
                linewidth = 0.75,
                fullrange = F, se = T) +
    geom_rug(aes(color = cluster)) +
    labs(x = 'Reciprocal of C1 opening stiffness (sec)',
         y = 'Time from C2 onset to C1 offset (ms)',
         col = 'Clusters',
         title = title) +
    theme_pubr()+
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=12),
          legend.position = "right",
          legend.direction = "vertical",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
}

plotDE <- scatterplot(DE, 'German')
plotSP <- scatterplot(SP, 'Spanish')
plotSG <- scatterplot(SG, 'Swiss German') 
figure <- ggarrange(plotDE, plotSP,plotSG, nrow = 1, ncol = 3)

annotate_figure(figure,
                bottom = text_grob("Reciprocal of C1 opening stiffness (sec)", 
                                   color = "black", face = "bold", size = 16),
                left = text_grob("Time from C2 onset to C1 offset (sec)", 
                                 color = "black", rot = 90,face = "bold", size = 16))

#Bottom panel: C2 closing sitffness
scatterplot <- function(data, title){
  plot1 <- ggplot(data, aes(x = recC2stiffness1, y = C2onset_C1offset)) +
    geom_point(aes(color = cluster), shape = 20, size = 2, alpha = 0.7) +
    geom_smooth(aes(color = cluster),
                method = "gam", formula = y ~ 0 + x,
                linewidth = 0.75,
                fullrange = F, se = T) +
    geom_rug(aes(color = cluster)) +
    labs(x = 'Reciprocal of C2 closing stiffness (sec)',
         y = 'Time from C2 onset to C1 offset (ms)',
         col = 'Clusters',
         title = title) +
    theme_pubr()+
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=12),
          legend.position = "right",
          legend.direction = "vertical",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
}

plotDE <- scatterplot(DE, 'German')
plotSP <- scatterplot(SP, 'Spanish')
plotSG <- scatterplot(SG, 'Swiss German') 
figure <- ggarrange(plotDE, plotSP,plotSG, nrow = 1, ncol = 3)

annotate_figure(figure,
                bottom = text_grob("Reciprocal of C2 closing stiffness (sec)", 
                                   color = "black", face = "bold", size = 16),
                left = text_grob("Time from C2 onset to C1 offset (sec)", 
                                 color = "black", rot = 90,face = "bold", size = 16))


##################################Table 3##################################
m1_DE_C1 = lmer(C2onset_C1offset~ recC1stiffness2+(1|subj), DE)
m1_DE_C2 = lmer(C2onset_C1offset~ recC2stiffness1+(1|subj), DE)
m1_SP_C1 = lmer(C2onset_C1offset~ recC1stiffness2+(1|subj), SP)
m1_SP_C2 = lmer(C2onset_C1offset~ recC2stiffness1+(1|subj), SP)
m1_SG_C1 = lmer(C2onset_C1offset~ recC1stiffness2+(1|subj), SG)
m1_SG_C2 = lmer(C2onset_C1offset~ recC2stiffness1+(1|subj), SG)

summary(m1_DE_C1) #resutls of linear-miexed models
get_gof(m1_DE_C1) #R-squred value


##################################Figure 7##################################
scatterplot <- function(data, title){
  plot1 <- ggplot(data, aes(x = recstiffmean, y = C2onset_C1offset)) +
    geom_point(aes(color = cluster), shape = 20, size = 2, alpha = 0.7) +
    geom_smooth(aes(color = cluster),
                method = "gam", formula = y ~ 0 + x,
                linewidth = 0.75,
                fullrange = F, se = T) +
    geom_rug(aes(color = cluster)) +
    labs(x = 'Reciprocal of the harmonic stiffness mean (sec)',
         y = 'Time from C2 onset to C1 offset (ms)',
         col = 'Clusters',
         title = title) +
    theme_pubr()+
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=12),
          legend.position = "right",
          legend.direction = "vertical",
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
}

plotDE <- scatterplot(DE, 'German')
plotSP <- scatterplot(SP, 'Spanish')
plotSG <- scatterplot(SG, 'Swiss German') 
figure <- ggarrange(plotDE, plotSP,plotSG, nrow = 1, ncol = 3)

annotate_figure(figure,
                bottom = text_grob("Reciprocal of the harmonic stiffness mean (sec)", 
                                   color = "black", face = "bold", size = 16),
                left = text_grob("Time from C2 onset to C1 offset (sec)", 
                                 color = "black", rot = 90,face = "bold", size = 16))


##################################Table 4##################################
#German
m2_DE_bl = lm(C2onset_C1offset~ 0+recstiffmean, DE[DE$cluster=='bl',])
m2_DE_pl = lm(C2onset_C1offset~ 0+recstiffmean, DE[DE$cluster=='pl',])
m2_DE_gl = lm(C2onset_C1offset~ 0+recstiffmean, DE[DE$cluster=='gl',])
m2_DE_kl = lm(C2onset_C1offset~ 0+recstiffmean, DE[DE$cluster=='kl',])
get_gof(m2_DE_bl)
#Spanish
m2_SP_bl = lm(C2onset_C1offset~ 0+recstiffmean, SP[SP$cluster=='bl',])
m2_SP_pl = lm(C2onset_C1offset~ 0+recstiffmean, SP[SP$cluster=='pl',])
m2_SP_gl = lm(C2onset_C1offset~ 0+recstiffmean, SP[SP$cluster=='gl',])
m2_SP_kl = lm(C2onset_C1offset~ 0+recstiffmean, SP[SP$cluster=='kl',])
get_gof(m2_SP_bl)
#Swiss German
m2_SG_bl = lm(C2onset_C1offset~ 0+recstiffmean, SG[SG$cluster=='bl',])
m2_SG_pl = lm(C2onset_C1offset~ 0+recstiffmean, SG[SG$cluster=='pl',])
m2_SG_gl = lm(C2onset_C1offset~ 0+recstiffmean, SG[SG$cluster=='gl',])
m2_SG_kl = lm(C2onset_C1offset~ 0+recstiffmean, SG[SG$cluster=='kl',])
get_gof(m2_SG_bl)


##################################Table 5##################################
data$voicing <- as.factor(data$voicing)
data$articulator <- as.factor(data$articulator)

contrasts(data$voicing) = contr.sum(2)
contrasts(data$articulator) = contr.sum(2)

DE = subset(data, language=="German")
SP = subset(data, language=="Spanish")
SG = subset(data, language=="Swiss German")

m3_DE = lmer(C2onset_C1offset~ 0+recstiffmean + recstiffmean:(voicing+articulator)
             + (0+recstiffmean|subj), DE)
m3_SP = lmer(C2onset_C1offset~ 0+recstiffmean + recstiffmean:(voicing+articulator)
             + (0+recstiffmean|subj), SP)
m3_SG = lmer(C2onset_C1offset~ 0+recstiffmean + recstiffmean:(voicing+articulator)
          + (0+recstiffmean|subj), SG)
summary(m3_DE)

##################################Figure 8##################################
em_art = emtrends(m2, ~recstiffmean:articulator, var='recstiffmean')
pairs(em_art)
em_voice = emtrends(m2, ~recstiffmean:voicing, var='recstiffmean')
pairs(em_voice)

stats_voicing = data.frame(language = c('German', 'German', 'Spanish', 'Spanish',
                                        'Swiss German', 'Swiss German'),
                           slope = c(2.438, 2.188, 2.247, 2.123, 2.421, 2.319),
                           se = c(0.016, 0.016, 0.013, 0.013, 0.014, 0.014),
                           grp = rep(c('voiced', 'voiceless'), 3))

stats_place = data.frame(language = c('German', 'German', 'Spanish', 'Spanish',
                                      'Swiss German', 'Swiss German'),
                         slope = c(2.322, 2.304, 2.059, 2.311, 2.513, 2.227),
                         se = c(0.016, 0.016, 0.017, 0.017, 0.015, 0.015),
                         grp = rep(c('velar', 'labial'), 3))

stats_voicing$language = factor(stats_voicing$language, levels = c('German', 'Spanish', 'Swiss German'))
stats_place$language = factor(stats_voicing$language, levels = c('German', 'Spanish', 'Swiss German'))


res_plot = function(data, legendtitle){
  ggplot(data, aes(x=language, y=slope, group=grp, color=grp)) +
    geom_pointrange(aes(ymin=slope-2*se, ymax=slope+2*se),
                    position = position_dodge2(width = 0.25))+
    labs(x = 'Languages',
         y = 'Slope',
         col = legendtitle,
         title = paste('Interactions with', legendtitle))+
    theme_light()+
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=11),
          axis.text.x = element_text(angle = 0),
          legend.direction = "vertical",
          legend.position = c(0.8,0.15),
          legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill='transparent'),
          plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'))
}

plot_voi = res_plot(stats_voicing, 'C1 voicing')
plot_pl = res_plot(stats_place, 'C1 place')+ theme(axis.title.y = element_blank())

ggarrange(plot_voi, plot_pl, nrow = 1, ncol = 2)
