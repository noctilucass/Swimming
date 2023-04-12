
# Packages used for the analysis
library(tidyverse)
library(lme4)
library(plyr)
library(lmerTest)
library(lme4)
library(car)
library(sjPlot)
library(effects)
library(cowplot)
library(emmeans)

#setting working directory and uploading data
setwd("~/Desktop/ECOLAB!/paper tesis/swimming R/NEW TRIAL/")
data<-read.table("full_mean.csv", sep=",", header=TRUE)

#Tidying data
data1<-data%>%group_by(Temp,Pote,Larva) %>%
  dplyr::summarise(Ist_Vel = mean(IV),
                   Tot_Time= sum(Time), Max_Time=max(Time))%>%rowid_to_column(var='ID')
data1$Temp<-as.factor(data1$Temp)


#Calculating median for each response variable
medIV <- ddply(data1, "Temp", summarise, grp.med=median(Ist_Vel))
head(medIV)
medTime <- ddply(data1, "Temp", summarise, grp.med=median(Tot_Time))
head(medTime)


#plotting as Kernel Densities
dplot_IV<-ggplot(data1, aes(Ist_Vel, color=Temp))+geom_density(kernel="gaussian", bw="")+
  geom_point(aes(x=5.30, y=0),colour="red",shape=17, size=6)+
  geom_point(aes(x=2.06, y=0),colour="orange",shape=17, size=6)+
  geom_point(aes(x=2.47, y=0),colour="blue",shape=17, size=6)+
  theme_classic() + labs(colour="Temperature", x="Instant Velocity (mm/s)", y="Kernel density")+ 
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.1))) +
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20))
dplot_IV

dplot_Time<-ggplot(data1, aes(Tot_Time, color=Temp))+geom_density(kernel="gaussian")+
  geom_point(aes(x=21.60, y=0),color="Red",shape=17, size=6)+
  geom_point(aes(x=51.25, y=0),color="orange",shape=17, size=6)+
  geom_point(aes(x=64.6, y=0),color="blue",shape=17, size=6)+
  theme_classic()+labs(colour="Temperature", x="Total Time Swimming (sec)", y="Kernel density")+ 
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.1))) +
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20))
dplot_Time

#Checking for normality
shapiro.test(data1$Ist_Vel)
hist(data1$Ist_Vel)
shapiro.test(data1$Tot_Time)
hist(data1$Tot_Time)

#Checking for heteroscedasticity
leveneTest(data1$Ist_Vel, data1$Temp)
leveneTest(data1$Tot_Time,data1$Temp)

#Linear mixed effects model with Log transformed data. Pote is a factor with random effects and larva is nested in Pote 
modIV<-lmer(log(Ist_Vel)~Temp + (1|Pote)+(Pote|Larva), data=data1)

#model diagnostic
plot(modIV,which = 1)
qqnorm(resid(modIV))
qqline(resid(modIV))

#model output
summary(modIV)
anova(modIV)
modTime<-lmer(log(Tot_Time)~Temp+(1|Pote) +(Pote|Larva), data=data1)

#model diagnostic
plot(modTime,which = 1)
qqnorm(resid(modTime))
qqline(resid(modTime))

#model output
summary(modTime)
anova(modTime)

#pairwise comparison using estimated means
emmeans(modIV, pairwise ~ Temp, adjust="bonferroni")
emmeans(modTime, pairwise ~ Temp, adjust="bonferroni")

#Plotting as boxplots with log transformed data
p1<-ggplot(data1, aes(x=Temp, y=log(Ist_Vel), fill=Temp))+geom_boxplot()+ 
  geom_jitter(color="black", size=0.6, alpha=0.3)+
  theme_classic()+labs(x="Temperature", y="IV")+theme(legend.position="none")
p1

p2<-ggplot(data1, aes(x=Temp, y=log(Tot_Time), fill=Temp))+geom_boxplot()+ 
  geom_jitter(color="black", size=0.6, alpha=0.3)+
  theme_classic() + theme(legend.position="none", axis.title = element_text(size = 10)) +
  labs(x="Temperature", y="TT")
p2

#plotting final kernel graphs with insets
plot.with.inset.IV <-
  ggdraw() +
  draw_plot(dplot_IV) +
  draw_plot(p1, x = .4, y = 0.5, width = .4, height = .4)
plot.with.inset.IV

plot.with.inset.Tot_Time <-
  ggdraw() +
  draw_plot(dplot_Time) +
  draw_plot(p2, x = .4, y = 0.5, width = .4, height = .4)
plot.with.inset.Tot_Time
