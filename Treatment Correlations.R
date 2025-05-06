#correlations

setwd('Year 4/Diss/Statistical Analysis/Data Collection')

Block1<-read.csv('Block1_Survival.csv',header=T)
str(Block1)
Block2<-read.csv('Block2_Survival.csv',header=T)
str(Block2)
Block3<-read.csv('Block3_Survival.csv',header=T)
str(Block3)
Block4<-read.csv('Block4_Survival.csv',header=T)
str(Block4)
Block5<-read.csv('Block5_Survival.csv',header=T)
str(Block5)

install.packages('dplyr')
library(dplyr)
install.packages('plyr')
library(plyr)
install.packages('readr')
library(readr)

#----# COMBINING CSV FILES #----#

survival<-list.files(pattern='*.csv',full.names=T)#creates an object holding all the blocks of survival
allBlocks<-survival%>%
  lapply(read_csv)%>%
  bind_rows()

#----# RESHAPING #----#

library(reshape)

reshape <- cbind(allBlocks[,c(1:5)],allBlocks[,c(6:26)]/allBlocks$d0)#turns numeric survival into proportional survival
reshape$new_treat <- paste(reshape$route,reshape$treatment)#combines route and treatments columns to create a new variable which defines our 4 treatments
str(reshape)
reshape<- melt(reshape, id = c("block","vial_id","species","treatment","route","new_treat"))#turns dataframe into long formatting - it leaves alone those columns defined under id but turns the remaining columns (in this case d0:d20) into a single column referred to as 'variable'
reshape[,c(3:6)] <- lapply(reshape[,c(3:6)],function(x)as.factor(x))#???

#----# FINDING MEANS #----#

filtered_mortality <- reshape %>%
  filter(!species %in% c("dsub", "dsuc", "dneb"))
  


mean_mortality<-aggregate(1-value~species+new_treat, data=filtered_mortality, FUN='mean')
colnames(mean_mortality)[3]<-'value'


SE<-function(x,n){
  sd(x)/sqrt(n)
} 

attempt <- reshape %>%group_by(new_treat,species) %>%
  dplyr::summarise(value = mean(value,na.rm=T),
    upr = mean(value)+SE(x=value,n=length(value)),
    lwr = mean(value)-SE(x=value,n=length(value)))
  )  


#----# PLOT CORRELATIONS #----#

oral_dcv<-subset(mean_mortality,mean_mortality$new_treat=='oral dcv')
oral_dcv$se<-paste()
oral_rin<-subset(mean_mortality,mean_mortality$new_treat=='oral ringers')
systemic_dcv<-subset(mean_mortality,mean_mortality$new_treat=='systemic dcv')
systemic_rin<-subset(mean_mortality,mean_mortality$new_treat=='systemic ringers')

#ORAL PATHWAY: Control Vs DCV 
plot(oral_rin$value,oral_dcv$value, pch=19, bty='l', col=alpha('grey',0.5), xlab='Oral Control Mean Proportion Dead', ylab='Oral DCV Mean Proportion Dead', xlim=c(0,1),ylim=c(0,1))
title('A',adj=0)
abline(lm(oral_dcv$value~oral_rin$value), col = "red", lwd = 1) #arguments = y~x
text(x=0.95, y=1.0, 'r = 0.42') #arguments = x, y
text(x=0.85, y=0.9, '95% CI = -0.13, 0.93')
?cor.test #arguments =x,y
cor.test(filtered_mortality$value[filtered_mortality$new_treat=='oral dcv'], filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], method='spearman')



#SYSTEMIC PATHWAY: Control Vs. DCV
plot(systemic_rin$value,systemic_dcv$value, pch=19, bty='l', col=alpha('grey',0.5), xlab='Injected Control Mean Proportion Dead', ylab='Injected DCV Mean Proportion Dead', xlim=c(0,1),ylim=c(0,1))
title('B',adj=0)
abline(lm(systemic_dcv$value~systemic_rin$value), col='red', lwd=1)
text(x=0.95, y=1.0, 'r = -0.04')
text(x=0.85,y=0.9, '95% CI = -0.66, 0.70')

#text(filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'], filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'],
     #labels=filtered_mortality$species, cex=0.5, pos=3)

cor.test(filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'], filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'], method='spearman')



#INFECTED: Oral Vs Systemic
plot(systemic_dcv$value,oral_dcv$value, pch=19, bty='l', col=alpha('grey',0.5), xlab='Injected DCV Mean Proportion Dead', ylab='Oral DCV Mean Proportion Dead', xlim=c(0,1),ylim=c(0,1))
title('C',adj=0)
abline(lm(oral_dcv$value~systemic_dcv$value), col='red', lwd=1)
text(x=0.95, y=1.0, 'r = 0.39')
text(x=0.85,y=0.9, '95% CI = -0.15, 0.86')

cor.test(filtered_mortality$value[filtered_mortality$new_treat=='oral dcv'], filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'], method= 'spearman')

text(mean_mortality$value[mean_mortality$new_treat=='systemic dcv'], mean_mortality$value[mean_mortality$new_treat=='oral dcv'],
    labels=mean_mortality$species, cex=0.5, pos=3)



#CONTROLS: Oral Vs Systemic
plot(systemic_rin$value,oral_rin$value, pch=19, bty='l', col=alpha('grey',0.5), xlab='Injected Control Mean Proportion Dead', ylab='Oral Control Mean Proportion Dead', xlim=c(0,1),ylim=c(0,1))
title('D',adj=0)
abline(lm(oral_rin$value~systemic_rin$value), col='red', lwd=1)
text(x=0.95, y=1.0, 'r = 0.61')
text(x=0.85,y=0.9, '95% CI = 0.04, 0.99')

cor.test(filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'], method= 'spearman')

#text(filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'],
 #    labels=filtered_mortality$species, cex=0.5, pos=3)



#RANDOM 1: Oral Control Vs Systemic DCV
plot(filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'],
     pch=19, bty='l', col='red', xlab='Oral Control Mean Proportion Dead', ylab='Systemic DCV Mean Proportion Dead', main= 'Spearmans Rho (p=0.7411, rho=0.07662338)',xlim=c(0,1), ylim=c(0,1))

cor.test(filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'], method= 'spearman')

abline(lm(filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv']~filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], col='black', lwd=2))

text(filtered_mortality$value[filtered_mortality$new_treat=='oral ringers'], filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'],
     labels=filtered_mortality$species, cex=0.5, pos=3)



#RANDOM 2: Systemic Control Vs Oral DCV
plot(filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'], filtered_mortality$value[filtered_mortality$new_treat=='oral dcv'],
     pch=19, bty='l', col='red', xlab='Systemic Control Mean Proportion Dead', ylab='Oral DCV Mean Proportion Dead', main= 'Spearmans Rho (p=0.01167, rho=0.5391361)',xlim=c(0,1), ylim=c(0,1))

cor.test(filtered_mortality$value[filtered_mortality$new_treat=='oral dcv'], filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'], method= 'spearman')

abline(lm(filtered_mortality$value[filtered_mortality$new_treat=='oral dcv']~filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'], col='black', lwd=2))

text(filtered_mortality$value[filtered_mortality$new_treat=='systemic ringers'], filtered_mortality$value[filtered_mortality$new_treat=='oral dcv'],
     labels=filtered_mortality$species, cex=0.5, pos=3)


#----# AGGREGATE FUNCTION LATER FOR PHYLO ANALYSIS 
aggregate(value~species+new_treat+block, data=reshape, FUN='mean') #pmm
#species column title needs to be animal for the pmm





