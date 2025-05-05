setwd("C:/Users/amber/Documents/Year 4/Diss/Statistical Analysis")

block1<-read.csv("Block1_Survival.csv",header=T)
library(ggplot2);library(reshape);library(dplyr)

##
head(block1)
block1 <- cbind(block1[,c(1:5)],block1[,c(6:26)]/block1$d0)
block1$new_treat <- paste(block1$route,block1$treatment)
str(block1)
block1 <- melt(block1, id = c("block","vial_id","species","treatment","route","new_treat"))
block1               
block1[,c(3:6)] <- lapply(block1[,c(3:6)],function(x)as.factor(x))

block1$value<-block1[1:1848,8]/block1[1:88,8]

(ggplot(data=na.omit(block1),aes(x = variable,y = value,colour= new_treat))+
    geom_point(size = 1)+
    geom_line(aes(group = new_treat))+
    facet_wrap(~species, labeller=labeller(species=
                                             c("daff"="D.affinis",
                                               "dana"="D.ananassae",
                                               "dari"="D.arizonae",
                                               "dbuz"="D.buzzatii",
                                               "dhyd"="D.hydei",
                                               "dimm"="D.immigrans",
                                               "dlac"="D.lacicola",
                                               "dlit"="D.littoralis",
                                               "dmel"="D.melanogaster",
                                               "dmic"="D.micromelanica",
                                               "dmon"="D.montana",
                                               "dneb"="D.nebulosa",
                                               "dpara"="D.paramelanica",
                                               "dpse"="D.pseudoobscura",
                                               "dsim"="D.simulans",
                                               "dsub"="D.subobscura",
                                               "dsuc"="D.sucinea",
                                               "dtei"="D.teisseri",
                                               "dvir"="D.virilis",
                                               "sleb"="S.labanonensis",
                                               "spat"="S.pattersoni",
                                               "ztub"="Z.tuberculatus")))+
    xlab("Days Post Infection")+
    ylab("Survival Proportion")+
    #scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
    scale_colour_manual(values=c("oral dcv"="red","oral ringers"="blue","systemic dcv"="darkred","systemic ringers"="darkblue"),
                        labels=c("oral dcv"="Oral DCV","oral ringers"="Oral Ringers","systemic dcv"="Injected DCV","systemic ringers"="Injected Ringers"))+
    labs(title ="Block 1",
         colour="Treatment",
         shape="Treatment")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))

xyplot(value~variable|factor(species), data=na.omit(block1),type="b", xlab=list(label="Day post-infection",cex=1.5), ylab=list(label="Proportion of flies alive",cex=1.5), strip=strip.custom(bg="white"), scales=list(x=list(at=(c(1,5,10,15,20)), labels=(c("","5","10","15","20")),cex=1.2), cex=1.2, alternating=1), col=c("red","blue"), par.strip.text=list(cex=1.2), pch=c(19,17), main="Block1", layout=c(4,1))

#################################################################################
install.packages("survival")
install.packages("survminer")

library(survminer)
library(survival)
Surv(new_treat ~ species) # fix this by looking it up

block1$variable <- as.numeric(gsub("d","",block1$variable))
str(block1)
ggsurvplot(block1,aes(value~variable))

####################################################################################

###SYSTEMIC INFECTIONS ONLY#####
block1systemic<-subset(block1,block1$route=="systemic")
(ggplot(data=na.omit(block1systemic),aes(x = variable,y = value,colour= new_treat))+
    geom_point(size = 1)+
    geom_line(aes(group = new_treat))+
    facet_wrap(~species, labeller=labeller(species=
                                             c("daff"="D.affinis",
                                               "dana"="D.ananassae",
                                               "dari"="D.arizonae",
                                               "dbuz"="D.buzzatii",
                                               "dhyd"="D.hydei",
                                               "dimm"="D.immigrans",
                                               "dlac"="D.lacicola",
                                               "dlit"="D.littoralis",
                                               "dmel"="D.melanogaster",
                                               "dmic"="D.micromelanica",
                                               "dmon"="D.montana",
                                               "dneb"="D.nebulosa",
                                               "dpara"="D.paramelanica",
                                               "dpse"="D.pseudoobscura",
                                               "dsim"="D.simulans",
                                               "dsub"="D.subobscura",
                                               "dsuc"="D.sucinea",
                                               "dtei"="D.teisseri",
                                               "dvir"="D.virilis",
                                               "sleb"="S.labanonensis",
                                               "spat"="S.pattersoni",
                                               "ztub"="Z.tuberculatus")))+
    xlab("Days Post Infection")+
    ylab("Survival Proportion")+
    #scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
    scale_colour_manual(values=c("systemic dcv"="darkred","systemic ringers"="darkblue"),
                        labels=c("systemic dcv"="DCV","systemic ringers"="Ringers"))+
    labs(title ="Block 1 - Injected Infection",
         colour="Treatment",
         shape="Treatment")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))
#####################################################################
###ORAL INFECTIONS ###

block1oral<-subset(block1,block1$route=="oral")
(ggplot(data=na.omit(block1oral),aes(x = variable,y = value,colour= new_treat))+
    geom_point(size = 1)+
    geom_line(aes(group = new_treat))+
    facet_wrap(~species, labeller=labeller(species=
                                             c("daff"="D.affinis",
                                               "dana"="D.ananassae",
                                               "dari"="D.arizonae",
                                               "dbuz"="D.buzzatii",
                                               "dhyd"="D.hydei",
                                               "dimm"="D.immigrans",
                                               "dlac"="D.lacicola",
                                               "dlit"="D.littoralis",
                                               "dmel"="D.melanogaster",
                                               "dmic"="D.micromelanica",
                                               "dmon"="D.montana",
                                               "dneb"="D.nebulosa",
                                               "dpara"="D.paramelanica",
                                               "dpse"="D.pseudoobscura",
                                               "dsim"="D.simulans",
                                               "dsub"="D.subobscura",
                                               "dsuc"="D.sucinea",
                                               "dtei"="D.teisseri",
                                               "dvir"="D.virilis",
                                               "sleb"="S.labanonensis",
                                               "spat"="S.pattersoni",
                                               "ztub"="Z.tuberculatus")))+
    xlab("Days Post Infection")+
    ylab("Survival Proportion")+
    #scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
    scale_colour_manual(values=c("oral dcv"="red","oral ringers"="blue"),
                        labels=c("oral dcv"="DCV","oral ringers"="Ringers"))+
    labs(title ="Block 1 - Oral Infection",
         colour="Treatment",
         shape="Treatment")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))
