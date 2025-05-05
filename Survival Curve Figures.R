setwd('Year 4/Diss/Statistical Analysis/Data Collection')

Block1<-read.csv('Block1_Survival.csv',header=T)
View(Block1)
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


#----# FINDING MEANS FOR EACH #----#

SE<-function(x,n){
  sd(x)/sqrt(n)
} 

survival_means<-reshape%>% group_by(new_treat,species,variable)%>% #could use na.omit(shape)%>% here instead of na.rm=T in the below line but creates resurrections
  dplyr::summarise(mean = mean(value,na.rm=T),#can use na.rm=T here but creates gaps in the survival plots
            upr = mean(value)+SE(x=value,n=length(value)),
            lwr = mean(value)-SE(x=value,n=length(value)))

filtered_data <- survival_means %>%
  filter(!species %in% c("dsub", "dsuc", "dneb"))

aggregate<-aggregate(value~species+new_treat+variable, data=reshape, FUN='mean')#COULD ALSO HAVE USED THIS CODE!!!
aggregate$mortality<-paste(1-aggregate$value)
str(aggregate)
aggregate$mortality<-as.numeric(aggregate$mortality)

#means
mean(aggregate$mortality[aggregate$variable=='d20'&aggregate$new_treat=='systemic dcv'])
mean(aggregate$mortality[aggregate$variable=='d20' & aggregate$new_treat=='oral dcv'])


#----# PLOTTING FIGURE #----#

library(ggplot2)

classic_colors <- c("#FF6363", "#00bbff",'#530000',"#0a11db") ## set colour
xlabels<-c('0','5','10','15','20')
figurecombine <- ggplot(data=na.omit(filtered_data), aes(x = variable, y = mean, colour = new_treat))+
  geom_point(aes(shape=new_treat), size = 1)+ ## picture parameter
  geom_line(aes(group = new_treat), lwd = 0.6)+
  scale_x_discrete(breaks = c('d0','d5','d10','d15','d20'),label=xlabels)+
  scale_color_manual(values = classic_colors) +
  labs(col='Treatment',shape='Treatment',)+
  facet_wrap(~species, ncol= 4, labeller=labeller(species=
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
                                           "dsal"="D.saltans",
                                           "dbif"="D.bifasciata",
                                           "dpara"="D.paramelanica",
                                           "dpse"="D.pseudoobscura",
                                           "dsim"="D.simulans",
                                           "dtei"="D.teissieri",
                                           "dvir"="D.virilis",
                                           "sleb"="S.lebanonensis",
                                           "spat"="S.pattersoni",
                                           "ztub"="Z.tuberculatus")))+
  ylab("Proportion of Flies Alive")+
  xlab("Days Post Infection")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 11),
        strip.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.text = element_text(size = 10,face = "italic"),
        legend.title = element_text(size = 11),
        legend.position = "right",
        legend.margin = margin(t = -10))

figurecombine


ggsave("SurvivalCurves.pdf", plot = figurecombine, height = 10, width = 12,dpi = 600 )

