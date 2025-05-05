setwd('Year 4/Diss/Statistical Analysis')

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
Block5$d2<-as.integer(Block5$d2)


install.packages('dplyr')
library(dplyr)
install.packages('plyr')
library(plyr)
install.packages('readr')
library(readr)

#--------------------#COMBINING CSV FILES#-------------------------#

survival<-list.files(pattern='*.csv',full.names=T)#creates an object holding all the blocks of survival
allBlocks<-survival%>%
  lapply(read_csv)%>%
  bind_rows()#combines files by columns 

#--------------------#MERGING AND MELTING FILES FOR PROPORTIONAL DATA#-------------------------#

library(reshape)

df_fig <- cbind(allBlocks[,c(1:5)],allBlocks[,c(6:26)]/allBlocks$d0)#turns numeric survival into proportional survival
df_fig$new_treat <- paste(df_fig$route,df_fig$treatment)#combines route and treatments columns to create a new variable which defines our 4 treatments
str(df_fig)
df_fig <- melt(df_fig, id = c("block","vial_id","species","treatment","route","new_treat"))#turns dataframe into long formatting - it leaves alone those columns defined under id but turns the remaining columns (in this case d0:d20) into a single column referred to as 'variable'
df_fig              
df_fig[,c(3:6)] <- lapply(df_fig[,c(3:6)],function(x)as.factor(x))#???

df_fig$value<-df_fig[1:9576,8]/df_fig[1:456,8]

#--------------------#CREATING FACET FIGURES#----------------------#

df_fig$block<-as.factor(df_fig$block)#otherwise recognises block number as a numeric 
str(df_fig)
#above has to be done before subsetting to ensure that subsets use block as a factor 

Con_Oral<-subset(df_fig,df_fig$treatment=='ringers'&df_fig$route=='oral')#oral control subsetted object
DCV_Oral<-subset(df_fig,df_fig$treatment=='dcv'&df_fig$route=='oral')#oral infected subsetted object
Con_Inj<-subset(df_fig,df_fig$treatment=='ringers'&df_fig$route=='systemic')#systemic control subsetted object
DCV_Inj<-subset(df_fig,df_fig$treatment=='dcv'&df_fig$route=='systemic')#systemic dcv subsetted object

library(ggplot2)

#PLOT FOR ORAL CONTROLS#
fig.Con_Oral<-(ggplot(data=na.omit(Con_Oral),aes(x = variable,y = value,colour= block))+
    geom_point(size = 1)+
    geom_line(aes(group = block))+
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
                                               "dsal"="D.saltans",
                                               "dbif"="D.bifasciata",
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
    scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
    scale_colour_manual(values=c("1"="orange","2"="blue","3"="green","4"="red","5"="pink"),
                        labels=c("1"="Block 1","2"="Block 2","3"="Block 3","4"="Block 4","5"="Block 5"))+
    labs(title ="Oral Control",
         colour="Block",
         shape="Block")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))

fig.Con_Oral
ggsave("Con_oral.pdf", plot = fig.Con_Oral, height = 10, width = 12,dpi = 600 )

#PLOT FOR ORAL DCV#
fig.DCV_Oral<-(ggplot(data=na.omit(DCV_Oral),aes(x = variable,y = value,colour= block))+
    geom_point(size = 1)+
    geom_line(aes(group = block))+
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
                                               "dsal"="D.saltans",
                                               "dbif"="D.bifasciata",
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
    scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
      scale_colour_manual(values=c("1"="orange","2"="blue","3"="green","4"="red","5"="pink"),
                          labels=c("1"="Block 1","2"="Block 2","3"="Block 3","4"="Block 4","5"="Block 5"))+
    labs(title ="Oral DCV",
         colour="Block",
         shape="Block")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))

fig.DCV_Oral
ggsave("DCV_oral.pdf", plot = fig.DCV_Oral, height = 10, width = 12,dpi = 600 )

#PLOT FOR INJECTED CONTROLS#
fig.Con_Inj<-(ggplot(data=na.omit(Con_Inj),aes(x = variable,y = value,colour= block))+
    geom_point(size = 1)+
    geom_line(aes(group = block))+
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
                                               "dsal"="D.saltans",
                                               "dbif"="D.bifasciata",
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
    scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
      scale_colour_manual(values=c("1"="orange","2"="blue","3"="green","4"="red","5"="pink"),
                          labels=c("1"="Block 1","2"="Block 2","3"="Block 3","4"="Block 4","5"="Block 5"))+
    labs(title ="Injected Control",
         colour="Block",
         shape="Block")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))

fig.Con_Inj
ggsave("Con_inj.pdf", plot = fig.Con_Inj, height = 10, width = 12,dpi = 600 )

#PLOT FOR INJECTED DCV#
fig.DCV_Inj<-(ggplot(data=na.omit(DCV_Inj),aes(x = variable,y = value,colour= block))+
    geom_point(size = 1)+
    geom_line(aes(group = block))+
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
                                               "dsal"="D.saltans",
                                               "dbif"="D.bifasciata",
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
    scale_x_discrete(breaks=c("d0","d5","d10","d15","d20"))+
      scale_colour_manual(values=c("1"="orange","2"="blue","3"="green","4"="red","5"="pink"),
                          labels=c("1"="Block 1","2"="Block 2","3"="Block 3","4"="Block 4","5"="Block 5"))+
    labs(title ="Injected DCV",
         colour="Block",
         shape="Block")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
    ))

fig.DCV_Inj
ggsave("DCV_inj.pdf", plot = fig.DCV_Inj, height = 10, width = 12,dpi = 600 )




