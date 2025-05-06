#----# PHYLOGENETIC MIXED MODEL #----#

#----# LIBRARIES #----#
library(tidyverse)
library(MCMCglmm)
library(ape)
install.packages('ggtree')
library(installr)
library(ggtree)
library(phytools)
library(geiger)
library(lme4)
library(nortest)
library(lmtest)
library(gt)
library(stringr)
library(dplyr)
library(plyr)
library(readr)
library(reshape)

#----# DIRECTORY #----#
setwd('Year 4/Diss/Statistical Analysis/Data Collection')

#----# DATA WRANGLING #----#
Block1<-read.csv('Block1_Survival.csv',header=T)
str(Block1)
Block2<-read.csv('Block2_Survival.csv',header=T)
str(Block2)
Block3<-read.csv('Block3_Survival.csv',header=T)
str(Block3)
Block3<-Block3[,-27]
Block4<-read.csv('Block4_Survival.csv',header=T)
str(Block4)
Block5<-read.csv('Block5_Survival.csv',header=T)
str(Block5)

## Combine Files
survival<-list.files(pattern='*.csv',full.names=T)#creates an object holding all the blocks of survival
allBlocks<-survival%>%
  lapply(read_csv)%>%
  bind_rows()
mean(allBlocks$d0,na.rm=T)

## Reshape
reshape <- cbind(allBlocks[,c(1:5)],allBlocks[,c(6:26)]/allBlocks$d0)#turns numeric survival into proportional survival
reshape$new_treat <- paste(reshape$route,reshape$treatment)#combines route and treatments columns to create a new variable which defines our 4 treatments
str(reshape)
reshape<- melt(reshape, id = c("block","vial_id","species","treatment","route","new_treat"))#turns dataframe into long formatting - it leaves alone those columns defined under id but turns the remaining columns (in this case d0:d20) into a single column referred to as 'variable'
reshape[,c(3:6)] <- lapply(reshape[,c(3:6)],function(x)as.factor(x))#???

## filter for species and columns of use
filtered_mortality <- reshape %>%
  filter(!species %in% c("dsub", "dsuc", "dneb"))

pmm1<-aggregate(value~species+new_treat, data=filtered_mortality, FUN='mean')
colnames(pmm1)[1]<-'species'

##Load Host Phylogeny
tree.host <- read.tree(file=paste("Tree_Host.nwk", sep = ''))
#tree.host <- read.nexus(file="tree.host.nexus")
plot(tree.host)
tree.host$tip.label <- gsub('\'', '' , tree.host$tip.label)
plot(tree.host)

##merge short names into full species name using second dataset
short.names<-read.csv('species_codes.csv',header=T)
colnames(short.names)[3]<-'species'
pmm1.1<-merge(pmm1,short.names) #creating 0
str(pmm1.1)
pmm1.1$animal <- as.factor(pmm1.1$animal)
pmm1.1$animal <- gsub(" ","", pmm1.1$animal) #gets rid of spaces from my tree animal names 

##Drop un-needed species
cat(sprintf("\nNo. of host species in phylogenetic tree: %s",length(tree.host$tip.label)))
cat(sprintf("\nNo. of host species in experimental data: %s",length(unique(pmm1.1$animal))))

my.tree <- drop.tip(tree.host, sort(tree.host$tip.label[
  as.character(unique(tree.host$tip.label)) %in% pmm1.1$animal == FALSE]))

cat(sprintf("\nNo. of host species in phylogenetic tree: %s",length(my.tree$tip.label)))
cat(sprintf("\nNo. of host species in experimental data: %s",length(unique(pmm1.1$animal))))

sort(unique(pmm1.1$animal)) #sanity check
sort(my.tree$tip.label) # check these are listed the same as above 
sort(my.tree$tip.label)==sort(unique(pmm1.1$animal)) # NEED TO LL BE TRUE 

plot(my.tree)

##CHECK THE SPECIES IN THE PHYLOGENY ARE CORRECT
cat(sprintf("\nHost species not present in phylogenetic tree: %s",
            paste(sort(unique(
              pmm1$animal)[as.character(unique(pmm1$animal) %in%
                                          my.tree$tip.label) == FALSE]), 
              collapse = ", "))) #should come back with nothing if all species in both datafiles are the same 

#----# PHYLOGENETIC MIXED MODEL #----#

##Define Priors
priors <- list(G = list(
  G1 = list(V = diag(4), nu = 4, alpha.mu = rep(0,4), alpha.V = diag(4) * 1000)),
  R = list(V = diag(4), nu = 0.002))
  
##Create/Run Model
is.ultrametric(my.tree)
pmm1.1$mortality <- 1-pmm1.1$value
pmm1.1$mortality <- as.numeric(pmm1.1$mortality)

phylo.model <-MCMCglmm(mortality ~ new_treat, random = ~us(new_treat):animal,
                                      rcov = ~idh(new_treat):units, pedigree = my.tree, prior = priors, data = pmm1.1,
                                      nitt = 130000, thin = 50, burnin = 30000, pr = TRUE)
?MCMCglmm
save(phylo.model,file = "phylo.model.RData")
load("phylo.model.RData")
##model summaries 
summary(phylo.model)
summary(phylo.model$Sol)
summary(phylo.model$VCV)

par(mfrow=c(2,2))
plot(phylo.model) #----# THESE LOOK ODD - THERE SEEMS TO BE A FEW WHICH DO NOT HAVE THE 'HAIRY CATEPILLAR' LOOK #----#
plot(phylo.model$VCV)


#----# CHECKS #----#
##autocorrelation
phylo.model.autocorr <- data.frame(t(autocorr.diag(phylo.model$VCV, lags = c(0:1000))))
print(phylo.model.autocorr)

phylo.model.autocorr$Comparison <- rownames(phylo.model.autocorr)
phylo.model.autocorr <- gather(phylo.model.autocorr, key = "Lag", value = "Autocorrelation", -Comparison)
phylo.model.autocorr$Lag <- as.numeric(str_split(phylo.model.autocorr$Lag, pattern = "\\.", simplify = TRUE)[,2])
phylo.model.autocorr$Comparison <- factor(phylo.model.autocorr$Comparison, levels = unique(phylo.model.autocorr$Comparison))

ggplot(data = phylo.model.autocorr) +
  ggtitle("Autocorrelation of Parameter Estimates") +
  geom_line(mapping = aes(x = Lag, y = Autocorrelation)) +
  geom_point(mapping = aes(x = Lag, y = Autocorrelation)) +
  geom_hline(yintercept = 0.1, color = "red") +
  facet_wrap(~Comparison, ncol = 4)+
  theme_bw()

##heteroscedasticity
phylo.model.residuals <- data.frame(
  predicted = predict(phylo.model, marginal=NULL),
  residuals = (pmm1.1$mortality) - predict(phylo.model, marginal=NULL))

ggplot(phylo.model.residuals, mapping = aes(x = predicted, y = residuals)) +
  ggtitle("Residuals VS Predicted (Heteroscedasticity)") +
  geom_smooth(method = lm, linetype=0) +
  geom_hline(yintercept = 0, color = "red") + 
  geom_point(color = "blue", alpha = 0.25)


phylo.model.residuals$bin <- cut(phylo.model.residuals$predicted, seq(from = min(phylo.model.residuals$predicted),
                                                                                                    to = max(phylo.model.residuals$predicted) + (abs(max(phylo.model.residuals$predicted)) + abs(min(phylo.model.residuals$predicted))/ 20),
                                                                                                    by = (abs(max(phylo.model.residuals$predicted)) + abs(min(phylo.model.residuals$predicted)))/ 20))

ggplot(data=subset(phylo.model.residuals, !is.na(bin)), mapping = aes(x = bin, y = residuals)) +
  ggtitle("Residuals VS Predicted: Boxplot") +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, color = "red") + 
  geom_point(position = position_jitter(w = 0.1, h = 0), color = "blue", alpha = 0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

##Normality
ggplot(phylo.model.residuals, mapping = aes(x = residuals)) +
  ggtitle("Distribution of residuals") +
  geom_histogram(aes(y =..density..)) +
  geom_vline(xintercept = 0, color = "red") + 
  stat_function(fun = dnorm, args = list(mean = mean(phylo.model.residuals$residuals), sd = sd(phylo.model.residuals$residuals)),
                linetype = "dotted", size = 1)

ggplot(phylo.model.residuals, mapping = aes(sample = residuals)) +
  ggtitle("QQ plot of residuals") + 
  stat_qq() + stat_qq_line()

#----# REPEATABILITY MANUAL #----#
##Repeatability

#variance/variance+residual
#both variances will be same column number - residual will be same treatment without interaction with animal

colnames(phylo.model$VCV)

OIrepeat<-phylo.model$VCV[,1]/(phylo.model$VCV[,1]+phylo.model$VCV[,17]) #----# what do these numbers mean #----#
mean(OIrepeat) #0.72
median(OIrepeat) # 0.76
HPDinterval(OIrepeat) #0.30 - 0.99
plot(OIrepeat)

OCrepeat<-phylo.model$VCV[,6]/(phylo.model$VCV[,6]+phylo.model$VCV[,18])
plot(OCrepeat)
mean(OCrepeat) #0.68
median(OCrepeat) # 0.74
HPDinterval(OCrepeat) #0.20 - 0.97

HPDinterval(OIrepeat-OCrepeat)

SIrepeat<-phylo.model$VCV[,11]/(phylo.model$VCV[,11]+phylo.model$VCV[,19])
plot(SIrepeat)
mean(SIrepeat) #0.65
median(SIrepeat)#0.70
HPDinterval(SIrepeat) #0.24 - 0.94

SCrepeat<-phylo.model$VCV[,16]/(phylo.model$VCV[,16]+phylo.model$VCV[,20])
plot(SCrepeat)
mean(SCrepeat) #0.53
median(SCrepeat)#0.56
HPDinterval(SCrepeat) #<0.001 - #0.88

#none cross zero but not very confident - repeatability is high 

#----# CORRELATIONS BTW TREATMENT #----#
colnames(phylo.model$VCV)
#correlations between Oral Pathway
OI_OC<-phylo.model$VCV[,5]/(sqrt(phylo.model$VCV[,1]*phylo.model$VCV[,6]))
plot(OI_OC)
mean(OI_OC) # 0.42
HPDinterval(OI_OC) # -0.13 - 0.93
  #no significant correlation: resluts of an infection are different in an oral pathway if infected with virus

#correlations between Systemic Pathway
SI_SC<-phylo.model$VCV[,12]/(sqrt(phylo.model$VCV[,11]*phylo.model$VCV[,16]))
mean(SI_SC) # -0.04
HPDinterval(SI_SC) # -0.66 - 0.70
  #not significant
   
#correlations between Infections
SI_OI<-phylo.model$VCV[,3]/(sqrt(phylo.model$VCV[,1]*phylo.model$VCV[,11]))
mean(SI_OI) # 0.39
HPDinterval(SI_OI) # -0.15 - 0.86
plot(SI_OI)
  #not significant

#correlations between Control
SC_OC<-phylo.model$VCV[,8]/(sqrt(phylo.model$VCV[,6]*phylo.model$VCV[,16]))
mean(SC_OC) # 0.61
HPDinterval(SC_OC) # 0.04 - 0.99
  #significant

#-----------------------------------------------------------------------------#
#----#  FIGURES #----# 
#-----------------------------------------------------------------------------#
#used to pull package 'ggtree' from Repositories - for futur use>
  #to see available packages go to:
                            #CRAN - https://cran.r-project.org/web/packages/available_packages_by_name.html
                            #CRAN(extras) - https://www.stats.ox.ac.uk/pub/RWin/bin/windows/contrib/
                            #Bioconductors (Bioc) - https://www.bioconductor.org/packages/release/BiocViews.html#___Software
                            #RForge - https://rforge.net/
                            #github (just shows how to download some) - https://github.com/search?l=R&q=R&type=Repositories

setRepositories()
ap<-available.packages()
View(ap)
'ggtree'%in% rownames(ap)

library(lattice)
library(scales)
library(gridExtra)
library(grid)
install.packages('ggtree')
library(ggtree)

#order species in phylogenetic tree - used later to attach bar and errors to correct species 
tree<- ggtree(my.tree, size = 0.25, right = T,ladderize = TRUE) + 
  geom_treescale()+
  geom_tiplab()+
  xlim(0,1) 
tree
plot(my.tree)

?get_taxa_name
?data.frame

phylo.order <- data.frame(c(1:21),get_taxa_name(tree))
names(phylo.order)[1]<-"ID"
names(phylo.order)[2]<- "animal"
phylo.order <- phylo.order[order(phylo.order$ID),]

#calculate mean mortality and standard error 
#set up standard error function
SE<-function(x,n){
  sd(x)/sqrt(n)
}

  #for Oral Control
 fig.OC <- filter(pmm1.1, new_treat == "oral ringers")
  mean(fig.OC$mortality) #6.08%
  fig.OC.means <- fig.OC %>% group_by(animal) %>% 
  dplyr::summarise(Mean = mean(mortality),
            std = SE(x = mortality,n=length(mortality)),
            upr = mean(mortality, na.rm = T)+std,
            lwr = mean(mortality, na.rm = T)-std,
  )
  str(fig.OC.means) #comes up with NAs
  mean(fig.OC.means$Mean) #sanity check - matches
  
  #for Oral DCV
  fig.OI <- filter(pmm1.1, new_treat == "oral dcv")
  mean(fig.OI$mortality) #10.40%
  fig.OI.means <- fig.OI %>% group_by(animal) %>% 
    dplyr::summarise(Mean = mean(mortality),
                     std = SE(x = mortality,n=length(mortality)),
                     upr = mean(mortality, na.rm = T)+std,
                     lwr = mean(mortality, na.rm = T)-std,
    )
  str(fig.OI.means) #comes up with NAs
  mean(fig.OI.means$Mean) #sanity check - matches
  
  #for Systemic Ringers
  fig.SC <- filter(pmm1.1, new_treat == "systemic ringers")
  mean(fig.SC$mortality) #8.41%
  fig.SC.means <- fig.SC %>% group_by(animal) %>% 
    dplyr::summarise(Mean = mean(mortality),
                     std = SE(x = mortality,n=length(mortality)),
                     upr = mean(mortality, na.rm = T)+std,
                     lwr = mean(mortality, na.rm = T)-std,
    )
  str(fig.SC.means) #comes up with NAs
  mean(fig.SC.means$Mean) #sanity check - matches
  
  #for Systemic DCV
  fig.SI <- filter(pmm1.1, new_treat == "systemic dcv")
  mean(fig.SI$mortality) #41.73
  fig.SI.means <- fig.SI %>% group_by(animal) %>% 
    dplyr::summarise(Mean = mean(mortality),
                     std = SE(x = mortality,n=length(mortality)),
                     upr = mean(mortality, na.rm = T)+std,
                     lwr = mean(mortality, na.rm = T)-std,
    )
  str(fig.SI.means) #comes up with NAs
  mean(fig.SI.means$Mean) #sanity check - matches

## fixing treatment groups onto host tree
  #Oral Control
  tree
  ?left_join
  ?left_join.data.frame
  OC.tree <- left_join(fig.OC.means,phylo.order,by = "animal")
  head(OC.tree)
  OC.tree <- OC.tree[order(OC.tree$ID),]
  head(OC.tree)
  OC.tree$animal <- factor(OC.tree$animal,levels = rev(OC.tree$animal))
  factor(OC.tree$animal)
  
  #Oral DCV
  OI.tree <- left_join(fig.OI.means,phylo.order,by = "animal")
  head(OC.tree)
  OI.tree <- OI.tree[order(OI.tree$ID),]
  head(OI.tree)
  OI.tree$animal <- factor(OI.tree$animal,levels = rev(OI.tree$animal))
  factor(OI.tree$animal)
  
  #Systemic control
  SC.tree <- left_join(fig.SC.means,phylo.order,by = "animal")
  head(SC.tree)
  SC.tree <- SC.tree[order(SC.tree$ID),]
  head(SC.tree)
  SC.tree$animal <- factor(SC.tree$animal,levels = rev(SC.tree$animal))
  factor(SC.tree$animal)
  
  #Systemic DCV
  SI.tree <- left_join(fig.SI.means,phylo.order,by = "animal")
  head(SI.tree)
  SI.tree <- SI.tree[order(SI.tree$ID),]
  head(SI.tree)
  SI.tree$animal <- factor(SI.tree$animal,levels = rev(SI.tree$animal))
  factor(SI.tree$animal)
  
  #----# CREATING BAR PLOT THEME #----#
  theme.my.own <- function(){
    theme_bw() + 
      theme(
        text = element_text(size = 10, color = "black"),
        strip.text = element_text(face = "italic", size = 10),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black"),
        panel.grid.minor = element_line(size = 0.5),
        legend.title.align = 0.5,
        legend.position = "bottom",
        legend.margin = margin(c(-20, 5, 5, 5)),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 9, color = "black")
      )
  }
  
  #----# SETTING UP BAR PLOT FOR EACH TREATMENT #----#
  
  #Oral Control
  OCPlot <- ggplot(OC.tree, aes(x = Mean, y = animal, fill = Mean)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(xmin = upr, xmax = lwr),size = 0.25,width = 0.3,color = "#000000",)+
    theme_bw()+
    scale_x_continuous(name = "", breaks = c(0,0.5,1))+
    scale_y_discrete(position = "right")+
    guides(fill=guide_colourbar(barwidth=4, barheight = 1, title.position = "bottom"))+
    facet_grid(cols = vars("Oral Control"))+
    scale_fill_gradientn(colours = c('#1e0848', '#300060',
                                     '#43006a', '#57096e', '#6b116f', '#81176d', '#991d69',
                                     '#b02363', '#ca2d59', '#e03b50', '#ee504a', '#f66b4d',
                                     '#fa8657', '#fca368', '#fcc17d', '#fcdf96', '#fbffb2'),
                         name = expression("Mean Mortality"),
                         breaks = c(0,0.5,1), limits = c(0, 1), oob = squish,
                         labels = c(0,0.5,1))+
    scale_x_continuous(name = "", breaks = c(0,0.5,1), limits = c(0,1), expand = c(0,0),
                       labels = c(0,0.5,1))+
    theme(text = element_text(size = 10, color = "black"),
          strip.text = element_text(face = "italic", size = 10),
          axis.text.y = element_text(face = "italic"),
          plot.margin = unit(c(0,0.1,0,0), "cm"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          panel.border=element_blank(),
          axis.line.y = element_line(color = "white"),
          panel.grid.minor = element_line(size = 0.5),
          legend.position = "bottom",
          legend.title.align=0.5,
          legend.margin=margin(c(-20,5,5,5)),
          legend.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 9, color = "black"),
    )
  
  OCPlot
  
  #Oral DCV
  OIPlot <- ggplot(OI.tree, aes(x = Mean, y = animal, fill = Mean))+
    geom_bar(stat='identity')+
    geom_errorbar(aes(xmin = upr, xmax = lwr),size = 0.25,width = 0.3,color = "#000000" )+
    scale_x_continuous(name = "", breaks = c(0,0.5,1))+
    guides(fill=guide_colourbar(barwidth=4, barheight = 1, title.position = "bottom"))+
    facet_grid(cols = vars("Oral DCV"))+
    scale_fill_gradientn(colours = c( '#1e0848', '#300060',
                                      '#43006a', '#57096e', '#6b116f', '#81176d', '#991d69',
                                      '#b02363', '#ca2d59', '#e03b50', '#ee504a', '#f66b4d',
                                      '#fa8657', '#fca368', '#fcc17d', '#fcdf96', '#fbffb2'),
                         name = expression("Mean Mortality"),
                         breaks = c(0,0.5,1), limits = c(0, 1), oob = squish,
                         labels = c(0,0.5,1))+
    scale_x_continuous(name = "", breaks = c(0,0.5,1), limits = c(0,1), expand = c(0,0),
                       labels = c(0,0.5,1))+
    theme.my.own()
  
  OIPlot
  
  ##Systemic Control
  SCPlot <- ggplot(SC.tree, aes(x = Mean, y = animal, fill = Mean)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(xmin = upr, xmax = lwr),size = 0.25,width = 0.3,color = "#000000" )+
    scale_x_continuous(name = "", breaks = c(0,0.5,1))+
    guides(fill=guide_colourbar(barwidth=4, barheight = 1, title.position = "bottom"))+
    facet_grid(cols = vars("Systemic Control"))+
    scale_fill_gradientn(colours = c( '#1e0848', '#300060',
                                      '#43006a', '#57096e', '#6b116f', '#81176d', '#991d69',
                                      '#b02363', '#ca2d59', '#e03b50', '#ee504a', '#f66b4d',
                                      '#fa8657', '#fca368', '#fcc17d', '#fcdf96', '#fbffb2'),
                         name = expression("Mean Mortality"),
                         breaks = c(0,0.5,1), limits = c(0, 1), oob = squish,
                         labels = c(0,0.5,1))+
    scale_x_continuous(name = "", breaks = c(0,0.5,1), limits = c(0,1), expand = c(0,0),
                       labels = c(0,0.5,1))+
    theme.my.own()
  
  SCPlot
  
  ##Systemic DCV
  SIPlot <- ggplot(SI.tree, aes(x = Mean, y = animal, fill = Mean)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(xmin = upr, xmax = lwr),size = 0.25,width = 0.3,color = "#000000" )+
    scale_x_continuous(name = "", breaks = c(0,0.5,1))+
    scale_y_discrete()+
    guides(fill=guide_colourbar(barwidth=4, barheight = 1, title.position = "bottom"))+
    facet_grid(cols = vars("Systemic DCV"))+
    scale_fill_gradientn(colours = c( '#1e0848', '#300060',
                                      '#43006a', '#57096e', '#6b116f', '#81176d', '#991d69',
                                      '#b02363', '#ca2d59', '#e03b50', '#ee504a', '#f66b4d',
                                      '#fa8657', '#fca368', '#fcc17d', '#fcdf96', '#fbffb2'),
                         name = expression("Mean mortality"),
                         breaks = c(0,0.5,1), limits = c(0, 1), oob = squish,
                         labels = c(0,0.5,1))+
    scale_x_continuous(name = "", breaks = c(0,0.5,1), limits = c(0,1), expand = c(0,0),
                       labels = c(0,0.5,1))+
    theme.my.own()
  
  
  SIPlot
  
  #----# COMBINE ALL MORTALITY PLOTS #----#
  
  windows(width = 8,height = 7.5)
  PhylogenyPlot <- ( SIPlot| OIPlot | SCPlot|OCPlot) #loads all 4 treatments in one image
  PhylogenyPlot
  
  #Combine phylogeny tree and barplot 
  tree<- ggtree(my.tree, size = 0.25, right = T,ladderize = TRUE,) + 
    geom_treescale()+
    theme(legend.position='bottom',legend.text = element_blank())
  tree #loads in host phylogeny
  
  windows(width = 8,height = 7.5)
  vplayout<- function(x,y) viewport(layout.pos.row = x,layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(400000,6))) #organises size of the image
  print(tree, vp = vplayout(23700:366450,1)) #pulls tree and puts into defined layout
  print(PhylogenyPlot, vp = vplayout(10000:400000,2:6)) #pulls treatment result combination and puts them in defined layout size
  
  citation()
  R.version
  