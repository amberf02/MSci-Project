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

## mean mortality on d20
mean(1-filtered_mortality$value[filtered_mortality$new_treat=='oral dcv'&filtered_mortality$variable=='d20'])
mean(1-filtered_mortality$value[filtered_mortality$new_treat=='systemic dcv'&filtered_mortality$variable=='d20'])

pmm1<-aggregate(value~species+new_treat+block, data=filtered_mortality, FUN='mean')
colnames(pmm1)[1]<-'species'
pmm1$mortality<-1-pmm1$value
pmm1<-pmm1[,-4]


##Load Host Phylogeny
tree.host <- read.tree(file=paste("mix_tree.newick", sep = ''))
#tree.host <- read.nexus(file="tree.host.nexus")
plot(tree.host)
tree.host$tip.label <- gsub('\'', '' , tree.host$tip.label)
plot(tree.host)

##merge short names into full species name using second dataset
short.names<-read.csv('species_codes.csv',header=T)
colnames(short.names)[3]<-'species'
pmm2<-merge(pmm1,short.names) #creating 0
str(pmm2)
pmm2$animal <- as.factor(pmm2$animal)
pmm2$animal <- gsub(" ","", pmm2$animal) #gets rid of spaces from my tree animal names 
pmm2$animal[pmm2$animal=='D.teisseri']<-'D.teissieri'

#---------------------------------------------------------------------------------------------#
#----# CORRECTING FOR HOMOSCEDASTICITY #----#

oral_dcv<-subset(pmm2,pmm2$new_treat=='oral dcv')
oral_ringers<-subset(pmm2,pmm2$new_treat=='oral ringers')
systemic_dcv<-subset(pmm2,pmm2$new_treat=='systemic dcv')
systemic_ringers<-subset(pmm2,pmm2$new_treat=='systemic ringers')

oral_route <- merge(oral_dcv, oral_ringers, by = c("species","speciesnumber","animal","block"), all = TRUE)
oral_route$mort.correct <- (oral_route$mortality.x-oral_route$mortality.y)/(1-oral_route$mortality.y)
oral_route$mort.correct <- ifelse(oral_route$mort.correct<0,0,oral_route$mort.correct)
oral_route$mort.correct <- ifelse(
  is.na(oral_route$mort.correct), 
  oral_route$mortality.x, 
  oral_route$mort.correct
)

oral_route<- oral_route[,-c(7,8)]
names(oral_route)[5]<-"treatment"

systemic_route <- merge(systemic_dcv, systemic_ringers, by = c("species","speciesnumber","animal","block"), all = TRUE)
systemic_route$mort.correct <- (systemic_route$mortality.x-systemic_route$mortality.y)/(1-systemic_route$mortality.y)
systemic_route$mort.correct <- ifelse(systemic_route$mort.correct<0,0,systemic_route$mort.correct) #not super necessary as corrected mrtalities are above 0
systemic_route$mort.correct <- ifelse(
  is.na(systemic_route$mort.correct), 
  systemic_route$mortality.x, 
  systemic_route$mort.correct
)

systemic_route <- systemic_route[,-c(7,8)]
names(systemic_route)[5]<-"treatment"
str(systemic_route)

corrected_mortality<-rbind(systemic_route,oral_route)
str(corrected_mortality)
unique(levels(corrected_mortality$treatment))
corrected_mortality <- subset(corrected_mortality,treatment == c("oral dcv","systemic dcv"))
corrected_mortality$treatment <- droplevels(corrected_mortality$treatment)
str(corrected_mortality)
names(corrected_mortality)[6] <- "prop.death"
str(corrected_mortality)

corrected_mortality[,c(1:5)] <- lapply(corrected_mortality[,c(1:5)], function(x)as.factor(x))
corrected_mortality[,c(6:7)] <- lapply(corrected_mortality[,c(6:7)], function(x)as.numeric(x))

corrected_mortality <- corrected_mortality[,-1] #removes short name species to replace for latin
corrected_mortality$species <- corrected_mortality$animal
unique(levels(corrected_mortality$animal))
unique(levels(corrected_mortality$treatment)) 
droplevels(corrected_mortality)
unique(levels(corrected_mortality$treatment)) 
unique(levels(corrected_mortality$animal))


#head(merged_data_noringer)
#str(merged_data_noringer)

#tree.host3 <- read.tree(file="tree_36hosts.nwk")


##Drop un-needed species
cat(sprintf("\nNo. of host species in phylogenetic tree: %s",length(tree.host$tip.label)))
cat(sprintf("\nNo. of host species in experimental data: %s",length(unique(corrected_mortality$animal))))

my.tree <- drop.tip(tree.host, sort(tree.host$tip.label[
  as.character(unique(tree.host$tip.label)) %in% corrected_mortality$animal == FALSE]))

cat(sprintf("\nNo. of host species in phylogenetic tree: %s",length(my.tree$tip.label)))
cat(sprintf("\nNo. of host species in experimental data: %s",length(unique(corrected_mortality$animal))))

sort(unique(corrected_mortality$animal)) #sanity check
sort(my.tree$tip.label) # check these are listed the same as above 
sort(my.tree$tip.label)==sort(unique(corrected_mortality$animal)) # NEED TO LL BE TRUE 

plot(my.tree)

##CHECK THE SPECIES IN THE PHYLOGENY ARE CORRECT
cat(sprintf("\nHost species not present in phylogenetic tree: %s",
            paste(sort(unique(
              corrected_mortality$animal)[as.character(unique(corrected_mortality$animal) %in%
                                                         my.tree$tip.label) == FALSE]), 
              collapse = ", "))) #should come back with nothing if all species in both datafiles are the same 


### prior with phylo only
prior.without.species1 <- list(
  G = list(
    G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000)),
  R = list(V = diag(2), nu = 0.002))

##Create/Run Model
is.ultrametric(my.tree)

corrected.phylo.model <-MCMCglmm(mort.correct ~ treatment, random = ~us(treatment):animal,
                                 rcov = ~idh(treatment):units, pedigree = my.tree, prior = prior.without.species1, data = corrected_mortality,
                                 nitt = 13000000, thin = 5000, burnin = 3000000, pr = TRUE)

##model summaries 
summary(corrected.phylo.model)
summary(corrected.phylo.model$Sol)
summary(corrected.phylo.model$VCV)

par(mfrow=c(2,2))
plot(corrected.phylo.model) 
plot(corrected.phylo.model$VCV)

#----# CHECKS #----#
##autocorrelation - all points fall below the red line
phylo.model.autocorr <- data.frame(t(autocorr.diag(corrected.phylo.model$VCV, lags = c(0:1000))))
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
  predicted = predict(corrected.phylo.model, marginal=NULL),
  residuals = (corrected_mortality$mort.correct - predict(corrected.phylo.model, marginal=NULL)))

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

colnames(corrected.phylo.model$VCV)

OIrepeat<-corrected.phylo.model$VCV[,1]/(corrected.phylo.model$VCV[,1]+corrected.phylo.model$VCV[,5]) 
mean(OIrepeat) #0.95
median(OIrepeat) # 0.96
HPDinterval(OIrepeat) #0.91 - 0.99
plot(OIrepeat)

SIrepeat<-corrected.phylo.model$VCV[,4]/(corrected.phylo.model$VCV[,4]+corrected.phylo.model$VCV[,6])
plot(SIrepeat)
mean(SIrepeat) #0.61
median(SIrepeat)#0.63
HPDinterval(SIrepeat) #0.33 - 0.87

#correlations between Infections
SI_OI<-corrected.phylo.model$VCV[,2]/(sqrt(corrected.phylo.model$VCV[,1]*corrected.phylo.model$VCV[,4]))
mean(SI_OI) # 0.55
HPDinterval(SI_OI) # 0.13 - 0.91
plot(SI_OI)
#significant

#-------------------------------------------------------------------------------------------------------#
#----# PLOTTING CORRELATION #----#
SE<-function(x,n){
  sd(x)/sqrt(n)
} 

plotting <- corrected_mortality %>%group_by(treatment,species) %>%
  dplyr::summarise(value = mean(mort.correct,na.rm=T),
                   upr = mean(mort.correct)+SE(x=mort.correct,n=length(mort.correct)),
                   lwr = mean(mort.correct)-SE(x=mort.correct,n=length(mort.correct)))
)


length(corrected_mortality$mort.correct[corrected_mortality$treatment=='systemic dcv'])
length(corrected_mortality$mort.correct[corrected_mortality$treatment=='oral dcv'])


head(plotting)
highlight_species <- c('D.affinis', 'D.pseudoobscura', 'D.simulans') #identifies species with influence
point_cols <- rep(alpha('grey',0.5), length(plotting$species)) #colour for all other species
point_cols[plotting$species %in% highlight_species] <- 'blue' #colour for influential species
plot(plotting$value[plotting$treatment == 'systemic dcv'],
     plotting$value[plotting$treatment == 'oral dcv'],
     pch=19, bty='l', col=point_cols,
     xlab='Injected DCV Mean Proportion Dead',
     ylab='Oral DCV Mean Proportion Dead',
     xlim=c(0,1), ylim=c(0,1),
     cex.axis=1.5,
     cex.lab=1.5,
     cex=1.5)

# regression line
abline(lm(plotting$value[plotting$treatment == 'oral dcv'] ~ plotting$value[plotting$treatment == 'systemic dcv']),
       col='red', lwd=1)

# text
text(x=0.95, y=1.0, 'r = 0.55,', cex=1.2)
text(x=0.86, y=0.9, '95% CI = 0.13, 0.91', cex=1.2)



#-----------------------------------------------------------------------------#
#----# MORTALITY BY PHYLOGENY FIGURE #----# 
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

#for Oral DCV
fig.OI <- filter(corrected_mortality, treatment == "oral dcv")
mean(fig.OI$mort.correct) #10.40%
fig.OI.means <- fig.OI %>% group_by(animal) %>% 
  dplyr::summarise(Mean = mean(mort.correct),
                   std = SE(x = mort.correct,n=length(mort.correct)),
                   upr = mean(mort.correct, na.rm = T)+std,
                   lwr = mean(mort.correct, na.rm = T)-std,
  )
str(fig.OI.means) #comes up with NAs
mean(fig.OI.means$Mean) #sanity check - matches


#for Systemic DCV
fig.SI <- filter(corrected_mortality, treatment == "systemic dcv")
mean(fig.SI$mort.correct) #41.73
fig.SI.means <- fig.SI %>% group_by(animal) %>% 
  dplyr::summarise(Mean = mean(mort.correct),
                   std = SE(x = mort.correct,n=length(mort.correct)),
                   upr = mean(mort.correct, na.rm = T)+std,
                   lwr = mean(mort.correct, na.rm = T)-std,
  )
str(fig.SI.means) #comes up with NAs
mean(fig.SI.means$Mean) #sanity check - matches

## fixing treatment groups onto host tree

#Oral DCV
OI.tree <- left_join(fig.OI.means,phylo.order,by = "animal")
head(OI.tree)
OI.tree <- OI.tree[order(OI.tree$ID),]
head(OI.tree)
OI.tree$animal <- factor(OI.tree$animal,levels = rev(OI.tree$animal))
factor(OI.tree$animal)


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


#Oral DCV
OIPlot <- ggplot(OI.tree, aes(x = Mean, y = animal, fill = Mean))+
  geom_bar(stat='identity') +
  geom_errorbar(aes(xmin = upr, xmax = lwr),size = 0.25,width = 0.3,color = "#000000",)+
  theme_bw()+
  scale_x_continuous(name = "", breaks = c(0,0.5,1))+
  scale_y_discrete(position = "right")+
  guides(fill=guide_colourbar(barwidth=4, barheight = 1, title.position = "bottom"))+
  facet_grid(cols = vars("Oral DCV"))+
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

OIPlot


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
PhylogenyPlot <- ( SIPlot| OIPlot ) #loads all 4 treatments in one image
PhylogenyPlot

#Combine phylogeny tree and barplot 
tree<- ggtree(my.tree, size = 0.25, right = T,ladderize = TRUE,) + 
  geom_treescale()+
  theme(legend.position='bottom',legend.text = element_blank())
tree #loads in host phylogeny

windows(width = 8,height = 7.5)
vplayout<- function(x,y) viewport(layout.pos.row = x,layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(415000,3))) #organises size of the image
print(tree, vp = vplayout(23700:371000,1)) #pulls tree and puts into defined layout
print(PhylogenyPlot, vp = vplayout(6000:415000,2:3)) #pulls treatment result combination and puts them in defined layout size

citation()
R.version



