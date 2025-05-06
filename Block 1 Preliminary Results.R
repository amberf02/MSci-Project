setwd("Year 4/Diss/Data Collection")

block1<-read.csv("Block1_Survival.csv",header=T)

?tapply

systemic<-subset(block1, block1$route=='systemic')
oral<-subset(block1, block1$route=='oral')
OC<-subset(block1, block1$treatment=='ringers')
OI<-subset(block1, block1$treatment=='dcv')
SC<-subset(block1, block1$treatment=='ringers')
SI<-subset(block1, block1$treatment=='dcv')

block1
str(block1)

########### CONTROL TREATMENT SUMMARY 
sum(OC$d0, na.rm=T) #801 tota flies in control treatment
mean(OC$d0, na.rm=T) #mean number per vial is 18.20455
min(tapply(OC$d0,OC$species, mean)) #9.5 - find out from stephen what this is coding for
max(tapply(OC$d0,OC$species, mean)) #20.5 - again find out from stephen what this is coding for 

########## INFECTED TREATMENT SUMMARY
sum(OI$d0, na.rm=T) #798 flies in infected treatment
mean(OI$d0, na.rm=T) #mean number per vial is 18.13636
min(tapply(OI$d0,OI$species, mean)) #10.5
max(tapply(OI$d0,OI$species, mean)) #20.5


###################ORAL INFECTIONS#######################################################

######### CALCULATING PROPORTION ALIVE
?dim #provides the dimensions of the dataframe eg. number of rows then columns
dim(oral) #44. 26
head(oral) 
names(oral)

oral_alive<-oral[,6:26]/oral$d0 #calculates for each specified column and row (here no set row but all columns between column 6 and 26) and divides this by the total number of flies in each row
oral_alive

oral_alive<-cbind(oral[,1:5],oral_alive) #adds our first 5 columns from the orginal dataset (block, vial_id, species, treatment and route) back into the dataset of proportion of living flies in oral treatments
head(oral_alive)

########## RESHAPE DATAFRAME
install.packages('reshape')
library(reshape)

melted_oral<-melt(oral_alive,id=c("vial_id","species","treatment", "route")) #seems to get rid of proportion alive?? Ask Stephen/Hongbo
melted_oral
head(melted_oral)

melted_oral$sp<-as.factor(melted_oral$species)
head(melted_oral)
str(melted_oral)
melted_OI<-subset(melted_oral, melted_oral$treatment=="dcv")
melted_OC<-subset(melted_oral, melted_oral$treatment=="ringers")


library(lattice)

########## BETWEEN VIAL VARIATION
xyplot(value~variable|species, data=melted_OC,type="b", groups=vial_id, main="CONTROL", xlab="day p.i.", ylab="prop alive", auto.key=TRUE, strip=strip.custom(bg="white"), ylim=c(0,1), na.remove=T)
xyplot(value~variable|species, data=melted_OI,type="b", groups=vial_id, main="Virus",  xlab="day p.i.", ylab="prop alive", auto.key=TRUE, strip=strip.custom(bg="white"), ylim=c(0,1))


######### PLOTTING THE CONTROL AND INFECTED GROUPS FOR EACH SPECIES
head(oral_alive)
dim(oral_alive)
#want to take the mean of data~species+trt
oral_alive_mean<-aggregate(oral_alive[,c(6:26)],by=list(species=oral_alive$species, treatment=oral_alive$treatment), FUN=mean)
oral_alive_mean

library(lattice)

########## ANOTHER RESHAPE
library(reshape)
mean.prop.alive_melted.oral<-melt(oral_alive_mean, id=c("species","treatment"))
mean.prop.alive_melted.oral


mean.prop.alive_melted.oral$species<-as.factor(mean.prop.alive_melted.oral$species)
head(mean.prop.alive_melted.oral)

OI2<-subset(mean.prop.alive_melted.oral, mean.prop.alive_melted.oral$treatment=="dcv")
OC2<-subset(mean.prop.alive_melted.oral, mean.prop.alive_melted.oral$treatment=="ringers")

### quick check - gaps where NAs 
xyplot(value~variable|species, data=OC2,type="b", main="Control", xlab="day p.i.", ylab="prop alive", strip=strip.custom(bg="white"), ylim=c(0,1))
xyplot(value~variable|species, data=OI2,type="b", main="Virus", xlab="day p.i.", ylab="prop alive", strip=strip.custom(bg="white"), ylim=c(0,1))

### reorder and create a column for each virus and control 

mean.prop.alive_melted.oral<-mean.prop.alive_melted.oral[with(mean.prop.alive_melted.oral, order(mean.prop.alive_melted.oral$species, mean.prop.alive_melted.oral$treatment, mean.prop.alive_melted.oral$variable)),]

OI3<-subset(mean.prop.alive_melted.oral, mean.prop.alive_melted.oral$treatment=="dcv")
OC3<-subset(mean.prop.alive_melted.oral, mean.prop.alive_melted.oral$treatment=="ringers")

both3<-merge(OI3,OC3,by = c("species","variable"))
head(both3)
xyplot((value.x)+(value.y)~variable|factor(species), data=both3,type="p",main="Virus", xlab="day p.i.", ylab="prop alive", strip=strip.custom(bg="white"))



both3<-both3[with(both3, order(variable)), ]
both3


mfrow=c(2,1)
head(both3)

### this gets rid of the "d" from variable col and turns it into numeric ###
both3$variable <- gsub("[A-Za-z]", "",   both3$variable)
both3$variable<-as.numeric(both3$variable)

### note na.omit used to remove NAs, as using day as numeric now bridges gap between points with line ###

xyplot((value.x)+(value.y)~variable|factor(species), data=na.omit(both3),type="b", xlab=list(label="Day post-infection",cex=1.5), ylab=list(label="Proportion of flies alive",cex=1.5), strip=strip.custom(bg="white"), scales=list(x=list(at=(c(1,5,10,15,20)), labels=(c("","5","10","15","20")),cex=1.2), cex=1.2, alternating=1), col=c("red","blue"), par.strip.text=list(cex=1.2), pch=c(19,17), main="Oral infection", layout=c(4,1)) 
  #####in 8 out of 22 oral infections the infected groups survived longer than the control - need to re-assess why this might be 




################################################################################################


###############SYSTEMIC INFECTIONS####################################


######### CALCULATING PROPORTION ALIVE
?dim #provides the dimensions of the dataframe eg. number of rows then columns
dim(systemic) #44. 26
head(systemic) 
names(systemic)

systemic_alive<-systemic[,6:26]/systemic$total_n #calculates for each specified column and row (here no set row but all columns between column 6 and 26) and divides this by the total number of flies in each row
systemic_alive

systemic_alive<-cbind(systemic[,1:5],systemic_alive) #adds our first 5 columns from the orginal dataset (block, vial_id, species, treatment and route) back into the dataset of proportion of living flies in oral treatments
head(systemic_alive)

########## RESHAPE DATAFRAME
install.packages('reshape')
library(reshape)

melted_systemic<-melt(systemic_alive,id=c("vial_id","species","treatment", "route")) #seems to get rid of proportion alive?? Ask Stephen/Hongbo
melted_systemic
head(melted_systemic)

melted_systemic$sp<-as.factor(melted_systemic$species)
head(melted_systemic)
str(melted_systemic)
melted_SI<-subset(melted_systemic, melted_systemic$treatment=="dvc")
melted_SC<-subset(melted_systemic, melted_systemic$treatment=="ringers")


library(lattice)

########## BETWEEN VIAL VARIATION
xyplot(value~variable|species, data=melted_SC,type="b", groups=vial_id, main="CONTROL", xlab="day p.i.", ylab="prop alive", auto.key=TRUE, strip=strip.custom(bg="white"), ylim=c(0,1), na.remove=T)
xyplot(value~variable|species, data=melted_SI,type="b", groups=vial_id, main="Virus",  xlab="day p.i.", ylab="prop alive", auto.key=TRUE, strip=strip.custom(bg="white"), ylim=c(0,1))


######### PLOTTING THE CONTROL AND INFECTED GROUPS FOR EACH SPECIES
head(systemic_alive)
dim(systemic_alive)
#want to take the mean of data~species+trt
systemic_alive_mean<-aggregate(systemic_alive[,c(6:26)],by=list(species=systemic_alive$species, treatment=systemic_alive$treatment), FUN=mean)
systemic_alive_mean

library(lattice)

########## ANOTHER RESHAPE
library(reshape)
mean.prop.alive_melted.systemic<-melt(systemic_alive_mean, id=c("species","treatment"))
mean.prop.alive_melted.systemic


mean.prop.alive_melted.systemic$species<-as.factor(mean.prop.alive_melted.systemic$species)
head(mean.prop.alive_melted.systemic)

SI2<-subset(mean.prop.alive_melted.systemic, mean.prop.alive_melted.systemic$treatment=="dvc")
SC2<-subset(mean.prop.alive_melted.systemic, mean.prop.alive_melted.systemic$treatment=="ringers")

### quick check - gaps where NAs 
xyplot(value~variable|species, data=SC2,type="b", main="Control", xlab="day p.i.", ylab="prop alive", strip=strip.custom(bg="white"), ylim=c(0,1))
xyplot(value~variable|species, data=SI2,type="b", main="Virus", xlab="day p.i.", ylab="prop alive", strip=strip.custom(bg="white"), ylim=c(0,1))

### reorder and create a column for each virus and control 

mean.prop.alive_melted.systemic<-mean.prop.alive_melted.systemic[with(mean.prop.alive_melted.systemic, order(mean.prop.alive_melted.systemic$species, mean.prop.alive_melted.systemic$treatment, mean.prop.alive_melted.systemic$variable)),]

SI3<-subset(mean.prop.alive_melted.systemic, mean.prop.alive_melted.systemic$treatment=="dvc")
SC3<-subset(mean.prop.alive_melted.systemic, mean.prop.alive_melted.systemic$treatment=="ringers")

both3systemic<-merge(SI3,SC3,by = c("species","variable"))
head(both3systemic)
xyplot((value.x)+(value.y)~variable|factor(species), data=both3systemic,type="p",main="Virus", xlab="day p.i.", ylab="prop alive", strip=strip.custom(bg="white"))



both3systemic<-both3systemic[with(both3systemic, order(variable)), ]
both3systemic


mfrow=c(1,1)
head(both3systemic)

### this gets rid of the "d" from variable col and turns it into numeric ###
both3systemic$variable <- gsub("[A-Za-z]", "",   both3systemic$variable)
both3systemic$variable<-as.numeric(both3systemic$variable)

### note na.omit used to remove NAs, as using day as numeric now bridges gap between points with line ###

xyplot((value.x)+(value.y)~variable|factor(species), data=na.omit(both3systemic),type="b", xlab=list(label="Day post-infection",cex=1.5), ylab=list(label="Proportion of flies alive",cex=1.5), strip=strip.custom(bg="white"), scales=list(x=list(at=(c(1,5,10,15,20)), labels=(c("","5","10","15","20")),cex=1.2), cex=1.2, alternating=1), col=c("red","blue"), par.strip.text=list(cex=1.2), pch=c(19,17), main="systemic infection", layout=c(4,1)) 
#####in 8 out of 22 oral infections the infected groups survived longer than the control - need to re-assess why this might be 

