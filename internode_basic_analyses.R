#internode_basic_analyses 


library(plyr)
library(ggplot2)
library(reshape2)
library(ggbiplot)


#initial visualization 
#some trait derivation 


#load clean internode dataset 
load("internode_data_clean.Rdata")



#trim dataset to relevant columns 
clean<-clean5[,c("Label","subplot_id","genotype","treatment","rep","tiller","sub","object","count_up","length","width","ratio")]

#order dataset 
clean<-clean[order(clean$subplot_id, clean$tiller, clean$sub, clean$count_up),]

#derive some new traits (on per observation basis)
#total stem height 
#total vegetative height 
#peduncle/total percent 
#peduncle:vegetative ratio 
#segment count 

#calculate summed lengths (total and vegetative)
sumtot<-aggregate(length~subplot_id+tiller+sub, data=clean, sum, na.rm=TRUE)
sumveg<-aggregate(length~subplot_id+tiller+sub, data=subset(clean, object=="internode"), sum, na.rm=TRUE)
colnames(sumtot)[4]<-"sum_total"
colnames(sumveg)[4]<-"sum_vegetative"
sumboth<-merge(sumtot, sumveg, by=c("subplot_id", "tiller", "sub"))

clean1<-merge(clean, sumboth, by=c("subplot_id","tiller","sub"))

#calculate peduncle specific traits 
ped<-subset(clean1, object=="peduncle")
ped$peduncle_total_percent<-ped$length/ped$sum_total
ped$peduncle_vegetative_ratio<-ped$length/ped$sum_vegetative 

colnames(ped)[9]<-"internode_count"
colnames(ped)[10]<-"peduncle_length"
colnames(ped)[11]<-"peduncle_width"
colnames(ped)[12]<-"peduncle_wl_ratio"




#Basal vs Distal split 
#break vegetative internodes into basal and distal portions using a ratio threshold (semi arbitrary at this point, based off my eye)
clean1$portion[clean1$ratio<0.05]<-"distal"
clean1$portion[is.na(clean1$portion)]<-"basal"

#Basal, distal trait calculation 
basal<-ddply(subset(clean1, portion=="basal"), c("subplot_id","tiller","sub"), summarise, basal_length_sum=sum(length), basal_count=length(count_up), basal_length_mean=mean(length), basal_length_max=max(length), basal_length_min=min(length), basal_width_mean=mean(width))
distal<-ddply(subset(clean1, portion=="distal"), c("subplot_id","tiller","sub"), summarise, distal_length_sum=sum(length), distal_count=length(count_up), distal_length_mean=mean(length), distal_length_max=max(length), distal_length_min=min(length), distal_width_mean=mean(width))


#join basal, distal, and peduncle data 
basal_distal<-merge(basal, distal, by=c("subplot_id","tiller","sub"), all=TRUE)
bdp<-merge(basal_distal, ped, by=c("subplot_id", "tiller", "sub"), all=TRUE)

bdp<-bdp[c(1:15,17,18,21:28)]
#derive some new traits 
bdp$basal_total_percent<-bdp$basal_length_sum/bdp$sum_total
bdp$distal_total_percent<-bdp$distal_length_sum/bdp$sum_total

#swing long to calculate genotype*treatment*tiller averages 
bdp_long<-melt(bdp, id.vars=c("subplot_id","tiller","sub","treatment","genotype"), variable.name="trait", value.name="data")
bdp_gtt<-ddply(bdp_long, c("genotype", "treatment", "tiller", "trait"), summarise, mean=mean(data))





#P:D:B visualization 
#absolute vs percent of total vs treatment response? 
#stacked bar chart? 




#internode trait PCA (trait reduction)
#wide by trait and still long by genotype, treatment, and tiller 
bdp_wide<-dcast(bdp_gtt, genotype+treatment+tiller~trait, value.var="mean")

#just look at culm for now
bdp_wide<-subset(bdp_wide, tiller=="culm")

#complete cases (in this iteration, ill lose ~80% of the data due to basal threshold)
bdp_wide1<-bdp_wide[complete.cases(bdp_wide),]


#pull out treatment column 
trt<-bdp_wide1[,2]



#join genotype and treatment
bdp_wide1$name<-paste(bdp_wide1$genotype, bdp_wide1$treatment, sep="_")
rownames(bdp_wide1)<-bdp_wide1$name
bdp_wide1$name<-NULL

bdp_wide2<-bdp_wide1[,c(4:25)]



pca_try<-prcomp(bdp_wide2, center=TRUE, scale.=TRUE)


pca_try$sdev
pca_try$rotation
pca_try$x



plot(pca_try, type="l")



plot(pca_try$x[,1:2])
biplot(pca_try)


pca_try_table<-rbind(pca_try$rotation, pca_try$sdev)


ggbiplot(pca_try, choices=c(1,2),obs.scale=1, var.scale=1, groups=trt, ellipse=TRUE)

ggbiplot(pca_try, choices=c(1,3),obs.scale=1, var.scale=1, groups=trt, ellipse=TRUE)

ggbiplot(pca_try, choices=c(2,3),obs.scale=1, var.scale=1, groups=trt, ellipse=TRUE)















#correlations 








