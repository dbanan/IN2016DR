load("internode_data_clean.Rdata")

#fit curve, keep order up or down in mind 

###
#try this again, but with percent_up instead to try to get more resolution at the distal internodes response to drought 
###


ggplot(clean5, aes(count_down, length, color=treatment, group=interaction(treatment, tiller)))+
  geom_smooth(aes(linetype=tiller), se=TRUE)+
  geom_point(aes(shape=object))+
  facet_wrap(~genotype)





#####


pdf(file="./loess_fit.pdf", onefile=TRUE)

genotypes<-unique(clean5$genotype)

for (i in seq_along(genotypes)){ 
  print(ggplot(subset(clean5, clean5$genotype==genotypes[i]), aes(count_up, length, color=treatment, group=treatment))+
          geom_smooth(se=TRUE)+
          geom_point(aes(shape=object))+
          facet_wrap(~tiller, ncol=2)+
          ggtitle(genotypes[i]))
} 

dev.off()




######

test_geno<-c("TB_0015", "TB_0047", "TB_0159", "TB_0281", "TB_0386", "TB_0418", "TB_0530", "TB_0581")

test_set<-subset(clean5, genotype=="TB_0281")

test_set$line_id<-paste(test_set$rep, test_set$treatment, test_set$sub, sep="_")


ggplot(test_set, 
       aes(count_up, length, color=treatment))+
  geom_line(aes(group=line_id))+
  geom_point()+
  facet_wrap(~tiller)
  
ggplot(test_set, 
       aes(count_up, ratio, color=treatment))+
  geom_line(aes(group=line_id))+
  #geom_smooth(se=FALSE)+
  geom_point()+ 
  facet_wrap(~tiller)

ggplot(clean6, 
       aes(ratio, count_up, color=treatment))+
  #geom_line(aes(group=line_id, alpha=0.2))+
  geom_smooth(se=FALSE)+
  geom_point()+ 
  facet_wrap(~tiller)




ggplot(test_set, 
       aes(count_up, length, color=treatment))+
  geom_point()+
  geom_line(aes(group=line_id))+
  facet_wrap(~tiller)



#####
#line plots for every observation of every genotype

clean6<-clean5
clean6$line_id<-paste(clean6$rep, clean6$treatment, clean6$sub, sep="_")


pdf(file="./lines_by_observation.pdf", onefile=TRUE)

genotypes<-unique(clean6$genotype)

for (i in seq_along(genotypes)){ 
  print(ggplot(subset(clean6, clean6$genotype==genotypes[i]), aes(count_up, length, color=treatment))+
          geom_line(aes(group=line_id))+
          geom_point()+
          facet_wrap(~tiller)+
          ggtitle(genotypes[i]))
} 

dev.off()




####
#break point peek look at all 
pdf(file="./break_point_all.pdf", onefile=TRUE)

genotypes<-unique(clean6$genotype)

for (i in seq_along(genotypes)){
print(ggplot(subset(clean6, clean6$genotype==genotypes[i]), 
       aes(count_up, ratio, color=treatment))+
  geom_line(aes(group=line_id))+
  #geom_smooth(se=FALSE)+
  geom_point()+ 
  facet_wrap(~tiller)+
  ggtitle(genotypes[i]))
}

dev.off()



















#curve fit 
#culms only for now 
culm<-subset(clean5, tiller=="culm")
#ignore B100 and BLANK for now 
culm<-culm[!(culm$genotype %in% c("B100", "BLANK")),]


#sort up
culm<-culm[order(culm$Label, culm$count_up),]

#sort down 
culm<-culm[order(culm$Label, culm$count_down),]


#break wet and dry apart 
#wet culm 
culm_wet<-subset(culm, treatment=="wet")
#dry culm
culm_dry<-subset(culm, treatment=="dry")

#perform loess fit on wet and dry separately 
fits_wet<-dlply(culm_wet, "genotype", function(df) loess(length~percent_down, data=df))
fits_dry<-dlply(culm_dry, "genotype", function(df) loess(length~percent_down, data=df))

###UP
fits_wet<-dlply(culm_wet, "genotype", function(df) loess(length~percent_up, data=df))
fits_dry<-dlply(culm_dry, "genotype", function(df) loess(length~percent_up, data=df))



#dummy positions to predict values on 
my.count<-seq(from=0, to=1.0, by=0.2)

####THIS IS THE WAY TO GO I THINK####
#predict by group using lapply 
#predict using model results by group 

#wet predictions 
predicted_wet<-lapply(fits_wet, function (x) predict(x, newdata=my.count))
#gather results and change to dataframe 
predicted_wet<-do.call(rbind, predicted_wet)
predicted_wet<-as.data.frame(predicted_wet)
#format a bit
colnames(predicted_wet)<-my.count
predicted_wet$genotype<-rownames(predicted_wet)
rownames(predicted_wet)<-NULL
predicted_wet$treatment<-"wet"

#dry predictions 
predicted_dry<-lapply(fits_dry, function (x) predict(x, newdata=my.count))
#gather results and change to dataframe 
predicted_dry<-do.call(rbind, predicted_dry)
predicted_dry<-as.data.frame(predicted_dry)
#format a bit
colnames(predicted_dry)<-my.count
predicted_dry$genotype<-rownames(predicted_dry)
rownames(predicted_dry)<-NULL
predicted_dry$treatment<-"dry"

#join wet and dry 
predicted<-rbind(predicted_wet, predicted_dry)

#melt long
predicted_long<-melt(predicted, id.vars=c("genotype", "treatment"), variable.name="percent_down", value.name="predicted_length")

###UP
predicted_long<-melt(predicted, id.vars=c("genotype", "treatment"), variable.name="percent_up", value.name="predicted_length")


#visualize 
ggplot()+geom_boxplot(data=predicted_long, aes(x=percent_up, y=predicted_length, fill=treatment))

#trim some outliers 
predicted_long$flag[predicted_long$predicted_length<(-50)]<-1
predicted_long$flag[predicted_long$predicted_length>(500)]<-1
predicted_long$flag[is.na(predicted_long$flag)]<-0
predicted_long1<-subset(predicted_long, flag==0)

ggplot()+geom_boxplot(data=predicted_long1, aes(x=percent_up, y=predicted_length, fill=treatment))

ggplot(data=predicted_long1, aes(x=percent_up, y=predicted_length, fill=treatment))+
  geom_boxplot()+geom_jitter()

lotsolines<-ggplot(data=predicted_long1, aes(x=percent_up, y=predicted_length, color=treatment, group=interaction(treatment,genotype)))+
  geom_line(alpha=0.5)
lotsolines
ggsave("/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/lotsolines.png",
       lotsolines,
       width=9.5,
       height=7)



percent_up_20<-ggplot()+geom_boxplot(data=predicted_long1, aes(x=percent_up, y=predicted_length, fill=treatment))
percent_up_20
ggsave("/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/percent_up_20.png",
       percent_up_20,
       width=9.5,
       height=7)

percent_up_10<-ggplot()+geom_boxplot(data=predicted_long1, aes(x=percent_up, y=predicted_length, fill=treatment))
percent_up_10
ggsave("/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/percent_up_10.png",
       percent_up_10,
       width=9.5,
       height=7)

percent_down_10<-ggplot()+geom_boxplot(data=predicted_long1, aes(x=percent_down, y=predicted_length, fill=treatment))
percent_down_10
ggsave("/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/percent_down_10.png",
       percent_down_10,
       width=9.5,
       height=7)

percent_down_20<-ggplot()+geom_boxplot(data=predicted_long1, aes(x=percent_down, y=predicted_length, fill=treatment))
percent_down_20
ggsave("/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/percent_down_20.png",
       percent_down_20,
       width=9.5,
       height=7)


#look at it with stacked bar charts 
test<-subset(predicted_long1, treatment=="dry")

ggplot()+geom_bar(data=test, aes(x=genotype, y=predicted_length, fill=percent_up), stat="identity")

ggplot()+geom_bar(data=test, aes(x=genotype, y=predicted_length, fill=percent_up), position=position_stack(reverse=TRUE), stat="identity")

###PLAY WITH SOME BIOMASS DATA###
biomass<-read.csv("/rsync/box/Darshi work/STORY internodes/16DR weights and traits Alonso.csv", header=T, stringsAsFactors=F, na.strings=".")
biomass1<-biomass[c(2,3,8)]
biomass2<-biomass1[!(biomass1$genotype %in% c("B100", "BLANK")),]
biomass2$per_plant_total_mass<-as.numeric(biomass2$per_plant_total_mass)

biomass_w<-dcast(biomass2, genotype~treatment, value.var="per_plant_total_mass")
biomass_w$abs_diff<-biomass_w$dry-biomass_w$wet
biomass_w$rel_diff<-biomass_w$abs_diff/biomass_w$wet

biomass_w1<-biomass_w[c(1,5)]
colnames(biomass_w1)<-c("genotype", "biomass_rel_diff")

#join with internode predictions 
join<-merge(test, biomass_w1)

join$biomass_rel_diff<-as.factor(join$biomass_rel_diff)

ggplot(data=join, aes(x=biomass_rel_diff, y=predicted_length))+
  geom_point(aes(color=factor(percent_up)))

###FIGURE
predicted_l_vs_BM<-ggplot(data=join, aes(x=biomass_rel_diff, y=predicted_length, color=percent_up))+
  geom_smooth(aes(group=percent_up), method="lm", se=FALSE)+
  geom_point(aes(color=factor(percent_up)))
predicted_l_vs_BM

png(file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/predicted_l_vs_BM.png")
dev.off()


ggsave("/rsync/box/Darshi work/LAB MEETINGS/BANAN 180207 internode loess fit percent/predicted_l_vs_BM.png",
       predicted_l_vs_BM,
       width=9.5,
       height=7)


#stack bar chart sorted by biomass relative difference 
#sort by TM trt response 


culm<-culm[order(culm$Label, culm$count_up),]

join1<-join[order(join$biomass_rel_diff),]
join1$genotype1<-as.factor(join1$genotype)
ggplot()+geom_bar(data=join1, aes(x=genotype1, y=predicted_length, fill=percent_up), position=position_stack(reverse=TRUE), stat="identity")



ggplot()+geom_bar(data=join, aes(x=genotype, y=predicted_length, fill=percent_down), position=position_stack(reverse=TRUE), stat="identity")







###ARCHIVED ATTEMPTS###
#This is the one that did it, using lapply
y.predicted<-lapply(fits_dry, function (x) predict(x, newdata = my.count))
#gather results and treat as data frame 
y.predicted<-do.call(rbind, y.predicted)
predicted_dry<-as.data.frame(y.predicted)
#format a bit 
colnames(predicted_dry)<-my.count
predicted_dry$genotype<-rownames(predicted_dry)
rownames(predicted_dry)<-NULL


#predict and plot by individual genotypes 
pred_dry<-predict(fits_dry$TB_0135, my.count, se=TRUE)
pred_wet<-predict(fits_wet$TB_0135, my.count, se=TRUE)

plot(pred_dry$fit ~ my.count, col="red")
points(pred_wet$fit ~ my.count, col="blue")

fits_dry$A10$x
fits_dry$A10$y


#a method heavy on pipes %>% 
dry_models<-culm_dry %>%
  group_by(genotype) %>%
  do(fit=loess(length~percent_down, data=.))

dry_simulated<-culm_dry %>%
  mutate(percent_down = percent_down+runif(length(percent_down)),
         length = length + rnorm(length(length)))

dry_simulated %>% group_by(genotype) %>% nest() %>% full_join(dry_models)

dry_simulated_predicted <- 
  dry_simulated %>% group_by(genotype) %>% nest() %>%
  full_join(dry_models) %>% 
  group_by(genotype) %>% 
  do(augment(.$fit[[1]], newdata = .$data[[1]]))

dry_simulated_predicted %>% 
  ggplot(aes(percent_down, length)) + 
  geom_point(shape=1) + 
  geom_ribbon(aes(percent_down, 
                  ymin=.fitted-1.96*.se.fit, 
                  ymax=.fitted+1.96*.se.fit),
              alpha=0.5, fill="black") + 
  geom_line(aes(percent_down, .fitted), size=1, color="red") + 
  facet_wrap(~genotype)


#loop help from will
listogenos <- unique(culm_dry$genotype)
#bad
for (i in listogenos){
  pred_will <- predict(fits_dry$i, my.count, se=TRUE)
  print(pred_will)
}

#work in progress
for (i in 1:210){
  a <- print(fits_dry[[i]])
  print(class(a))
}

#playground (what will left)
for (i in 1:210){
  pred_will <- predict(fits_dry[[i]], my.count, se=TRUE)
  print(pred_will$fit)
  print(class(pred_will$fit))
}

#playground (where i try to get what I want)

for (i in 1:210){
  pred_will <- predict(fits_dry[[i]], my.count, se=TRUE)
  pred_y<-print(pred_will$fit)
}



