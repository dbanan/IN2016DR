#internode analysis 180104.R 
#darshi
#1/4/2018
#full data cleanup and analysis 

#packages 
library(tidyr)
library(ggplot2)
library(zoo)
library(segmented)

library(plyr)
library(dplyr)
library(reshape2)
library(corrplot)

install.packages("data.table")
library(data.table)

install.packages("broom")
library(broom)



#INFILES
#design and scan files 
design_16dr<-read.csv("/rsync/box/Setaria/2016 Setaria/16DR_subplot_genotype_list_FINAL.csv", header=T, stringsAsFactors=F, na.strings=".")
scan_16dr<-read.csv("/rsync/box/Darshi work/STORY internodes/16dr harvest scan.csv", header=T, stringsAsFactors=F, na.strings=".")
#format
colnames(design_16dr)<-c("subplot_id","rep","awning","treatment","genotype")
scan_16dr1<-scan_16dr[,c(4,9)]
colnames(scan_16dr1)<-c("subplot_id","scan")

#image measurement output file (raw data)
output<-read.csv("/rsync/box/Darshi work/STORY internodes/raw data/Internodes All Results 171003.csv", header=T, stringsAsFactor=F, na.strings=".")
colnames(output)[1]<-"Obs"

#MISLABELS
#change some mislabeled scan entries ("-" when it should be "_")
output$Label[output$Label=="$aStem_biomass16DR44301-C1.JPG"]<-"$aStem_biomass16DR44301_C1.JPG"
output$Label[output$Label=="$aStem_biomass16DR10201-T2.JPG"]<-"$aStem_biomass16DR10201_T2.JPG"

#MISPLACED 
#the below images have the entire internode set misplaced on a diagonal, subset them out for now 
slant<-c("$aStem_biomass16DR70301_C2.JPG",
         "$aStem_biomass16DR46801_C3.JPG",
         "$aStem_biomass16DR78001_C1.JPG",
         "$aStem_biomass16DR78001_C3.JPG",
         "$aStem_biomass16DR81501_C3.JPG",
         "$aStem_biomass16DR33501_C3.JPG",
         "$aStem_biomass16DR33501_C2.JPG",
         "$aStem_biomass16DR58401_C2.JPG"
)

slant1<-output[output$Label %in% slant,]
output1<-output[!(output$Label %in% slant),]

#all images in 11901 are from early attempts on a blue background
#11901_C1A is the one use instead of 11901_C1 (improperly placed segments)
output1<-output1[!output1$Label=="$aStem_biomass16DR11901_C1.JPG",]
#rename 11901_C1A for later ease of formatting 
output1$Label[output1$Label=="$aStem_biomass16DR11901_C1A.JPG"]<-"$aStem_biomass16DR11901_C1.JPG"

#ARTIFACTS 
#other specific images with artifacts to remove for now (clippers, smears, more paper) and process individually 
later<-c("$aStem_biomass16DR80101_C2.JPG",
         "$aStem_biomass16DR80101_T1.JPG",
         "$aStem_biomass16DR50601_T1.JPG"
)

later1<-output1[output1$Label %in% later,]
output2<-output1[!(output1$Label %in% later),]

#paper labels (these are large and can use an size threshold to id)
output2$artifact[output2$Area>1200]<-1
#"crumbs" (these are bits of stem left on the stage after cleanup, can use a size and position threshold to id)
output2$artifact[output2$X<200&output2$Area<28]<-1

output2$artifact[is.na(output2$artifact)]<-0
artifact<-subset(output2, artifact==1)
output3<-subset(output2, artifact==0)

#OBSERVATION COUNT 
#add a raw observation count for each label 
output3$label_obs<-as.numeric(ave(output3$Obs, output3$Label, FUN=seq_along))

#REFERENCE OBJECTS 
#these are the dimes used to establish pixel:mm 
output3$reference[output3$Feret>17&output3$Feret<20&output3$MinFeret>17&output3$MinFeret<20]<-1
#15501_C3 and 92201_C3's reference objects were not identified because of artifacts next to the dime, the first object is in fact the reference object 
output3$reference[output3$Label=="$aStem_biomass16DR15501_C3.JPG"&output3$label_obs==1]<-1
output3$reference[output3$Label=="$aStem_biomass16DR92201_C3.JPG"&output3$label_obs==1]<-1

#DIAGONALS
#some internode segments were diagonally placed on the camera stage and may not be assessed in the correct order by ImageJ
#use Feret and FeretAngle to ID 
output3$diagonal[output3$Feret>400&output3$FeretAngle>20&output3$FeretAngle<151]<-1
#these are peduncles
diagonal<-subset(output3, diagonal==1)
diagonal1<-diagonal[,c(2)]

#SUBSET data based on number of REFERENCE objects to sum broken peduncles 
#count number of reference objects per Label  
ref_count<-ddply(output3, c("Label"), summarise, ref_N=sum(!is.na(reference)))
#join reference object counts with data
output4<-merge(output3, ref_count, by="Label")

#subset data into two groups: 1 ref vs 2 ref (this ignores the image set with no reference objects)
ref_zero<-subset(output4, ref_N==0)
ref_one<-subset(output4, ref_N==1)
ref_two<-subset(output4, ref_N==2)

#REF_ZERO
#ref_zero...either (1) truly without reference object or (2) reference object not properly identified, didn't meet size threshold criteria 
#swing wide to look at segment length trends 
ref_zero_w<-dcast(ref_zero, Label~label_obs, value.var="Feret")
#all images in 11901 are from early days on a blue background, the bag label artifacts have already been purged from dataset
#however, the images came out fuzzy so artifacts made it through ImageJ, need to flag them out
ref_zero$artifact[ref_zero$label_obs>1 & ref_zero$Feret<7]<-1
ref_zero1<-subset(ref_zero, artifact==0)
#id objects 
#in this set, peduncle should be the max observation for each Label 
ref_one_pedz<-aggregate(label_obs ~ Label, ref_zero1, max)
colnames(ref_one_pedz)<-c("Label", "max")
ref_zero2<-merge(ref_zero1, ref_one_pedz, by=c("Label"))
ref_zero2$object[ref_zero2$max==ref_zero2$label_obs]<-"peduncle"
ref_zero2$object[is.na(ref_zero2$object)]<-"internode"

ref_zero_ready<-subset(ref_zero2, select=c(Label, Obs, Area, Feret, MinFeret, label_obs, object))


#REF_ONE
#either (1) correct, (2) still impaired by artifacts that slipped through threshold, (3) impaired by diagonals 
#assess diagonals first (which are peduncles)
ref_one_p<-ref_one[ref_one$Label %in% diagonal1,]
ref_one_p$object[ref_one_p$diagonal==1]<-"peduncle"
ref_one_p$object[ref_one_p$reference==1]<-"dime"
ref_one_p$object[is.na(ref_one_p$object)]<-"internode"

#assess correct vs artifact 
ref_one_m<-ref_one[!(ref_one$Label %in% diagonal1),]
#separate "ref_one" into Labels in correct order vs not (not includes misplaced and artifact impairment)
ref_one_m$misplaced[ref_one_m$reference!=ref_one_m$label_obs]<-1

#list of Labels with misplaced reference objects
misplaced<-subset(ref_one_m, misplaced==1)
misplaced1<-misplaced[,c(1)]

#subset of ref_one with misplaced reference objects 
ref_one_m1<-ref_one_m[ref_one_m$Label %in% misplaced1,]
ref_one_w<-dcast(ref_one_m1, Label~label_obs, value.var="Feret")
#the first obs before the reference objects in these images appear to all be errant bits of stem, flag as artifacts 
ref_one_m1$artifact[ref_one_m1$label_obs==1]<-1
ref_one_m2<-ref_one_m1[!(ref_one_m1$artifact==1),]

#subset of ref_one with correctly placed reference objects 
ref_one_c<-ref_one_m[!(ref_one_m$Label %in% misplaced1),]
plot(ref_one_c$label_obs, ref_one_c$Feret)
#coarse plot looks OK 

#stack correct images with fixed misplaced 
ref_one_s<-rbind(ref_one_c, ref_one_m2)
#id objects
ref_one_s$object[ref_one_s$reference==1]<-"dime"
#in this set, peduncle should be the max observation for each Label 
ref_one_ped<-aggregate(label_obs ~ Label, ref_one_s, max)
colnames(ref_one_ped)<-c("Label", "max")
ref_one_s1<-merge(ref_one_s, ref_one_ped, by=c("Label"))
ref_one_s1$object[ref_one_s1$max==ref_one_s1$label_obs]<-"peduncle"
ref_one_s1$object[is.na(ref_one_s1$object)]<-"internode"

#stack with cleaned diagonals 
ref_one_s2<-subset(ref_one_s1, select=-c(max, misplaced))
ref_one1<-rbind(ref_one_s2, ref_one_p)
#remove reference objects 
ref_one1<-ref_one1[!ref_one1$object=="dime",]
ref_one_ready<-subset(ref_one1, select=c(Label, Obs, Area, Feret, MinFeret, label_obs, object))


#REF_TWO
#ref_two...second reference object signifies peduncle segments need to be summed
ref_two_w<-dcast(ref_two, Label~label_obs, value.var="Feret")
#a couple images that need special treatment 
special<-c("$aStem_biomass16DR64801_C3.JPG", #peduncle assessed before 2nd ref object
           "$aStem_biomass16DR87701_C3.JPG" #stem scrap next to segment 
)
special1<-ref_two[ref_two$Label %in% special,]
ref_two1<-ref_two[!(ref_two$Label %in% special),]

#for the rest, peduncle sums 
#id reference objects that are not first 
ref_two1$reference[ref_two1$reference!=ref_two1$label_obs]<-2
#id peduncle first step (apply flag to na)
ref_two2<-ddply(ref_two1, .(Label), function(x) replace(x, TRUE, lapply(x, na.locf, na.rm = FALSE)))
#re id reference objects 
ref_two2$object[ref_two2$Feret>17&ref_two2$Feret<20&ref_two2$MinFeret>17&ref_two2$MinFeret<20]<-"dime"
#keep non-reference objects 
ref_two3<-ref_two2[is.na(ref_two2$object),]
#distinguish vegetative internodes from peduncles 
ref_two3$object[ref_two3$reference==1]<-"internode"
ref_two3$object[ref_two3$reference==2]<-"peduncle"
#trim away a few columns containing NA (aggregate cannot sum NA)
ref_two4<-subset(ref_two3, select = -c(diagonal))
#sum by Label and object type (head's up, this results in peduncle label_obs values that will be way high and incorrect, need to address down the line)
ref_two_sum<-aggregate(.~Label+object, ref_two4, sum)
#subset summed peduncles 
ref_two_ped<-subset(ref_two_sum, object=="peduncle")
#remove peduncles from ref_two2 and keep vegetative internodes 
ref_two5<-subset(ref_two4, object=="internode")
#rejoin summed peduncles with vegetative internodes 
ref_two6<-rbind(ref_two5, ref_two_ped)

ref_two_ready<-subset(ref_two6, select=c(Label, Obs, Area, Feret, MinFeret, label_obs, object))


#MEGA-STACK
#stack the cleaned up output parts 
clean<-rbind(ref_zero_ready, ref_one_ready, ref_two_ready)
#now that everything is together....
#break into internodes and peduncles for counts and calculations 
internode<-subset(clean, object=="internode")
peduncle<-subset(clean, object=="peduncle")

#internodes
#order based on "Label" and "label_obs" 
internode<-internode[order(internode$Label, internode$label_obs),]
#re-count observation UP 
internode$count_up<-ave(internode$label_obs, internode$Label, FUN=seq_along)
#reverse order and count DOWN
internode<-internode[order(internode$Label, -internode$count_up),]
internode$count_down<-ave(internode$count_up, internode$Label, FUN=seq_along)

#identify max count 
max<-aggregate(count_up ~ Label, internode, max)
colnames(max)<-c("Label", "max")
internode1<-merge(internode, max, by=c("Label"))
#this reveals that 62 images have only one recognized internode. Many of these are due to basal internodes not meeting the size thresholds used in ImageJ macro
#will have to address this later 

#peduncles 
peduncle1<-merge(peduncle, max, by=c("Label"))
peduncle1$count_up<-peduncle1$max+1
peduncle1$count_down<-peduncle1$max+1

#rejoin peduncles and internodes 
clean1<-rbind(peduncle1, internode1)

#express segment as a percent of max observations, up and down 
clean1$percent_up<-clean1$count_up/clean1$max
clean1$percent_down<-clean1$count_down/clean1$max

#break Label, add experimental design, calculate traits 
#split label into useful information to merge and analyze on (need the \\ to distinguish . as character rather than operator)
clean1$label<-clean1$Label
clean2<-separate(clean1, label, into=c("organ", "campaign", "stem"), sep="_")
clean3<-separate(clean2, stem, into=c("twig", "file"), sep="\\.")
clean3$scan<-paste(clean3$organ, clean3$campaign, sep="_")

#merge to scan and design files 
clean4<-merge(clean3, scan_16dr1, by="scan")
clean5<-merge(clean4, design_16dr, by="subplot_id")

#identify subsampling 
clean5$sub[clean5$twig=="C1"]<-1
clean5$sub[clean5$twig=="C2"]<-2
clean5$sub[clean5$twig=="C3"]<-3
clean5$sub[clean5$twig=="T1"]<-1
clean5$sub[clean5$twig=="T2"]<-2
clean5$sub[clean5$twig=="T3"]<-3

#distinguish culm and tiller  
clean5$tiller[clean5$twig=="C1"]<-"culm"
clean5$tiller[clean5$twig=="C2"]<-"culm"
clean5$tiller[clean5$twig=="C3"]<-"culm"
clean5$tiller[clean5$twig=="T1"]<-"tiller"
clean5$tiller[clean5$twig=="T2"]<-"tiller"
clean5$tiller[clean5$twig=="T3"]<-"tiller"

#calculate basic traits of interest 
clean5$length<-clean5$Feret
clean5$width<-clean5$Area/clean5$Feret
clean5$ratio<-clean5$width/clean5$length

save(clean5, file="internode_data_clean.Rdata")






