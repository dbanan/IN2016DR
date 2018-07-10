#internode_geno_geo

#work with full size map (trimmed down to 2016 genotypes)

geno<-unique(design_16dr$Genotype)


snp_map<-read.csv("~/Desktop/setaria_10k_phylo.snp.csv", header=T, stringsAsFactors=F)
snp_map<-read.csv("~/Desktop/2016_setaria_mRNA_hapmap.csv", header=T, stringsAsFactors=F)



snp_sub<-snp_map[,colnames(snp_map)%in%geno]

d_snp<-dist(snp_sub)



d_snp<-dist(snp_map)
hc_snp<-hclust(d_snp)

dd <- dist(scale(USArrests), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")


dd<-dist(scale(snp_map), method="euclidean")






#work with filtered map 
snp_012<-read.table("~/Desktop/snp_100k/snp_100k_phylo.012", header=F, sep="\t", stringsAsFactors=FALSE)
snp_012<-snp_012[-1]

snp_pos<-read.table("~/Desktop/snp_100k/snp_100k_phylo.012.pos", header=F, sep="\t", stringsAsFactors=FALSE)
snp_pos$marker<-paste(snp_pos$V1, snp_pos$V2, sep="_")
snp_pos<-as.vector(snp_pos[,3])

snp_indv<-read.table("~/Desktop/snp_100k/snp_100k_phylo.012.indv", header=F, sep="\t", stringsAsFactors=FALSE)
  
colnames(snp_012)<-snp_pos  
rownames(snp_012)<-snp_indv$V1  



snp_map<-t(snp_012)

snp_map<-as.data.frame(snp_map)


d_snp<-dist(snp_012)
hc_snp<-hclust(d_snp)

plot(hc_snp)



