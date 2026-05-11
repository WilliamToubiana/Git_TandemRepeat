#################################################################
## TR proportion from per-base coverage of indG Illumina reads ##
#################################################################

### Proportion along chromosomes
################################

#######################################
### TR Proportion along chromosomes ###
#######################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRprop<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation. I then assign them a value of 0.
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW



### Proportion between chromosomes
##################################

sum_cov_GW<-aggregate(mergeTRprop$sum_cov_GW~mergeTRprop$chr, FUN=sum)
sum_cov_TR<-aggregate(mergeTRprop$sum_cov_TR~mergeTRprop$chr, FUN=sum)
sum_cov<-cbind(sum_cov_GW, sum_cov_TR)
colnames(sum_cov)<-c("chr", "sum_cov_GW", "chr2", "sum_cov_TR")

sum_cov$propTR<-sum_cov$sum_cov_TR/sum_cov$sum_cov_GW

sum_cov_subset<-subset(sum_cov, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                         | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                         | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

sum_cov_subset$chr <- factor(sum_cov_subset$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))

ggplot(data=sum_cov_subset, aes(x=chr, y=propTR)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(propTR, digits = 3)), vjust=1.6, size=3)+ ylab("proportion of tandem repeats")+
  theme_classic()



#Correlation with chromosome size 
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
Tps_LG_length<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")

sum_cov_chr<-merge(sum_cov, Tps_LG_length, by.x="chr", by.y="scaffolds")

sum_cov_chr$LGs<-gsub("Tps_LRv5b_scf", "chr", sum_cov_chr$chr)

ggplot(data=sum_cov_chr, aes(x=length, y=propTR)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()


#############################################################
## TR proportion from TRF estimates on the genome assembly ##
#############################################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tps$array_length<-(Tps$end-Tps$start)+1

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tps$chr,ranges=IRanges(start=Tps$start,end=Tps$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind_sum<-aggregate(nonoverlapWind$width~nonoverlapWind$seqnames, FUN=sum)


nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tps_LRv5b_scf1" | seqnames=="Tps_LRv5b_scf2" | seqnames=="Tps_LRv5b_scf3" | seqnames=="Tps_LRv5b_scf4" | seqnames=="Tps_LRv5b_scf5"
                           | seqnames=="Tps_LRv5b_scf6" | seqnames=="Tps_LRv5b_scf7" | seqnames=="Tps_LRv5b_scf8" | seqnames=="Tps_LRv5b_scf9" | seqnames=="Tps_LRv5b_scf10"
                           | seqnames=="Tps_LRv5b_scf11" | seqnames=="Tps_LRv5b_scf12")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion

nonoverlapWind_chr_sum$LG_length <- c(253718452, 178265049, 132376913, 102503173, 89066744, 87859127, 83130533, 76885976, 76156695, 72695627, 67869710, 42522306) #chromosome length estimates from the assembly (file.fasta.fai)

colnames (nonoverlapWind_chr_sum) <- c("LGs", "TR_length", "LG_length")
nonoverlapWind_chr_sum$Prop_repeated_region<-nonoverlapWind_chr_sum$TR_length/nonoverlapWind_chr_sum$LG_length

ggplot(data=nonoverlapWind_chr_sum, aes(x=LG_length, y=Prop_repeated_region)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()+ylim(0,0.18)



#########################################################################
## Correlation TR proportion from per-base coverage and  TRF estimates ##
#########################################################################

merge<-merge(nonoverlapWind_chr_sum, sum_cov_subset, by.x="LGs", by.y="chr", all.x =T)
merge$chr<-c("chr1", "chr10", "chr11", "chr12", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9")

ggplot(data=merge, aes(x=Prop_repeated_region, y=propTR)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=chr), vjust=1.6, size=3)+ ylab("TR proportion TRF estimates")+ xlab("TR proportion coverage")+
  geom_smooth(method=lm) +
  theme_classic()

cor.test(merge$Prop_repeated_region, merge$propTR, method = "spearman")
summary(lm(merge$propTR~merge$Prop_repeated_region))
