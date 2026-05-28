#################################################################
## TR proportion from per-base coverage of indG Illumina reads ##
#################################################################

### Proportion between chromosomes
##################################

##Tps
#####

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt

Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$V1, Tps_cov_genome$V3, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$V1, Tps_cov_TR$V3, sep="-")

merge<-merge(Tps_cov_genome, Tps_cov_TR[4:5], by="chm_pos", all.x=T)
merge[is.na(merge)] <- 0 #NA are retrieved when windows are found without TR annotation

#subset only chromosomes
merge2<-subset(merge, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | 
                 V1=="Tps_LRv5b_scf4" |V1=="Tps_LRv5b_scf5" | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | 
                 V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10" | V1=="Tps_LRv5b_scf11" | 
                 V1=="Tps_LRv5b_scf12")

colnames(merge2)=c("chm_pos", "chr", "wind_start", "wind_end", "sum_cov_GW", "sum_cov_TR")


#TR proportion across chromosomes and plot
sum_cov_GW<-aggregate(merge2$sum_cov_GW~merge2$chr, FUN=sum)
sum_cov_TR<-aggregate(merge2$sum_cov_TR~merge2$chr, FUN=sum)
sum_cov<-cbind(sum_cov_GW, sum_cov_TR)
colnames(sum_cov)<-c("chr", "sum_cov_GW", "chr2", "sum_cov_TR")

sum_cov$propTR<-sum_cov$sum_cov_TR/sum_cov$sum_cov_GW
sum_cov$chr <- factor(sum_cov$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))

ggplot(data=sum_cov, aes(x=chr, y=propTR)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(propTR, digits = 3)), vjust=1.6, size=3)+ ylab("proportion of tandem repeats")+
  theme_classic()


#Correlation with chromosome size and plot
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/chromosome_lengths")
Tps_LG_length<-read.table("Tps_chm_size_mtDNAv350.txt", header=F, sep="", quote = "")
colnames(Tps_LG_length)=c("scaffold","length")
Tps_LG_length$chr<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", Tps_LG_length$scaffold)))

sum_cov_chr<-merge(sum_cov, Tps_LG_length, by.x="chr", by.y="scaffold")

ggplot(data=sum_cov_chr, aes(x=length, y=propTR)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=chr), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()





#plot   
# I need to scale the second variable so it fits on the same plot
# Example: divide effective population size to bring it into same range
sum_cov_chr$sum_cov_noTR<-sum_cov_chr$sum_cov_GW-sum_cov_chr$sum_cov_TR
scaling_factor <- max(sum_cov_chr$sum_cov_noTR) / max(sum_cov_chr$sum_cov_TR)

ggplot(sum_cov_chr, aes(x = length)) +
  geom_point(aes(y = sum_cov_noTR), color = "mediumvioletred", size = 3) +
  geom_line(aes(y = sum_cov_noTR), color = "mediumvioletred") +
  geom_point(aes(y = sum_cov_TR * scaling_factor), color = "darkgreen", size = 3, shape = 17) +
  geom_line(aes(y = sum_cov_TR * scaling_factor), color = "darkgreen", linetype = "dashed") +
  scale_y_continuous(
    name = "Genome Coverage without TRs",
    sec.axis = sec_axis(~ . / scaling_factor, name = "TR Coverage")
  ) +
  xlab("Chromosome length") +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "mediumvioletred"),
    axis.text.y.left = element_text(color = "mediumvioletred"),
    axis.title.y.right = element_text(color = "darkgreen"),
    axis.text.y.right = element_text(color = "darkgreen"),
  )




##Tcm
#####

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_other_Timema")

Tcm_cov_genome<-read.table("Tcm-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tcm-filtered_indG_minimap2_GW_coverage.txt genomes/Tcm_chm_size_mtDNAv350_w250000.bed > Tcm-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
Tcm_cov_TR<-read.table("Tcm-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tcm-filtered_indG_minimap2_TR_coverage.txt genomes/Tcm_chm_size_mtDNAv350_w250000.bed > Tcm-filtered_indG_minimap2_TR_sum_coverage_250kb.txt

Tcm_cov_genome$chm_pos<-paste(Tcm_cov_genome$V1, Tcm_cov_genome$V3, sep="-")
Tcm_cov_TR$chm_pos<-paste(Tcm_cov_TR$V1, Tcm_cov_TR$V3, sep="-")

merge<-merge(Tcm_cov_genome, Tcm_cov_TR[4:5], by="chm_pos", all.x=T)
merge[is.na(merge)] <- 0 #NA are retrieved when windows are found without TR annotation

merge$scaff_number<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", merge$V1)))

#subset only chromosomes
merge2<-subset(merge, scaff_number=="1" | scaff_number=="2" | scaff_number=="3" | 
                 scaff_number=="4" |scaff_number=="5" | scaff_number=="6" | scaff_number=="7" | 
                 scaff_number=="8" | scaff_number=="9" | scaff_number=="10" | scaff_number=="11" | 
                 scaff_number=="12")

colnames(merge2)=c("chm_pos", "scaffold", "wind_start", "wind_end", "sum_cov_GW", "sum_cov_TR", "chr")


#TR proportion across chromosomes and plot
sum_cov_GW<-aggregate(merge2$sum_cov_GW~merge2$chr, FUN=sum)
sum_cov_TR<-aggregate(merge2$sum_cov_TR~merge2$chr, FUN=sum)
sum_cov<-cbind(sum_cov_GW, sum_cov_TR)
colnames(sum_cov)<-c("chr", "sum_cov_GW", "chr2", "sum_cov_TR")

sum_cov$propTR<-sum_cov$sum_cov_TR/sum_cov$sum_cov_GW
sum_cov$chr <- as.character(sum_cov$chr)
sum_cov$chr <- factor(sum_cov$chr,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

ggplot(data=sum_cov, aes(x=chr, y=propTR)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(propTR, digits = 3)), vjust=1.6, size=3)+ ylab("proportion of tandem repeats")+
  theme_classic()


#Correlation with chromosome size and plot
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/chromosome_lengths")
Tcm_LG_length<-read.table("Tcm_chm_size_mtDNAv350.txt", header=F, sep="", quote = "")
colnames(Tcm_LG_length)=c("scaffold","length")
Tcm_LG_length$chr<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", Tcm_LG_length$scaffold)))
Tcm_LG_length2<-aggregate(Tcm_LG_length$length~Tcm_LG_length$chr, FUN=sum)
colnames(Tcm_LG_length2)=c("chr","length")

sum_cov_chr<-merge(sum_cov, Tcm_LG_length2, by="chr")

ggplot(data=sum_cov_chr, aes(x=length, y=propTR)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=chr), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()




##Tpa
#####

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_other_Timema")

Tpa_cov_genome<-read.table("Tpa-filtered_indG.1_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tpa-filtered_indG.1_minimap2_GW_coverage.txt genomes/Tpa_chm_size_mtDNAv350_w250000.bed > Tpa-filtered_indG.1_minimap2_GW_sum_coverage_250kb.txt
Tpa_cov_TR<-read.table("Tpa-filtered_indG.1_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tpa-filtered_indG.1_minimap2_TR_coverage.txt genomes/Tpa_chm_size_mtDNAv350_w250000.bed > Tpa-filtered_indG.1_minimap2_TR_sum_coverage_250kb.txt

Tpa_cov_genome$chm_pos<-paste(Tpa_cov_genome$V1, Tpa_cov_genome$V3, sep="-")
Tpa_cov_TR$chm_pos<-paste(Tpa_cov_TR$V1, Tpa_cov_TR$V3, sep="-")

merge<-merge(Tpa_cov_genome, Tpa_cov_TR[4:5], by="chm_pos", all.x=T)
merge[is.na(merge)] <- 0 #NA are retrieved when windows are found without TR annotation

merge$scaff_number<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", merge$V1)))

#subset only chromosomes
merge2<-subset(merge, scaff_number=="1" | scaff_number=="2" | scaff_number=="3" | 
                 scaff_number=="4" |scaff_number=="5" | scaff_number=="6" | scaff_number=="7" | 
                 scaff_number=="8" | scaff_number=="9" | scaff_number=="10" | scaff_number=="11" | 
                 scaff_number=="12" | scaff_number=="13" | scaff_number=="14")

colnames(merge2)=c("chm_pos", "scaffold", "wind_start", "wind_end", "sum_cov_GW", "sum_cov_TR", "chr")


#TR proportion across chromosomes and plot
sum_cov_GW<-aggregate(merge2$sum_cov_GW~merge2$chr, FUN=sum)
sum_cov_TR<-aggregate(merge2$sum_cov_TR~merge2$chr, FUN=sum)
sum_cov<-cbind(sum_cov_GW, sum_cov_TR)
colnames(sum_cov)<-c("chr", "sum_cov_GW", "chr2", "sum_cov_TR")

sum_cov$propTR<-sum_cov$sum_cov_TR/sum_cov$sum_cov_GW
sum_cov$chr <- as.character(sum_cov$chr)
sum_cov$chr <- factor(sum_cov$chr,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12","13","14"))

ggplot(data=sum_cov, aes(x=chr, y=propTR)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(propTR, digits = 3)), vjust=1.6, size=3)+ ylab("proportion of tandem repeats")+
  theme_classic()


#Correlation with chromosome size and plot
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/chromosome_lengths")
Tpa_LG_length<-read.table("Tpa_chm_size_mtDNAv350.txt", header=F, sep="", quote = "")
colnames(Tpa_LG_length)=c("scaffold","length")
Tpa_LG_length$chr<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", Tpa_LG_length$scaffold)))

sum_cov_chr<-merge(sum_cov, Tpa_LG_length, by="chr")

ggplot(data=sum_cov_chr, aes(x=length, y=propTR)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=chr), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()

