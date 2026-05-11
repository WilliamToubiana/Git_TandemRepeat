#################################################################
## TR proportion from per-base coverage of indG Illumina reads ##
#################################################################

### Proportion between chromosomes
##################################

library(ggplot2)

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
sum_cov_chr$chr<-c("chr1","chr10","chr11","chr12","chr2","chrX","chr4","chr5","chr6","chr7","chr8","chr9")

ggplot(data=sum_cov_chr, aes(x=length, y=propTR)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=chr), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  scale_y_continuous(limits = c(0, NA))+
  geom_smooth(method='lm')+
  theme_classic()
  


##Stats

cor.test(sum_cov_chr$propTR, sum_cov_chr$length, method = "spearman")
cor.test(sum_cov_chr$propTR, sum_cov_chr$length, method = "pearson")

model<-lm(sum_cov_chr$propTR~sum_cov_chr$length)
summary(model)

sum_cov_auto<-subset(sum_cov_chr, chr!="chrX")
model<-lm(sum_cov_auto$propTR~sum_cov_auto$length)
summary(model)


#########################################################################################
### TR CN, sequence length, density and diversity comparisons between X and autosomes ###
#########################################################################################

## T.poppense TR annotation
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
Tps_cleanedRegs_main <- readRDS("Tps_cleaned_withMainRepresentant.rds")

Tps_cleanedRegs_main$chr_start<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$start, sep = "_")
Tps_cleanedRegs_main$chr_end<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$end, sep = "_")


library(dplyr)
df_filtered <- Tps_cleanedRegs_main %>%
  group_by(chr_start) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering


df_filtered2 <- df_filtered %>%
  group_by(chr_end) %>%              # Group by the end position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering


#Recalculate copy number based on array size and length of representative TR sequences 
df_filtered2$array_size<-abs(df_filtered2$end-df_filtered2$start)+1
df_filtered2$copy_nb2<-df_filtered2$array_size/df_filtered2$L_motif_representant


#Summary per chromosome
window_summary <- df_filtered2 %>%
  group_by(chr) %>%
  summarise(
    total_repeats = n(),
    avg_motif = mean(L_motif_representant),
    median_motif = median(L_motif_representant),
    median_CN = median(copy_nb2),
    mean_CN = mean(copy_nb2),
    motif_diversity = n_distinct(motif_representant),
    .groups = "drop",
    short_motif=sum(L_motif_representant < 10),
  )

summary(window_summary)



window_summary_chr<-subset(window_summary, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                           | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                           | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

window_summary_chr$LG_length <- c(253718452, 72695627, 67869710, 42522306, 178265049, 132376913, 102503173, 89066744, 87859127, 83130533, 76885976, 76156695)
window_summary_chr$LGs<-c("chr1", "ch10","chr11","chr12","chr2","chrX","chr4","chr5","chr6","chr7","chr8","chr9")
window_summary_chr$propTR<-sum_cov$propTR
window_summary_chr$chr <- factor(window_summary_chr$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3",
                                                             "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6",
                                                             "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9",
                                                             "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))
window_summary_chr$repeat_density<-window_summary_chr$total_repeats/window_summary_chr$LG_length

window_summary_autosomes<-subset(window_summary_chr, chr!="Tps_LRv5b_scf3")

window_summary_small<-subset(window_summary_chr, chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                           | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                           | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")


##Plots
ggplot(data=window_summary_chr, aes(x=propTR, y=median_CN)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("median copy number")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')+ylim(0,14)

model<-lm(window_summary_chr$median_CN~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$median_CN~window_summary_autosomes$propTR)
summary(model)


ggplot(data=window_summary_chr, aes(x=propTR, y=mean_CN)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("mean copy number")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')+ylim(0,32.5)

model<-lm(window_summary_chr$mean_CN~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$mean_CN~window_summary_autosomes$propTR)
summary(model)


ggplot(data=window_summary_chr, aes(x=propTR, y=median_motif)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("median TR sequence length")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')

model<-lm(window_summary_chr$median_motif~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$median_motif~window_summary_autosomes$propTR)
summary(model)


ggplot(data=window_summary_chr, aes(x=propTR, y=avg_motif)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1, size=5.5)+ ylab("mean TR sequence length")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')+ylim(0,30)

model<-lm(window_summary_chr$avg_motif~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$avg_motif~window_summary_autosomes$propTR)
summary(model)


ggplot(data=window_summary_chr, aes(x=propTR, y=total_repeats)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("total TR sequences")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')+ylim(0,45000)

model<-lm(window_summary_chr$total_repeats~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$total_repeats~window_summary_autosomes$propTR)
summary(model)


ggplot(data=window_summary_chr, aes(x=propTR, y=repeat_density)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("TR array density")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')+ylim(0,0.00025)

model<-lm(window_summary_chr$repeat_density~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$repeat_density~window_summary_autosomes$propTR)
summary(model)


window_summary_chr$motif_diversity_corrected<-window_summary_chr$motif_diversity/window_summary_chr$total_repeats
ggplot(data=window_summary_chr, aes(x=propTR, y=motif_diversity_corrected)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("TR sequence diversity / total TR sequences")+ xlab("TR proportion")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')+ylim(0,0.5)

model<-lm(window_summary_chr$motif_diversity_corrected~window_summary_chr$propTR)
summary(model)
model<-lm(window_summary_autosomes$motif_diversity_corrected~window_summary_autosomes$propTR)
summary(model)

