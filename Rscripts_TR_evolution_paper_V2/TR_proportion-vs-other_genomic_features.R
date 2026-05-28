library("ggVennDiagram")
library(ggplot2)
library(GenomicRanges)
library(dplyr)

##################################
# Correlation TR - recombination #
##################################

##Tandem repeats
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/100000)] += $3; next} {print $1, $2, $3, a[$1, int($2/100000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w100000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_100kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/100000)] += $3; next} {print $1, $2, $3, a[$1, int($2/100000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w100000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_100kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRprop<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW


merge2TRprop<-subset(mergeTRprop, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | 
                       chr=="Tps_LRv5b_scf4" |chr=="Tps_LRv5b_scf5" | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | 
                       chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10" | chr=="Tps_LRv5b_scf11" | 
                       chr=="Tps_LRv5b_scf12")

merge2TRprop$scaff_number<-as.numeric(as.character(gsub("^.*scf","", merge2TRprop$chr)))

merge3TRprop<-merge2TRprop[order(merge2TRprop$scaff_number, merge2TRprop$start2),]
merge3TRprop$cumul_start<-seq(0, by = 250000, length.out = nrow(merge3TRprop))
merge3TRprop$chrom<-ifelse(merge3TRprop$scaff_number=="3", "sex-chrom", "autosome")
merge4TRprop<-merge3TRprop[merge3TRprop$scaff_number=="3",]



TR_array<-ggplot(merge3TRprop, aes(x=cumul_start, y=TRprop, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge3TRprop$scaff_number)))/2))+
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))




##recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
  group_by(scaffold,start)%>%
  summarize(median_rho=median(rho))

rho_g$start2<-rho_g$start+1
rho_g$chm_pos<-paste(rho_g$scaffold, rho_g$start2, sep="-")


##Merge TR and recomb datasets
merge<-merge(merge3TRprop, rho_g, by="chm_pos")

rho<-ggplot(merge, aes(x=cumul_start, y=median_rho, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge$scaff_number)))/2))+
  ylab("rho")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+ylim(0,0.0015)+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


#filtering extreme
merge_sub <- merge %>%
  group_by(scaffold)%>%
  mutate(quantile=quantile(median_rho,0.99))%>%
  filter(median_rho<quantile)


###Extract and compare values of the Generalized Additive Models per chromosome
library(mgcv)

merge_sub2<-subset(merge_sub, chr=="Tps_LRv5b_scf1")

gamrho<-gam(median_rho ~ s(cumul_start, bs = "cs", k=30), data=merge_sub2)#per chromosomes
gamTR<-gam(TRprop ~ s(cumul_start, bs = "cs", k=30), data=merge_sub2)#per chromosomes

predicted_values_gamrho <- data.frame(cumul_start = merge_sub2$cumul_start, 
                                      fitted = predict(gamrho))
predicted_values_gamTR <- data.frame(cumul_start = merge_sub2$cumul_start, 
                                     fitted = predict(gamTR))


summary(lm(predicted_values_gamTR$fitted~predicted_values_gamrho$fitted))
cor.test(predicted_values_gamrho$fitted, predicted_values_gamTR$fitted, method = "spearman")


ggplot(merge_sub2, aes(x = cumul_start, y = median_rho)) +
  geom_point() +
  geom_line(data = predicted_values_gamrho, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))

ggplot(merge_sub2, aes(x = cumul_start, y = TRprop)) +
  geom_point() +
  geom_line(data = predicted_values_gamTR, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))




###############################
# Correlation TR - GC content #
###############################
##Tandem repeats
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/100000)] += $3; next} {print $1, $2, $3, a[$1, int($2/100000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w100000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_100kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/100000)] += $3; next} {print $1, $2, $3, a[$1, int($2/100000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w100000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_100kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRprop<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW


merge2TRprop<-subset(mergeTRprop, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | 
                       chr=="Tps_LRv5b_scf4" |chr=="Tps_LRv5b_scf5" | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | 
                       chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10" | chr=="Tps_LRv5b_scf11" | 
                       chr=="Tps_LRv5b_scf12")

merge2TRprop$scaff_number<-as.numeric(as.character(gsub("^.*scf","", merge2TRprop$chr)))

merge3TRprop<-merge2TRprop[order(merge2TRprop$scaff_number, merge2TRprop$start2),]
merge3TRprop$cumul_start<-seq(0, by = 250000, length.out = nrow(merge3TRprop))
merge3TRprop$chrom<-ifelse(merge3TRprop$scaff_number=="3", "sex-chrom", "autosome")
merge4TRprop<-merge3TRprop[merge3TRprop$scaff_number=="3",]



TR_array<-ggplot(merge3TRprop, aes(x=cumul_start, y=TRprop, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge3TRprop$scaff_number)))/2))+
  ylab("Proportion of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



##GC content
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("GC_content_Tps_250kb.txt", header=FALSE) #bedtools nuc -fi Tps_LRv5b_mtDNAv350.fasta -bed Tps_chm_size_mtDNAv350_w100000.bed  | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"($7+$8)/($6+$7+$8+$9+1)}' > GC_content_Tps_100kb.txt
colnames(data_GC)<-c("chromosome", "start", "end", "perc_GC")

data_GC$start2<-data_GC$start+1
data_GC$chm_pos<-paste(data_GC$chromosome, data_GC$start2, sep="-")


#chromosomes only
chms<-subset(data_GC, chromosome=="Tps_LRv5b_scf1" | chromosome=="Tps_LRv5b_scf2" | chromosome=="Tps_LRv5b_scf3" | chromosome=="Tps_LRv5b_scf4" | chromosome=="Tps_LRv5b_scf5" 
             | chromosome=="Tps_LRv5b_scf6" | chromosome=="Tps_LRv5b_scf7" | chromosome=="Tps_LRv5b_scf8" | chromosome=="Tps_LRv5b_scf9" | chromosome=="Tps_LRv5b_scf10"
             | chromosome=="Tps_LRv5b_scf11" | chromosome=="Tps_LRv5b_scf12")

summary(chms)

##Merge TR and recomb datasets
merge2<-merge(merge3TRprop, chms, by="chm_pos")

GC_content<-ggplot(merge2, aes(x=cumul_start, y=perc_GC, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge2$scaff_number)))/2))+
  ylab("Percent of GC")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


###Extract and compare values of the Generalized Additive Models

merge_sub<-merge2
merge_sub<-subset(merge2, chr=="Tps_LRv5b_scf1")

library(mgcv)
gamGC<-gam(perc_GC ~ s(cumul_start, bs = "cs", k=30), data=merge_sub)#per chromosomes
gamTR<-gam(TRprop ~ s(cumul_start, bs = "cs", k=30), data=merge_sub)#per chromosomes

predicted_values_gamGC <- data.frame(cumul_start = merge_sub$cumul_start, 
                                     fitted = predict(gamGC))
predicted_values_gamTR <- data.frame(cumul_start = merge_sub$cumul_start, 
                                     fitted = predict(gamTR))

summary(lm(predicted_values_gamTR$fitted~predicted_values_gamGC$fitted))
cor.test(predicted_values_gamGC$fitted, predicted_values_gamTR$fitted, method = "spearman")




ggplot(merge_sub, aes(x = cumul_start, y = perc_GC)) +
  geom_point() +
  geom_line(data = predicted_values_gamGC, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))

ggplot(merge_sub, aes(x = cumul_start, y = TRprop)) +
  geom_point() +
  geom_line(data = predicted_values_gamTR, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))





# Plots TRprop vs recomb vs GC content
library(gridExtra)
grid.arrange(TR_array, rho, GC_content, nrow=3, ncol=1)


library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    TR_array + theme(legend.position = "none")
    , rho
    , GC_content + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(rho)
    , ggplot()
    , get_legend(GC_content)
    , ncol =1)
  , rel_widths = c(12,0)
)




###############################
# Correlation TR - centromere #
###############################




###############################################################
###Correlation recombination and TR load between chromosomes###
###############################################################

###Tandem repeat proportion between chromosomes Timema
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tps_subset<-subset(Tps, copy_nb >= 5)
Tps_subset$array_length<-(Tps_subset$end-Tps_subset$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tps_subset$chr,ranges=IRanges(start=Tps_subset$start,end=Tps_subset$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tps_LRv5b_scf1" | seqnames=="Tps_LRv5b_scf2" | seqnames=="Tps_LRv5b_scf3" | 
                             seqnames=="Tps_LRv5b_scf4" |seqnames=="Tps_LRv5b_scf5" | seqnames=="Tps_LRv5b_scf6" | seqnames=="Tps_LRv5b_scf7" | 
                             seqnames=="Tps_LRv5b_scf8" | seqnames=="Tps_LRv5b_scf9" | seqnames=="Tps_LRv5b_scf10" | seqnames=="Tps_LRv5b_scf11" | 
                             seqnames=="Tps_LRv5b_scf12")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)

#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
Tps_LG_length<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")

nonoverlapWind_chr_sum$LG_length <- Tps_LG_length$length
nonoverlapWind_chr_sum$Prop_repeated_region<-nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`/nonoverlapWind_chr_sum$LG_length

colnames(nonoverlapWind_chr_sum) <- c("scaffold", "TR_length", "LG_length", "Prop_repeated_region")


###recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
  group_by(scaffold)%>%
  summarize(median_rho=median(rho))


merge<-merge(nonoverlapWind_chr_sum, rho_g, by="scaffold")
merge$scaffold<-factor(merge$scaffold, levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))
merge_sub<-subset(merge, scaffold!="Tps_LRv5b_scf3")

cor.test(merge$median_rho, merge$Prop_repeated_region, method = "spearman")#chrX included
summary(lm(merge$Prop_repeated_region~merge$median_rho))#chrX included

cor.test(merge_sub$median_rho, merge_sub$Prop_repeated_region, method = "spearman")#chrX excluded
summary(lm(merge_sub$Prop_repeated_region~merge_sub$median_rho))#chrX excluded

cor.test(merge$median_rho, merge$Prop_repeated_region, method = "spearman")#chrX included
summary(lm(merge$Prop_repeated_region~merge$LG_length))#chrX included





##########################################
### TR distributions along chromosomes ###
##########################################
##Tandem repeats
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/100000)] += $3; next} {print $1, $2, $3, a[$1, int($2/100000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w100000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_100kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/100000)] += $3; next} {print $1, $2, $3, a[$1, int($2/100000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w100000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_100kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRprop<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW


merge2TRprop<-subset(mergeTRprop, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | 
                       chr=="Tps_LRv5b_scf4" |chr=="Tps_LRv5b_scf5" | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | 
                       chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10" | chr=="Tps_LRv5b_scf11" | 
                       chr=="Tps_LRv5b_scf12")

merge2TRprop$scaff_number<-as.numeric(as.character(gsub("^.*scf","", merge2TRprop$chr)))

merge3TRprop<-merge2TRprop[order(merge2TRprop$scaff_number, merge2TRprop$start2),]
merge3TRprop$cumul_start<-seq(0, by = 250000, length.out = nrow(merge3TRprop))
merge3TRprop$chrom<-ifelse(merge3TRprop$scaff_number=="3", "sex-chrom", "autosome")
merge4TRprop<-merge3TRprop[merge3TRprop$scaff_number=="3",]


ggplot(chms1, aes(x=start, y=percTR))+
  geom_point(size=0.5) +
  ylab("Percent of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")+theme(axis.text.y = element_text(angle = 90))

TR_array<-ggplot(merge3TRprop, aes(x=cumul_start, y=TRprop, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge3TRprop$scaff_number)))/2))+
  geom_point(data=merge4TRprop, aes(x=merge4TRprop$cumul_start, y=merge4TRprop$TRprop), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")



##recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
       group_by(scaffold,start)%>%
       summarize(median_rho=median(rho))

#filtering extreme
rho_gf <- rho_g %>%
       group_by(scaffold)%>%
       mutate(quantile=quantile(median_rho,0.99))%>%
       filter(median_rho<quantile)

#keeping only big scaf
scaffolds <- paste("Tps_LRv5b_scf",1:12,sep="")

rho_gf$scaff_number<-as.numeric(as.character(gsub("^.*scf","", rho_gf$scaffold)))
rho_gf<-rho_gf %>% 
       filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms3<-rho_gf[order(rho_gf$scaff_number, rho_gf$start),]
chms3$cumul_start<-seq(0, by = 10000, length.out = nrow(chms3))
chms3$chrom<-ifelse(chms3$scaff_number=="3", "sex-chrom", "autosome")
chms4<-chms3[chms3$scaff_number=="3",]


recomb<-ggplot(chms3, aes(x=cumul_start, y=median_rho, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms3$scaff_number)))/2))+
  geom_point(data=chms4, aes(x=chms4$cumul_start, y=chms4$median_rho), stat="identity", colour="maroon", size=0.5) +
  ylab("rho")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+ylim(0,0.0015)



##GC content
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("GC_content_Tps_250kb.txt", header=FALSE) #bedtools nuc -fi Tps_LRv5b_mtDNAv350.fasta -bed Tps_chm_size_mtDNAv350_w100000.bed  | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"($7+$8)/($6+$7+$8+$9+1)}' > GC_content_Tps_100kb.txt
colnames(data_GC)<-c("chromosome", "start", "end", "perc_GC")

data_GC$start2<-data_GC$start+1
data_GC$chm_pos<-paste(data_GC$chromosome, data_GC$start2, sep="-")



#Plot

chms<-subset(data_GC, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms3<-chms[order(chms$scaff_number, chms$V2),]
chms3$cumul_start<-seq(0, by = 10000, length.out = nrow(chms3))
chms3$chrom<-ifelse(chms3$scaff_number=="3", "sex-chrom", "autosome")
chms4<-chms3[chms3$scaff_number=="3",]

GC_content<-ggplot(chms3, aes(x=cumul_start, y=V5, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms3$scaff_number)))/2))+
  geom_point(data=chms4, aes(x=chms4$cumul_start, y=chms4$V5), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of GC")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")






##ChIP cenh3 testes

setwd("/Users/wtoubian/Desktop/Chip-seq/species-comparison/tracks/Tps_testes_cenh3_R1")
cenh3<-read.table("Tps_testes_cenh3_1_coverage_DR_10kb.txt",
                  header=FALSE)

setwd("/Users/wtoubian/Desktop/Chip-seq/species-comparison/tracks/Tps_testes_input_R1")
input1<-read.table("Tps_testes_input1_coverage_DR_10kb.txt",
                   header=FALSE)

data1<-cbind(cenh3,input1[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1$coverage_cenh3_norm<-data1$coverage_cenh3/66295120 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratiocenH3_norm<-data1$coverage_cenh3_norm/data1$coverage_input_norm
data1$log2cenH3_norm<-log2(data1$ratiocenH3_norm)

data2<-subset(data1, coverage_cenh3>1 | coverage_input>1)
data2$start_chr<-paste(data2$start, data2$Scaffold_name, sep="-")

mergeCenh3_wind<-merge(windows, data2, by="start_chr", all.x=T)
mergeCenh3_wind[is.na(mergeCenh3_wind)] <- 0

##Plots

chms<-subset(mergeCenh3_wind, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms9<-chms[order(chms$scaff_number, chms$start),]
chms9$cumul_start<-seq(0, by = 10000, length.out = nrow(chms9))
chms9$chrom<-ifelse(chms9$scaff_number=="3", "sex-chrom", "autosome")
chms10<-chms9[chms9$scaff_number=="3",]


cenh3_testes<-ggplot(chms9, aes(x=cumul_start, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms9$scaff_number)))/2))+
  scale_y_continuous(limits=c(-1,7)) +
  geom_point(data=chms10, aes(x=chms10$cumul_start, y=chms10$log2cenH3_norm), stat="identity", colour="maroon", size=0.5) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept=c(2,0,0), linetype=c("dashed","solid", "dashed"), color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white", method="gam", formula = y~s(x, bs="cs", k=20))


###Extract coordinates smooth line intersecting with X-axis
library(mgcv)
chr<-subset(chms9, V1=="Tps_LRv5b_scf12") #example for one chromosome
gam<-gam(log2cenH3_norm ~ s(V2, bs = "cs", k=30), data=chr)
chr$predict<-predict(gam, newdata = chr)
chr$predict_abs<-abs(chr$predict)

ggplot(chr, aes(x=V2, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_y_continuous(limits=c(-1,7)) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept=c(2,0,0), linetype=c("dashed","solid", "dashed"), color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))+geom_vline(xintercept =c(46640000, 60030000))


#All Plots
library(gridExtra)
grid.arrange(cenh3_testes, TR_array, h3k9_testes, nrow=3, ncol=1)


library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    cenh3_testes + theme(legend.position = "none")
    , TR_array
    , h3k9_testes + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(h3k9_testes)
    , ncol =1)
  , rel_widths = c(12,3)
)



plot_grid(
  plot_grid(
    GC_content + theme(legend.position = "none")
    , TR_array
    , GC_content_noTR + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(GC_content)
    , ggplot()
    , get_legend(GC_content_noTR)
    , ncol =1)
  , rel_widths = c(12,3)
)


plot_grid(
  plot_grid(
    cenh3_testes + theme(legend.position = "none")
    , TR_array
    , recomb + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(recomb)
    , ncol =1)
  , rel_widths = c(12,0)
)


