#######################################
### TR Proportion along chromosomes ###
#######################################

###Tps
######
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRpropTps<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRpropTps)
mergeTRpropTps[is.na(mergeTRpropTps)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRpropTps$TRprop<-mergeTRpropTps$sum_cov_TR/mergeTRpropTps$sum_cov_GW


###Tcm
######
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tcm")

Tcm_cov_genome<-read.table("Tcm-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
colnames(Tcm_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tcm_cov_TR<-read.table("Tcm-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt
colnames(Tcm_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tcm_cov_genome$start2<-Tcm_cov_genome$start+1
Tcm_cov_TR$start2<-Tcm_cov_TR$start+1
Tcm_cov_genome$chm_pos<-paste(Tcm_cov_genome$chr, Tcm_cov_genome$start2, sep="-")
Tcm_cov_TR$chm_pos<-paste(Tcm_cov_TR$chr, Tcm_cov_TR$start2, sep="-")

mergeTRpropTcm<-cbind(Tcm_cov_genome, Tcm_cov_TR[4])
summary(mergeTRpropTcm)
mergeTRpropTcm[is.na(mergeTRpropTcm)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRpropTcm$TRprop<-mergeTRpropTcm$sum_cov_TR/mergeTRpropTcm$sum_cov_GW





#########################
### ChIP cenh3 testes ###
#########################

###Tps
######
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
cenh3Tps<-read.table("Tps_cenh3_testes_1_GW_coverage_DR_250kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_coverage_DR_100kb.txt

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
input1Tps<-read.table("Tps_input_testes_1_GW_coverage_DR_250kb.txt",
                   header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data1Tps<-cbind(cenh3Tps,input1Tps[4])
colnames(data1Tps)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1Tps$coverage_cenh3_norm<-data1Tps$coverage_cenh3/66295120 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1Tps$coverage_input_norm<-data1Tps$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1Tps$ratiocenH3_norm<-data1Tps$coverage_cenh3_norm/data1Tps$coverage_input_norm
data1Tps$log2cenH3_norm<-log2(data1Tps$ratiocenH3_norm)
data1Tps$start2<-data1Tps$start+1
data1Tps$chm_pos<-paste(data1Tps$Scaffold_name, data1Tps$start2, sep="-")
data1Tps$log2cenH3_norm[!is.finite(data1Tps$log2cenH3_norm)] <- 0



###Tcm
######
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
cenh3Tcm<-read.table("Tcm_testes_M-Fbodies_cenh3_1_coverage_DR_250kb.txt",
                     header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_coverage_DR_100kb.txt

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
input1Tcm<-read.table("Tcm_testes_M-Fbodies_input_coverage_DR_250kb.txt",
                      header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data1Tcm<-cbind(cenh3Tcm,input1Tcm[4])
colnames(data1Tcm)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1Tcm$coverage_cenh3_norm<-data1Tcm$coverage_cenh3/53128790 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1Tcm$coverage_input_norm<-data1Tcm$coverage_input/80611894 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1Tcm$ratiocenH3_norm<-data1Tcm$coverage_cenh3_norm/data1Tcm$coverage_input_norm
data1Tcm$log2cenH3_norm<-log2(data1Tcm$ratiocenH3_norm)
data1Tcm$start2<-data1Tcm$start+1
data1Tcm$chm_pos<-paste(data1Tcm$Scaffold_name, data1Tcm$start2, sep="-")
data1Tcm$log2cenH3_norm[!is.finite(data1Tcm$log2cenH3_norm)] <- 0





###################################################
### TR annotation with representative sequences ###
###################################################

###Tps
######

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
Tps_cleanedRegs_main <- readRDS("Tps_cleaned_withMainRepresentant.rds")

Tps_cleanedRegs_main$chr_start<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$start, sep = "_")
Tps_cleanedRegs_main$chr_end<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$end, sep = "_")

library(dplyr)
df_filtered_tps <- Tps_cleanedRegs_main %>%
  group_by(chr_start) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering


df_filtered2_tps <- df_filtered_tps %>%
  group_by(chr_end) %>%              # Group by the end position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering



###Tcm
######

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
Tcm_cleanedRegs_main <- readRDS("Tcm_cleaned_withMainRepresentant.rds")

Tcm_cleanedRegs_main$chr_start<-paste(Tcm_cleanedRegs_main$chr, Tcm_cleanedRegs_main$start, sep = "_")
Tcm_cleanedRegs_main$chr_end<-paste(Tcm_cleanedRegs_main$chr, Tcm_cleanedRegs_main$end, sep = "_")

library(dplyr)
df_filtered_tcm <- Tcm_cleanedRegs_main %>%
  group_by(chr_start) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering


df_filtered2_tcm <- df_filtered_tcm %>%
  group_by(chr_end) %>%              # Group by the end position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering





##################################################################
### Assign to each repeat array a non-overlapping 250kb window ###
##################################################################

###Tps
######

df_filtered2_tps$window_start <- (floor(df_filtered2_tps$start / 250000) * 250000)+1
df_filtered2_tps$window_end <- df_filtered2_tps$window_start + 249999  # End of the 100kb window
df_filtered2_tps$chm_pos<-paste(df_filtered2_tps$chr, df_filtered2_tps$window_start, sep="-")


###Tcm
######

df_filtered2_tcm$window_start <- (floor(df_filtered2_tcm$start / 250000) * 250000)+1
df_filtered2_tcm$window_end <- df_filtered2_tcm$window_start + 249999  # End of the 100kb window
df_filtered2_tcm$chm_pos<-paste(df_filtered2_tcm$chr, df_filtered2_tcm$window_start, sep="-")






#######################################################
### Merge TR prop, annotation and ChIP and datasets ###
#######################################################

###Tps
######

mergeTps<-merge(mergeTRpropTps, data1Tps[c(9,11)], by="chm_pos")
merge2Tps<-merge(mergeTps, df_filtered2_tps[c(10,11,19)], by="chm_pos")

chmsTps<-subset(mergeTps, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

chmsTps$scaff_number<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", chmsTps$chr)))

chmsTps<-chmsTps[order(chmsTps$scaff_number, chmsTps$start2),]
chmsTps$cumul_start<-seq(0, by = 250000, length.out = nrow(chmsTps))

###Tcm
######

mergeTcm<-merge(mergeTRpropTcm, data1Tcm[c(9,11)], by="chm_pos")
merge2Tcm<-merge(mergeTcm, df_filtered2_tcm[c(10,11,19)], by="chm_pos")

chmsTcm<-subset(mergeTcm, chr=="Tcm_LRv5a_scf1.1" | chr=="Tcm_LRv5a_scf1.2" | chr=="Tcm_LRv5a_scf2" | chr=="Tcm_LRv5a_scf3" | chr=="Tcm_LRv5a_scf4"
                | chr=="Tcm_LRv5a_scf5.1"| chr=="Tcm_LRv5a_scf5.2" | chr=="Tcm_LRv5a_scf6.1" | chr=="Tcm_LRv5a_scf6.2" 
                | chr=="Tcm_LRv5a_scf7" | chr=="Tcm_LRv5a_scf8.1" | chr=="Tcm_LRv5a_scf8.2" | chr=="Tcm_LRv5a_scf8.3" 
                | chr=="Tcm_LRv5a_scf9" | chr=="Tcm_LRv5a_scf10" | chr=="Tcm_LRv5a_scf11.1" | chr=="Tcm_LRv5a_scf11.2" | chr=="Tcm_LRv5a_scf12")

chmsTcm$scaff_number<-as.numeric(as.character(gsub(".*scf(\\d+).*", "\\1", chmsTcm$chr)))

chmsTcm<-chmsTcm[order(chmsTcm$scaff_number, chmsTcm$start2),]
chmsTcm$cumul_start<-seq(0, by = 250000, length.out = nrow(chmsTcm))


TR_prop<-ggplot(chmsTcm, aes(x=cumul_start, y=TRprop, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chmsTcm$scaff_number)))/2))+
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))






#############################################################
### Percent of shared TR sequences per chromosomal regions ###
##############################################################


###Tps
######

chmTps<-subset(chmsTps, scaff_number=="1")
chmTps_subset<-subset(chmTps, TRprop>0.15)

ggplot(chmTps_subset, aes(x=start2, y=TRprop))+
  geom_point(size=0.5) +
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )



###Tcm
######

chmTcm<-subset(chmsTcm, scaff_number=="1")
chmTcm_subset<-subset(chmTcm, TRprop>0.15)

ggplot(chmTcm_subset, aes(x=start2, y=TRprop))+
  geom_point(size=0.5) +
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )


#chr1: Tcm_LRv5a_scf1.1-1-Tcm_LRv5a_scf1.1-13000001 & Tcm_LRv5a_scf1.2-211250001-Tcm_LRv5a_scf1.2-226161097

telo1<-merge2Tcm[1:1079]






