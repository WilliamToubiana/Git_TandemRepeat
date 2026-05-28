#################################################################
## TR proportion from per-base coverage of indG Illumina reads ##
#################################################################

### Proportion along chromosomes
################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("sum_coverage_genome_100kb.txt", header=F, sep="\t", quote="") #This was generated from the perl script "script_coverage_100kbW_2.pl" loacted in "/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/scripts"
Tps_cov_TR<-read.table("sum_coverage_repeat_100kb.txt", header=F, sep="\t", quote="") #This was generated from the perl script "script_coverage_100kbW_2.pl" loacted in "/work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/scripts"

Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$V1, Tps_cov_genome$V3, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$V1, Tps_cov_TR$V3, sep="-")

merge<-merge(Tps_cov_genome, Tps_cov_TR[4:5], by="chm_pos", all.x=T)
merge[is.na(merge)] <- 0 #NA are retrieved when windows are found without TR annotation
merge$TRprop<-merge$V4.y/merge$V4.x

merge2<-subset(merge, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | 
                 V1=="Tps_LRv5b_scf4" |V1=="Tps_LRv5b_scf5" | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | 
                 V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10" | V1=="Tps_LRv5b_scf11" | 
                 V1=="Tps_LRv5b_scf12")

colnames(merge2)=c("chm_pos", "chr", "wind_start", "wind_end", "sum_cov_GW", "sum_cov_TR", "TRprop")

merge2$scaff_number<-as.numeric(as.character(gsub("^.*scf","", merge2$chr)))

merge3<-merge2[order(merge2$scaff_number, merge2$wind_start),]
merge3$cumul_start<-seq(0, by = 100000, length.out = nrow(merge3))
merge3$chrom<-ifelse(merge3$scaff_number=="3", "sex-chrom", "autosome")
merge4<-merge3[merge3$scaff_number=="3",]


TR<-ggplot(merge3, aes(x=cumul_start, y=TRprop, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge3$scaff_number)))/2))+
  scale_y_continuous(limits=c(0,1)) +
  geom_point(data=merge4, aes(x=merge4$cumul_start, y=merge4$TRprop), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))

  
  


### Proportion between chromosomes
##################################


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


### Proportion along chromosomes
################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
TR_overlap_windows<-read.table("TR_nonoverlapWind_chr_intersection100kbWind.txt", header=F, sep="\t", quote="")
#The table "TR_nonoverlapWind_chr.txt" was generated from "Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt" to remove overlaps between TR arrays
##The table "TR_nonoverlapWind_chr.txt" was generated to intersect non-overlapping TR positions with 100kb windows across every chromosomes.
#This was computed as: bedtools intersect -a TR_annotation_timema/minimal_rotations/TR_nonoverlapWind_chr.txt -b genomes/Tps_chm_size_mtDNAv350_w100000.bed > TR_annotation_timema/minimal_rotations/TR_nonoverlapWind_chr_intersection100kbWind.txt


colnames(TR_overlap_windows) <- c("chr", "start", "end")
TR_overlap_windows$array<-(TR_overlap_windows$end-TR_overlap_windows$start)

TR_overlap_windows$window_start <- floor(TR_overlap_windows$start / 100000) * 100000
TR_overlap_windows$window_end <- TR_overlap_windows$window_start + 99999  # End of the 100kb window
TR_overlap_windows$chm_pos<-paste(TR_overlap_windows$chr, TR_overlap_windows$window_end, sep="-")


setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
chromosome_ranges<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")
colnames(chromosome_ranges)<-c("chr", "max_position")

windows <- chromosome_ranges %>%
  rowwise() %>%
  do({
    chrom <- .$chr
    max_pos <- .$max_position
    start_pos <- seq(0, max_pos, by = 100000)
    end_pos <- pmin(start_pos + 99999, max_pos)
    data.frame(chr = chrom, window_start = start_pos, window_end = end_pos)
  })

windows$chm_pos<-paste(windows$chr, windows$window_end, sep="-")


merge3<-merge(windows, TR_overlap_windows[c(4,7)], by="chm_pos", all.x=T)
merge3[is.na(merge3)] <- 0 #NA are retrieved when windows are found without TR annotation


merge4<-aggregate(merge3$array~merge3$chm_pos, FUN=sum)
colnames(merge4)<-c("chm_pos", "array")
merge4$TRprop2<-merge4$array/100000

library(tidyr)
merge4<-merge3 %>%
  separate_wider_delim(chm_pos, "-", names=c("chr", "start", "end"))


merge_cov_TRF<-merge(merge2, merge4[c(1,3)], by="chm_pos", all.x=T)


ggplot(merge_cov_TRF, aes(x=TRprop, y=TRprop2)) + geom_point()




setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
Tps_cleanedRegs_main <- readRDS("Tps_cleaned_withMainRepresentant.rds")

Tps_cleanedRegs_main$chr_start<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$start, sep = "_")
Tps_cleanedRegs_main$chr_end<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$end, sep = "_")

df_filtered <- Tps_cleanedRegs_main %>%
  group_by(chr_start) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the longest motif
  ungroup()                                 # Ungroup after filtering


df_filtered2 <- df_filtered %>%
  group_by(chr_end) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the longest motif
  ungroup()                                 # Ungroup after filtering

df_filtered2$start_window <- floor(df_filtered2$start / 100000) * 100000
df_filtered2$end_window <- df_filtered2$start_window + 99999  # End of the 10kb window

df_filtered2$chr_end_window<-paste(df_filtered2$chr, df_filtered2$end_window, sep = "_")

mean_motif_length<-aggregate(df_filtered2$L_motif_representant~df_filtered2$chr_end_window, FUN=mean)
colnames(mean_motif_length)<-c("chr_end_window", "mean_motif_length")

mean_motif_length2<-mean_motif_length %>%
  separate_wider_delim(chr_end_window, "_", names=c("species", "ref", "scf", "end_window"))


mean_motif_length2$scaff_number<-as.numeric(as.character(gsub("^.*scf","", mean_motif_length2$scf)))
mean_motif_length2$end_window<-as.numeric(mean_motif_length2$end_window)
mean_motif_length2$start_window<-mean_motif_length2$end_window-99999

mean_motif_length3<-mean_motif_length2[order(mean_motif_length2$scaff_number, mean_motif_length2$start_window),]

mean_motif_length4<-subset(mean_motif_length3, scaff_number=="1" | scaff_number=="2" | scaff_number=="3" | 
                             scaff_number=="4" |scaff_number=="5" | scaff_number=="6" | scaff_number=="7" | 
                             scaff_number=="8" | scaff_number=="9" | scaff_number=="10" | scaff_number=="11" | 
                             scaff_number=="12")

mean_motif_length4$cumul_start<-seq(0, by = 100000, length.out = nrow(mean_motif_length4))

ggplot(mean_motif_length4, aes(x=cumul_start, y=mean_motif_length, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(mean_motif_length4$scaff_number)))/2))+
  ylab("mean motif length")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+ylim(0,200)






#####################################################
## Distribution genomic features along chromosomes ##
#####################################################


### ChIP H3K9me3 somatic males
##############################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/h3k9")

h3k9<-read.table("Tps_male_organs_h3k9A_coverage_DR_100kb.txt",
                 header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_male_organs_h3k9_R1/Tps_male_organs_h3k9A_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_male_organs_h3k9_R1/Tps_male_organs_h3k9A_coverage_DR_100kb.txt

input<-read.table("Tps_male_organs_input3_coverage_DR_100kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_male_organs_input_R3/Tps_male_organs_input3_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_male_organs_input_R3/Tps_male_organs_input3_coverage_DR_100kb.txt

data1<-cbind(h3k9,input[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_h3k9", "coverage_input")

data1$coverage_h3k9_norm<-data1$coverage_h3k9/62681606 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/59528909 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratioH3k9_norm<-data1$coverage_h3k9_norm/data1$coverage_input_norm
data1$log2H3k9_norm<-log2(data1$ratioH3k9_norm)

#data2<-subset(data1, coverage_h3k9>1 | coverage_input>1)
#data2$start_chr<-paste(data2$start, data2$Scaffold_name, sep="-")

#mergeH3k9_wind<-merge(windows, data2, by="start_chr", all.x=T)
#mergeH3k9_wind[is.na(mergeH3k9_wind)] <- 0


##Plots

chms<-subset(data1, Scaffold_name=="Tps_LRv5b_scf1" | Scaffold_name=="Tps_LRv5b_scf2" | Scaffold_name=="Tps_LRv5b_scf3" | Scaffold_name=="Tps_LRv5b_scf4" | Scaffold_name=="Tps_LRv5b_scf5" 
             | Scaffold_name=="Tps_LRv5b_scf6" | Scaffold_name=="Tps_LRv5b_scf7" | Scaffold_name=="Tps_LRv5b_scf8" | Scaffold_name=="Tps_LRv5b_scf9" | Scaffold_name=="Tps_LRv5b_scf10"
             | Scaffold_name=="Tps_LRv5b_scf11" | Scaffold_name=="Tps_LRv5b_scf12")

#chms<-subset(mergeH3k9_wind, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
#             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
#             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")

chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$Scaffold_name)))
#chms<-chms %>% 
#  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms7<-chms[order(chms$scaff_number, chms$start),]
chms7$cumul_start<-seq(0, by = 100000, length.out = nrow(chms7))
chms7$chrom<-ifelse(chms7$scaff_number=="3", "sex-chrom", "autosome")
chms8<-chms7[chms7$scaff_number=="3",]


h3k9_soma<-ggplot(chms7, aes(x=cumul_start, y=log2H3k9_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms7$scaff_number)))/2))+
  scale_y_continuous(limits=c(-2,2)) +
  geom_point(data=chms8, aes(x=chms8$cumul_start, y=chms8$log2H3k9_norm), stat="identity", colour="maroon", size=0.5) +
  ylab("log2(H3K9me3/Input)")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))




### ChIP H3K9me3 testes
#######################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/h3k9")

h3k9<-read.table("Tps_testes_h3k9_coverage_DR_100kb.txt",
                 header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_h3k9/Tps_testes_h3k9_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_h3k9/Tps_testes_h3k9_coverage_DR_100kb.txt

input<-read.table("Tps_testes_input1_coverage_DR_100kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data1<-cbind(h3k9,input[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_h3k9", "coverage_input")

data1$coverage_h3k9_norm<-data1$coverage_h3k9/57113559 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratioH3k9_norm<-data1$coverage_h3k9_norm/data1$coverage_input_norm
data1$log2H3k9_norm<-log2(data1$ratioH3k9_norm)

#data2<-subset(data1, coverage_h3k9>1 | coverage_input>1)
#data2$start_chr<-paste(data2$start, data2$Scaffold_name, sep="-")

#mergeH3k9_wind<-merge(windows, data2, by="start_chr", all.x=T)
#mergeH3k9_wind[is.na(mergeH3k9_wind)] <- 0


##Plots

chms<-subset(data1, Scaffold_name=="Tps_LRv5b_scf1" | Scaffold_name=="Tps_LRv5b_scf2" | Scaffold_name=="Tps_LRv5b_scf3" | Scaffold_name=="Tps_LRv5b_scf4" | Scaffold_name=="Tps_LRv5b_scf5" 
             | Scaffold_name=="Tps_LRv5b_scf6" | Scaffold_name=="Tps_LRv5b_scf7" | Scaffold_name=="Tps_LRv5b_scf8" | Scaffold_name=="Tps_LRv5b_scf9" | Scaffold_name=="Tps_LRv5b_scf10"
             | Scaffold_name=="Tps_LRv5b_scf11" | Scaffold_name=="Tps_LRv5b_scf12")

#chms<-subset(mergeH3k9_wind, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
#             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
#             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")

chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$Scaffold_name)))
#chms<-chms %>% 
#  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms7<-chms[order(chms$scaff_number, chms$start),]
chms7$cumul_start<-seq(0, by = 100000, length.out = nrow(chms7))
chms7$chrom<-ifelse(chms7$scaff_number=="3", "sex-chrom", "autosome")
chms8<-chms7[chms7$scaff_number=="3",]


h3k9_testes<-ggplot(chms7, aes(x=cumul_start, y=log2H3k9_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms7$scaff_number)))/2))+
  scale_y_continuous(limits=c(-1,1)) +
  geom_point(data=chms8, aes(x=chms8$cumul_start, y=chms8$log2H3k9_norm), stat="identity", colour="maroon", size=0.5) +
  ylab("log2(H3K9me3/Input)")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))




### ChIP cenh3 testes
#####################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
cenh3<-read.table("Tps_testes_cenh3_1_coverage_DR_100kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_coverage_DR_100kb.txt

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
input1<-read.table("Tps_testes_input1_coverage_DR_100kb.txt",
                   header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w100000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data1<-cbind(cenh3,input1[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1$coverage_cenh3_norm<-data1$coverage_cenh3/66295120 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratiocenH3_norm<-data1$coverage_cenh3_norm/data1$coverage_input_norm
data1$log2cenH3_norm<-log2(data1$ratiocenH3_norm)

#data2<-subset(data1, coverage_h3k9>1 | coverage_input>1)
#data2$start_chr<-paste(data2$start, data2$Scaffold_name, sep="-")

#mergeH3k9_wind<-merge(windows, data2, by="start_chr", all.x=T)
#mergeH3k9_wind[is.na(mergeH3k9_wind)] <- 0


##Plots

chms<-subset(data1, Scaffold_name=="Tps_LRv5b_scf1" | Scaffold_name=="Tps_LRv5b_scf2" | Scaffold_name=="Tps_LRv5b_scf3" | Scaffold_name=="Tps_LRv5b_scf4" | Scaffold_name=="Tps_LRv5b_scf5" 
             | Scaffold_name=="Tps_LRv5b_scf6" | Scaffold_name=="Tps_LRv5b_scf7" | Scaffold_name=="Tps_LRv5b_scf8" | Scaffold_name=="Tps_LRv5b_scf9" | Scaffold_name=="Tps_LRv5b_scf10"
             | Scaffold_name=="Tps_LRv5b_scf11" | Scaffold_name=="Tps_LRv5b_scf12")

#chms<-subset(mergeH3k9_wind, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
#             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
#             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")

chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$Scaffold_name)))
#chms<-chms %>% 
#  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms9<-chms[order(chms$scaff_number, chms$start),]
chms9$cumul_start<-seq(0, by = 100000, length.out = nrow(chms9))
chms9$chrom<-ifelse(chms9$scaff_number=="3", "sex-chrom", "autosome")
chms10<-chms9[chms9$scaff_number=="3",]


cenh3_testes<-ggplot(chms9, aes(x=cumul_start, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms9$scaff_number)))/2))+
  scale_y_continuous(limits=c(-1,3)) +
  geom_point(data=chms10, aes(x=chms10$cumul_start, y=chms10$log2cenH3_norm), stat="identity", colour="maroon", size=0.5) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))




### GC content
##############

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("GC_content_Tps_100kb.txt", header=FALSE) #bedtools nuc -fi Tps_LRv5b_mtDNAv350.fasta -bed Tps_chm_size_mtDNAv350_w100000.bed  | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"($7+$8)/($6+$7+$8+$9+1)}' > GC_content_Tps_100kb.txt

#Plot

chms<-subset(data_GC, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
#chms<-chms %>% 
#  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms3<-chms[order(chms$scaff_number, chms$V2),]
chms3$cumul_start<-seq(0, by = 100000, length.out = nrow(chms3))
chms3$chrom<-ifelse(chms3$scaff_number=="3", "sex-chrom", "autosome")
chms4<-chms3[chms3$scaff_number=="3",]

GC_content<-ggplot(chms3, aes(x=cumul_start, y=V4, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms3$scaff_number)))/2))+
  scale_y_continuous(limits=c(0.3,0.5)) +
  geom_point(data=chms4, aes(x=chms4$cumul_start, y=chms4$V4), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of GC")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))





### GC content without TR
#########################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("GC_content_Tps_100kb_withoutTR.txt", header=FALSE) #bedtools nuc -fi Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000.mask -bed Tps_chm_size_mtDNAv350_mask_w100000.bed  | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"($7+$8)/($6+$7+$8+$9+1)}' > GC_content_Tps_100kb_withoutTR.txt


#Plot

chms<-subset(data_GC, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
#chms<-chms %>% 
#  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms5<-chms[order(chms$scaff_number, chms$V2),]
chms5$cumul_start<-seq(0, by = 100000, length.out = nrow(chms5))
chms5$chrom<-ifelse(chms5$scaff_number=="3", "sex-chrom", "autosome")
chms6<-chms5[chms5$scaff_number=="3",]

GC_content_noTR<-ggplot(chms5, aes(x=cumul_start, y=V4, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms5$scaff_number)))/2))+
  scale_y_continuous(limits=c(0.3,0.5)) +
  geom_point(data=chms6, aes(x=chms6$cumul_start, y=chms6$V4), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of GC")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))





### recombination
#################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection100kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
  group_by(scaffold,start)%>%
  summarize(median_rho=median(rho))

#filtering extreme
#rho_gf <- rho_g %>%
#  group_by(scaffold)%>%
#  mutate(quantile=quantile(median_rho,0.99))%>%
#  filter(median_rho<quantile)

#keeping only big scaf
#scaffolds <- paste("Tps_LRv5b_scf",1:12,sep="")

#rho_gf$scaff_number<-as.numeric(as.character(gsub("^.*scf","", rho_gf$scaffold)))
rho_g$scaff_number<-as.numeric(as.character(gsub("^.*scf","", rho_g$scaffold)))

#rho_gf<-rho_gf %>% 
#  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


#chms1<-rho_gf[order(rho_gf$scaff_number, rho_gf$start),]
chms1<-rho_g[order(rho_g$scaff_number, rho_g$start),]
chms1$cumul_start<-seq(0, by = 100000, length.out = nrow(chms1))
chms1$chrom<-ifelse(chms1$scaff_number=="3", "sex-chrom", "autosome")
chms2<-chms1[chms1$scaff_number=="3",]


recomb<-ggplot(chms1, aes(x=cumul_start, y=median_rho, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=chms2$cumul_start, y=chms2$median_rho), stat="identity", colour="maroon", size=0.5) +
  ylab("rho")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))+
  ylim(0,0.005)




###################################################
## Distribution genomic features all chromosomes ##
###################################################

library(cowplot)
plot_grid(
  plot_grid(
    TR
    , cenh3_testes + theme(legend.position = "none")
    , recomb + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(recomb)
    , ggplot()
    , get_legend(GC_content)
    , ncol =1)
  , rel_widths = c(12,-0.5)
)


plot_grid(
  plot_grid(
    TR
    , GC_content + theme(legend.position = "none")
    , GC_content_noTR + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(GC_content)
    , ggplot()
    , get_legend(GC_content_noTR)
    , ncol =1)
  , rel_widths = c(12,-0.5)
)


plot_grid(
  plot_grid(
    TR
    , cenh3_testes + theme(legend.position = "none")
    , h3k9_soma + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(h3k9_soma)
    , ncol =1)
  , rel_widths = c(12,-0.5)
)


##################################################
## Distribution genomic features per chromosome ##
##################################################

TR_chr<-subset(merge3, chr=="Tps_LRv5b_scf1")
cenh3_chr<-subset(chms9, Scaffold_name=="Tps_LRv5b_scf1")
h3k9_chr<-subset(chms7, Scaffold_name=="Tps_LRv5b_scf1")
GC_chr<-subset(chms3, V1=="Tps_LRv5b_scf1")
recomb_chr<-subset(chms1, scaffold=="Tps_LRv5b_scf1")


recomb<-ggplot(recomb_chr, aes(x=start, y=median_rho))+
  geom_point(size=0.5) +
  ylab("rho")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))+
  ylim(0,0.005)


TR<-ggplot(TR_chr, aes(x=wind_start, y=TRprop))+
  geom_point(size=0.5) +
  ylab("Proportion tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


h3k9<-ggplot(h3k9_chr, aes(x=start, y=log2H3k9_norm))+
  geom_point(size=0.5) +
  ylab("log2(H3K9me3/Input)")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


cenh3<-ggplot(cenh3_chr, aes(x=start, y=log2cenH3_norm))+
  geom_point(size=0.5) +
  ylab("log2(cenh3/Input)")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


GC<-ggplot(GC_chr, aes(x=V2, y=V4))+
  geom_point(size=0.5) +
  ylab("%GC")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))




library(cowplot)

theme_set(theme_minimal())


plot_grid(
  plot_grid(
   TR
    , recomb + theme(legend.position = "none")
    , GC + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(recomb)
    , ggplot()
    , get_legend(GC)
    , ncol =1)
  , rel_widths = c(12,0)
)

plot_grid(
  plot_grid(
    TR
    , cenh3 + theme(legend.position = "none")
    , h3k9 + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3)
    , ggplot()
    , get_legend(h3k9)
    , ncol =1)
  , rel_widths = c(12,0)
)


plot_grid(
  plot_grid(
    GC + theme(legend.position = "none")
    , recomb + theme(legend.position = "none")
    , TR
    , cenh3 + theme(legend.position = "none")
    , h3k9 + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(GC)
    , ggplot()
    , get_legend(recomb)
    , ggplot()
    , get_legend(cenh3)
    , ggplot()
    , get_legend(h3k9)
    , ncol =1)
  , rel_widths = c(12,0)
)



