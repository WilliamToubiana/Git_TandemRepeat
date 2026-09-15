##################################################
### Generate a presence (1) absence (0) matrix ###
##################################################
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")

Tps<-read.table("Tps_cleaned.txt", header=T, sep="\t", quote="")
Tcm<-read.table("Tcm_cleaned.txt", header=T, sep="\t", quote="")
Tce<-read.table("Tce_cleaned.txt", header=T, sep="\t", quote="")
Tpa<-read.table("Tpa_cleaned.txt", header=T, sep="\t", quote="")
Tbi<-read.table("Tbi_cleaned.txt", header=T, sep="\t", quote="")


# Combine all motifs from the 6 dataframes and get unique motifs
all_motifs <- unique(c(Tps$aomfi, Tcm$aomfi, Tce$aomfi, Tpa$aomfi, Tbi$aomfi))

# Sort the motifs
all_motifs <- sort(all_motifs)

# Create an empty matrix for presence-absence (motifs as rows, species as columns)
presence_absence_matrix <- matrix(0, nrow = length(all_motifs), ncol = 5)

# Assign row names (motifs) and column names (species)
rownames(presence_absence_matrix) <- all_motifs
colnames(presence_absence_matrix) <- c("Tps", "Tcm", "Tce", "Tpa", "Tbi")

# Fill the matrix with 1s for presence of motifs in each species
presence_absence_matrix[all_motifs %in% Tps$aomfi, 1] <- 1
presence_absence_matrix[all_motifs %in% Tcm$aomfi, 2] <- 1
presence_absence_matrix[all_motifs %in% Tce$aomfi, 3] <- 1
presence_absence_matrix[all_motifs %in% Tpa$aomfi, 4] <- 1
presence_absence_matrix[all_motifs %in% Tbi$aomfi, 5] <- 1

head(presence_absence_matrix)


#convert to dataframe 
presence_absence_df<-data.frame(presence_absence_matrix)
presence_absence_df$aomfi<-rownames(presence_absence_df)

#select only motifs present in a single species and sum rows
presence_absence_df_Tps<-subset(presence_absence_df, Tps=="1") #118243 motifs in total
presence_absence_df_Tps$shared<-rowSums( presence_absence_df_Tps[,c(1,2,3,4,5)])#shared across Timema (T.douglasi excluded)


#Upset plot
library(UpSetR)
presence_absence_df_Tps2<-presence_absence_df_Tps[1:5]

upset(
  presence_absence_df_Tps2,
  sets = colnames(presence_absence_df_Tps2),
  keep.order = TRUE,
  order.by = "freq")


#merge with T.poppense TR annotation
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


mergeTps<-merge(df_filtered2, presence_absence_df_Tps, by.x="motif_representant", by.y="aomfi", all.x =T)
summary(mergeTps)


#Recalculate copy number based on array size and length of representative TR sequences 
mergeTps$array_size<-abs(mergeTps$end-mergeTps$start)+1
mergeTps$copy_nb2<-mergeTps$array_size/mergeTps$L_motif_representant


#Assign to each repeat array a non-overlapping 250kb window
mergeTps$window_start <- (floor(mergeTps$start / 250000) * 250000)+1
mergeTps$window_end <- mergeTps$window_start + 249999  # End of the 250kb window
mergeTps$chm_pos<-paste(mergeTps$chr, mergeTps$window_start, sep="-")


#Calculate the number and proportion of shared repeats among species per 250kb windows
#change "mergeTps" by "mergeTps2" for the analysis without AACCT motifs

subset_chroms <- c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4",
                   "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8",
                   "Tps_LRv5b_scf9", "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12")


# First, compute for each sequence how many unique chromosomes (in the subset) it appears in
seq_chr_counts <- mergeTps %>%
  filter(chr %in% subset_chroms) %>%
  group_by(motif_representant) %>%
  summarise(chromosome_count_subset = n_distinct(chr), .groups = "drop")


# Now join that info back into the main summarization
window_summary <- mergeTps %>%
  left_join(seq_chr_counts, by = "motif_representant") %>%
  group_by(chr, window_start) %>%
  summarise(
    total_repeats = n(),
    shared_repeats = sum(shared >= 2),
    non_shared_repeats = sum(shared == "1"),
    avg_shared = mean(shared),
    median_shared = median(shared),
    proportion_shared = shared_repeats / total_repeats,
    proportion_non_shared = non_shared_repeats / total_repeats,
    avg_array_size = mean(array_size),
    median_array_size = median(array_size),
    avg_motif = mean(L_motif_representant),
    median_motif = median(L_motif_representant),
    median_CN = median(copy_nb2),
    mean_CN = mean(copy_nb2),
    motif_diversity = n_distinct(motif_representant),
    short_motif = sum(L_motif_representant < 10),
    # new column: average chromosome count among sequences in this window
    avg_chromosome_count = mean(chromosome_count_subset, na.rm = TRUE),
    median_chromosome_count = median(chromosome_count_subset, na.rm = TRUE),
    sum_shared_chr = sum(chromosome_count_subset == "12"),
    prop_shared_chr = sum(chromosome_count_subset == "12")/n(),
    .groups = "drop"
  )


summary(window_summary) #NA values correspond to regions with non shared sequences (between chromosomes) and that are placed on unanchored scaffolds

window_summary$chm_pos<-paste(window_summary$chr, window_summary$window_start, sep="-")



#########################
### ChIP cenh3 testes ###
#########################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
cenh3<-read.table("Tps_cenh3_testes_1_GW_coverage_DR_250kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_coverage_DR_100kb.txt

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
input1<-read.table("Tps_input_testes_1_GW_coverage_DR_250kb.txt",
                   header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data1<-cbind(cenh3,input1[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1$coverage_cenh3_norm<-data1$coverage_cenh3/66295120 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratiocenH3_norm<-data1$coverage_cenh3_norm/data1$coverage_input_norm
data1$log2cenH3_norm<-log2(data1$ratiocenH3_norm)
data1$start2<-data1$start+1
data1$chm_pos<-paste(data1$Scaffold_name, data1$start2, sep="-")
summary(data1)
data1$log2cenH3_norm[!is.finite(data1$log2cenH3_norm)] <- 0
data1$region_size<-abs(data1$stop-data1$start)



#######################################
### TR Proportion along chromosomes ###
#######################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW") #NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRprop<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation. I then assign them a value of 0.
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW
summary(mergeTRprop)



#############################################################
### Tcm vs Tps Illumina reads alignment onto Tps assembly ###
#############################################################

##TR coverage
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")
#Tps_cov<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="")
Tps_cov<-read.table("Tps-filtered_indG_to_Tps_genome_minimap2_TR_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tps_cov)<-c("chr", "start", "end", "sum_cov_tps")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tps_cov$sum_cov_tps<-as.numeric(Tps_cov$sum_cov_tps)
summary(Tps_cov)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/Tcm_reads_against_Tps_genome")
#Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="")
Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_TR_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tcm_cov)<-c("chr", "start", "end", "sum_cov_tcm")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tcm_cov$sum_cov_tcm<-as.numeric(Tcm_cov$sum_cov_tcm)
summary(Tcm_cov)


##Genome-wide coverage
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")
#Tps_cov<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") 
Tps_cov<-read.table("Tps-filtered_indG_to_Tps_genome_minimap2_GW_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tps_cov)<-c("chr", "start", "end", "sum_cov_tps")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tps_cov$sum_cov_tps<-as.numeric(Tps_cov$sum_cov_tps)
summary(Tps_cov)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/Tcm_reads_against_Tps_genome")
#Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="")
Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_GW_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tcm_cov)<-c("chr", "start", "end", "sum_cov_tcm")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tcm_cov$sum_cov_tcm<-as.numeric(Tcm_cov$sum_cov_tcm)
summary(Tcm_cov)


Tps_cov$start2<-Tps_cov$start+1
Tcm_cov$start2<-Tcm_cov$start+1
Tps_cov$chm_pos<-paste(Tps_cov$chr, Tps_cov$start2, sep="-")
Tcm_cov$chm_pos<-paste(Tcm_cov$chr, Tcm_cov$start2, sep="-")

mergeTpsTcm<-cbind(Tps_cov, Tcm_cov[4])
summary(mergeTpsTcm)
mergeTpsTcm[is.na(mergeTpsTcm)] <- 0 #NAs are retrieved when no reads are mapped on the entire scaffold.
mergeTpsTcm$sum_cov_tps<-mergeTpsTcm$sum_cov_tps+1 #I added +1 to all coverage values to prevent inf values on the ratios below. 
mergeTpsTcm$sum_cov_tcm<-mergeTpsTcm$sum_cov_tcm+1 #I added +1 to all coverage values to prevent inf values on the ratios below. 
summary(mergeTpsTcm)


mergeTpsTcm$sum_cov_tps_norm<-mergeTpsTcm$sum_cov_tps/341373875 #normalized by the number of mapped reads
mergeTpsTcm$sum_cov_tcm_norm<-mergeTpsTcm$sum_cov_tcm/547587902 #normalized by the number of mapped reads

mergeTpsTcm$ratio_norm<-mergeTpsTcm$sum_cov_tcm_norm/mergeTpsTcm$sum_cov_tps_norm
summary(mergeTpsTcm)

mergeTpsTcm$log2ratio_norm<-log2(mergeTpsTcm$ratio_norm)
summary(mergeTpsTcm)
mergeTpsTcm$start2<-mergeTpsTcm$start+1
mergeTpsTcm$chm_pos<-paste(mergeTpsTcm$chr, mergeTpsTcm$start2, sep="-")
#mergeTpsTcm$log2ratio_norm[!is.finite(mergeTpsTcm$log2ratio_norm)] <- 0
summary(mergeTpsTcm)




#########################################
### Prop TR aligned along chromosomes ###
#########################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_alignments")

TR_align_Tps<-read.table("proportion_tr_aligned_per_window.sorted.bed", header=T, sep="\t", quote="") #

TR_align_Tps$start2<-TR_align_Tps$start+1
TR_align_Tps$chm_pos<-paste(TR_align_Tps$chrom, TR_align_Tps$start2, sep="-")



#######################################
### Recombination along chromosomes ###
#######################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#filtering extreme
rho_sub <- rho %>%
  mutate(quantile=quantile(rho_weighted,0.99))%>%
  filter(rho_weighted<quantile)

#weighted rho per genetic distance
rho_sub$overlap<-rho_sub$right_SNP-rho_sub$left_SNP
rho_sub$rho_weighted<-rho_sub$rho*rho_sub$overlap

#Correction per window
rho_g <- rho_sub %>%
  group_by(scaffold,start)%>%
  summarize(corr_rho=sum(rho_weighted)/sum(overlap))


rho_g$start2<-rho_g$start+1
rho_g$chm_pos<-paste(rho_g$scaffold, rho_g$start2, sep="-")



##################
### GC content ###
##################


data_GC<-data_GC[c(1,2,3,13)]
colnames(data_GC)<-c("chromosome", "start", "end", "perc_GC")

data_GC$start2<-data_GC$start+1
data_GC$chm_pos<-paste(data_GC$chromosome, data_GC$start2, sep="-")

summary(data_GC)

##GC content
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("GC_content_Tps_250kb.txt", header=FALSE) #bedtools nuc -fi Tps_LRv5b_mtDNAv350.fasta -bed Tps_chm_size_mtDNAv350_w250000.bed  | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"($7+$8)/($6+$7+$8+$9+1)}' > GC_content_Tps_250kb.txt

#Sum per chr
data_GC2 <- data_GC %>%
  group_by(V1,V2)%>%
  summarize(V6=sum(V6),
            V7=sum(V7),
            V8=sum(V8),
            V9=sum(V9))

colnames(data_GC2)<-c("chromosome", "start", "A", "C", "G","T")
summary(data_GC2)

data_GC2$GC_perc<-(data_GC2$C+data_GC2$G)/(data_GC2$A+data_GC2$C+data_GC2$G+data_GC2$T)
summary(data_GC2)

data_GC2$start2<-data_GC2$start+1
data_GC2$chm_pos<-paste(data_GC2$chromosome, data_GC$start2, sep="-")
summary(data_GC2)



#################################################################
### TR sequence turnover centromere vs non-centromere regions ###
#################################################################

##Combine TR prop and shared repeat datasets
combined<-merge(mergeTRprop, window_summary[3:22], by="chm_pos", all.x=T)
summary(combined)
  
##Subdivide windows into centromere vs non-centromere based on log2 ratios
data1$centromere<-ifelse(data1$log2cenH3_norm>0.25, "centromere", "non-centromere")

##Combine TR prop, centromere and shared repeat datasets
combined0<-merge(combined, data1[c(9,11,12,13)], by="chm_pos", all.x=T)
summary(combined2)

##Combine TR prop, centromere, shared repeat and TR aligned datasets
combined1<-merge(combined0, TR_align_Tps[c(4,5,6,8)], by="chm_pos", all.x = T)
summary(combined3)

##Combine TR prop, centromere, shared repeat, TR aligned and Tcm mapping reads datasets
combined2<-merge(combined1, mergeTpsTcm[c(6,10,11)], by="chm_pos", all.x=T)
summary(combined4)

##Combine TR prop, centromere, shared repeat and recombination datasets
combined3<-merge(combined2, rho_g[c(3,5)], by="chm_pos", all.x = T)
summary(combined4)

##Combine TR prop, centromere, shared repeat, recombination and GC estimates datasets
combined4<-merge(combined3, data_GC2[c(7,9)], by="chm_pos", all.x = T)
summary(combined6)



##Scatterplot and density plots with TR proportion
##################################################
combined4_subset<-subset(combined4, chr!="Tps_mtDNA_v350")
summary(combined4_subset)

combined4_subset2<-subset(combined4, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                          | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                          | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")
summary(combined4_subset2)


mycols <- c("centromere" = "red",
  "non-centromere" = "black")

stats1 <- combined4_subset2 %>%
  group_by(centromere) %>%
  summarise(
    mean = mean(TRprop, na.rm = TRUE),
    median = median(TRprop, na.rm = TRUE))

stats2 <- combined4_subset2 %>%
  group_by(centromere) %>%
  summarise(
    mean = mean(log2ratio_norm, na.rm = TRUE),
    median = median(log2ratio_norm, na.rm = TRUE))

scatterPlot <- ggplot() +
  geom_point(data = combined4_subset2 %>% filter(centromere == "non-centromere"),
             aes(TRprop, log2ratio_norm),
             color = "black", alpha = 0.4) +
  geom_point(data = combined4_subset2 %>% filter(centromere == "centromere"),
             aes(TRprop, log2ratio_norm),
             color = "red", alpha = 0.6) +
  geom_smooth(data = combined4_subset2,
              aes(TRprop, log2ratio_norm, color = centromere),
              method = "lm", se = TRUE) +
  scale_color_manual(values = c("centromere" = "red",
                                "non-centromere" = "black")) +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(-3.5,-0.5)) +
  ylab("log2(Tcm/Tps)")+  xlab("TR proportion")


# Marginal density plot of x (top panel)
xdensity <- ggplot(combined4_subset2,
                   aes(TRprop, fill = centromere)) +
  
  geom_density(alpha = 0.5) +
  
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  
  ## Mean lines
  geom_vline(data = stats1,
             aes(xintercept = mean,
                 color = centromere),
             linetype = "dashed",
             linewidth = 0.5) +
  
  ## Median lines
  geom_vline(data = stats1,
             aes(xintercept = median,
                 color = centromere),
             linetype = "solid",
             linewidth = 0.5) +
  
  theme_bw() +
  xlab("TR proportion")+
  theme(legend.position = "none")


# Marginal density plot of y (right panel)
ydensity <- ggplot(combined4_subset2,
                   aes(log2ratio_norm, fill = centromere)) +
  
  geom_density(alpha = 0.5) +
  
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  
  ## Mean lines
  geom_vline(data = stats2,
             aes(xintercept = mean,
                 color = centromere),
             linetype = "dashed",
             linewidth = 0.5) +
  
  ## Median lines
  geom_vline(data = stats2,
             aes(xintercept = median,
                 color = centromere),
             linetype = "solid",
             linewidth = 0.5) +
  
  theme_bw() + coord_flip()+
  xlab("log2(Tcm/Tps)")+
  theme(legend.position = "none")


blankPlot <- ggplot() + 
  geom_blank(aes(1, 1)) +
  theme_void()


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))




##Stats
#coverage ratios

combined4_subset2$centromere <- as.factor(combined4_subset2$centromere)
combined4_subset2$centromere <- relevel(
  combined4_subset2$centromere,
  ref = "non-centromere")

library(lmerTest)
model<-lmer(log2ratio_norm~centromere*TRprop+(1|chr),data=combined4_subset2)
summary(model)


#perc aligned TRs
combined4_subset2$centromere <- as.factor(combined4_subset2$centromere)
combined4_subset2$centromere <- relevel(
  combined4_subset2$centromere,
  ref = "non-centromere")

n <- nrow(combined4_subset2)
combined4_subset2$prop2 <- (combined4_subset2$prop * (n - 1) + 0.5) / n #shifts "0" into tiny positive values and "1" into slightly below 1 values 

library(glmmTMB)
model <- glmmTMB(prop2 ~ centromere * TRprop + (1|chr),
     family = beta_family(),
     data = combined4_subset2)

summary(model)


#proportion shared TR sequences (Tps-Tcm)
combined4_subset2$centromere <- as.factor(combined4_subset2$centromere)
combined4_subset2$centromere <- relevel(
  combined4_subset2$centromere,
  ref = "non-centromere")

n <- nrow(combined6_subset2)
combined4_subset2$proportion_shared2 <- (combined4_subset2$proportion_shared * (n - 1) + 0.5) / n #shifts "0" into tiny positive values and "1" into slightly below 1 values 

library(glmmTMB)
model <- glmmTMB(proportion_shared2 ~ centromere * TRprop + (1|chr),
                 family = beta_family(),
                 data = combined4_subset2)

summary(model)


model <- glmmTMB(total_repeats ~ centromere * TRprop + (1|chr),
                 family = nbinom2(),
                 data = combined4_subset2)

summary(model)


#nb TR arrays and shared sequences (Tps-Tcm)
combined4_subset2$centromere <- as.factor(combined4_subset2$centromere)
combined4_subset2$centromere <- relevel(
  combined4_subset2$centromere,
  ref = "non-centromere")

library(glmmTMB)
model <- glmmTMB(total_repeats ~ centromere * TRprop + (1|chr),
                 family = nbinom2(),
                 data = combined4_subset2)

summary(model)

model <- glmmTMB(shared_repeats ~ centromere * TRprop + (1|chr),
                 family = nbinom2(),
                 data = combined4_subset2)

summary(model)


##Subdividing TR prop into categories

combined4_subset2$TRcat <- cut(
  combined4_subset2$TRprop,
  breaks = c(0, 1/3, 2/3, 1),
  labels = c("low", "moderate", "high"),
  include.lowest = TRUE)

combined4_subset3<-combined4_subset2[!is.na(combined4_subset2$TRcat),]

ggplot(combined4_subset3,
       aes(x = TRcat,
           y = log2ratio_norm,
           fill = centromere)) +
  geom_boxplot(outlier.alpha = 0.2,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c('black','red'))+
  labs(x = "Categories of TR proportion",
       y = "log2(Tcm/Tps)",
       fill = "Region") +
  theme_bw()


ggplot(combined4_subset3,
       aes(x = TRcat,
           y = avg_array_size,
           fill = centromere))+ 
  geom_boxplot(width = 0.5,
               outlier.alpha = 0.2,
               position = position_dodge()) +
  scale_fill_manual(values = c("black", "red")) +
  labs(x = "Categories of TR proportion",
       y = "mean TR array size (bp)",
       fill = "Region") +
  theme_bw()


##Median and means
library(dplyr)

combined4_subset3 %>%
  group_by(TRcat, centromere) %>%
  summarise(
    median_log2 = median(log2ratio_norm, na.rm = TRUE),
    mean_log2 = mean(log2ratio_norm, na.rm = TRUE),
    n = n()
  )


#########################################################
## Distribution TR sequence turnover along chromosomes ##
#########################################################

combined4_subset2<-subset(combined4, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                                | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                                | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined4_subset2$scaff_number<-as.numeric(as.character(gsub("^.*scf","", combined4_subset2$chr)))

combined4_subset2<-combined4_subset2[order(combined4_subset2$scaff_number, combined4_subset2$start2),]
combined4_subset2$cumul_start<-seq(0, by = 250000, length.out = nrow(combined4_subset2))

# Calculate scaffold boundaries (start and end positions)
library(dplyr)
scaffold_bounds <- combined4_subset2 %>%
  group_by(scaff_number) %>%
  summarise(min_pos = min(cumul_start),
            max_pos = max(cumul_start)) %>%
  mutate(fill_col = rep(c("white", "grey90"), length.out = n()))



TR_turnover<-ggplot(combined4_subset2, aes(x=cumul_start, y=log2ratio_norm, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size=0.5, alpha = 0.5) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "black"))+
  ylab("log2(Tcm/Tps)")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.6, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



cenh3_testes<-ggplot(combined4_subset2, 
                     aes(x = cumul_start, y = log2cenH3_norm, color = centromere)) +
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "black")) +
  geom_smooth(aes(group = scaff_number), size = 0.5, color = "white",
              method = "gam", formula = y ~ s(x, bs = "cs", k = 30)) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ylab("log2(CenH3/Input)") +
  xlab("Genomic coordinates") +
  ylim(-1, 3)


TR_array<-ggplot(combined4_subset2, aes(x=cumul_start, y=TRprop, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.4) +
  scale_fill_identity() +
  geom_point(size=0.5, alpha = 0.5) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "black")) +
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


# Plots TRprop vs cenh3 vs avg_shared
library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    TR_array + theme(legend.position = "none")
    , cenh3_testes
    , TR_turnover + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(TR_turnover)
    , ncol =1)
  , rel_widths = c(12,0)
)



##################################################################################################################
### Iterations and Subsampling of "non-centromere" windows to match "centromere" distribution of TR proportion ###
##################################################################################################################

library(dplyr)
library(broom)
library(ggplot2)

combined4_subset2<-subset(combined4, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                          | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                          | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")


#set.seed(42)    # reproducible overall

# Parameters
n_iter <- 1000        # number of random subsamples (change to 1000 if you want)
n_bins  <- 12        # bins for TRprop quantile-matching (as before)

# Split data
cen    <- combined4_subset2 %>% filter(centromere == "centromere")
noncen <- combined4_subset2 %>% filter(centromere == "non-centromere")

# Precompute quantile breaks (same bins used every iteration)
breaks <- quantile(combined4_subset2$TRprop, probs = seq(0, 1, length.out = n_bins + 1),
                   na.rm = TRUE, type = 7)

# Assign bins to both sets (keeps binning constant across iterations)
cen <- cen %>% mutate(bin = cut(TRprop, breaks = breaks, include.lowest = TRUE))
noncen <- noncen %>% mutate(bin = cut(TRprop, breaks = breaks, include.lowest = TRUE))

# How many to sample per bin (target = counts in centromere)
target_per_bin <- table(cen$bin)

# quick helper to sample from noncen for a given bin counts (handles small bins)
sample_noncen_for_bin <- function(bin_name, n_target) {
  sub <- noncen %>% filter(bin == bin_name)
  n_avail <- nrow(sub)
  if (n_avail == 0) {
    # no available non-centromere in this bin: return 0 rows
    return(sub[0, ])
  }
  # if fewer available than needed, sample with replacement and warn once (but do not stop)
  replace_flag <- ifelse(n_avail < n_target, TRUE, FALSE)
  sub %>% sample_n(size = n_target, replace = replace_flag)
}

# Precompute centromere model (doesn't change each iteration)
#to modify depending on the metric (here log2 coverage ratio is shown)
lm_cen <- lm(log2ratio_norm ~ TRprop, data = cen)
lm_cen_sum <- summary(lm_cen)
cen_stats <- list(
  slope = coef(lm_cen)["TRprop"],
  intercept = coef(lm_cen)["(Intercept)"],
  r2 = lm_cen_sum$r.squared,
  mean_seq_id = mean(cen$log2ratio_norm, na.rm = TRUE),
  n = nrow(cen)
)

# Container for iteration results
results <- vector("list", n_iter)

# Main loop: iterate subsampling + fit model on sampled noncen
for (i in seq_len(n_iter)) {
  # build sampled non-centromere dataframe by sampling within each bin
  sampled_list <- mapply(
    sample_noncen_for_bin,
    names(target_per_bin),
    as.integer(target_per_bin),
    SIMPLIFY = FALSE
  )
  noncen_sub <- bind_rows(sampled_list)
  
  # In case some bins had zero available noncen, the total may be < nrow(cen).
  # You can choose to skip this iteration or proceed — here we proceed and record the actual n.
  # Fit lm on sampled noncentromeres (if enough points)
  if (nrow(noncen_sub) >= 3) {
    lm_non <- lm(log2ratio_norm ~ TRprop, data = noncen_sub)
    s <- summary(lm_non)
    non_stats <- list(
      slope = coef(lm_non)["TRprop"],
      intercept = coef(lm_non)["(Intercept)"],
      r2 = s$r.squared,
      mean_log2 = mean(noncen_sub$log2ratio_norm, na.rm = TRUE),
      n = nrow(noncen_sub)
    )
  } else {
    non_stats <- list(slope = NA, intercept = NA, r2 = NA, mean_log2 = NA, n = nrow(noncen_sub))
  }
  
  results[[i]] <- tibble(
    iter = i,
    cen_slope = cen_stats$slope,
    cen_intercept = cen_stats$intercept,
    cen_r2 = cen_stats$r2,
    cen_mean_log2 = cen_stats$mean_log2,
    cen_n = cen_stats$n,
    non_slope = non_stats$slope,
    non_intercept = non_stats$intercept,
    non_r2 = non_stats$r2,
    non_mean_log2 = non_stats$mean_log2,
    non_n = non_stats$n,
    slope_diff = non_stats$slope - cen_stats$slope
  )
}

res_df <- bind_rows(results)  
summary(res_df)

# Summaries
summary_stats <- res_df %>%
  summarize(
    non_slope_mean = mean(non_slope, na.rm = TRUE),
    non_slope_sd = sd(non_slope, na.rm = TRUE),
    slope_diff_mean = mean(slope_diff, na.rm = TRUE),
    slope_diff_sd = sd(slope_diff, na.rm = TRUE),
    prop_non_less_than_cen = mean(non_slope < cen_slope, na.rm = TRUE),
    median_non_n = median(non_n, na.rm = TRUE)
  )

print(summary_stats)

# Quick diagnostics plots
p1 <- ggplot(res_df, aes(x = non_slope)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = cen_stats$slope, color = "red", linetype = "dashed") +
  ggtitle("Distribution of non-centromere slopes across iterations\n(vertical dashed = centromere slope)") +
  xlab("slope (log2(Tcm/Tps) ~ TRprop)")

p2 <- ggplot(res_df, aes(x = non_intercept)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = cen_stats$intercept, color = "red", linetype = "dashed") +
  ggtitle("Distribution of non-centromere intercepts across iterations\n(vertical dashed = centromere intercept") +
  xlab("intercept (log2(Tcm/Tps) ~ TRprop)")

# Print plots
print(p1); print(p2)
 



#############################
## Relationship TR metrics ##
#############################

combined4_subset2<-subset(combined4, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                          | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                          | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined4_subset2$region<-combined6_subset2$end-combined4_subset2$start
combined4_subset2$array_density<-combined4_subset2$total_repeats/combined4_subset2$region

summary(combined4_subset2)

#change "median_rho" by "perc_GC" to get the association with GC content
ggplot(combined4_subset2,aes(log10(corr_rho), avg_motif)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  scale_y_continuous(limits = c(0, NA))+
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length")+
  xlab("rho values (log)")

model<-lm(combined4_subset2$avg_motif~log10(combined4_subset2$corr_rho))
summary(model)


ggplot(combined4_subset2,aes(log10(corr_rho), median_motif)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  scale_y_continuous(limits = c(0, NA))+
  theme_classic(base_size = 17)+
  ylab("median TR sequence length")+
  xlab("rho values (log)")

model<-lm(combined4_subset2$median_motif~log10(combined4_subset2$corr_rho))
summary(model)


ggplot(combined4_subset2,aes(log10(corr_rho), mean_CN)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  scale_y_continuous(limits = c(0, NA))+
  ylab("mean copy number")+
  xlab("rho values (log)")

model<-lm(combined4_subset2$mean_CN~log10(combined4_subset2$corr_rho))
summary(model)


ggplot(combined4_subset2,aes(log10(corr_rho), log10(median_CN))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  scale_y_continuous(limits = c(0, NA))+
  ylab("median copy number")+
  xlab("rho values (log)")

model<-lm(log10(combined4_subset2$median_CN)~log10(combined4_subset2$corr_rho))
summary(model)


ggplot(combined4_subset2,aes(log10(corr_rho), array_density)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  scale_y_continuous(limits = c(0, NA))+
  theme_classic(base_size = 17)+ 
  ylab("TR array density")+
  xlab("rho values (log)")

model<-lm(log10(combined4_subset2$array_density)~log10(combined4_subset2$corr_rho))
summary(model)


combined4_subset2$motif_diversity_corrected<-combined4_subset2$motif_diversity/combined4_subset2$total_repeats
ggplot(combined4_subset2,aes(log10(corr_rho), log10(motif_diversity_corrected))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR sequence diversity / total TR sequence")+
  xlab("TR proportion")

model<-lm(log10(combined4_subset2$motif_diversity_corrected)~log10(combined4_subset2$corr_rho))
summary(model)

