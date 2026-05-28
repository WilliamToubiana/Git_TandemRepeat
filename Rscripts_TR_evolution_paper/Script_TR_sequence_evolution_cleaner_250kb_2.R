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
presence_absence_df_Tps$shared<-rowSums( presence_absence_df_Tps[,c(1,2)])#shared with T.californicum only


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


# Now join that info back into your main summarization
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




#######################################
### Recombination along chromosomes ###
#######################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho","chr", "window_start", "window_end")

#Assign to each pair of SNP a non-overlapping 250kb window
rho$window_start <- (floor(rho$left_SNP / 250000) * 250000)+1
rho$window_end <- rho$window_start + 249999  # End of the 250kb window

#Median per window
rho_g <- rho %>%
  group_by(scaffold,window_start)%>%
  summarize(median_rho=median(rho))

rho_g$chm_pos<-paste(rho_g$scaffold, rho_g$window_start, sep="-")




##################
### GC content ###
##################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("GC_content_Tps_250kb.txt", header=FALSE) #bedtools nuc -fi Tps_LRv5b_mtDNAv350.fasta -bed Tps_chm_size_mtDNAv350_w250000.bed  | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"($7+$8)/($6+$7+$8+$9+1)}' > GC_content_Tps_250kb.txt
colnames(data_GC)<-c("chromosome", "start", "end", "perc_GC")

data_GC$start2<-data_GC$start+1
data_GC$chm_pos<-paste(data_GC$chromosome, data_GC$start2, sep="-")

summary(data_GC)



#############################################################
### Tcm vs Tps Illumina reads alignment onto Tps assembly ###
#############################################################

##TR coverage
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")
Tps_cov<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="")
Tps_cov<-read.table("Tps-filtered_indG_to_Tps_genome_minimap2_TR_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tps_cov)<-c("chr", "start", "end", "sum_cov_tps")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tps_cov$sum_cov_tps<-as.numeric(Tps_cov$sum_cov_tps)
summary(Tps_cov)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/Tcm_reads_against_Tps_genome")
Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="")
Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_TR_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tcm_cov)<-c("chr", "start", "end", "sum_cov_tcm")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tcm_cov$sum_cov_tcm<-as.numeric(Tcm_cov$sum_cov_tcm)
summary(Tcm_cov)


##Genome-wide coverage
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")
Tps_cov<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") 
Tps_cov<-read.table("Tps-filtered_indG_to_Tps_genome_minimap2_GW_coverage_mean250kb.bed", header=F, sep="\t", quote="")
colnames(Tps_cov)<-c("chr", "start", "end", "sum_cov_tps")#NAs are retrieved when no reads are mapped on the entire scaffold. 0 are retrieved when no reads are mapped in a window but reads are mapped on another window from the scaffold.
Tps_cov$sum_cov_tps<-as.numeric(Tps_cov$sum_cov_tps)
summary(Tps_cov)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/Tcm_reads_against_Tps_genome")
Tcm_cov<-read.table("Tcm-filtered_indG_to_Tps_genome_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="")
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


#######################
### Gene annotation ###
#######################
setwd("/Users/wtoubian/Desktop/Gene_position/new annotations")

gene_annotation<-read.table("Tps.v1.gtf", header=F, sep="\t", quote="")
colnames(gene_annotation)<-c("chr", "annotation", "category", "start", "end", "dot", "plus", "number", "category2")

genes<-subset(gene_annotation, category=="gene")

genes$window_start <- (floor(genes$start / 250000) * 250000)+1
genes$window_end <- genes$window_start + 249999  # End of the 250kb window
genes$count<-1

test <- genes %>%
  group_by(chr, window_start) %>%
  summarise(
    total_genes = n(),
    density     = sum(count) / 250000,
    .groups = "drop"
  ) 

test$chm_pos<-paste(test$chr, test$window_start, sep="-")

test3<-merge(combined6, test[2:5], by="chm_pos", all.x = T)

test4<-subset(test3, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
              | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
              | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

summary(test4)

test4[is.na(test4)] <- 0 #NAs are retrieved when windows have no gene annotation

test4$scaff_number<-as.numeric(as.character(gsub("^.*scf","", test4$chr)))
test4<-test4[order(test4$scaff_number, test4$start2),]
test4$cumul_start<-seq(0, by = 250000, length.out = nrow(test4))

ggplot(test4, aes(x=cumul_start, y=sum_counts, color = centromere))+
  scale_fill_identity() +
  geom_point(size=0.5, alpha = 0.6) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey"))+
  ylab("log2(Tcm/Tps)")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


#########################################
### Prop TR aligned along chromosomes ###
#########################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_alignments")

TR_align_Tps<-read.table("proportion_tr_aligned_per_window.sorted.bed", header=T, sep="\t", quote="") #

TR_align_Tps$start2<-TR_align_Tps$start+1
TR_align_Tps$chm_pos<-paste(TR_align_Tps$chrom, TR_align_Tps$start2, sep="-")





#################################################################
### TR sequence turnover centromere vs non-centromere regions ###
#################################################################

##Combine TR prop and shared repeat datasets
combined<-merge(mergeTRprop, window_summary[3:20], by="chm_pos", all.x=T)
summary(combined)
  
##Subdivide windows into centromere vs non-centromere based on log2 ratios
data1$centromere<-ifelse(data1$log2cenH3_norm>0.25, "centromere", "non-centromere")

##Combine TR prop, centromere and shared repeat datasets
combined2<-merge(combined, data1[c(9,11,12,13)], by="chm_pos", all.x=T)
summary(combined2)

##Combine TR prop, centromere, shared repeat and recombination datasets
combined3<-merge(combined2, rho_g[c(3,4)], by="chm_pos", all.x = T)
summary(combined3)

##Combine TR prop, centromere, shared repeat, recombination and TR aligned datasets
combined4<-merge(combined3, TR_align_Tps[c(4,5,6,8)], by="chm_pos", all.x = T)
summary(combined4)

##Combine TR prop, centromere, shared repeat, recombination, TR aligned and Tcm mapping reads datasets
combined5<-merge(combined4, mergeTpsTcm[c(6,11)], by="chm_pos", all.x=T)
summary(combined5)

##Combine TR prop, centromere, shared repeat, recombination , TR aligned, Tcm mapping reads and GC estimates datasets
combined6<-merge(combined5, data_GC[c(4,6)], by="chm_pos", all.x = T)
summary(combined6)

##remove windows without recombination estimates
combined7<-na.omit(combined5)
summary(combined7)





##Scatterplot and density plots with TR proportion
##################################################
combined6_subset<-subset(combined6, chr!="Tps_mtDNA_v350")
summary(combined6_subset)
combined6_subset2<-subset(combined6_subset, region_size>=10000)
combined6_subset2<-combined6_subset2[!is.na(combined6_subset2$TRprop),]
summary(combined6_subset2)

combined6_subset2<-subset(combined6, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                          | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                          | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")
summary(combined6_subset2)

scatterPlot <- ggplot() +
  geom_point(data = combined6_subset2 %>% filter(centromere == "non-centromere"),
             aes(TRprop, log2ratio_norm),
             color = "black", alpha = 0.4) +
  geom_point(data = combined6_subset2 %>% filter(centromere == "centromere"),
             aes(TRprop, log2ratio_norm),
             color = "red", alpha = 0.6) +
  geom_smooth(data = combined6_subset2,
              aes(TRprop, log2ratio_norm, color = centromere),
              method = "lm", se = TRUE) +
  scale_color_manual(values = c("centromere" = "red",
                                "non-centromere" = "black")) +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(-5.5,-1)) +
  ylab("log2(Tcm/Tps)")

# Marginal density plot of x (top panel)
xdensity <- ggplot(combined6_subset2, aes(TRprop, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('red','black')) +
  theme_bw() +
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined6_subset2, aes(log2ratio_norm, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('red','black')) + 
  theme_bw() +
  theme(legend.position = "none")+coord_flip()+
  xlab("log2(Tcm/Tps)")


blankPlot <- ggplot() + 
  geom_blank(aes(1, 1)) +
  theme_void()


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))



##Stats
#coverage ratios
library("car")
qqPlot(combined6_subset2$TRprop)
qqPlot(combined6_subset2$log2ratio_norm)

library(npANCOVA)
combined6_subset2[is.na(combined6_subset2)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
results<-McSweeny_Porter(combined6_subset2, log2ratio_norm~TRprop+centromere)
print(results$regression_equation_interaction)
print(results$group_effect)
print(results$interaction_effect)

combined6_subset2$centromere <- as.factor(combined6_subset2$centromere)
combined6_subset2$centromere <- relevel(
  combined6_subset2$centromere,
  ref = "non-centromere")

library(lme4)
model<-lmer(log2ratio_norm~centromere*TRprop+(1|chr),data=combined6_subset2)
summary(model)


library(effectsize)
eta_squared(model, partial = TRUE)


#perc aligned TRs
library("car")
qqPlot(combined6_subset2$TRprop)
qqPlot(combined6_subset2$prop)

library(npANCOVA)
results<-McSweeny_Porter(combined6_subset2, prop~TRprop+centromere)
print(results$regression_equation_covariate_factor)
print(results$group_effect)
print(results$interaction_effect)
print(results$regression_equation_interaction)

model2<-lm(combined6_subset2$prop~combined6_subset2$TRprop*combined6_subset2$centromere)
summary(model2)
anova(model2)





#########################################################
## Distribution TR sequence turnover along chromosomes ##
#########################################################

combined6_subset3<-subset(combined6, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                                | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                                | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined6_subset3$scaff_number<-as.numeric(as.character(gsub("^.*scf","", combined6_subset3$chr)))

combined6_subset3<-combined6_subset3[order(combined6_subset3$scaff_number, combined6_subset3$start2),]
combined6_subset3$cumul_start<-seq(0, by = 250000, length.out = nrow(combined6_subset3))

# Calculate scaffold boundaries (start and end positions)
library(dplyr)
scaffold_bounds <- combined6_subset3 %>%
  group_by(scaff_number) %>%
  summarise(min_pos = min(cumul_start),
            max_pos = max(cumul_start)) %>%
  mutate(fill_col = rep(c("white", "grey90"), length.out = n()))



TR_turnover<-ggplot(combined6_subset3, aes(x=cumul_start, y=log2ratio_norm, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size=0.5, alpha = 0.6) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey"))+
  ylab("log2(Tcm/Tps)")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



cenh3_testes<-ggplot(combined6_subset3, 
                     aes(x = cumul_start, y = log2cenH3_norm, color = centromere)) +
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey")) +
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


TR_array<-ggplot(combined6_subset3, aes(x=cumul_start, y=TRprop, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.4) +
  scale_fill_identity() +
  geom_point(size=0.5) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey")) +
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

#set.seed(42)    # reproducible overall

# Parameters
n_iter <- 1000        # number of random subsamples (change to 1000 if you want)
n_bins  <- 12        # bins for TRprop quantile-matching (as before)

# Split data
cen    <- combined6_subset2 %>% filter(centromere == "centromere")
noncen <- combined6_subset2 %>% filter(centromere == "non-centromere")

# Precompute quantile breaks (same bins used every iteration)
breaks <- quantile(combined6_subset2$TRprop, probs = seq(0, 1, length.out = n_bins + 1),
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
  xlab("slope (log2cenH3_norm ~ TRprop)")

p2 <- ggplot(res_df, aes(x = non_mean_log2)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = cen_stats$mean_seq_id, color = "red", linetype = "dashed") +
  ggtitle("Distribution of mean(log2cenH3_norm) for sampled non-centromere sets") +
  xlab("mean(log2cenH3_norm)")

# Print plots
print(p1); print(p2)




########################################################
### Choose representative iteration from subsampling ###
########################################################

library(dplyr)

#set.seed(123)

# Split data
cen  <- combined6_subset2 %>% filter(centromere == "centromere")
noncen <- combined6_subset2 %>% filter(centromere == "non-centromere")

# Create TRprop bins
n_bins <- 12   # adjust for smoother matching
breaks <- quantile(combined6_subset2$TRprop, probs = seq(0, 1, length.out = n_bins + 1))

cen$bin <- cut(cen$TRprop, breaks = breaks, include.lowest = TRUE)
noncen$bin <- cut(noncen$TRprop, breaks = breaks, include.lowest = TRUE)

# Count how many centromere points in each bin
cen_bin_counts <- table(cen$bin)

# Sample non-centromere points to match counts
noncen_sub <- noncen %>%
  group_by(bin) %>%
  sample_n(size = cen_bin_counts[bin][1], replace = FALSE) %>%
  ungroup()

# Combine subsampled non-centromeres with all centromeres
matched <- bind_rows(cen, noncen_sub)


#Plots
scatterPlot <- ggplot(matched, aes(TRprop, log2ratio_norm, color = centromere)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c('red','grey')) +
  geom_smooth(method = lm) +
  theme_bw() +
  theme(legend.position = c(0,1), legend.justification = c(-5.5, -1)) +
  ylab("log2(Tcm/Tps)")

xdensity <- ggplot(matched, aes(TRprop, fill = centromere)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c('red','black')) +
  theme_bw() +
  theme(legend.position = "none")

ydensity <- ggplot(matched, aes(log2ratio_norm, fill = centromere)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c('red','black')) +
  theme_bw() +
  theme(legend.position = "none") +
  coord_flip() +
  xlab("log2(Tcm/Tps)")

grid.arrange(
  xdensity, blankPlot,
  scatterPlot, ydensity,
  ncol = 2, nrow = 2,
  widths = c(4, 1.4),
  heights = c(1.4, 4)
)




#############################
## Relationship TR metrics ##
#############################
#change the x-axis by TRprop and perc_GC

#combined6_subset$region<-combined6_subset$end-combined6_subset$start
#combined6_subset$array_density<-combined6_subset$total_repeats/combined6_subset$region
combined6_subset<-subset(combined6, chr!="Tps_mtDNA_v350")
summary(combined6_subset)

combined6_subset2<-subset(combined6_subset, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                          | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                          | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")
summary(combined6_subset2)

combined6_subset2<- subset(combined6_subset2, TRprop > 0.15)

ggplot(combined6_subset2,aes(log10(TRprop), log10(avg_motif))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  scale_y_continuous(limits = c(0, NA))+
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length")+
  xlab("TR proportion")

model<-lm(log10(combined6_subset2$avg_motif)~log10(combined6_subset2$TRprop))
summary(model)


ggplot(combined6_subset2,aes(log10(TRprop), log10(median_motif))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  scale_y_continuous(limits = c(0, NA))+
  theme_classic(base_size = 17)+
  ylab("median TR sequence length")+
  xlab("TR proportion")

model<-lm(log10(combined6_subset2$median_motif)~log10(combined6_subset2$TRprop))
summary(model)


ggplot(combined6_subset2,aes(log10(TRprop), log10(mean_CN))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  scale_y_continuous(limits = c(0, NA))+
  ylab("mean copy number")+
  xlab("TR proportion")

model<-lm(log10(combined6_subset2$mean_CN)~log10(combined6_subset2$TRprop))
summary(model)


ggplot(combined6_subset2,aes(log10(TRprop), log10(median_CN))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  scale_y_continuous(limits = c(0, NA))+
  ylab("median copy number")+
  xlab("TR proportion")

model<-lm(log10(combined6_subset2$median_CN)~log10(combined6_subset2$TRprop))
summary(model)


ggplot(combined6_subset2,aes(log10(TRprop), log10(total_repeats))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  scale_y_continuous(limits = c(0, NA))+
  theme_classic(base_size = 17)+ 
  ylab("total TR sequences")+
  xlab("TR proportion")

model<-lm(log10(combined6_subset2$total_repeats)~log10(combined6_subset2$TRprop))
summary(model)


combined6_subset2$motif_diversity_corrected<-combined6_subset2$motif_diversity/combined6_subset2$total_repeats
ggplot(combined6_subset2,aes(log10(TRprop), log10(motif_diversity_corrected))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR sequence diversity / total TR sequence")+
  xlab("TR proportion")

model<-lm(log10(combined6_subset2$motif_diversity_corrected)~log10(combined6_subset2$TRprop))
summary(model)




##############################################################################
## Distribution of between-chromosome TR sequence sharing along chromosomes ##
##############################################################################
combined6_subset2<-subset(combined6, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                          | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                          | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined6_subset2$chr2 <- paste0("chr", sub(".*scf", "", combined6_subset2$chr))

combined6_subset2$chr2<-factor(combined6_subset2$chr2,levels = c("chr1", "chr2", "chr3",
                                                                           "chr4", "chr5", "chr6",
                                                                           "chr7", "chr8", "chr9",
                                                                           "chr10", "chr11", "chr12"))


ggplot(combined6_subset2, aes(TRprop, prop_shared_chr, color = avg_motif)) +
  geom_point(alpha = 0.6, size = 1.8) +
  facet_wrap(~ chr2, scales = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(face = "bold")
  ) +
  scale_color_viridis_c(option = "H") +
  labs(
    x = "TR proportion",
    y = "Proportion of shared TRs",
    color = "Avg. motif length",
  )


ggplot(combined6_subset2, aes(TRprop, sum_shared_chr, color = avg_motif)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    color = "darkgrey",
    linewidth = 0.8
  ) +
  facet_wrap(~ chr2, scales = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(face = "bold")
  ) +
  scale_color_viridis_c(option = "H") +
  labs(
    x = "TR proportion",
    y = "Number of shared TRs",
    color = "Avg. motif length",
  )
##############################################################################


















##Scatterplot and density plots with recombination
##################################################
combined7_subset<-combined7
summary(combined7_subset)
combined7_subset<-subset(combined7, chr=="Tps_LRv5b_scf1")
combined7_subset<-subset(combined7, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                         | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                         | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined7_subset2<-subset(combined7_subset, avg_shared>0)
combined7_subset2<-subset(combined7_subset2, TRprop>0.1)

#filtering extreme
combined7_subset3 <- combined7_subset2 %>%
  mutate(quantile=quantile(median_rho,0.95))%>%
  filter(median_rho<quantile)

combined3_subset3<-subset(combined3_subset3, TRprop>0.1)

scatterPlot <- ggplot(combined3_subset3,aes(median_rho, avg_shared, color=centromere)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) +
  geom_smooth(method=lm) +
  theme(legend.position=c(0,1), legend.justification=c(2,1))+
  ylab("Avg. Number of Species Sharing TR Sequences")


# Marginal density plot of x (top panel)
xdensity <- ggplot(combined3_subset3, aes(median_rho, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined3_subset3, aes(avg_shared, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+coord_flip()+
  ylab("Avg. Number of Species Sharing TR Sequences")




blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))



#Stats
model<-lm(combined3_subset2$proportion_species_specific~combined3_subset2$median_rho*combined3_subset2$centromere)
anova(model)

model<-lm(combined3_subset2$proportion_species_specific~combined3_subset2$median_rho+combined3_subset2$centromere)
anova(model)



##Correlations TR metrics
#########################
#change the x-axis by TRprop and perc_GC

ggplot(combined7_subset3,aes(median_rho, total_repeats)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("total TR sequences")+
  xlab("rho")

model<-lm(combined7_subset3$total_repeats~combined7_subset3$median_rho)
summary(model)


ggplot(combined7_subset3,aes(median_rho, avg_motif)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length")+
  xlab("rho")

model<-lm(combined7_subset3$avg_motif~combined7_subset3$median_rho)
summary(model)


ggplot(combined7_subset3,aes(median_rho, median_motif)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("median TR sequence length")+
  xlab("rho")

model<-lm(combined7_subset3$median_motif~combined7_subset3$median_rho)
summary(model)


ggplot(combined7_subset3,aes(median_rho, mean_CN)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ ylim(0,300)+
  ylab("mean copy number")+
  xlab("rho")

model<-lm(combined7_subset3$mean_CN~combined7_subset3$median_rho)
summary(model)


ggplot(combined7_subset3,aes(median_rho, median_CN)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("median copy number")+
  xlab("rho")

model<-lm(combined7_subset3$median_CN~combined7_subset3$median_rho)
summary(model)


combined7_subset3$motif_diversity_corrected<-combined7_subset3$motif_diversity/combined7_subset3$total_repeats
ggplot(combined7_subset3,aes(median_rho, motif_diversity_corrected)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR sequence diversity / total TR sequence")+
  xlab("rho")

model<-lm(combined7_subset3$motif_diversity_corrected~combined7_subset3$median_rho)
summary(model)
















##Scatterplot and density plots with GC content
###############################################
combined4_subset<-combined4
combined4_subset2<-subset(combined4_subset, avg_shared>0)
combined4_subset2<-subset(combined4_subset2, TRprop>0.1)
combined4_subset3<-subset(combined4_subset2, perc_GC>0.3 & perc_GC<0.5)

#filtering extreme
combined4_subset3 <- combined4_subset2 %>%
  mutate(quantile=quantile(perc_GC,0.99))%>%
  mutate(quantile2=quantile(perc_GC,0.01))%>%
  filter(perc_GC<quantile)%>%
  filter(perc_GC>quantile2)
  

scatterPlot <- ggplot(combined4_subset3,aes(perc_GC, avg_shared, color=centromere)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) +
  geom_smooth(method=lm) +
  theme(legend.position=c(0,0), legend.justification=c(2,1))+ ylim(1, 5)

# Marginal density plot of x (top panel)
xdensity <- ggplot(combined4_subset3, aes(perc_GC, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined4_subset3, aes(avg_shared, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+coord_flip()




blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))



#Stats
model<-lm(combined3_subset2$proportion_species_specific~combined3_subset2$median_rho*combined3_subset2$centromere)
anova(model)

model<-lm(combined3_subset2$proportion_species_specific~combined3_subset2$median_rho+combined3_subset2$centromere)
anova(model)

ggplot(combined4_subset, aes(x=centromere, y=proportion_species_specific)) + 
  geom_boxplot()+
  ylab("Proportion unshared repeats")+
  theme_minimal()+
  theme(axis.text=element_text(size=15),
         axis.title=element_text(size=15))
















































































combined4_subset_subset2<-subset(combined4_subset_subset, avg_shared>0)
combined4_subset_subset2<-subset(combined4_subset_subset2, TRprop>0.1)

scatterPlot <- ggplot(combined4_subset_subset2,aes(TRprop, avg_shared, color=centromere)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  geom_smooth(method=lm) +
  theme(legend.position=c(0,1), legend.justification=c(2,1))+
  ylab("Avg. Number of Species Sharing TR Sequences")

# Marginal density plot of x (top panel)
xdensity <- ggplot(combined4_subset_subset2, aes(TRprop, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined4_subset_subset2, aes(avg_shared, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+coord_flip()+
  xlab("Avg. Number of Species Sharing TR Sequences")


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))



scatterPlot <- ggplot(combined4_subset_subset2,aes(perc_GC, avg_shared, color=centromere)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  geom_smooth(method=lm) +
  theme(legend.position=c(0,1), legend.justification=c(2,1))+
  ylab("Avg. Number of Species Sharing TR Sequences")

# Marginal density plot of x (top panel)
xdensity <- ggplot(combined4_subset_subset2, aes(perc_GC, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined4_subset_subset2, aes(avg_shared, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+coord_flip()+
  xlab("Avg. Number of Species Sharing TR Sequences")


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))





################################################
## Distribution shared repeats per chromosome ##
################################################

chms1_sub<-subset(chms1, scaff_number=="12")#Shared vs Tps_specific #change chromosome number
chms9_sub<-subset(chms9, scaff_number=="12")#CenH3 #change chromosome number
merge3_sub<-subset(merge3TRprop, scaff_number=="12")#GW TR proportion #change chromosome number


cenh3_testes<-ggplot(chms9_sub, aes(x=start, y=log2cenH3_norm))+
  geom_point(size=0.3) +
  scale_y_continuous(limits=c(-1,3)) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept= 0, linetype="solid", color = "black") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


TR<-ggplot(merge3_sub, aes(x=wind_start, y=TRprop))+
  geom_point(size=0.3) +
  ylab("Proportion of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


prop_s_ns<-ggplot() +
  # Add non-shared motifs
  geom_point(data = chms1_sub, 
             aes(x = window_start, y = proportion_non_shared, color = "Non-shared"), 
             size = 0.2) +
  geom_smooth(data = chms1_sub, 
              aes(x = window_start, y = proportion_non_shared, group = scaff_number, color = "Non-shared"), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 20), size=0.5) +
  
  # Add shared motifs
  geom_point(data = chms1_sub, 
             aes(x = window_start, y = proportion_shared, color = "Shared"), 
             size = 0.2) +
  geom_smooth(data = chms1_sub, 
              aes(x = window_start, y = proportion_shared, group = scaff_number, color = "Shared"), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 20), size=0.5) +
  
  # Customize labels and theme
  ylab("Proportion of motifs") +
  xlab("Genomic coordinates") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    legend.title = element_blank()
  ) +
  # Set colors manually
  scale_color_manual(
    name = "Motif Type",
    values = c("Non-shared" = "red", "Shared" = "blue"), labels=NULL
  )



count_ns<-ggplot(chms1_sub, aes(x=window_start, y=non_shared_repeats))+
  geom_point(size=0.3) +
  ylab("number of Tps-specific repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    legend.position = "none"
  )+ ylim(0,50)+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


count_s<-ggplot(chms1_sub, aes(x=window_start, y=shared_repeats))+
  geom_point(size=0.3) +
  ylab("number of shared repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    legend.position = "none"
  )+ ylim(0,50)+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))



#prop_ns<-ggplot(chms1_sub, aes(x=window_start, y=proportion_non_shared))+
#  geom_point(size=0.5) +
#  ylab("Proportion of Tps-specific repeats")+
#  xlab("Genomic coordinates")+
#  theme(
#    panel.background = element_blank(),
#    axis.line = element_line(colour = "black"),
#    legend.position = "none"
#  )+
#  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


#prop_s<-ggplot(chms1_sub, aes(x=window_start, y=proportion_shared))+
#  geom_point(size=0.5) +
#  ylab("Proportion of shared repeats")+
#  xlab("Genomic coordinates")+
#  theme(
#    panel.background = element_blank(),
#    axis.line = element_line(colour = "black"),
#    legend.position = "none"
#  )+
#  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))


#Plots
library(cowplot)
plot_grid(
  plot_grid(
    cenh3_testes
    , TR + theme(legend.position = "none")
    , count_s + theme(legend.position = "none")
    , count_ns + theme(legend.position = "none")
    , prop_s_ns + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(TR)
    , ggplot()
    , get_legend(count_s)
    , ggplot()
    , get_legend(count_ns)
    , ggplot()
    , get_legend(prop_s_ns)
    , ncol =1)
  , rel_widths = c(12,-0.35)
)




plot_grid(
  plot_grid(
    cenh3_testes
    , TR + theme(legend.position = "none")
    , prop_s + theme(legend.position = "none")
    , prop_ns + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(TR)
    , ggplot()
    , get_legend(prop_s)
    , ggplot()
    , get_legend(prop_ns)
    , ncol =1)
  , rel_widths = c(12,-0.35)
)




### All chromosomes comparisons
##############################

#CenH3
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


#TR proportion
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


#Shared vs Tps-specific

prop_ns<-ggplot(chms1, aes(x=cumul_start, y=proportion_non_shared, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=cumul_start, y=proportion_non_shared), stat="identity", colour="maroon", size=0.5) +
  ylab("Proportion of Tps_specific repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))

prop_s<-ggplot(chms1, aes(x=cumul_start, y=proportion_shared, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=cumul_start, y=proportion_shared), stat="identity", colour="maroon", size=0.5) +
  ylab("Proportion of shared repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))

count_ns<-ggplot(chms1, aes(x=cumul_start, y=non_shared_repeats, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=cumul_start, y=non_shared_repeats), stat="identity", colour="maroon", size=0.5) +
  ylab("number of Tps_specific repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))

count_s<-ggplot(chms1, aes(x=cumul_start, y=shared_repeats, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=cumul_start, y=shared_repeats), stat="identity", colour="maroon", size=0.5) +
  ylab("number of shared repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+ ylim(0,50)+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))













































##################################################################
### Plot proportion of Tps specific motifs between chromosomes ###
##################################################################

#Calculate the number and proportion of shared repeats per chromosome
window_summary2 <- chms_ord_condensed2 %>%
  group_by(chr) %>%
  summarise(
    total_repeats = n(),
    non_shared_repeats = sum(shared == "1"),
    shared_repeats = total_repeats-non_shared_repeats,
    proportion_shared = shared_repeats / total_repeats,
    proportion_non_shared = non_shared_repeats / total_repeats
  )

window_summary2 <- chms_ord_condensed2 %>%
  group_by(chr) %>%
  summarise(
    total_repeats = n(),
    shared_repeats = sum(shared > 3),
    non_shared_repeats = sum(shared == "1"),
    proportion_shared = shared_repeats / total_repeats,
    proportion_non_shared = non_shared_repeats / total_repeats
  )


window_summary2$chr <- factor(window_summary2$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3",
                                                             "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6",
                                                             "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9",
                                                             "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))

ggplot(data=window_summary2, aes(x=chr, y=proportion_non_shared)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(proportion_non_shared, digits = 3)), vjust=1.6, size=3)+ ylab("proportion of unique tandem repeats")+
  theme_classic()

ggplot(data=window_summary2, aes(x=chr, y=shared_repeats)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(shared_repeats, digits = 3)), vjust=1.6, size=3)+ ylab("Number of shared repeats")+
  theme_classic()

ggplot(data=window_summary2, aes(x=chr, y=non_shared_repeats)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(non_shared_repeats, digits = 3)), vjust=1.6, size=3)+ ylab("Number of Tps-specific repeats")+
  theme_classic()




################################################
### Motifs shared between chromosomes of Tps ###
################################################
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")

Tps<-read.table("Tps_cleaned.txt", header=T, sep="\t", quote="")

#merge with Tps dataframe of all motifs
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
Tps_cleanedRegs_main <- readRDS("Tps_cleaned_withMainRepresentant.rds")

Tps_cleanedRegs_main$chr_start<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$start, sep = "_")
Tps_cleanedRegs_main$chr_end<-paste(Tps_cleanedRegs_main$chr, Tps_cleanedRegs_main$end, sep = "_")

df_filtered <- Tps_cleanedRegs_main %>%
  group_by(chr_start) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering


df_filtered2 <- df_filtered %>%
  group_by(chr_end) %>%              # Group by the start position
  filter(motif_length == min(motif_length)) %>%  # Keep only the rows with the shortest motif
  ungroup()                                 # Ungroup after filtering


chms<-subset(df_filtered2, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
             | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
             | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")


chms_ord<-chms[order(chms$chr, chms$start),]


chms_ord <- chms_ord %>%
  group_by(motif_representant) %>%
  mutate(shared_nb_chr = n_distinct(chr)) %>%
  ungroup()
#shared repeats among different numbers of chromosomes


chms_ord <- chms_ord %>%
  group_by(motif_representant) %>%
  mutate(shared = ifelse(n() > 1, 1, 0)) %>%
  ungroup()
#shared within and between repeats

chms_ord <- chms_ord %>%
  group_by(motif_representant, chr) %>%
  mutate(shared_within_chr = ifelse(n() > 1, 1, 0)) %>%
  ungroup() # Remove grouping after operation
#shared repeats within the same chromosome

chms_ord <- chms_ord %>%
  group_by(motif_representant) %>%
  mutate(shared_between_chr = ifelse(n_distinct(chr) > 1, 1, 0)) %>%
  ungroup() # Remove grouping after operation
#shared repeats between chromosomes


## Number of motifs shared among chromosomes vs genome representation 

chms_ord$representation<-as.numeric(chms_ord$nMotifsRepresented)
chms_ord$count<-1
chms_ord$shared_nb_chr<-as.character(chms_ord$shared_nb_chr)
chms_ord$shared_nb_chr <- factor(chms_ord$shared_nb_chr,levels = c("1", "2", "3",
                                                             "4", "5", "6",
                                                             "7", "8", "9",
                                                             "10", "11", "12"))


#Number of motifs shared among chromosomes
ggplot(data=chms_ord, aes(x=shared_nb_chr, y=count))+ 
  geom_bar(stat="identity")+
  xlab("number of chromsomes")+
  ggtitle("shared tandem repeat sequences")+
  theme_minimal()


ggplot(chms_ord, aes(x=shared_nb_chr, y=representation)) +
  geom_point(size=2, shape=23)+ ylab("number of representation in the genome")+ 
  xlab("number of chromsomes")+
  ggtitle("shared tandem repeat sequences")+
  theme_minimal()


chms_ord_sub<-subset(chms_ord, shared_nb_chr=="12")

chms_ord_sub2 <- chms_ord %>%
  group_by(shared_nb_chr, motif_representant, chr) %>%
  summarise(nb_occurence = length(chr))

chms_ord_sub2$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms_ord_sub2$chr)))

chms_ord_sub3<-subset(chms_ord_sub2, shared_nb_chr=="12" | shared_nb_chr=="1" | shared_nb_chr=="11")

ggplot(chms_ord_sub2, aes(x=as.factor(scaff_number), y=log10(nb_occurence))) +
  geom_boxplot()+
  facet_wrap(.~shared_nb_chr)+
  xlab("number of chromsomes")+
  ylab("number of representations per chromosome")+ 
  ggtitle("shared tandem repeat sequences")+
  theme_minimal()

ggplot(chms_ord_sub3, aes(x=as.factor(scaff_number), y=log10(nb_occurence))) +
  geom_boxplot()+
  facet_wrap(.~shared_nb_chr)+
  xlab("number of chromsomes")+
  ylab("number of representation per chromosome (log10)")+ 
  ggtitle("shared tandem repeat sequences", subtitle = "1 point = 1 sequence")+
  theme_minimal()


ggplot(chms_ord_sub3, aes(x=shared_nb_chr, y=log10(size))) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  xlab("number of chromosomes")+
  ylab("sequence length (log)")+ 
  ggtitle("shared tandem repeat sequences")+
  theme_minimal()

geom_point(size=2, shape=23)

hist(chms_ord_sub$shared_nb_chr, breaks=12)


ggplot(chms_ord_sub, aes(x=shared_nb_chr, y=representation)) +
  geom_point(size=2, shape=23)+ ylab("number of representation in the genome")+ 
  xlab("number of chromosomes")+
  ggtitle("shared motifs")+
  theme()


window_summary <- chms_ord %>%
  group_by(motif_representant) %>%
  summarise(
    total_repeats = n(),
    shared_repeats_both = sum(shared == "1"),
    shared_repeats_within = sum(shared_within_chr == "1"),
    shared_repeats_between = sum(shared_between_chr == "1"),
    proportion_shared_repeats_both = shared_repeats_both / total_repeats,
    proportion_shared_repeats_within = shared_repeats_within / total_repeats,
    proportion_shared_repeats_between = shared_repeats_between / total_repeats
  )


chms_ord$shared_nb_chr<-as.numeric(chms_ord$shared_nb_chr)

window_summary <- chms_ord %>%
  group_by(motif_representant, shared_nb_chr) %>%
  summarise(
    counts = count(chms_ord, shared_nb_chr),
  )

window_summary2<-count(window_summary, shared_nb_chr)

ggplot(data=window_summary2, aes(x=shared_nb_chr, y=n))+ 
  geom_bar(stat="identity")+
  theme_minimal()




#Assign to each repeat array a non-overlapping 100kb window
chms_ord$window_start <- (floor(chms_ord$start / 100000) * 100000)+1
chms_ord$window_end <- chms_ord$window_start + 99999  # End of the 100kb window
chms_ord$chm_pos<-paste(chms_ord$chr, chms_ord$window_start, sep="-")

#Calculate the number and proportion of repeats shared among Tps chromosomes per 100kb windows
window_summary <- chms_ord %>%
  group_by(chr, window_start) %>%
  summarise(
    total_repeats = n(),
    shared_repeats_both = sum(shared == "1"),
    shared_repeats_within = sum(shared_within_chr == "1"),
    shared_repeats_between = sum(shared_between_chr == "1"),
    proportion_shared_repeats_both = shared_repeats_both / total_repeats,
    proportion_shared_repeats_within = shared_repeats_within / total_repeats,
    proportion_shared_repeats_between = shared_repeats_between / total_repeats
  )


summary(window_summary)

windows$chm_pos<-paste(windows$chr, windows$window_start, sep="-")
window_summary$chm_pos<-paste(window_summary$chr, window_summary$window_start, sep="-")


#Combine non-overlapping 100kb windows assigned per array with non-overlapping 100kb windows chromosome-wide (i.e., include windows without TR arrays)
merge<-merge(windows, window_summary[3:10], by="chm_pos", all.x=T)
summary(merge)
merge[is.na(merge)] <- 0 #NA are retrieved when windows are found without TR annotation 

merge$scaff_number<-as.numeric(as.character(gsub("^.*scf","", merge$chr)))
chms1<-merge[order(merge$scaff_number, merge$window_start),]
chms1$cumul_start<-seq(1, by = 100000, length.out = nrow(chms1))
chms1$chrom<-ifelse(chms1$scaff_number=="3", "sex-chrom", "autosome")
chms2<-chms1[chms1$scaff_number=="3",]


ggplot(chms1, aes(x=cumul_start, y=proportion_shared_repeats_between, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=cumul_start, y=proportion_shared_repeats_between), stat="identity", colour="maroon", size=0.5) +
  ylab("Proportion of tandem repeats NOT shared")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "gray61", method="gam", formula = y~s(x, bs="cs", k=20))




window_summary2 <- chms_ord_condensed2 %>%
  group_by(chr) %>%
  summarise(
    total_repeats = n(),
    shared_repeats_both = sum(shared_all_chr == "1"),
    shared_repeats_within = sum(shared_within_chr == "1"),
    shared_repeats_between = sum(shared_between_chr == "1"),
    proportion_shared_repeats_both = shared_repeats_both / total_repeats,
    proportion_shared_repeats_within = shared_repeats_within / total_repeats,
    proportion_shared_repeats_between = shared_repeats_between / total_repeats
  )


window_summary2$chr <- factor(window_summary2$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3",
                                                             "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6",
                                                             "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9",
                                                             "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))

ggplot(data=window_summary2, aes(x=chr, y=proportion_shared_repeats_between)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(proportion_shared_repeats_between, digits = 3)), vjust=1.6, size=3)+ ylab("proportion of tandem repeats NOT shared")+
  theme_classic()






## Proportion of Repeat motifs shared between pairs of chromosomes
##################################################################

# Step 1: Filter for shared motifs
shared_data <- chms_ord_condensed2 %>%
  filter(shared_between_chr == 1)

# Get all unique chromosomes
library(stringr)
chromosomes <- unique(chms_ord_condensed2$chr)
chromosome_order <- chromosomes[order(as.numeric(str_extract(chromosomes, "\\d+$")))]

# Create all possible chromosome pairs (including same chromosomes)
library(tidyr)
all_pairs <- expand.grid(chromosome1 = chromosome_order, chromosome2 = chromosome_order)

# Map motifs to chromosomes
motif_chromosomes <- shared_data %>%
  group_by(motif_representant) %>%
  summarize(chromosomes = list(unique(chr))) %>%
  unnest(chromosomes)

# Create all pairs of chromosomes sharing a motif
chromosome_pairs <- motif_chromosomes %>%
  rename(chromosome1 = chromosomes) %>%
  inner_join(motif_chromosomes, by = "motif_representant", relationship = "many-to-many") %>%
  filter(chromosome1 < chromosomes) %>%
  rename(chromosome2 = chromosomes)

# Calculate the number of unique motifs for each chromosome
chromosome_totals <- chms_ord_condensed2 %>%
  group_by(chr) %>%
  summarize(total_unique_by_chromosome = length(unique(motif_representant))) %>%
  rename(chromosome = chr)

# Calculate the proportion of shared motifs between pairs
pair_proportions <- chromosome_pairs %>%
  count(chromosome1, chromosome2) %>%
  left_join(chromosome_totals, by = c("chromosome1" = "chromosome")) %>%
  rename(total1 = total_unique_by_chromosome) %>%
  left_join(chromosome_totals, by = c("chromosome2" = "chromosome")) %>%
  rename(total2 = total_unique_by_chromosome) %>%
  mutate(proportion = n / (total1 + total2))


## Count shared motifs for each pair
#pair_counts <- chromosome_pairs %>%
#  count(chromosome1, chromosome2)

## Create a matrix of shared motifs
#chromosome_matrix <- reshape2::acast(
#  pair_counts, 
#  chromosome1 ~ chromosome2, 
#  value.var = "n", 
#  fill = 0
#)


# Merge counts with the full set of pairs, default to 0 for missing pairs
all_pairs <- left_join(all_pairs, pair_proportions, by = c("chromosome1", "chromosome2"))

## Create a matrix of shared motif proportions
chromosome_matrix <- reshape2::acast(
  all_pairs, 
  chromosome1 ~ chromosome2, 
  value.var = "proportion", 
  fill = 0
)

## Make the matrix symmetrical
sym_matrix <- (chromosome_matrix + t(chromosome_matrix))


# Visualize with a heatmap
library(ggplot2)
library(reshape2)
heatmap_data <- melt(sym_matrix, na.rm = TRUE)

heatmap_data$Var1 <- factor(heatmap_data$Var1, levels = chromosome_order)
heatmap_data$Var2 <- factor(heatmap_data$Var2, levels = chromosome_order)


ggplot(heatmap_data, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "Proportion Shared Motifs") +
  theme_minimal() +
  labs(x = "Chromosome 1", y = "Chromosome 2", title = "Proportion Shared Motif Sequences Between Chromosomes")





##################################################################################
## Similar Motifs between chromosomes (allowing mutations based on Levenshtein) ##
##################################################################################

##List of unique set of motifs in Tps
Tps$Seq <- paste0(">seq_", seq_len(nrow(Tps)))


##Pairs of unique motifs with at least 80% similarities
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/Levenstein_Marion")
table<-read.table("above80_Oct2024.txt", header=T, sep="\t", quote="")
colnames(table)<-c("block", "Seq1", "Seq2", "Levenstein_distance", "perc_divergence", "motif_length_seq1", "motif_length_seq2")

table_sub<-subset(table, perc_divergence<=10) #subset pairs with at least 90% similarities
table_sub2<-subset(table_sub, motif_length_seq1>2 & motif_length_seq2>2)
table_sub3<-table_sub2[c(2,3,5)]

# Flatten table_sub3 to a single vector of unique DNA motifs that have at least another motif with similar sequence
table_sub3_flat<-unique(c(table_sub3$Seq1, table_sub3$Seq2))

##Combine List of unique set of motifs with Pair of similar motifs
merge<-merge(table_sub3, Tps[c(7,10)], by.x="Seq1", by.y="Seq", all.x=T)
merge2<-merge(merge, Tps[c(7,10)], by.x="Seq2", by.y="Seq", all.x=T)
colnames(merge2)<-c("Seq1", "Seq2", "perc_divergence", "motif_representant_1", "motif_representant_2")
merge2$levenstein<-"similar_to"

# Flatten merge2 to a single vector of unique DNA motifs that have at least another motif with similar sequence
merge2_flat <- unique(c(merge2$motif_representant_1, merge2$motif_representant_2))

venn <- list(Tps_unique=unique, related_motifs=table_sub3_flat)
ggVennDiagram(venn) #2845/72286=0.03935755

# Add a new column to chms_ord_condensed2 based on whether the sequence exists in merge2_flat
chms_ord <- chms_ord %>%
  mutate(similarity_status = ifelse(motif_representant %in% merge2_flat, "similar_to", "no_similarity"))


#Calculate the number and proportion of repeats shared among Tps chromosomes per 100kb windows
window_summary <- chms_ord %>%
  group_by(chr, window_start) %>%
  summarise(
    total_repeats = n(),
    similar_repeats = sum(similarity_status == "similar_to"),
    not_similar_repeats = sum(similarity_status == "no_similarity"),
    proportion_similar_repeats = similar_repeats / total_repeats,
    proportion_not_similar_repeats = not_similar_repeats / total_repeats,
  )





##Reduce chms_ord_condensed2

similarity_table<-chms_ord_condensed2[c(3,4,17,18,19,20,21,22)]

# Remove duplicate rows
similarity_table_unique <- similarity_table %>% distinct()

#Count the number of DNA motifs with similar sequences based on their chromosomal representation
result <- similarity_table_unique %>%
  group_by(nb_chromosomes = shared_chr, similarity = similarity_status) %>%
  summarise(count = n(), .groups = "drop")



















### Test for the Library Hypothesis in Tps
##########################################
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")

Tps<-read.table("Tps_cleaned.txt", header=T, sep="\t", quote="")
Tdi<-read.table("Tdi_cleaned.txt", header=T, sep="\t", quote="")
Tcm<-read.table("Tcm_cleaned.txt", header=T, sep="\t", quote="")
Tce<-read.table("Tce_cleaned.txt", header=T, sep="\t", quote="")


# Combine all motifs from the 6 dataframes and get unique motifs
all_motifs <- unique(c(Tps$aomfi, Tdi$aomfi, Tcm$aomfi, Tce$aomfi))

# Sort the motifs
all_motifs <- sort(all_motifs)

# Create an empty matrix for presence-absence (motifs as rows, species as columns)
presence_absence_matrix <- matrix(0, nrow = length(all_motifs), ncol = 4)
#presence_absence_matrix <- matrix(0, nrow = length(all_motifs), ncol = 7)

# Assign row names (motifs) and column names (species)
rownames(presence_absence_matrix) <- all_motifs
colnames(presence_absence_matrix) <- c("Tps", "Tdi", "Tcm", "Tce")
#colnames(presence_absence_matrix) <- c("Tps", "Tdi", "Tcm", "Tsi", "Tce", "Tpa", "Tbi")

# Fill the matrix with 1s for presence of motifs in each species
presence_absence_matrix[all_motifs %in% Tps$aomfi, 1] <- 1
presence_absence_matrix[all_motifs %in% Tdi$aomfi, 2] <- 1
presence_absence_matrix[all_motifs %in% Tcm$aomfi, 3] <- 1
presence_absence_matrix[all_motifs %in% Tce$aomfi, 4] <- 1

head(presence_absence_matrix)


#convert to dataframe and sum rows
presence_absence_df<-data.frame(presence_absence_matrix)
presence_absence_df$aomfi<-rownames(presence_absence_df)
#presence_absence_df$shared <- rowSums( presence_absence_df[,c(1,4)]) #depending on the species comparisons (e.g., for Tps vs Tdi --> presence_absence_df[,1:2])
#presence_absence_df$sum <- rowSums( presence_absence_df[,1:7] )

#select only motifs present in specific species
presence_absence_df_Tps<-subset(presence_absence_df, Tps=="1")
presence_absence_df_Tps$shared<-rowSums( presence_absence_df_Tps[,c(1:2)])

Tps_unique<-subset(presence_absence_df_Tps, shared=="1")
Tps_unique2<-subset(Tps_unique, Tcm=="1") #subset repeats that are shared between Tps and Tcm but not between Tps and Tdi
Tps_unique3<-subset(Tps_unique, Tce=="1")
#merge with Tps dataframe of all motifs
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
Tps_cleanedRegs_main <- readRDS("Tps_cleaned_withMainRepresentant.rds")

merge_Tps_unique<-merge(Tps_unique, Tps_cleanedRegs_main, by="aomfi", all.x = T)




Tdi_full_dataset<-read.table("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt", header=T, sep="\t", quote="")

venn <- list(Tps_unique=Tps_unique$aomfi, Tdi_full_dataset=Tdi_full_dataset$aomfi)
ggVennDiagram(venn) #2845/72286=0.03935755

venn <- list(Tps_unique=Tps_unique2$aomfi, Tdi_full_dataset=Tdi_full_dataset$aomfi)
ggVennDiagram(venn) #533/4282=0.1244745 


venn <- list(Tps=Tps_full_dataset$aomfi, Tdi=Tdi_full_dataset$aomfi)
ggVennDiagram(venn) #2845/72286=0.03935755


shared<-subset(merge_Tps_unique, (aomfi %in% Tdi_full_dataset$aomfi))





#Levenstein between chromosomes

LevensteinTps<-read.table("above80_Oct2024.txt", header=T, sep="\t", quote="")


















# Create 10kb windows for each chromosome
create_windows <- function(chrom, max_position, window_size = 10000) {
  start <- seq(1, max_position, by = window_size)
  end <- pmin(start + window_size - 1, max_position)
  return(data.frame(chr = chrom, window_start = start, window_end = end))
}


# Find the maximum chromosome end position (this depends on your actual data)
chrom_max_positions <- chms %>%
  group_by(chr) %>%
  summarise(max_position = max(end))

# Create windows for all chromosomes
windows <- chrom_max_positions %>%
  rowwise() %>%
  do(create_windows(.$chr, .$max_position))

# Classify sequences into shared (num_species_shared > 1) or non-shared (num_species_shared == 1)
chms <- chms %>%
  mutate(shared = ifelse(sum > 1, TRUE, FALSE))

# Create a GRanges object for sequences
gr_sequences <- GRanges(seqnames = chms$chr,
                        ranges = IRanges(start = chms$start, end = chms$end),
                        shared = chms$shared)

# Create a GRanges object for windows
gr_windows <- GRanges(seqnames = windows$chr,
                      ranges = IRanges(start = windows$window_start, end = windows$window_end))

# Find overlaps between sequences and windows
overlaps <- findOverlaps(gr_sequences, gr_windows)


# Create a data frame to store the overlap information
overlap_df <- data.frame(
  sequence_idx = queryHits(overlaps),
  window_idx = subjectHits(overlaps)
)

chms <- chms %>%
  mutate(sequence_idx = row_number())


# Join the overlap data with the sequences to get shared/non-shared info
overlap_df <- overlap_df %>%
  left_join(chms, by = "sequence_idx")

# Summarize the proportion of shared vs non-shared sequences for each window
window_summary <- overlap_df %>%
  left_join(chms[overlap_df$sequence_idx, ], by = "sequence_idx") %>%
  group_by(window_idx) %>%
  summarise(
    total_sequences = n(),
    shared_sequences = sum(shared),
    non_shared_sequences = sum(!shared),
    proportion_shared = shared_sequences / total_sequences,
    proportion_non_shared = non_shared_sequences / total_sequences
  )


# Step 7: Add a unique identifier (row number) to the windows dataframe for joining purposes
windows <- windows %>%
  mutate(window_idx = row_number())

# Step 8: Combine the window summary with the window information using window_idx for joining
window_summary_final <- window_summary %>%
  left_join(windows, by = "window_idx")









###recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs.bed", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")



#only keep assembled scaffolds
chms<-subset(mergeTps, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
             | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
             | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")


chms <- chms %>%
  mutate(shared = ifelse(sum > 1, "shared", "non_shared"))


setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
chromosome_ranges<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")
colnames(chromosome_ranges)<-c("chr", "max_position")

#chromosome_ranges <- chms %>%
#  group_by(chr) %>%
#  summarise(max_position = max(end))


windows <- chromosome_ranges %>%
  rowwise() %>%
  do({
    chrom <- .$chr
    max_pos <- .$max_position
    start_pos <- seq(1, max_pos, by = 100000)
    end_pos <- pmin(start_pos + 99999, max_pos)
    data.frame(chr = chrom, window_start = start_pos, window_end = end_pos)
  })


gr_repeats <- GRanges(
  seqnames = chms$chr,
  ranges = IRanges(start = chms$start, end = chms$end)
)


gr_repeats <- GRanges(
  seqnames = chms$chr,
  ranges = IRanges(start = chms$start, end = chms$end),
  shared = chms$shared
)


gr_windows <- GRanges(
  seqnames = windows$chr,
  ranges = IRanges(start = windows$window_start, end = windows$window_end)
)

overlaps <- findOverlaps(gr_repeats, gr_windows)

overlap_df <- data.frame(
  repeat_idx = queryHits(overlaps),
  window_idx = subjectHits(overlaps)
)


overlap_df <- overlap_df %>%
  mutate(window_start = windows$window_start[window_idx],
         window_end = windows$window_end[window_idx])


chms <- chms %>%
  mutate(repeat_idx = row_number()) %>%
  left_join(overlap_df, by = "repeat_idx")



window_summary <- chms %>%
  group_by(chr, window_start, window_end) %>%
  summarise(
    total_repeats = n(),
    median_motif_length = median(motif_length),
  )


window_summary <- chms %>%
  group_by(chr, window_start, window_end) %>%
  summarise(
    total_repeats = n(),
    shared_repeats = sum(shared == "shared"),
    non_shared_repeats = sum(shared == "non_shared"),
    proportion_shared = shared_repeats / total_repeats,
    proportion_non_shared = non_shared_repeats / total_repeats
  )



window_summary$scaff_number<-as.numeric(as.character(gsub("^.*scf","", window_summary$chr)))
chms1<-window_summary[order(window_summary$scaff_number, window_summary$window_start),]
chms1$cumul_start<-seq(0, by = 100000, length.out = nrow(chms1))
chms1$chrom<-ifelse(chms1$scaff_number=="3", "sex-chrom", "autosome")
chms2<-chms1[chms1$scaff_number=="3",]

ggplot(chms1, aes(x=cumul_start, y=shared_repeats, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=chms2$cumul_start, y=chms2$shared_repeats), stat="identity", colour="maroon", size=0.5) +
  ylab("Number of shared tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")


ggplot(chms1, aes(x=cumul_start, y=proportion_shared, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=chms2$cumul_start, y=chms2$proportion_shared), stat="identity", colour="maroon", size=0.5) +
  ylab("Proportion of shared tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")

#filter out motifs that are present in all species
universal_motifs <- rownames(presence_absence_matrix)[rowSums(presence_absence_matrix) == ncol(presence_absence_matrix)]
filtered_motif_matrix <- presence_absence_matrix[!rownames(presence_absence_matrix) %in% universal_motifs, ]




# Perform ancestral state reconstruction using ace (from the ape package)
ancestral_states <- list()

# Reconstruct each motif across the phylogeny using the ace function
for (i in 1:nrow(filtered_motif_matrix)) {
  motif <- filtered_motif_matrix[i, ]
  
  # Reconstruct using ace for binary characters
  ancestral_states[[rownames(filtered_motif_matrix)[i]]] <- ace(motif, tree, type = "discrete")
}




# Loop through each motif and reconstruct its ancestral state
for (i in 1:nrow(filtered_motif_matrix)) {
  # Extract the binary presence/absence data for one motif
  motif <- filtered_motif_matrix[i, ]
  
  # Convert motif data to a named vector (names correspond to species)
  motif_data <- setNames(as.numeric(motif), colnames(filtered_motif_matrix))
  
  # Ensure that the motif data is a vector with correct species names
  print(motif_data)
  
  # Reconstruct using ace for binary characters
  ancestral_states[[rownames(filtered_motif_matrix)[i]]] <- ace(motif, tree, type = "discrete")
}






root_ancestral_states <- sapply(ancestral_states, function(x) x$lik.anc[1, ])
#this is what is expected from chatGPT
#root_ancestral_states
#       ATGC       GGTA       TATA       etc
#        1          0          1          0

ancestral_motifs_per_species <- colSums(test & (root_ancestral_states == 1)) / colSums(motif_binary_matrix)
ancestral_motifs_per_species <- colSums(test & (root_ancestral_states > 0.5)) / colSums(motif_binary_matrix)

ancestral_motifs_per_species <- colSums(motif_binary_matrix & (root_ancestral_states == 1)) / colSums(motif_binary_matrix)



# Plot tree with ancestral states for one motif (example)
plot(tree, show.tip.label = TRUE)
nodelabels(pie = ancestral_states[["ATGC"]]$lik.anc, piecol = c("white", "black"), cex = 0.6)








setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")
rho<-read.table("RHO_Tpop_mean_5_runs.bed", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho")












