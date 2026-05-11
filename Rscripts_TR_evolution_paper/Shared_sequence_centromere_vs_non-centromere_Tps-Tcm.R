##################################################
### Generate a presence (1) absence (0) matrix ###
##################################################
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")

Tps<-read.table("Tps_cleaned.txt", header=T, sep="\t", quote="")
Tcm<-read.table("Tcm_cleaned.txt", header=T, sep="\t", quote="")

# Combine all motifs from the 6 dataframes and get unique motifs
all_motifs <- unique(c(Tps$aomfi, Tcm$aomfi))

# Sort the motifs
all_motifs <- sort(all_motifs)

# Create an empty matrix for presence-absence (motifs as rows, species as columns)
presence_absence_matrix <- matrix(0, nrow = length(all_motifs), ncol = 2)

# Assign row names (motifs) and column names (species)
rownames(presence_absence_matrix) <- all_motifs
colnames(presence_absence_matrix) <- c("Tps", "Tcm")

# Fill the matrix with 1s for presence of motifs in each species
presence_absence_matrix[all_motifs %in% Tps$aomfi, 1] <- 1
presence_absence_matrix[all_motifs %in% Tcm$aomfi, 2] <- 1

head(presence_absence_matrix)


#convert to dataframe 
presence_absence_df<-data.frame(presence_absence_matrix)
presence_absence_df$aomfi<-rownames(presence_absence_df)

presence_absence_df$shared<-rowSums( presence_absence_df[,c(1,2)])#shared between T.poppense and T.californicum

#merge with T.poppense TR annotation
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


mergeTps<-merge(df_filtered2_tps, presence_absence_df, by.x="motif_representant", by.y="aomfi", all.x =T)
summary(mergeTps)


#Assign to each repeat array a non-overlapping 250kb window
mergeTps$window_start <- (floor(mergeTps$start / 250000) * 250000)+1
mergeTps$window_end <- mergeTps$window_start + 249999  # End of the 250kb window
mergeTps$chm_pos<-paste(mergeTps$chr, mergeTps$window_start, sep="-")



#merge with T.californicum TR annotation
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


mergeTcm<-merge(df_filtered2_tcm, presence_absence_df, by.x="motif_representant", by.y="aomfi", all.x =T)
summary(mergeTcm)


#Assign to each repeat array a non-overlapping 250kb window
mergeTcm$window_start <- (floor(mergeTcm$start / 250000) * 250000)+1
mergeTcm$window_end <- mergeTcm$window_start + 249999  # End of the 250kb window
mergeTcm$chm_pos<-paste(mergeTcm$chr, mergeTcm$window_start, sep="-")



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



#############################
### ChIP cenh3 testes Tps ###
#############################

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
data1$log2cenH3_norm[!is.finite(data1$log2cenH3_norm)] <- 0
data1$region_size<-abs(data1$stop-data1$start)

data1$centromere<-ifelse(data1$log2cenH3_norm>0.25, "centromere", "non-centromere") #~8% of the genome



#############################
### ChIP cenh3 testes Tcm ###
#############################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
cenh3_tcm<-read.table("Tcm_testes_M-Fbodies_cenh3_1_coverage_DR_250kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_coverage_DR_100kb.txt

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
input1_tcm<-read.table("Tcm_testes_M-Fbodies_input_coverage_DR_250kb.txt",
                   header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data2<-cbind(cenh3_tcm,input1_tcm[4])
colnames(data2)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data2$coverage_cenh3_norm<-data2$coverage_cenh3/53128790 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data2$coverage_input_norm<-data2$coverage_input/80611894 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data2$ratiocenH3_norm<-data2$coverage_cenh3_norm/data2$coverage_input_norm
data2$log2cenH3_norm<-log2(data2$ratiocenH3_norm)
data2$start2<-data2$start+1
data2$chm_pos<-paste(data2$Scaffold_name, data2$start2, sep="-")
data2$log2cenH3_norm[!is.finite(data2$log2cenH3_norm)] <- 0

data2$centromere<-ifelse(data2$log2cenH3_norm>0, "centromere", "non-centromere") #~8% of the genome 
#matching the threshold to Tps (i.e., log2cenH3_norm>0.25), makes only the "centromere" intercept to vary and leads to the same interpretation.

mergeTps2<-merge(mergeTps, mergeTRprop[c(6,8)], by="chm_pos", all.x=T)
mergeTps2<-merge(mergeTps2, data1[c(9,11,12,13)], by="chm_pos", all.x=T)
mergeTps2$chr2 <- sub("^.*_", "", mergeTps2$chr)

mergeTcm2<-merge(mergeTcm, data2[11:12], by="chm_pos", all.x=T)
mergeTcm2$chr <- sub("\\..*$", "", mergeTcm2$chr)
mergeTcm2$chr2 <- sub("^.*_", "", mergeTcm2$chr)


#Add matching column to Tps df based on sequence and centromere position matches
key_tbl <- mergeTcm2 %>%
  distinct(motif_representant, centromere, chr2)

mergeTps2 <- mergeTps2 %>%
  mutate(match = if_else(
    paste(motif_representant, centromere, chr2) %in%
      paste(key_tbl$motif_representant, key_tbl$centromere, key_tbl$chr2),
    1L, 0L
  ))


#key_tbl <- mergeTcm2 %>%
#  distinct(motif_representant, centromere)

#mergeTps2 <- mergeTps2 %>%
#  mutate(match = if_else(
#    paste(motif_representant, centromere) %in%
#      paste(key_tbl$motif_representant, key_tbl$centromere),
#    1L, 0L
#  ))

summary(mergeTps2)

# Main summarization
window_summary <- mergeTps2 %>%
  group_by(chr, window_start, window_end) %>%
  summarise(
    total_repeats = n(),
    avg_motif = mean(L_motif_representant),
    median_motif = median(L_motif_representant),
    median_CN = median(copy_nb),
    mean_CN = mean(copy_nb),
    motif_diversity = n_distinct(motif_representant),
    prop_shared = sum(match == "1")/n(),
    TRprop=unique(TRprop),
    region_size=unique(region_size),
    centromere=unique(centromere),
    log2cenH3_norm=unique(log2cenH3_norm),
    .groups = "drop"
  )

summary(window_summary) #NA values correspond to regions with non shared sequences (between chromosomes) and that are placed on unanchored scaffolds



##PLots
ggplot(window_summary,aes(TRprop, prop_shared, color=centromere)) + 
  geom_point() + 
  scale_color_manual(values = c('red','grey')) +
  geom_smooth(method=lm) +
  theme(legend.position=c(0,1), legend.justification=c(2,1))+
  ylab("shared")

#stats
window_summary<-subset(window_summary, chr!="Tps_mtDNA_v350")
window_summary<-subset(window_summary, region_size>=10000)
window_summary<-window_summary[!is.na(window_summary$TRprop),]

results <- McSweeny_Porter(window_summary, prop_shared~TRprop+centromere)
print(results$regression_equation_interaction)
print(results)
print(results$group_effect)
print(results$interaction_effect)

library(effectsize)
eta_squared(model, partial = TRUE)


##chromosomes only
window_summary2<-subset(window_summary, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                                | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                                | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

summary(window_summary2)

library(ggplot2)

scatterPlot <- ggplot() +
  geom_point(data = window_summary2 %>% filter(centromere == "non-centromere"),
             aes(TRprop, prop_shared),
             color = "grey", alpha = 0.4) +
  geom_point(data = window_summary2 %>% filter(centromere == "centromere"),
             aes(TRprop, prop_shared),
             color = "red", alpha = 0.6) +
  geom_smooth(data = window_summary2,
              aes(TRprop, prop_shared, color = centromere),
              method = "lm", se = TRUE) +
  scale_color_manual(values = c("centromere" = "red",
                                "non-centromere" = "grey")) +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(-5.5,-1)) +
  ylab("Proportion shared")

# Marginal density plot of x (top panel)
xdensity <- ggplot(window_summary2, aes(TRprop, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('red','black')) + 
  theme_bw() +
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(window_summary2, aes(prop_shared, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('red','black')) + 
  theme_bw() +
  theme(legend.position = "none")+coord_flip()+
  xlab("Proportion shared")


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



#stats
library(npANCOVA)
results <- McSweeny_Porter(window_summary2, prop_shared~TRprop+centromere)
print(results$regression_equation_interaction)
print(results)
print(results$group_effect)
print(results$interaction_effect)


model<-lm(window_summary2$prop_shared~window_summary2$TRprop*window_summary2$centromere)
summary(model)
anova(model)

library(effectsize)
eta_squared(model, partial = TRUE)





#########################################################
## Distribution TR sequence turnover along chromosomes ##
#########################################################

window_summary2$scaff_number<-as.numeric(as.character(gsub("^.*scf","", window_summary2$chr)))

window_summary2<-window_summary2[order(window_summary2$scaff_number, window_summary2$window_start),]
window_summary2$cumul_start<-seq(0, by = 250000, length.out = nrow(window_summary2))

# Calculate scaffold boundaries (start and end positions)
library(dplyr)
scaffold_bounds <- window_summary2 %>%
  group_by(scaff_number) %>%
  summarise(min_pos = min(cumul_start),
            max_pos = max(cumul_start)) %>%
  mutate(fill_col = rep(c("white", "grey90"), length.out = n()))



TR_turnover<-ggplot(window_summary2, aes(x=cumul_start, y=prop_shared, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size=0.5, alpha = 0.6) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey"))+
  ylab("Proportion shared")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



cenh3_testes<-ggplot(window_summary2, 
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


TR_array<-ggplot(window_summary2, aes(x=cumul_start, y=TRprop, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.4) +
  scale_fill_identity() +
  geom_point(size=0.5) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey")) +
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept = 0.125, linetype = "dashed", color = "black") +
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
cen    <- window_summary2 %>% filter(centromere == "centromere")
noncen <- window_summary2 %>% filter(centromere == "non-centromere")

# Precompute quantile breaks (same bins used every iteration)
breaks <- quantile(window_summary2$TRprop, probs = seq(0, 1, length.out = n_bins + 1),
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
lm_cen <- lm(prop_shared ~ TRprop, data = cen)
lm_cen_sum <- summary(lm_cen)
cen_stats <- list(
  slope = coef(lm_cen)["TRprop"],
  intercept = coef(lm_cen)["(Intercept)"],
  r2 = lm_cen_sum$r.squared,
  mean_seq_id = mean(cen$prop_shared, na.rm = TRUE),
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
    lm_non <- lm(prop_shared ~ TRprop, data = noncen_sub)
    s <- summary(lm_non)
    non_stats <- list(
      slope = coef(lm_non)["TRprop"],
      intercept = coef(lm_non)["(Intercept)"],
      r2 = s$r.squared,
      mean_prop = mean(noncen_sub$prop_shared, na.rm = TRUE),
      n = nrow(noncen_sub)
    )
  } else {
    non_stats <- list(slope = NA, intercept = NA, r2 = NA, mean_prop = NA, n = nrow(noncen_sub))
  }
  
  results[[i]] <- tibble(
    iter = i,
    cen_slope = cen_stats$slope,
    cen_intercept = cen_stats$intercept,
    cen_r2 = cen_stats$r2,
    cen_mean_prop = cen_stats$mean_prop,
    cen_n = cen_stats$n,
    non_slope = non_stats$slope,
    non_intercept = non_stats$intercept,
    non_r2 = non_stats$r2,
    non_mean_prop = non_stats$mean_prop,
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
  xlab("slope (Proportion shared ~ TRprop)")

p2 <- ggplot(res_df, aes(x = non_mean_prop)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = cen_stats$mean_seq_id, color = "red", linetype = "dashed") +
  ggtitle("Distribution of mean(log2cenH3_norm) for sampled non-centromere sets") +
  xlab("mean(Proportion shared)")

# Print plots
print(p1); print(p2)



########################################################
### Choose representative iteration from subsampling ###
########################################################

library(dplyr)

#set.seed(123)

# Split data
cen  <- window_summary2 %>% filter(centromere == "centromere")
noncen <- window_summary2 %>% filter(centromere == "non-centromere")

# Create TRprop bins
n_bins <- 12   # adjust for smoother matching
breaks <- quantile(window_summary2$TRprop, probs = seq(0, 1, length.out = n_bins + 1))

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
scatterPlot <- ggplot(matched, aes(TRprop, prop_shared, color = centromere)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c('red','grey')) +
  geom_smooth(method = lm) +
  theme_bw() +
  theme(legend.position = c(0,1), legend.justification = c(-5.5, -1)) +
  ylab("Proportion shared")

xdensity <- ggplot(matched, aes(TRprop, fill = centromere)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c('red','black')) +
  theme_bw() +
  theme(legend.position = "none")

ydensity <- ggplot(matched, aes(prop_shared, fill = centromere)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c('red','black')) +
  theme_bw() +
  theme(legend.position = "none") +
  coord_flip() +
  xlab("Proportion shared")

grid.arrange(
  xdensity, blankPlot,
  scatterPlot, ydensity,
  ncol = 2, nrow = 2,
  widths = c(4, 1.4),
  heights = c(1.4, 4)
)

