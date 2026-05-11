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
 
df_filtered2$window_start <- (floor(df_filtered2$start / 250000) * 250000)+1
df_filtered2$window_end <- df_filtered2$window_start + 249999  # End of the 250kb window
df_filtered2$chm_pos<-paste(df_filtered2$chr, df_filtered2$window_start, sep="-")


## T.poppense centromere
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
data1$centromere<-ifelse(data1$log2cenH3_norm>0.25, "centromere", "non-centromere")


library(dplyr)
df2_annotated <- merge(
       df_filtered2,
       data1[, c("chm_pos", "centromere")],
       by = "chm_pos",
       all.x = TRUE)

summary(df2_annotated)



## Plots
stats <- df2_annotated %>%
  group_by(centromere) %>%
  summarise(
    mean_size = mean(L_motif_representant, na.rm = TRUE),
    median_size = median(L_motif_representant, na.rm = TRUE)
  )

custom_colors <- c("non-centromere" = "#E69F00", "centromere" = "#999999")

ggplot(df2_annotated, aes(x = L_motif_representant, fill = centromere, color = centromere)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  geom_vline(
    data = stats,
    aes(xintercept = mean_size, color = centromere),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    data = stats,
    aes(xintercept = median_size, color = centromere),
    linetype = "solid",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    subtitle = "Dashed line = Mean | Solid line = Median",
    x = "TR sequence length (bp)",
    y = "Density",
    fill = "Region type",
    color = "Region type"
  ) +
  xlim(0,500)+
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )




stats2 <- df2_annotated %>%
  group_by(centromere) %>%
  summarise(
    mean_size = mean(copy_nb, na.rm = TRUE),
    median_size = median(copy_nb, na.rm = TRUE)
  )

ggplot(df2_annotated, aes(x = copy_nb, fill = centromere, color = centromere)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  geom_vline(
    data = stats2,
    aes(xintercept = mean_size, color = centromere),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    data = stats2,
    aes(xintercept = median_size, color = centromere),
    linetype = "solid",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    subtitle = "Dashed line = Mean | Solid line = Median",
    x = "Copy number",
    y = "Density",
    fill = "Region type",
    color = "Region type"
  ) +
  xlim(0,500)+
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )
