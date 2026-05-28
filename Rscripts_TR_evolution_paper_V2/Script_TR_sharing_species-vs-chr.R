##################################################
### Generate a presence (1) absence (0) matrix ###
##################################################
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")

Tps<-read.table("Tps_cleaned.txt", header=T, sep="\t", quote="")
Tcm<-read.table("Tcm_cleaned.txt", header=T, sep="\t", quote="")
Tce<-read.table("Tce_cleaned.txt", header=T, sep="\t", quote="")
Tpa<-read.table("Tpa_cleaned.txt", header=T, sep="\t", quote="")
Tbi<-read.table("Tbi_cleaned.txt", header=T, sep="\t", quote="")
Tdi<-read.table("Tsi_cleaned.txt", header=T, sep="\t", quote="")


# Combine all motifs from the 6 dataframes and get unique motifs
all_motifs <- unique(c(Tps$aomfi, Tcm$aomfi, Tce$aomfi, Tpa$aomfi, Tbi$aomfi, Tdi$aomfi))

# Sort the motifs
all_motifs <- sort(all_motifs)

# Create an empty matrix for presence-absence (motifs as rows, species as columns)
presence_absence_matrix <- matrix(0, nrow = length(all_motifs), ncol = 6)

# Assign row names (motifs) and column names (species)
rownames(presence_absence_matrix) <- all_motifs
colnames(presence_absence_matrix) <- c("Tps", "Tcm", "Tce", "Tpa", "Tbi", "Tdi")

# Fill the matrix with 1s for presence of motifs in each species
presence_absence_matrix[all_motifs %in% Tps$aomfi, 1] <- 1
presence_absence_matrix[all_motifs %in% Tcm$aomfi, 2] <- 1
presence_absence_matrix[all_motifs %in% Tce$aomfi, 3] <- 1
presence_absence_matrix[all_motifs %in% Tpa$aomfi, 4] <- 1
presence_absence_matrix[all_motifs %in% Tbi$aomfi, 5] <- 1
presence_absence_matrix[all_motifs %in% Tdi$aomfi, 6] <- 1

head(presence_absence_matrix)


#convert to dataframe 
presence_absence_df<-data.frame(presence_absence_matrix)
presence_absence_df$aomfi<-rownames(presence_absence_df)

#select only motifs present in a single species and sum rows
presence_absence_df_Tps<-subset(presence_absence_df, Tps=="1") #118243 motifs in total
presence_absence_df_Tps$shared<-rowSums( presence_absence_df_Tps[,c(1,2,3,4,5,6)])#shared across Timema (T.douglasi excluded)


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


## Distribution of motif sizes based on chromosome sharing

mergeTps_subset<-subset(mergeTps, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                        | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                        | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

motif_chrom_counts <- mergeTps_subset %>%
  distinct(motif_representant, chr) %>%  # Remove duplicate motif-chromosome pairs
  group_by(motif_representant) %>%
  summarise(num_chromosomes = n()) # Count unique chromosomes per motif

motif_chrom_counts$motif_length<-nchar(motif_chrom_counts$motif_representant)

ggplot(motif_chrom_counts, aes(x = motif_length, 
                               color = as.factor(num_chromosomes), 
                               fill = as.factor(num_chromosomes))) + 
  geom_density(alpha = 0.3) +
  labs(
    x = "Motif Size (bp)",
    y = "Density",
    color = "Chromosome sharing",
    fill = "Chromosome sharing",
    title = "Density of Motif Sizes based on chromosome sharing"
  ) +
  theme_minimal()+xlim(0,150)


motif_chrom_counts_subset1<-subset(motif_chrom_counts, num_chromosomes=="1")
motif_chrom_counts_subset11<-subset(motif_chrom_counts, num_chromosomes=="11")
motif_chrom_counts_subset12<-subset(motif_chrom_counts, num_chromosomes=="12")

shared1<-ggplot(motif_chrom_counts_subset1, aes(x = motif_length, 
                                                color = as.factor(num_chromosomes), 
                                                fill = as.factor(num_chromosomes))) + 
  geom_density(alpha = 0.3) +
  labs(
    x = "Motif Size (bp)",
    y = "Density",
  ) +
  theme(legend.position="none")+xlim(0,150)

shared11<-ggplot(motif_chrom_counts_subset11, aes(x = motif_length, 
                                                  color = as.factor(num_chromosomes), 
                                                  fill = as.factor(num_chromosomes))) + 
  geom_density(alpha = 0.3) +
  labs(
    x = "Motif Size (bp)",
    y = "Density",
  ) +
  theme(legend.position="none")+xlim(0,150)

shared12<-ggplot(motif_chrom_counts_subset12, aes(x = motif_length, 
                                                  color = as.factor(num_chromosomes), 
                                                  fill = as.factor(num_chromosomes))) + 
  geom_density(alpha = 0.3) +
  labs(
    x = "Motif Size (bp)",
    y = "Density",
  ) +
  theme(legend.position="none")+xlim(0,150)


library("gridExtra")
grid.arrange(shared1, shared11, shared12, 
             ncol = 3, nrow = 1)

## Number of motifs shared among species vs genome representation 

ggplot(mergeTps, aes(x=shared, y=nMotifsRepresented)) +
  geom_point(size=2, shape=23)+ ylab("number of representation in the genome")+ 
  xlab("number of species")+
  ylim(0,61000)+
  ggtitle("shared motifs")+
  theme()



