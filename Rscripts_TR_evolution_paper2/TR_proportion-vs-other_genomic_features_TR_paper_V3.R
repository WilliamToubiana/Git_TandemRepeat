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


##Merge TR and recomb datasets
merge<-merge(merge3TRprop, rho_g, by="chm_pos")


rho_plot<-ggplot(merge, aes(x=cumul_start, y=corr_rho, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(merge$scaff_number)))/2))+
  ylab("rho")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


###Linear Model and correlation raw data

model<-lm(TRprop ~ log10(corr_rho) * chr, data = merge)
summary(model)

cor.test(log10(merge$corr_rho), merge$TRprop, method = "spearman")



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
data_GC2<-data_GC[c(1,2,6,7,8,9)]
colnames(data_GC2)<-c("chromosome", "start", "A", "C", "G","T")
summary(data_GC2)

data_GC2$GC_perc<-(data_GC2$C+data_GC2$G)/(data_GC2$A+data_GC2$C+data_GC2$G+data_GC2$T)
summary(data_GC2)

data_GC2$start2<-data_GC2$start+1
data_GC2$chm_pos<-paste(data_GC2$chromosome, data_GC2$start2, sep="-")


#chromosomes only
chms<-subset(data_GC2, chromosome=="Tps_LRv5b_scf1" | chromosome=="Tps_LRv5b_scf2" | chromosome=="Tps_LRv5b_scf3" | chromosome=="Tps_LRv5b_scf4" | chromosome=="Tps_LRv5b_scf5" 
             | chromosome=="Tps_LRv5b_scf6" | chromosome=="Tps_LRv5b_scf7" | chromosome=="Tps_LRv5b_scf8" | chromosome=="Tps_LRv5b_scf9" | chromosome=="Tps_LRv5b_scf10"
             | chromosome=="Tps_LRv5b_scf11" | chromosome=="Tps_LRv5b_scf12")

summary(chms)


data_GC3 <- data_GC %>%
  group_by(V1)%>%
  summarize(V6=sum(V6),
            V7=sum(V7),
            V8=sum(V8),
            V9=sum(V9))

colnames(data_GC3)<-c("chromosome", "A", "C", "G","T")
summary(data_GC3)

data_GC3$GC_perc<-(data_GC3$C+data_GC3$G)/(data_GC3$A+data_GC3$C+data_GC3$G+data_GC3$T)
summary(data_GC3)

chms2<-subset(data_GC3, chromosome=="Tps_LRv5b_scf1" | chromosome=="Tps_LRv5b_scf2" | chromosome=="Tps_LRv5b_scf3" | chromosome=="Tps_LRv5b_scf4" | chromosome=="Tps_LRv5b_scf5" 
             | chromosome=="Tps_LRv5b_scf6" | chromosome=="Tps_LRv5b_scf7" | chromosome=="Tps_LRv5b_scf8" | chromosome=="Tps_LRv5b_scf9" | chromosome=="Tps_LRv5b_scf10"
             | chromosome=="Tps_LRv5b_scf11" | chromosome=="Tps_LRv5b_scf12")

##Merge TR and recomb datasets
merge2<-merge(merge3TRprop, chms, by="chm_pos")

GC_content<-ggplot(merge2, aes(x=cumul_start, y=GC_perc, color = as.factor(scaff_number)))+
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


###Linear Model and correlation raw data

model<-lm(TRprop ~ perc_GC * chr, data = merge2)
summary(model)

cor.test(merge2$perc_GC, merge2$TRprop, method = "spearman")



########################################
# Plots TRprop vs recomb vs GC content #
########################################

library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    TR_array + theme(legend.position = "none")
    , rho_plot
    , GC_content + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(rho_plot)
    , ggplot()
    , get_legend(GC_content)
    , ncol =1)
  , rel_widths = c(12,0)
)





#############################
## Relationship TR metrics ##
#############################


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
df_filtered2$array_size<-abs(df_filtered2$end-df_filtered2$start)+1
df_filtered2$copy_nb2<-df_filtered2$array_size/df_filtered2$L_motif_representant

#Assign to each repeat array a non-overlapping 250kb window
df_filtered2$window_start <- (floor(df_filtered2$start / 250000) * 250000)+1
df_filtered2$window_end <- df_filtered2$window_start + 249999  # End of the 250kb window
df_filtered2$chm_pos<-paste(df_filtered2$chr, df_filtered2$window_start, sep="-")


merge3<-merge(merge2, rho_g, by="chm_pos")
merge4<-merge(df_filtered2, merge3[c(1,8,17,21,22)], by="chm_pos")#all.x=T recovers TR sequences associated with NA values for rho,TRprop and GC.
summary(merge4)

#per window
window_summary <- merge4 %>%
  group_by(chr, window_start) %>%
  summarise(
    window_end = unique(window_end),
    total_repeats = n(),
    avg_array_size = mean(array_size),
    median_array_size = median(array_size),
    avg_motif = mean(L_motif_representant),
    median_motif = median(L_motif_representant),
    median_CN = median(copy_nb2),
    mean_CN = mean(copy_nb2),
    motif_diversity = n_distinct(motif_representant),
    median_rho = unique(median_rho),
    GC_perc = unique(GC_perc),
    .groups = "drop"
  )

#per chr
window_summary <- merge4 %>%
  group_by(chr) %>%
  summarise(
    total_repeats = n(),
    avg_array_size = mean(array_size),
    median_array_size = median(array_size),
    avg_motif = mean(L_motif_representant),
    median_motif = median(L_motif_representant),
    median_CN = median(copy_nb2),
    mean_CN = mean(copy_nb2),
    motif_diversity = n_distinct(motif_representant),
    median_rho = median(median_rho),
    mean_rho = mean(median_rho),
    GC_perc = mean(GC_perc),
    .groups = "drop"
  )

summary(window_summary)
window_summary$region<-window_summary$window_end-window_summary$window_start
window_summary$array_density<-window_summary$total_repeats/window_summary$region

window_summary$LG_length <- c(253718452, 72695627, 67869710, 42522306, 178265049, 132376913, 102503173, 89066744, 87859127, 83130533, 76885976, 76156695)
window_summary$array_density<-window_summary$total_repeats/window_summary$LG_length

#Relationships with rho
ggplot(window_summary,aes(log10(median_rho), log10(avg_motif))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length (log)")+
  xlab("rho values (log)")

model<-lm(log10(window_summary$avg_motif)~log10(window_summary$median_rho))
summary(model)


ggplot(window_summary,aes(log10(median_rho), log10(mean_CN))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean copy number (log)")+
  xlab("rho values (log)")

model<-lm(log10(window_summary$mean_CN)~log10(window_summary$median_rho))
summary(model)


ggplot(window_summary,aes(log10(median_rho), log10(array_density))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR array density (log)")+
  xlab("rho values (log)")

model<-lm(log10(window_summary$array_density)~log10(window_summary$median_rho))
summary(model)



#Relationships with GC
ggplot(window_summary,aes(GC_perc, log10(avg_motif))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length (log)")+
  xlab("GC percent")

model<-lm(log10(window_summary$avg_motif)~window_summary$GC_perc)
summary(model)


ggplot(window_summary,aes(GC_perc, log10(mean_CN))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean copy number (log)")+
  xlab("GC percent")

model<-lm(log10(window_summary$mean_CN)~window_summary$GC_perc)
summary(model)


ggplot(window_summary,aes(GC_perc, log10(array_density))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR array density (log)")+
  xlab("GC percent")

model<-lm(log10(window_summary$array_density)~window_summary$GC_perc)
summary(model)


##scaled values
window_summary$avg_motif_scaled<-(window_summary$avg_motif-mean(window_summary$avg_motif))/sd(window_summary$avg_motif)
window_summary$mean_CN_scaled<-(window_summary$mean_CN-mean(window_summary$mean_CN))/sd(window_summary$mean_CN)
window_summary$array_density_scaled<-(window_summary$array_density-mean(window_summary$array_density))/sd(window_summary$array_density)

ggplot(window_summary,aes(log10(median_rho), avg_motif_scaled)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length (log)")+
  xlab("rho values (log)")

model<-lm((window_summary$avg_motif_scaled)~log10(window_summary$median_rho))
summary(model)


ggplot(window_summary,aes(log10(median_rho), (mean_CN_scaled))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean copy number (log)")+
  xlab("rho values (log)")

model<-lm((window_summary$mean_CN_scaled)~log10(window_summary$median_rho))
summary(model)


ggplot(window_summary,aes(log10(median_rho), (array_density_scaled))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR array density (log)")+
  xlab("rho values (log)")

model<-lm((window_summary$array_density_scaled)~log10(window_summary$median_rho))
summary(model)


#Relationships with GC
ggplot(window_summary,aes(GC_perc, (avg_motif_scaled))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean TR sequence length (log)")+
  xlab("GC percent")

model<-lm((window_summary$avg_motif_scaled)~window_summary$GC_perc)
summary(model)


ggplot(window_summary,aes(GC_perc, (mean_CN_scaled))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+
  ylab("mean copy number (log)")+
  xlab("GC percent")

model<-lm((window_summary$mean_CN_scaled)~window_summary$GC_perc)
summary(model)


ggplot(window_summary,aes(GC_perc, (array_density_scaled))) + 
  geom_point() + 
  geom_smooth(method=lm) +
  theme_classic(base_size = 17)+ 
  ylab("TR array density (log)")+
  xlab("GC percent")

model<-lm((window_summary$array_density_scaled)~window_summary$GC_perc)
summary(model)
