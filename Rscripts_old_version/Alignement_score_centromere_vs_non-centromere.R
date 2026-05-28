#############################
# Mapping 250kb Tps windows #
#############################

##Tandem repeats
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/whole_genome_alignment")

data<-read.table("Tps-genome_250kb_to_Tcm_genome_minimap2_s40.paf", header=F, sep="\t", quote="", fill = TRUE)
#data<-read.table("Tps-genome_250kb_to_Tcm_genome_minimap2_s2500.paf", header=F, sep="\t", quote="", fill = TRUE)
colnames(data) <- c("query", "query_length", "query_start", "query_end", "strand", "target", "target_length", "target_start", "target_end",
                    "number_matches", "alignment_block_length","mapping_quality","edit_distance","matching_score","alignment_score","number_ambiguous_bases",
                    "type_alignment","V18","V19","V20","V21","V22","V23","V24")

data$cov<-abs(data$query_start-data$query_end)

library(stringr)
data[c('AS', 'i', 'alignment_score2')] <- str_split_fixed(data$alignment_score, ':', 3)
data[c('NM', 'i', 'edit_distance2')] <- str_split_fixed(data$edit_distance, ':', 3)
data[c('ms', 'i', 'matching_score2')] <- str_split_fixed(data$matching_score, ':', 3)
summary(data)

data$alignment_score2<-as.numeric(data$alignment_score2)
data$edit_distance2<-as.numeric(data$edit_distance2)
data$matching_score2<-as.numeric(data$matching_score2)
summary(data)

library(dplyr)
window_summary <- data %>%
  group_by(query) %>%
  summarise(
    alignment_score_avg = mean(alignment_score2),
    edit_distance_avg = mean(edit_distance2),
    matching_score_avg = mean(matching_score2),
  )

summary(window_summary)

window_summary[c('chr', 'start_end')] <- str_split_fixed(window_summary$query, ':', 2)
window_summary[c('window_start', 'window_end')] <- str_split_fixed(window_summary$start_end, '-', 2)
window_summary$window_start<-as.numeric(window_summary$window_start)
window_summary$window_start2<-window_summary$window_start+1
window_summary$chm_pos<-paste(window_summary$chr, window_summary$window_start2, sep="-")



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
data1$log2cenH3_norm[!is.finite(data1$log2cenH3_norm)] <- 0




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
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
summary(mergeTRprop)



######################################################
### Alignment centromere vs non-centromere regions ###
######################################################

##Combine TR prop and shared repeat datasets
combined<-merge(mergeTRprop, window_summary[c(2,3,4,10)], by="chm_pos", all.x=T)
summary(combined)

##Subdivide windows into centromere vs non-centromere based on log2 ratios
data1$centromere<-ifelse(data1$log2cenH3_norm>0.25, "centromere", "non-centromere")

##Combine TR prop, centromere and shared repeat datasets
combined2<-merge(combined, data1[c(9,11,12)], by="chm_pos", all.x=T)
summary(combined2)





##Scatterplot and density plots with TR proportion
##################################################
combined2<-subset(combined2, TRprop>0)
combined2<-subset(combined2, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                  | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                  | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

scatterPlot <- ggplot(combined2,aes(TRprop, alignment_score_avg, color=centromere)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  geom_smooth(method=lm) +
  theme(legend.position=c(0,1), legend.justification=c(-5.5,-1))+
  ylab("Alignment score")

# Marginal density plot of x (top panel)
xdensity <- ggplot(combined2, aes(TRprop, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined2, aes(alignment_score_avg, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")+coord_flip()+
  xlab("Alignment score")




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
model<-lm(combined2$alignment_score_avg~combined2$TRprop*combined2$centromere)
anova(model)




#####################################################
## Distribution alignment scores along chromosomes ##
#####################################################

combined3<-subset(combined2, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                                | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                                | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined3$scaff_number<-as.numeric(as.character(gsub("^.*scf","", combined3$chr)))

combined3<-combined3[order(combined3$scaff_number, combined3$start2),]
combined3$cumul_start<-seq(0, by = 250000, length.out = nrow(combined3))



library(ggplot2)
alignment<-ggplot(combined3, aes(x=cumul_start, y=alignment_score_avg, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(combined3$scaff_number)))/2))+
  ylab("Alignment score")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



cenh3_testes<-ggplot(combined3, aes(x=cumul_start, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(combined3$scaff_number)))/2))+
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept= 0.25, linetype="dashed", color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


TR_array<-ggplot(combined3, aes(x=cumul_start, y=TRprop, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(combined3$scaff_number)))/2))+
  ylab("Proportion of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



# Plots TRprop vs recomb vs GC content
library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    TR_array + theme(legend.position = "none")
    , cenh3_testes
    , alignment + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(alignment)
    , ncol =1)
  , rel_widths = c(12,0)
)











