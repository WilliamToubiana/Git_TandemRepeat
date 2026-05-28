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




########################################
# Plots TRprop vs recomb vs GC content #
########################################

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




