library(ggplot2)
library(ggpubr)

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed")
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution")

data <- read.table("Reseq_TR_prop_minimap2_TRpaper.txt", sep="\t", 
                   header=TRUE, fill=TRUE)
data$sp_ind<-paste(data$sp,data$ind)
#data$prop<-data$covTR_minimap/data$covGW_minimap
data$prop5<-data$covTR_minimap_5copies/data$covGW_minimap
head(data)


###############################################################
### Relationships between TR proportion, genome size and Ne ### 
###############################################################

## TR proportion between sexual species
data$species=factor(data$sp, levels=c("Tps", "Tcm", "Tce", "Tpa", "Tbi"))

TRprop<-ggplot(data, aes(x=species, y=prop5)) +
  geom_boxplot()+ ylab("TR proportion")+
  theme_minimal()+ylim(0.025,0.2)


## Genome size estimates between sexual species

data$estimated_genome_size<-data$covGW_minimap/data$mean_coverage_GW #mean_coverage_GW calculated as:cut -f3 /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt | awk '{sum+=$1} END {print sum/NR}' 
#data$estimated_genome_size<-data$covGW_minimap/data$mode_coverage_GW #mode_coverage_GW calculated as:awk '{print $3}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tbi-filtered_indG_minimap2_GW_coverage.txt | sort | uniq -c | sort -nr | head -n 1

GenomeSex<-ggplot(data, aes(x=species, y=estimated_genome_size)) +
  geom_boxplot()+ ylab("Estimated Genome size")+
  theme_minimal()+ylim(1.05e+09,1.4e+09)


## Estimated effective population size between species
setwd("~/Desktop/Presentation/genomic_days")
data2<-read.table("Estimated_pop_size.txt", header=T, sep="\t", quote="")

data2$species<-factor(data2$species, levels=c("T.poppensis", "T.californicum", "T.cristinae", "T.podura", "T.bartmani"))
Ne<-ggplot(data=data2, aes(x=species, y=Ne)) +
  geom_bar(stat="identity", color="black", fill="white")+ 
  ylab("Effective population size")


## Combined TRprop, genomeSize and Ne plots
library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    TRprop + theme(legend.position = "none")
    , GenomeSex
    , Ne+ theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(GenomeSex)
    , ggplot()
    , get_legend(Ne)
    , ncol =1)
  , rel_widths = c(12,0)
)


## Estimated mean autosome size (as proxy for recombination rate)

table<-data.frame(
  species = c("T.poppense", "T.californicum", "T.cristinae (old)", "T.cristinae (new)", "T.podura"),
  genome_size_noX = c(1207460021,1158934791,1087191189,1108158110,1017359884),
  autosome_size = c(1130673392, 1090713397, NA, 1078799922, 976679347),
  chromosome_number = c(12, 12, 13, 13, 14)
)

table$mean_autosome_size_genome<-table$genome_size_noX/table$chromosome_number
table$mean_autosome_size_autosome<-table$autosome_size/table$chromosome_number

table$species<-factor(table$species, levels=c("T.poppense", "T.californicum", "T.cristinae (old)", "T.cristinae (new)", "T.podura"))

ggplot(table, aes(x=species, y=mean_autosome_size_genome)) +
  geom_boxplot()+ ylab("Mean autosome size")+
  theme_minimal()

ggplot(table, aes(x=species, y=mean_autosome_size_autosome)) +
  geom_boxplot()+ ylab("Mean autosome size")+
  theme_minimal()

## Correlation between TRprop and genomeSize/Ne

data$estimated_TR_size<-data$covTR_minimap_5copies/data$mean_coverage_TR #mean_coverage_TR calculated as:cut -f3 /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tbi-filtered_indG_minimap2_TR_coverage.txt | awk '{sum+=$1} END {print sum/NR}' 
data$estimated_nnTR_size<-data$estimated_genome_size-data$estimated_TR_size

#Summarize data
summary_df <- data %>%
  group_by(species) %>%
  summarise(
    genome_mean = mean(estimated_genome_size, na.rm = TRUE),
    genome_sd = sd(estimated_genome_size, na.rm = TRUE),
    tr_mean = mean(estimated_TR_size, na.rm = TRUE),
    tr_sd = sd(estimated_TR_size, na.rm = TRUE),
    notr_mean = mean(estimated_nnTR_size, na.rm = TRUE),
    notr_sd = sd(estimated_nnTR_size, na.rm = TRUE),
    tr_prop = mean(prop5, na.rm = TRUE),
    .groups = "drop"
  )
    
summary_df$Ne<-data2$Ne
summary_df <- as.data.frame(summary_df)
rownames(summary_df) <- summary_df$species


   
# I need to scale the second variable so it fits on the same plot
# Example: divide effective population size to bring it into same range
#scaling_factor <- max(summary_df$genome_mean) / max(summary_df$Ne)

#ggplot(summary_df, aes(x = tr_prop)) +
#  geom_point(aes(y = genome_mean), color = "steelblue", size = 3) +
#  geom_line(aes(y = genome_mean), color = "steelblue") +
#  geom_point(aes(y = Ne * scaling_factor), color = "darkred", size = 3, shape = 17) +
#  geom_line(aes(y = Ne * scaling_factor), color = "darkred", linetype = "dashed") +
#  scale_y_continuous(
#    name = "Genome Size",
#    sec.axis = sec_axis(~ . / scaling_factor, name = "Effective Population Size")
#  ) +
#  xlab("TR proportion") +
#  theme_minimal() +
#  theme(
#    axis.title.y.left = element_text(color = "steelblue"),
#    axis.text.y.left = element_text(color = "steelblue"),
#    axis.title.y.right = element_text(color = "darkred"),
#    axis.text.y.right = element_text(color = "darkred"),
#  )




###########################################################################
### Ne and Genome size variation explained by variation in TR abundance ###
###########################################################################
#Plot
ggplot(summary_df, aes(x = Ne, y = tr_prop, label = species)) +
  geom_point(size = 0.05) +
  geom_text(vjust = -1, size = 3) +
  xlab("Effective Population Size") +
  ylab("TR proportion") +
  theme_minimal() +
  geom_smooth(method='lm',formula=y~x) + ylim(0,0.18)

cor.test(summary_df$tr_prop, summary_df$Ne, method = "pearson")

model<-lm(summary_df$tr_prop~summary_df$Ne)
summary(model)



#PLot
ggplot(summary_df, aes(x = genome_mean, y = tr_prop, label = species)) +
  geom_point(size = 0.05) +
  geom_text(vjust = -1, size = 3) +
  xlab("Genome Size Estimate") +
  ylab("TR proportion") +
  theme_minimal() +
  geom_smooth(method='lm',formula=y~x)+ ylim(0,0.18)

#Linear model to estimate the variation in genome size explained by TR size variation
model <- lm(genome_mean ~ tr_prop, data = summary_df)
summary(model)

cor.test(summary_df$tr_prop, summary_df$genome_mean, method = "pearson")






###Stats with phylogenetic correction
library(ape)
library(phytools)
library(caper)

newick_str <- "(Tbi:0.00174220447502769168,(((Tms:0.00358987165276624613,Tce:0.00304432455403352408)100:0.01271016184978783839,((Tdi:0.00335602954628538971,Tps:0.00392474709480761241)100:0.00576515566828786873,(Tcm:0.00310888700338952072,Tsi:0.00267156539116077075)100:0.00705381430340124016)100:0.01046875224079147140)100:0.02522943290742237984,(Tge:0.00531692311750064234,Tpa:0.00583613471885804924)100:0.00376148450849715177)100:0.00636176950690947769,Tte:0.00197612024312248573);"

tree <- read.tree(text = newick_str)

keep_species <- c("Tbi", "Tpa", "Tce", "Tcm", "Tps") #keep only sexual species
pruned_phy <- drop.tip(tree, setdiff(tree$tip.label, keep_species))
plot(pruned_phy, show.node.label = TRUE, cex = 0.9)
edgelabels(round(pruned_phy$edge.length, 4), cex = 0.6)

# Comparative data object
pruned_phy$node.label <- NULL #Remove internal node labels (support values)

comp_data <- comparative.data(
  phy = pruned_phy,
  data = summary_df,
  names.col = "species",
  vcv = TRUE
)

# Run PGLS for TRprop vs Ne (Phylogenetic Generalised Least Square)
pgls_model <- pgls(tr_prop ~ Ne, data = comp_data)
summary(pgls_model)

# Run PGLS for genome size vs TR abundance estimates (Phylogenetic Generalised Least Square) 
pgls_model <- pgls(genome_mean ~ tr_prop, data = comp_data)
summary(pgls_model)


#corrected estimates
genome_pic <- pic(summary_df$genome_mean, pruned_phy) #corrected values for genome size estimates
#tr_pic <- pic(summary_df$tr_mean, pruned_phy)#corrected values for TR size estimates
tr_pic <- pic(summary_df$tr_prop, pruned_phy)#corrected values for TR proportion estimates
ne_pic <- pic(summary_df$Ne, pruned_phy)#corrected values for Ne estimates

#Correlation on corrected estimates
cor.test(tr_pic, genome_pic, method = "pearson")
cor.test(tr_pic, ne_pic, method = "pearson")



#################################################
### Correlation Reseq vs TRF TRprop estimates ###
#################################################

##Tps
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tps$array_length<-(Tps$end-Tps$start)+1

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tps$chr,ranges=IRanges(start=Tps$start,end=Tps$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind_sum<-aggregate(nonoverlapWind$width~nonoverlapWind$seqnames, FUN=sum)
Tps<-sum(nonoverlapWind_sum$`nonoverlapWind$width`)/1362169457 #genome size estimated from the genome assembly (grep -v ">" file.fasta | wc -c)


##Tcm
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tcm<-read.table("Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tcm$array_length<-(Tcm$end-Tcm$start)+1

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tcm$chr,ranges=IRanges(start=Tcm$start,end=Tcm$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind_sum<-aggregate(nonoverlapWind$width~nonoverlapWind$seqnames, FUN=sum)
Tcm<-sum(nonoverlapWind_sum$`nonoverlapWind$width`)/1317820820 #genome size estimated from the genome assembly (grep -v ">" file.fasta | wc -c)


##Tce
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tce<-read.table("Tcr_CEN4280_h2.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tce$array_length<-(Tce$end-Tce$start)+1

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tce$chr,ranges=IRanges(start=Tce$start,end=Tce$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind_sum<-aggregate(nonoverlapWind$width~nonoverlapWind$seqnames, FUN=sum)
Tce<-sum(nonoverlapWind_sum$`nonoverlapWind$width`)/1257119929 #genome size estimated from the genome assembly (grep -v ">" file.fasta | wc -c)


##Tpa
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tpa<-read.table("Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tpa$array_length<-(Tpa$end-Tpa$start)+1

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tpa$chr,ranges=IRanges(start=Tpa$start,end=Tpa$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind_sum<-aggregate(nonoverlapWind$width~nonoverlapWind$seqnames, FUN=sum)
Tpa<-sum(nonoverlapWind_sum$`nonoverlapWind$width`)/1164797593 #genome size estimated from the genome assembly (grep -v ">" file.fasta | wc -c)


##Tbi
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tbi<-read.table("Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tbi$array_length<-(Tbi$end-Tbi$start)+1

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tbi$chr,ranges=IRanges(start=Tbi$start,end=Tbi$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind_sum<-aggregate(nonoverlapWind$width~nonoverlapWind$seqnames, FUN=sum)
Tbi<-sum(nonoverlapWind_sum$`nonoverlapWind$width`)/1228449545 #genome size estimated from the genome assembly (grep -v ">" file.fasta | wc -c)

summary<-data.frame(
  species=c("Tps","Tcm","Tce","Tpa","Tbi"),
  TRprop=c(Tps,Tcm,Tce,Tpa,Tbi))
summary$species<-factor(summary$species, levels=c("Tps", "Tcm", "Tce", "Tpa", "Tbi"))

ggplot(summary, aes(x=species, y=TRprop)) +
  geom_boxplot()+ ylab("TR proportion")+
  theme_minimal()+ylim(0.025,0.2)
