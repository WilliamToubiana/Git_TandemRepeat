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
  geom_boxplot(width=0.4)+ ylab("TR proportion")+
  geom_point(position=position_dodge(width=0.75), size=1)+
  theme_minimal()+ylim(0.025,0.2)


## Genome size estimates between sexual species

data$estimated_genome_size<-data$covGW_minimap/data$mean_coverage_GW #mean_coverage_GW calculated as:cut -f3 /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt | awk '{sum+=$1} END {print sum/NR}' 
#data$estimated_genome_size<-data$covGW_minimap/data$mode_coverage_GW #mode_coverage_GW calculated as:awk '{print $3}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tbi-filtered_indG_minimap2_GW_coverage.txt | sort | uniq -c | sort -nr | head -n 1

GenomeSex<-ggplot(data, aes(x=species, y=estimated_genome_size)) +
  geom_boxplot(width=0.4)+ 
  geom_point(position=position_dodge(width=0.75), size=1)+
  ylab("Estimated Genome size")+
  theme_minimal()+ylim(1.05e+09,1.4e+09)


## Estimated effective population size between species
setwd("~/Desktop/Presentation/genomic_days")
data2<-read.table("Estimated_pop_size.txt", header=T, sep="\t", quote="")

data2$species<-factor(data2$species, levels=c("T.poppensis", "T.californicum", "T.cristinae", "T.podura", "T.bartmani"))
Ne<-ggplot(data=data2, aes(x=species, y=Ne)) +
  geom_bar(stat="identity", color="black", fill="white", width = 0.4)+ 
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
tree.mid <- midpoint.root(tree)

keep_species <- c("Tbi", "Tpa", "Tce", "Tcm", "Tps") #keep only sexual species
pruned_phy <- drop.tip(tree.mid, setdiff(tree.mid$tip.label, keep_species))
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

