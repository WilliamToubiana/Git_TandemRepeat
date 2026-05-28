library("ggVennDiagram")
library(ggplot2)
library(GenomicRanges)
library(dplyr)


#############################
# TR characteristics in Tps #
#############################

### Chromosome specificity of TR motifs

##Data preparation GW
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

Tps_uniq<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse2_5copies_CMwRC.txt", header=T, sep="\t", quote="")
head(Tps_uniq)

Tps_uniq$total_length<-Tps_uniq$Counts_s*Tps_uniq$monomer_length
Tps_uniq$proportion<-(Tps_uniq$total_length/1339836934)*100


#Removing duplicate words in "Scaffold_name"
Tps_uniq$Scaffold_name = as.data.frame(sapply(Tps_uniq$Scaffold_name, function(x) gsub("\"", "", x))) #Removing quotes

reduce_row = function(i) {
  split = strsplit(i, split=",")[[1]]
  paste(unique(split), collapse =
          ",") 
}

Tps_uniq$Scaffold_name = apply(Tps_uniq, 1, reduce_row) #Removing duplicates
Tps_uniq<-Tps_uniq[order(Tps_uniq$Counts_s, decreasing = TRUE),]

#data<-subset(data, monomer_length>100)
#data<-subset(data, monomer_length>10 & monomer_length<101)

#Count number of anchored scaffolds (from Tps_LRv5_scf1 to Tps_LRv5_scf12 = LGs) per monomer

Tps_uniq$Count_LGs <- rowSums(sapply(c("\\bTps_LRv5b_scf1\\b", "\\bTps_LRv5b_scf2\\b", "\\bTps_LRv5b_scf3\\b", "\\bTps_LRv5b_scf4\\b", "\\bTps_LRv5b_scf5\\b", "\\bTps_LRv5b_scf6\\b", "\\bTps_LRv5b_scf7\\b", "\\bTps_LRv5b_scf8\\b", "\\bTps_LRv5b_scf9\\b", "\\bTps_LRv5b_scf10\\b", "\\bTps_LRv5b_scf11\\b", "\\bTps_LRv5b_scf12\\b"),
                                     function(x) grepl(x, Tps_uniq$Scaffold_name)))

Tps_uniq$Count_LGs=factor(Tps_uniq$Count_LGs, levels=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))


##Data preparation centromere
setwd("/Users/wtoubian/Desktop/Chip-seq/species-comparison/results_Marion_24-05-2024/enriched_motifs/levenstein/Tps_testes_cenh3_R1")
table<-read.table("Network_database.txt", header=T, sep="\t", quote="")


# Function to check presence of centromere repeats in the GW dataset (from the dataset of enriched centromere repeats)
check_presence <- function(row, df2) {
  # Split the comma-separated string into a vector
  strings_vector <- unlist(strsplit(row, ","))
  # Check if any string from df2 is in the strings_vector
  if (any(table$motif %in% strings_vector)) {
    return("centromere repeat")
  } else {
    return("no centromere repeat")
  }
}


# Apply the function to each row of df1
Tps_uniq$feature <- sapply(Tps_uniq$monomer_seqs, check_presence, df2 = table)


#subset dataset
centromere<-subset(Tps_uniq, feature=="centromere repeat")
non_centromere<-subset(Tps_uniq, feature=="no centromere repeat")

#Counts
centromere_c<-count(centromere, Count_LGs)
non_centromere_c<-count(non_centromere, Count_LGs)
  

##Merge the two datasets and plots

summary<-merge(non_centromere_c, centromere_c, by="Count_LGs", all.x=T)
summary[is.na(summary)] <- 0
colnames(summary)<-c("Count_LGs", "non centromere", "centromere")
summary$feature1<-"non centromere"
summary$feature2<-"centromere"

Count_LGs_repeated <- rep(summary$Count_LGs, times = 2)

# Combine the second and third columns using rbind
combined <- c(summary$`non centromere`, summary$centromere)
combined2 <- c(summary$feature1, summary$feature2)

# Create the new data frame
data <- data.frame(
  Count_LGs = Count_LGs_repeated,
  counts = combined,
  feature = combined2
)



ggplot(data=data, aes(x=Count_LGs, y=counts))+ 
  geom_bar(stat="identity", aes(fill = feature))+
  theme_minimal()

ggplot(data=data, aes(x=Count_LGs, y=counts))+ 
  geom_bar(stat="identity", aes(fill = feature))+
  scale_y_continuous(limits=c(0, 1300))+
  theme_minimal()




###Tandem repeat proportion between chromosomes Timema
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tps_subset<-subset(Tps, copy_nb >= 5)
Tps_subset$array_length<-(Tps_subset$end-Tps_subset$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tps_subset$chr,ranges=IRanges(start=Tps_subset$start,end=Tps_subset$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tps_LRv5b_scf1" | seqnames=="Tps_LRv5b_scf2" | seqnames=="Tps_LRv5b_scf3" | 
                             seqnames=="Tps_LRv5b_scf4" |seqnames=="Tps_LRv5b_scf5" | seqnames=="Tps_LRv5b_scf6" | seqnames=="Tps_LRv5b_scf7" | 
                             seqnames=="Tps_LRv5b_scf8" | seqnames=="Tps_LRv5b_scf9" | seqnames=="Tps_LRv5b_scf10" | seqnames=="Tps_LRv5b_scf11" | 
                             seqnames=="Tps_LRv5b_scf12")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)

#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
Tps_LG_length<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")

table=data.frame(c("chr1", "chr2", "chrX", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chrX", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))
table$LG_length <- Tps_LG_length$length
table$TR_length <- nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`
table$Prop_repeated_region<-table$TR_length/table$LG_length

ggplot(data=table, aes(x=LG_length, y=Prop_repeated_region)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()


ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()

#stats all chromosomes
model<-lm(table$Prop_repeated_region~table$LG_length)
summary(model)

#stats without the X
table2<-subset(table, LGs!="chr3")
model<-lm(table2$Prop_repeated_region~table2$LG_length)
summary(model)



#number of motifs between chromosomes

Tps_subset_chr<- subset(Tps_subset, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | 
                          chr=="Tps_LRv5b_scf4" |chr=="Tps_LRv5b_scf5" | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | 
                          chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10" | chr=="Tps_LRv5b_scf11" | 
                          chr=="Tps_LRv5b_scf12")

nonoverlapWind_chr$count<-1
nonoverlapWind_chr_sum2<-aggregate(nonoverlapWind_chr$count~nonoverlapWind_chr$seqnames, FUN=sum)

Tps_subset_chr$count<-1
Tps_subset_chr_sum<-aggregate(Tps_subset_chr$count~Tps_subset_chr$chr, FUN=sum)


table$motif_count<-Tps_subset_chr_sum$`Tps_subset_chr$count`
table$motif_count<-nonoverlapWind_chr_sum2$`nonoverlapWind_chr$count`


ggplot(data=table, aes(x=LGs, y=motif_count)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(motif_count, digits = 3)), vjust=1.6, size=3)+ ylab("number of tandem repear arrays")+
  theme_classic()




#array/motif length and copy number between chromosomes

Tps_subset_chr$chr <- factor(Tps_subset_chr$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", 
                                                           "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10",
                                                           "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))


nonoverlapWind_array_mean<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=mean)
table$array_mean<-nonoverlapWind_array_mean$`nonoverlapWind_chr$width`

ggplot(data=table, aes(x=LGs, y=array_mean)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(array_mean, digits = 3)), vjust=1.6, size=3)+ ylab("number of tandem repear arrays")+
  theme_classic()


Tps_subset_chr_motif<-aggregate(Tps_subset_chr$motif_length~Tps_subset_chr$chr, FUN=mean)
colnames(Tps_subset_chr_motif)<-c("LGs", "motif_length")
ggplot(data=Tps_subset_chr_motif, aes(x=LGs, y=motif_length)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(motif_length, digits = 3)), vjust=1.6, size=3)+ ylab("number of tandem repear arrays")+
  theme_classic()


Tps_subset_chr_copies<-aggregate(Tps_subset_chr$copy_nb~Tps_subset_chr$chr, FUN=mean)
colnames(Tps_subset_chr_copies)<-c("LGs", "copy_nb")
ggplot(data=Tps_subset_chr_copies, aes(x=LGs, y=copy_nb)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(copy_nb, digits = 3)), vjust=1.6, size=3)+ ylab("number of tandem repear arrays")+
  theme_classic()





ggplot(Tps_subset_chr, aes(x=chr, y=array_length)) + ylim(0,2500) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.05) + theme_classic()

ggplot(nonoverlapWind_chr, aes(x=seqnames, y=width)) + ylim(0,2500) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.05) + theme_classic()

ggplot(Tps_subset_chr, aes(x=chr, y=motif_length)) + ylim(0,100) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.05) + theme_classic()

ggplot(Tps_subset_chr, aes(x=chr, y=copy_nb)) + ylim(0,60) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.05) + theme_classic()




###############################################################
###Correlation recombination and TR load between chromosomes###
###############################################################

###Tandem repeat proportion between chromosomes Timema
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tps_subset<-subset(Tps, copy_nb >= 5)
Tps_subset$array_length<-(Tps_subset$end-Tps_subset$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tps_subset$chr,ranges=IRanges(start=Tps_subset$start,end=Tps_subset$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tps_LRv5b_scf1" | seqnames=="Tps_LRv5b_scf2" | seqnames=="Tps_LRv5b_scf3" | 
                             seqnames=="Tps_LRv5b_scf4" |seqnames=="Tps_LRv5b_scf5" | seqnames=="Tps_LRv5b_scf6" | seqnames=="Tps_LRv5b_scf7" | 
                             seqnames=="Tps_LRv5b_scf8" | seqnames=="Tps_LRv5b_scf9" | seqnames=="Tps_LRv5b_scf10" | seqnames=="Tps_LRv5b_scf11" | 
                             seqnames=="Tps_LRv5b_scf12")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)

#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")
Tps_LG_length<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")

nonoverlapWind_chr_sum$LG_length <- Tps_LG_length$length
nonoverlapWind_chr_sum$Prop_repeated_region<-nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`/nonoverlapWind_chr_sum$LG_length

colnames(nonoverlapWind_chr_sum) <- c("scaffold", "TR_length", "LG_length", "Prop_repeated_region")


###recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
  group_by(scaffold)%>%
  summarize(median_rho=median(rho))


merge<-merge(nonoverlapWind_chr_sum, rho_g, by="scaffold")
merge$scaffold<-factor(merge$scaffold, levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10", "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))
merge_sub<-subset(merge, scaffold!="Tps_LRv5b_scf3")

cor.test(merge$median_rho, merge$Prop_repeated_region, method = "spearman")#chrX included
summary(lm(merge$Prop_repeated_region~merge$median_rho))#chrX included

cor.test(merge_sub$median_rho, merge_sub$Prop_repeated_region, method = "spearman")#chrX excluded
summary(lm(merge_sub$Prop_repeated_region~merge_sub$median_rho))#chrX excluded

cor.test(merge$median_rho, merge$Prop_repeated_region, method = "spearman")#chrX included
summary(lm(merge$Prop_repeated_region~merge$LG_length))#chrX included

###TR distributions along chromosomes

##Tandem repeats
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tps_subset<-subset(Tps, copy_nb >= 5)
Tps_subset$array_length<-(Tps_subset$end-Tps_subset$start)+1


#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tps_subset$chr,ranges=IRanges(start=Tps_subset$start,end=Tps_subset$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tps_LRv5b_scf1" | seqnames=="Tps_LRv5b_scf2" | seqnames=="Tps_LRv5b_scf3" | 
                             seqnames=="Tps_LRv5b_scf4" |seqnames=="Tps_LRv5b_scf5" | seqnames=="Tps_LRv5b_scf6" | seqnames=="Tps_LRv5b_scf7" | 
                             seqnames=="Tps_LRv5b_scf8" | seqnames=="Tps_LRv5b_scf9" | seqnames=="Tps_LRv5b_scf10" | seqnames=="Tps_LRv5b_scf11" | 
                             seqnames=="Tps_LRv5b_scf12")
nonoverlapWind_chr$scaff_number<-as.numeric(as.character(gsub("^.*scf","", nonoverlapWind_chr$seqnames)))
nonoverlapWind_chr<-nonoverlapWind_chr %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values

write.table(nonoverlapWind_chr, file = "TR_nonoverlapWind_chr.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#The table "TR_nonoverlapWind_chr.txt" was generated to intersect non-overlapping TR positions with 10kb windows across every chromosomes.
#This was computed as: bedtools intersect -a TR_annotation_timema/minimal_rotations/TR_nonoverlapWind_chr.txt -b genomes/Tps_chm_size_mtDNAv350_w10000.bed > TR_annotation_timema/minimal_rotations/TR_nonoverlapWind_chr_intersection10kbWind.txt
#Then 10kb windows was associated with each intersected positions by computing: bedtools intersect -a TR_annotation_timema/minimal_rotations/TR_nonoverlapWind_chr_intersection10kbWind.txt -b genomes/Tps_chm_size_mtDNAv350_w10000.bed -wa -wb > TR_annotation_timema/minimal_rotations/TR_nonoverlapWind_chr_intersection10kbWind_10kbWindows.txt
#The final outcome was used to calculate the percentage of TR per 10kb windows (see below)


setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

TR_overlap_windows<-read.table("TR_nonoverlapWind_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(TR_overlap_windows)<-c("chromosome", "start_TR", "end_TR", "TR_array", "V5", "chrm_number", "chromosome_2", "start_wind", "end_wind")

TR_overlap_windows$TR_length<-(TR_overlap_windows$end_TR-TR_overlap_windows$start_TR)
TR_overlap_windows$start_wind_chr<-paste(TR_overlap_windows$start_wind, TR_overlap_windows$chromosome, sep="-")
TR_overlap_windows_sum<-aggregate(TR_overlap_windows$TR_length~TR_overlap_windows$start_wind_chr, FUN=sum)

library(stringr)
colnames(TR_overlap_windows_sum)[1]<-"start_chr"

#merge with all 10kb windows
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution")
windows<-read.table("Tps_chm_size_mtDNAv350_w250000.bed", header=F, sep="\t", quote="")
windows$start_chr<-paste(windows$V2, windows$V1, sep="-")

mergeTR_wind<-merge(windows, TR_overlap_windows_sum, by="start_chr", all.x=T)
mergeTR_wind[is.na(mergeTR_wind)] <- 0
colnames(mergeTR_wind)<-c("start_chr", "chromosome", "start", "end", "TR_length")
mergeTR_wind$percTR<-mergeTR_wind$TR_length/10000

chms<-subset(mergeTR_wind, chromosome=="Tps_LRv5b_scf1" | chromosome=="Tps_LRv5b_scf2" | chromosome=="Tps_LRv5b_scf3" | chromosome=="Tps_LRv5b_scf4" | chromosome=="Tps_LRv5b_scf5" 
             | chromosome=="Tps_LRv5b_scf6" | chromosome=="Tps_LRv5b_scf7" | chromosome=="Tps_LRv5b_scf8" | chromosome=="Tps_LRv5b_scf9" | chromosome=="Tps_LRv5b_scf10"
             | chromosome=="Tps_LRv5b_scf11" | chromosome=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$chromosome)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms1<-chms[order(chms$scaff_number, chms$start),]
chms1$cumul_start<-seq(0, by = 10000, length.out = nrow(chms1))
chms1$chrom<-ifelse(chms1$scaff_number=="3", "sex-chrom", "autosome")
chms2<-chms1[chms1$scaff_number=="3",]


ggplot(chms1, aes(x=start, y=percTR))+
  geom_point(size=0.5) +
  ylab("Percent of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")+theme(axis.text.y = element_text(angle = 90))

TR_array<-ggplot(chms1, aes(x=cumul_start, y=percTR, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms1$scaff_number)))/2))+
  geom_point(data=chms2, aes(x=chms2$cumul_start, y=chms2$percTR), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of tandem repeats")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")



##recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
       group_by(scaffold,start)%>%
       summarize(median_rho=median(rho))

#filtering extreme
rho_gf <- rho_g %>%
       group_by(scaffold)%>%
       mutate(quantile=quantile(median_rho,0.99))%>%
       filter(median_rho<quantile)

#keeping only big scaf
scaffolds <- paste("Tps_LRv5b_scf",1:12,sep="")

rho_gf$scaff_number<-as.numeric(as.character(gsub("^.*scf","", rho_gf$scaffold)))
rho_gf<-rho_gf %>% 
       filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms3<-rho_gf[order(rho_gf$scaff_number, rho_gf$start),]
chms3$cumul_start<-seq(0, by = 10000, length.out = nrow(chms3))
chms3$chrom<-ifelse(chms3$scaff_number=="3", "sex-chrom", "autosome")
chms4<-chms3[chms3$scaff_number=="3",]


recomb<-ggplot(chms3, aes(x=cumul_start, y=median_rho, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms3$scaff_number)))/2))+
  geom_point(data=chms4, aes(x=chms4$cumul_start, y=chms4$median_rho), stat="identity", colour="maroon", size=0.5) +
  ylab("rho")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+ylim(0,0.0015)



##GC content
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("Tps_chm_size_mtDNAv350_w10000_GC_content.txt", header=FALSE)


#Plot

chms<-subset(data_GC, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms3<-chms[order(chms$scaff_number, chms$V2),]
chms3$cumul_start<-seq(0, by = 10000, length.out = nrow(chms3))
chms3$chrom<-ifelse(chms3$scaff_number=="3", "sex-chrom", "autosome")
chms4<-chms3[chms3$scaff_number=="3",]

GC_content<-ggplot(chms3, aes(x=cumul_start, y=V5, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms3$scaff_number)))/2))+
  geom_point(data=chms4, aes(x=chms4$cumul_start, y=chms4$V5), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of GC")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")


##GC content without TR
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/GC_content")

data_GC<-read.table("Tps_chm_size_mtDNAv350_w10000_GC_content_withoutTR.txt", header=FALSE)


#Plot

chms<-subset(data_GC, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms5<-chms[order(chms$scaff_number, chms$V2),]
chms5$cumul_start<-seq(0, by = 10000, length.out = nrow(chms5))
chms5$chrom<-ifelse(chms5$scaff_number=="3", "sex-chrom", "autosome")
chms6<-chms5[chms5$scaff_number=="3",]

GC_content_noTR<-ggplot(chms5, aes(x=cumul_start, y=V5, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms5$scaff_number)))/2))+
  geom_point(data=chms6, aes(x=chms6$cumul_start, y=chms6$V5), stat="identity", colour="maroon", size=0.5) +
  ylab("Percent of GC without TR")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white")






##ChIP cenh3 testes

setwd("/Users/wtoubian/Desktop/Chip-seq/species-comparison/tracks/Tps_testes_cenh3_R1")
cenh3<-read.table("Tps_testes_cenh3_1_coverage_DR_10kb.txt",
                  header=FALSE)

setwd("/Users/wtoubian/Desktop/Chip-seq/species-comparison/tracks/Tps_testes_input_R1")
input1<-read.table("Tps_testes_input1_coverage_DR_10kb.txt",
                   header=FALSE)

data1<-cbind(cenh3,input1[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1$coverage_cenh3_norm<-data1$coverage_cenh3/66295120 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratiocenH3_norm<-data1$coverage_cenh3_norm/data1$coverage_input_norm
data1$log2cenH3_norm<-log2(data1$ratiocenH3_norm)

data2<-subset(data1, coverage_cenh3>1 | coverage_input>1)
data2$start_chr<-paste(data2$start, data2$Scaffold_name, sep="-")

mergeCenh3_wind<-merge(windows, data2, by="start_chr", all.x=T)
mergeCenh3_wind[is.na(mergeCenh3_wind)] <- 0

##Plots

chms<-subset(mergeCenh3_wind, V1=="Tps_LRv5b_scf1" | V1=="Tps_LRv5b_scf2" | V1=="Tps_LRv5b_scf3" | V1=="Tps_LRv5b_scf4" | V1=="Tps_LRv5b_scf5" 
             | V1=="Tps_LRv5b_scf6" | V1=="Tps_LRv5b_scf7" | V1=="Tps_LRv5b_scf8" | V1=="Tps_LRv5b_scf9" | V1=="Tps_LRv5b_scf10"
             | V1=="Tps_LRv5b_scf11" | V1=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$V1)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values


chms9<-chms[order(chms$scaff_number, chms$start),]
chms9$cumul_start<-seq(0, by = 10000, length.out = nrow(chms9))
chms9$chrom<-ifelse(chms9$scaff_number=="3", "sex-chrom", "autosome")
chms10<-chms9[chms9$scaff_number=="3",]


cenh3_testes<-ggplot(chms9, aes(x=cumul_start, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_color_manual(values = rep(c("black", "grey"), length(levels(as.factor(chms9$scaff_number)))/2))+
  scale_y_continuous(limits=c(-1,7)) +
  geom_point(data=chms10, aes(x=chms10$cumul_start, y=chms10$log2cenH3_norm), stat="identity", colour="maroon", size=0.5) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept=c(2,0,0), linetype=c("dashed","solid", "dashed"), color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white", method="gam", formula = y~s(x, bs="cs", k=20))


###Extract coordinates smooth line intersecting with X-axis
library(mgcv)
chr<-subset(chms9, V1=="Tps_LRv5b_scf12") #example for one chromosome
gam<-gam(log2cenH3_norm ~ s(V2, bs = "cs", k=30), data=chr)
chr$predict<-predict(gam, newdata = chr)
chr$predict_abs<-abs(chr$predict)

ggplot(chr, aes(x=V2, y=log2cenH3_norm, color = as.factor(scaff_number)))+
  geom_point(size=0.5) +
  scale_y_continuous(limits=c(-1,7)) +
  ylab("log2(CenH3/Input)")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept=c(2,0,0), linetype=c("dashed","solid", "dashed"), color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))+geom_vline(xintercept =c(46640000, 60030000))


#All Plots
library(gridExtra)
grid.arrange(cenh3_testes, TR_array, h3k9_testes, nrow=3, ncol=1)


library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    cenh3_testes + theme(legend.position = "none")
    , TR_array
    , h3k9_testes + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(h3k9_testes)
    , ncol =1)
  , rel_widths = c(12,3)
)



plot_grid(
  plot_grid(
    GC_content + theme(legend.position = "none")
    , TR_array
    , GC_content_noTR + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(GC_content)
    , ggplot()
    , get_legend(GC_content_noTR)
    , ncol =1)
  , rel_widths = c(12,3)
)


plot_grid(
  plot_grid(
    cenh3_testes + theme(legend.position = "none")
    , TR_array
    , recomb + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(recomb)
    , ncol =1)
  , rel_widths = c(12,0)
)


##################################
# Correlation TR - recombination #
##################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

TR_overlap_windows<-read.table("TR_nonoverlapWind_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(TR_overlap_windows)<-c("chromosome", "start_TR", "end_TR", "TR_array", "V5", "chrm_number", "chromosome_2", "start_wind", "end_wind")

TR_overlap_windows$TR_length<-(TR_overlap_windows$end_TR-TR_overlap_windows$start_TR)
TR_overlap_windows$start_wind_chr<-paste(TR_overlap_windows$start_wind, TR_overlap_windows$chromosome, sep="-")
TR_overlap_windows_sum<-aggregate(TR_overlap_windows$TR_length~TR_overlap_windows$start_wind_chr, FUN=sum)

library(stringr)
colnames(TR_overlap_windows_sum)[1]<-"start_chr"

#merge with all 250kb windows
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution")
windows<-read.table("Tps_chm_size_mtDNAv350_w250000.bed", header=F, sep="\t", quote="")
windows$start_chr<-paste(windows$V2, windows$V1, sep="-")

mergeTR_wind<-merge(windows, TR_overlap_windows_sum, by="start_chr", all.x=T)
mergeTR_wind[is.na(mergeTR_wind)] <- 0
colnames(mergeTR_wind)<-c("start_chr", "chromosome", "start", "end", "TR_length")
mergeTR_wind$percTR<-mergeTR_wind$TR_length/10000

chms<-subset(mergeTR_wind, chromosome=="Tps_LRv5b_scf1" | chromosome=="Tps_LRv5b_scf2" | chromosome=="Tps_LRv5b_scf3" | chromosome=="Tps_LRv5b_scf4" | chromosome=="Tps_LRv5b_scf5" 
             | chromosome=="Tps_LRv5b_scf6" | chromosome=="Tps_LRv5b_scf7" | chromosome=="Tps_LRv5b_scf8" | chromosome=="Tps_LRv5b_scf9" | chromosome=="Tps_LRv5b_scf10"
             | chromosome=="Tps_LRv5b_scf11" | chromosome=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$chromosome)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values



##recombination
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/recombination_map")

rho<-read.table("RHO_Tpop_mean_5_runs_chr_intersection100kbWind.txt", header=F, sep="\t", quote="")
colnames(rho) <- c("scaffold", "left_SNP", "right_SNP", "rho", "scaffold_bis", "start", "end")

#Median per window
rho_g <- rho %>%
  group_by(scaffold,start)%>%
  summarize(median_rho=median(rho))

rho_g$start_chr<-paste(rho_g$start, rho_g$scaffold, sep="-")

merge<-merge(chms, rho_g, by="start_chr")
merge1<-merge[order(merge$scaff_number, merge$start.x),]
merge1$cumul_start<-seq(0, by = 250000, length.out = nrow(merge1))

merge_sub<-subset(merge1, chromosome=="Tps_LRv5b_scf12")

#filtering extreme
merge_sub2 <- merge_sub %>%
  group_by(scaffold)%>%
  mutate(quantile=quantile(median_rho,0.99))%>%
  filter(median_rho<quantile)


#summary(lm(merge_sub2$percTR~merge_sub2$median_rho))
#cor.test(merge_sub2$median_rho, merge_sub2$percTR, method = "spearman")


###Extract and compare values of the Generalized Additive Models
library(mgcv)
gamrho<-gam(median_rho ~ s(cumul_start, bs = "cs", k=20), data=merge_sub2)#per chromosomes
gamTR<-gam(percTR ~ s(cumul_start, bs = "cs", k=20), data=merge_sub2)#per chromosomes

predicted_values_gamrho <- data.frame(cumul_start = merge_sub2$cumul_start, 
                                       fitted = predict(gamrho))
predicted_values_gamTR <- data.frame(cumul_start = merge_sub2$cumul_start, 
                                     fitted = predict(gamTR))



ggplot(merge_sub2, aes(x = cumul_start, y = median_rho)) +
  geom_point() +
  geom_line(data = predicted_values_gamrho, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))

ggplot(merge_sub2, aes(x = cumul_start, y = percTR)) +
  geom_point() +
  geom_line(data = predicted_values_gamTR, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))


summary(lm(predicted_values_gamTR$fitted~predicted_values_gamrho$fitted))
cor.test(predicted_values_gamrho$fitted, predicted_values_gamTR$fitted, method = "spearman")




####################################
# Correlation TR - heterochromatin #
####################################

##windows
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution")
windows<-read.table("Tps_chm_size_mtDNAv350_w250000.bed", header=F, sep="\t", quote="")
windows$start_chr<-paste(windows$V2, windows$V1, sep="-")


##TRs
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

TR_overlap_windows<-read.table("TR_nonoverlapWind_chr_intersection250kbWind.txt", header=F, sep="\t", quote="")
colnames(TR_overlap_windows)<-c("chromosome", "start_TR", "end_TR", "TR_array", "V5", "chrm_number", "chromosome_2", "start_wind", "end_wind")

TR_overlap_windows$TR_length<-(TR_overlap_windows$end_TR-TR_overlap_windows$start_TR)
TR_overlap_windows$start_wind_chr<-paste(TR_overlap_windows$start_wind, TR_overlap_windows$chromosome, sep="-")
TR_overlap_windows_sum<-aggregate(TR_overlap_windows$TR_length~TR_overlap_windows$start_wind_chr, FUN=sum)

library(stringr)
colnames(TR_overlap_windows_sum)[1]<-"start_chr"


#merge with all 10kb windows
mergeTR_wind<-merge(windows, TR_overlap_windows_sum, by="start_chr", all.x=T)
mergeTR_wind[is.na(mergeTR_wind)] <- 0
colnames(mergeTR_wind)<-c("start_chr", "chromosome", "start", "end", "TR_length")
mergeTR_wind$percTR<-mergeTR_wind$TR_length/250000



##h3k9
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/h3k9")

h3k9<-read.table("Tps_male_organs_h3k9A_coverage_DR_250kb.txt",
                 header=FALSE)

input<-read.table("Tps_male_organs_input3_coverage_DR_250kb.txt",
                  header=FALSE)

data1<-cbind(h3k9,input[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_h3k9", "coverage_input")

data1$coverage_h3k9_norm<-data1$coverage_h3k9/62681606 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/59528909 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratioH3k9_norm<-data1$coverage_h3k9_norm/data1$coverage_input_norm
data1$log2H3k9_norm<-log2(data1$ratioH3k9_norm)

data2<-subset(data1, coverage_h3k9>1 | coverage_input>1)
data2$start_chr<-paste(data2$start, data2$Scaffold_name, sep="-")

mergeH3k9_wind<-merge(windows, data2, by="start_chr", all.x=T)
mergeH3k9_wind[is.na(mergeH3k9_wind)] <- 0

merge<-merge(mergeTR_wind, mergeH3k9_wind, by="start_chr")


chms<-subset(merge, chromosome=="Tps_LRv5b_scf1" | chromosome=="Tps_LRv5b_scf2" | chromosome=="Tps_LRv5b_scf3" | chromosome=="Tps_LRv5b_scf4" | chromosome=="Tps_LRv5b_scf5" 
             | chromosome=="Tps_LRv5b_scf6" | chromosome=="Tps_LRv5b_scf7" | chromosome=="Tps_LRv5b_scf8" | chromosome=="Tps_LRv5b_scf9" | chromosome=="Tps_LRv5b_scf10"
             | chromosome=="Tps_LRv5b_scf11" | chromosome=="Tps_LRv5b_scf12")
chms$scaff_number<-as.numeric(as.character(gsub("^.*scf","", chms$chromosome)))
chms<-chms %>% 
  filter_all(all_vars(!is.infinite(.))) #filters out infinite values

chms1<-chms[order(chms$scaff_number, chms$start.x),]
chms1$cumul_start<-seq(0, by = 250000, length.out = nrow(chms1))

chms1_sub<-subset(chms1, chromosome=="Tps_LRv5b_scf10")

#summary(lm(chms1_sub$percTR~chms1_sub$log2H3k9_norm))
#cor.test(chms1$log2H3k9_norm, chms1$percTR, method = "spearman")
#cor.test(chms1_sub$log2H3k9_norm, chms1_sub$percTR, method = "spearman")



###Extract and compare values of the Generalized Additive Models
library(mgcv)
gamh3k9<-gam(log2H3k9_norm ~ s(cumul_start, bs = "cs", k=20), data=chms1_sub)#per chromosomes
gamTR<-gam(percTR ~ s(cumul_start, bs = "cs", k=20), data=chms1_sub)#per chromosomes

predicted_values_gamh3k9 <- data.frame(cumul_start = chms1_sub$cumul_start, 
                                   fitted = predict(gamh3k9))
predicted_values_gamTR <- data.frame(cumul_start = chms1_sub$cumul_start, 
                                       fitted = predict(gamTR))


ggplot(chms1_sub, aes(x = cumul_start, y = log2H3k9_norm)) +
  geom_point() +
  geom_line(data = predicted_values_gamh3k9, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))

ggplot(chms1_sub, aes(x = cumul_start, y = percTR)) +
  geom_point() +
  geom_line(data = predicted_values_gamTR, aes(x = cumul_start, y = fitted), color = "blue") +
  geom_smooth(colour = "red", method="gam", formula = y~s(x, bs="cs", k=30))


summary(lm(predicted_values_gamTR$fitted~predicted_values_gamh3k9$fitted))
cor.test(predicted_values_gamh3k9$fitted, predicted_values_gamTR$fitted, method = "spearman")




############################################
# TR proportion comparison between species #
############################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

##Tps
Tps<-read.table("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tps_subset<-subset(Tps, copy_nb >= 5)
Tps_subset$genome<-1339836934
Tps_subset2<-subset(Tps_subset, motif_length <= 10)
#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tps_subset$chr,ranges=IRanges(start=Tps_subset$start,end=Tps_subset$end))
nonoverlapWind_Tps<-reduce(myranges)
nonoverlapWind_Tps<-as.data.frame(nonoverlapWind_Tps)
nonoverlapWind_Tps$species<-"Tps"



##Tdi
Tdi<-read.table("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt", header=T, sep="\t", quote="")
Tdi_subset<-subset(Tdi, copy_nb >= 5)
Tdi_subset$genome<-1292596824

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tdi_subset$chr,ranges=IRanges(start=Tdi_subset$start,end=Tdi_subset$end))
nonoverlapWind_Tdi<-reduce(myranges)
nonoverlapWind_Tdi<-as.data.frame(nonoverlapWind_Tdi)
nonoverlapWind_Tdi$species<-"Tdi"



##Tcm
Tcm<-read.table("Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tcm_subset<-subset(Tcm, copy_nb >= 5)
Tcm_subset$genome<-1296216039
Tcm_subset2<-subset(Tcm_subset, motif_length <= 10)
#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tcm_subset$chr,ranges=IRanges(start=Tcm_subset$start,end=Tcm_subset$end))
nonoverlapWind_Tcm<-reduce(myranges)
nonoverlapWind_Tcm<-as.data.frame(nonoverlapWind_Tcm)
nonoverlapWind_Tcm$species<-"Tcm"


##Tsi
Tsi<-read.table("Tsi_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tsi_subset<-subset(Tsi, copy_nb >= 5)
Tsi_subset$genome<-1324360681

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tsi_subset$chr,ranges=IRanges(start=Tsi_subset$start,end=Tsi_subset$end))
nonoverlapWind_Tsi<-reduce(myranges)
nonoverlapWind_Tsi<-as.data.frame(nonoverlapWind_Tsi)
nonoverlapWind_Tsi$species<-"Tsi"


##Tce
Tce<-read.table("Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tce_subset<-subset(Tce, copy_nb >= 5)
Tce_subset$genome<-1211435093#1313940951
Tce_subset2<-subset(Tce_subset, motif_length <= 10)
#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tce_subset$chr,ranges=IRanges(start=Tce_subset$start,end=Tce_subset$end))
nonoverlapWind_Tce<-reduce(myranges)
nonoverlapWind_Tce<-as.data.frame(nonoverlapWind_Tce)
nonoverlapWind_Tce$species<-"Tce"


##Tms
Tms<-read.table("Tms_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tms_subset<-subset(Tms, copy_nb >= 5)
Tms_subset$genome<-1168105536

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tms_subset$chr,ranges=IRanges(start=Tms_subset$start,end=Tms_subset$end))
nonoverlapWind_Tms<-reduce(myranges)
nonoverlapWind_Tms<-as.data.frame(nonoverlapWind_Tms)
nonoverlapWind_Tms$species<-"Tms"



##Tpa
Tpa<-read.table("Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tpa_subset<-subset(Tpa, copy_nb >= 5)
Tpa_subset$genome<-1145701453
Tpa_subset2<-subset(Tpa_subset, motif_length <= 10)
#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tpa_subset$chr,ranges=IRanges(start=Tpa_subset$start,end=Tpa_subset$end))
nonoverlapWind_Tpa<-reduce(myranges)
nonoverlapWind_Tpa<-as.data.frame(nonoverlapWind_Tpa)
nonoverlapWind_Tpa$species<-"Tpa"


##TgeA
TgeA<-read.table("TgeA_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
TgeA_subset<-subset(TgeA, copy_nb >= 5)
TgeA_subset$genome<-1121025499

#Remove overlaps between repeats
myranges<-GRanges(seqnames=TgeA_subset$chr,ranges=IRanges(start=TgeA_subset$start,end=TgeA_subset$end))
nonoverlapWind_TgeA<-reduce(myranges)
nonoverlapWind_TgeA<-as.data.frame(nonoverlapWind_TgeA)
nonoverlapWind_TgeA$species<-"TgeA"


##TgeH20
TgeH20<-read.table("TgeH20_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
TgeH20_subset<-subset(TgeH20, copy_nb >= 5)
TgeH20_subset$genome<-1117223061

#Remove overlaps between repeats
myranges<-GRanges(seqnames=TgeH20_subset$chr,ranges=IRanges(start=TgeH20_subset$start,end=TgeH20_subset$end))
nonoverlapWind_TgeH20<-reduce(myranges)
nonoverlapWind_TgeH20<-as.data.frame(nonoverlapWind_TgeH20)
nonoverlapWind_TgeH20$species<-"TgeH20"


##Tbi
Tbi<-read.table("Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tbi_subset<-subset(Tbi, copy_nb >= 5)
Tbi_subset$genome<-1208310017
Tbi_subset2<-subset(Tbi_subset, motif_length <= 10)
#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tbi_subset$chr,ranges=IRanges(start=Tbi_subset$start,end=Tbi_subset$end))
nonoverlapWind_Tbi<-reduce(myranges)
nonoverlapWind_Tbi<-as.data.frame(nonoverlapWind_Tbi)
nonoverlapWind_Tbi$species<-"Tbi"


##Tte
Tte<-read.table("Tte_LRv1_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotations.txt", header=T, sep="\t", quote="")
Tte_subset<-subset(Tte, copy_nb >= 5)
Tte_subset$genome<-1270759584

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tte_subset$chr,ranges=IRanges(start=Tte_subset$start,end=Tte_subset$end))
nonoverlapWind_Tte<-reduce(myranges)
nonoverlapWind_Tte<-as.data.frame(nonoverlapWind_Tte)
nonoverlapWind_Tte$species<-"Tte"


###Store datasets (example for one species)
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")

write.table(TgeA_subset, file = "TgeA_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TgeA_subset[c(1,2,3,6)], file = "TgeA_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)



###
data2<-rbind(do.call("rbind", list(nonoverlapWind_Tps, nonoverlapWind_Tdi, nonoverlapWind_Tcm, nonoverlapWind_Tsi, nonoverlapWind_Tce, nonoverlapWind_Tms, 
                                   nonoverlapWind_Tpa, nonoverlapWind_TgeA, nonoverlapWind_TgeH20, nonoverlapWind_Tbi, nonoverlapWind_Tte)))

data2$species <- factor(data2$species,levels = c("Tps", "Tdi", "Tcm", "Tsi", "Tce", "Tms", "Tpa", "TgeA", "TgeH20", "Tbi", "Tte"))
data2$count<-1
data2_TR_sum<-aggregate(data2$count~data2$species, FUN=sum)
colnames(data2_TR_sum)<-c("species", "TR_number")
data2_TR_sum_sex<-subset(data2_TR_sum, species=="Tps"| species=="Tcm"| species=="Tce"| species=="Tpa"| species=="Tbi")

data2_array_mean<-aggregate(data2$width~data2$species, FUN=mean)
colnames(data2_array_mean)<-c("species", "TR_array_mean")
data2_array_mean_sex<-subset(data2_array_mean, species=="Tps"| species=="Tcm"| species=="Tce"| species=="Tpa"| species=="Tbi")

ggplot(data=data2_TR_sum_sex, aes(x=species, y=TR_number)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(TR_number, digits = 3)), vjust=1.6, size=3)+ ylab("number of tandem repeat arrays")+
  theme_classic()

ggplot(data=data2_array_mean_sex, aes(x=species, y=TR_array_mean)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(TR_array_mean, digits = 3)), vjust=1.6, size=3)+ ylab("average tandem repeat array lengths")+
  theme_classic()









data3<-rbind(do.call("rbind", list(Tps_subset, Tdi_subset, Tcm_subset, Tsi_subset, Tce_subset, Tms_subset, 
                                   Tpa_subset, TgeA_subset, TgeH20_subset, Tbi_subset, Tte_subset)))
data3$species <- factor(data3$species,levels = c("Tps", "Tdi", "Tcm", "Tsi", "Tce", "Tms", "Tpa", "TgeA", "TgeH20", "Tbi", "Tte"))

data3_sex<-subset(data3, species=="Tps"| species=="Tcm"| species=="Tce"| species=="Tpa"| species=="Tbi")


library(ggridges)
ggplot(
  data3_sex, 
  aes(x = `motif_length`, y = `species`, fill = stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "motif_length", option = "C") +
  labs(title = 'Density plot motif length') + xlim(0,200)


data3_motif_length<-aggregate(data3$motif_length~data3$species, FUN=mean)
colnames(data3_motif_length)<-c("species", "motif_length_mean")
data3_motif_length_sex<-subset(data3_motif_length, species=="Tps"| species=="Tcm"| species=="Tce"| species=="Tpa"| species=="Tbi")

ggplot(data=data3_motif_length_sex, aes(x=species, y=motif_length_mean)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(motif_length_mean, digits = 3)), vjust=1.6, size=3)+ ylab("average motif lengths")+
  theme_classic()


ggplot(
  data2, 
  aes(x = `copy_nb`, y = `species`, fill = stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "motif_length", option = "C") +
  labs(title = 'Density plot motif length') + xlim(0,100)



#################################
## Tbi contigs to Tpa assembly ##
#################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/Tbi-Tte_contigs_to_Tpa_assembly")
data<-read.table("TBI_TPA.50.csv", header=F, sep="", quote="")
autosomes<-subset(data, V2=="Tpa_LRv5a_scf2" | V2=="Tpa_LRv5a_scf3" | V2=="Tpa_LRv5a_scf4" | V2=="Tpa_LRv5a_scf5" | V2=="Tpa_LRv5a_scf6" | V2=="Tpa_LRv5a_scf7" | V2=="Tpa_LRv5a_scf8" | V2=="Tpa_LRv5a_scf8" | V2=="Tpa_LRv5a_scf10" | V2=="Tpa_LRv5a_scf11" | V2=="Tpa_LRv5a_scf12" | V2=="Tpa_LRv5a_scf13" | V2=="Tpa_LRv5a_scf14" )
chrX<-subset(data, V2=="Tpa_LRv5a_scf1")

data<-read.table("TTE_TPA.50.csv", header=F, sep="", quote="")
autosomes<-subset(data, V2=="Tpa_LRv5a_scf2" | V2=="Tpa_LRv5a_scf3" | V2=="Tpa_LRv5a_scf4" | V2=="Tpa_LRv5a_scf5" | V2=="Tpa_LRv5a_scf6" | V2=="Tpa_LRv5a_scf7" | V2=="Tpa_LRv5a_scf8" | V2=="Tpa_LRv5a_scf8" | V2=="Tpa_LRv5a_scf10" | V2=="Tpa_LRv5a_scf11" | V2=="Tpa_LRv5a_scf12" | V2=="Tpa_LRv5a_scf13" | V2=="Tpa_LRv5a_scf14" )

chrX<-subset(data, V2=="Tpa_LRv5a_scf1")







data2_sum<-aggregate(data2$copy_nb~data2$species, FUN=sum)
data2_sum<-aggregate(data2$region_size~data2$species, FUN=sum)
data2_sum<-aggregate(data2$propTR~data2$species, FUN=sum)

colnames(data2_sum)<-c("species", "copy_nb")

ggplot(data=data2_sum, aes(x=species, y=copy_nb)) +
  geom_bar(stat="identity", color="red", fill="white")+
  geom_text(aes(label=round(copy_nb, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()


#Count nb of monomers

Tps_subset_chr<- subset(Tps_subset, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | 
                          chr=="Tps_LRv5b_scf4" |chr=="Tps_LRv5b_scf5" | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | 
                          chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10" | chr=="Tps_LRv5b_scf11" | 
                          chr=="Tps_LRv5b_scf12")

Tps_subset_chr$count<-1
Tps_subset_chr_sum<-aggregate(Tps_subset_chr$count~Tps_subset_chr$chr, FUN=sum)

Tps_LG_length<-read.table("Tps_scaffold_length.txt", header=T, sep="", quote = "")

table=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
table$LG_length <- Tps_LG_length$length
table$motif_count<-Tps_subset_chr_sum$`Tps_subset_chr$count`


ggplot(data=table, aes(x=LGs, y=motif_count)) +
  geom_bar(stat="identity", color="red", fill="white")+
  geom_text(aes(label=round(motif_count, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()


#motif length across chromosomes

Tps_subset_chr$chr <- factor(Tps_subset_chr$chr,levels = c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", 
                                                  "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10",
                                                  "Tps_LRv5b_scf11", "Tps_LRv5b_scf12"))


ggplot(Tps_subset_chr, aes(x=chr, y=motif_length)) + ylim(0,50) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.05)

ggplot(Tps_subset_chr, aes(x=chr, y=copy_nb)) + ylim(0,50) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.05)


ggplot(
  Tps_subset_chr, 
  aes(x = `motif_length`, y = `chr`, fill = stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "motif_length", option = "C") +
  labs(title = 'Density plot motif length') + xlim(0,200)


ggplot(
  Tps_subset_chr, 
  aes(x = `copy_nb`, y = `chr`, fill = stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "motif_length", option = "C") +
  labs(title = 'Density plot motif length') + xlim(0,200)

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tps_subset$chr,ranges=IRanges(start=Tps_subset$start,end=Tps_subset$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tps_LRv5b_scf1" | seqnames=="Tps_LRv5b_scf2" | seqnames=="Tps_LRv5b_scf3" | 
                   seqnames=="Tps_LRv5b_scf4" |seqnames=="Tps_LRv5b_scf5" | seqnames=="Tps_LRv5b_scf6" | seqnames=="Tps_LRv5b_scf7" | 
                   seqnames=="Tps_LRv5b_scf8" | seqnames=="Tps_LRv5b_scf9" | seqnames=="Tps_LRv5b_scf10" | seqnames=="Tps_LRv5b_scf11" | 
                   seqnames=="Tps_LRv5b_scf12")

nonoverlapWind_chr_sum<-aggregate(data_chr$width~data_chr$seqnames, FUN=sum)


#Calculate proportion
table$TR_length <- nonoverlapWind_chr_sum$`data_chr$width`
table$Prop_repeated_region<-table$TR_length/table$LG_length


ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="red", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()



table$nb_monomer <- lapply(df.list, function(x) nrow(x)) #count number of repeat per scaffold

res=lapply(df.list, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table$Repeated_region_size <- lapply(df.list, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table$Counts_median <- lapply(df.list, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table$Counts_sum <- lapply(df.list, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table$nb_monomer=as.numeric(table$nb_monomer)
table$Repeated_region_size=as.numeric(table$Repeated_region_size)
table$Counts_sum=as.numeric(table$Counts_sum)
table$Counts_median=as.numeric(table$Counts_median)
table$Prop_repeated_region= table$Repeated_region_size/table$LG_length
table$LGs <- factor(table$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
table$species<-"Tps"
head(table)
summary(table)





#Tps
Tdi<-read.table("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_minimal_rotation.txt", header=T, sep="\t", quote="")
Tdi_LG_length<-read.table("Tdi_scaffold_length.txt", header=T, sep="", quote = "")
Tdi_subset<-subset(Tdi, copy_nb >= 5)

#Remove overlaps between repeats
library(GenomicRanges)
myranges<-GRanges(seqnames=Tdi_subset$chr,ranges=IRanges(start=Tdi_subset$start,end=Tdi_subset$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)
nonoverlapWind$Repeated_region_size<- nonoverlapWind$end-nonoverlapWind$start

Tdi_LG1<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf1")
Tdi_LG2<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf2")
Tdi_LG3<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf3")
Tdi_LG4<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf4")
Tdi_LG5<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf5")
Tdi_LG6<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf6")
Tdi_LG7<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf7")
Tdi_LG8<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf8")
Tdi_LG9<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf9")
Tdi_LG10<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf10")
Tdi_LG11<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf11")
Tdi_LG12<- subset(nonoverlapWind, seqnames=="Tdi_LRv5a_scf12")

df.list2= list(Tdi_LG1, Tdi_LG2, Tdi_LG3, Tdi_LG4, Tdi_LG5, Tdi_LG6, Tdi_LG7, Tdi_LG8, Tdi_LG9, Tdi_LG10, Tdi_LG11, Tdi_LG12)

table2=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
colnames (table2) [1] <- "LGs"
table2$LG_length <- Tdi_LG_length$length

table2$Repeated_region_size <- lapply(df.list, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold


table2$nb_monomer <- lapply(df.list2, function(x) nrow(x)) #count number of repeat per scaffold

res2=lapply(df.list2, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table2$Repeated_region_size <- lapply(res2, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table2$Counts_median <- lapply(df.list2, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table2$Counts_sum <- lapply(df.list2, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table2$nb_monomer=as.numeric(table2$nb_monomer)
table2$Repeated_region_size=as.numeric(table2$Repeated_region_size)
table2$Counts_sum=as.numeric(table2$Counts_sum)
table2$Counts_median=as.numeric(table2$Counts_median)
table2$Prop_repeated_region= table2$Repeated_region_size/table2$LG_length
table2$LGs <- factor(table2$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
table2$species<-"Tdi"
head(table2)
summary(table2)

ggplot(data=table2, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="blue", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()


Merge1=rbind(table, table2)
Merge1$species <- factor(Merge1$species,levels = c("Tps", "Tdi"))


ggplot(data=Merge1, aes(x=LGs, y=Prop_repeated_region, fill=species)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()


Tps$species<-"Tps"
Tdi$species<-"Tdi"
bindTpsTdi<-rbind(Tps, Tdi)

sharedTps<-merge(Tps, Tdi[,c(1,6)], by="X_remonomsqes", all.x = F, all.y = F)
sharedTdi<-merge(Tdi, Tps[,c(1,6)], by="X_remonomsqes", all.x = F, all.y = F)
bindsharedTpsTdi<-rbind(sharedTps, sharedTdi)

ggplot(bindsharedTpsTdi, aes(x=species, y=Counts_s)) + 
  ylim(0,500)+
  geom_violin()



## Tcm vs Tsi
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tcm")

Tcm<-read.table("Tcm_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")

Tcm_LG1.1<-read.table("Tcm_LRv5a_scf1.1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG1.2<-read.table("Tcm_LRv5a_scf1.2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG1<-rbind(Tcm_LG1.1, Tcm_LG1.2)
Tcm_LG2<-read.table("Tcm_LRv5a_scf2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG3<-read.table("Tcm_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG4<-read.table("Tcm_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG5.1<-read.table("Tcm_LRv5a_scf5.1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG5.2<-read.table("Tcm_LRv5a_scf5.2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG5<-rbind(Tcm_LG5.1, Tcm_LG5.2)
Tcm_LG6.1<-read.table("Tcm_LRv5a_scf6.1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG6.2<-read.table("Tcm_LRv5a_scf6.2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG6<-rbind(Tcm_LG6.1, Tcm_LG6.2)
Tcm_LG7<-read.table("Tcm_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG8.1<-read.table("Tcm_LRv5a_scf8.1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG8.2<-read.table("Tcm_LRv5a_scf8.2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG8.3<-read.table("Tcm_LRv5a_scf8.3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG8<-rbind(Tcm_LG8.1, Tcm_LG8.2, Tcm_LG8.3)
Tcm_LG9<-read.table("Tcm_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG10<-read.table("Tcm_LRv5a_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG11.1<-read.table("Tcm_LRv5a_scf11.1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG11.2<-read.table("Tcm_LRv5a_scf11.2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG11<-rbind(Tcm_LG11.1, Tcm_LG11.2)
Tcm_LG12<-read.table("Tcm_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tcm_LG_length<-read.table("Tcm_scaffold_length.txt", header=T, sep="", quote = "")
Tcm_LG_length<-Tcm_LG_length[1:12,]
#Create summary table of LG length, number of distinct repeats and overall repeated region size per scaffold

df.list= list(Tcm_LG1, Tcm_LG2, Tcm_LG3, Tcm_LG4, Tcm_LG5, Tcm_LG6, Tcm_LG7, Tcm_LG8, Tcm_LG9, Tcm_LG10, Tcm_LG11, Tcm_LG12)

table=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
colnames (table) [1] <- "LGs"
table$LG_length <- Tcm_LG_length$length

table$nb_monomer <- lapply(df.list, function(x) nrow(x)) #count number of repeat per scaffold

res=lapply(df.list, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table$Repeated_region_size <- lapply(res, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table$Counts_median <- lapply(df.list, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table$Counts_sum <- lapply(df.list, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table$nb_monomer=as.numeric(table$nb_monomer)
table$Repeated_region_size=as.numeric(table$Repeated_region_size)
table$Counts_sum=as.numeric(table$Counts_sum)
table$Counts_median=as.numeric(table$Counts_median)
table$Prop_repeated_region= table$Repeated_region_size/table$LG_length
table$LGs <- factor(table$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
table$species<-"Tcm"
head(table)
summary(table)

ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="red", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3.5)+ ylab("Proportion tandem repeats")+
  theme_classic()

ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="red", fill="white")+
  ylab("Proportion tandem repeats")+
  theme_classic()

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")

Tsi<-read.table("Tsi_LRv5b.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")

Tsi_LG1<-read.table("Tsi_LRv5b_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG2<-read.table("Tsi_LRv5b_scf2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG3<-read.table("Tsi_LRv5b_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG4<-read.table("Tsi_LRv5b_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG5<-read.table("Tsi_LRv5b_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG6<-read.table("Tsi_LRv5b_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG7<-read.table("Tsi_LRv5b_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG8<-read.table("Tsi_LRv5b_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG9<-read.table("Tsi_LRv5b_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG10<-read.table("Tsi_LRv5b_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG11<-read.table("Tsi_LRv5b_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG12<-read.table("Tsi_LRv5b_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tsi_LG_length<-read.table("Tsi_scaffold_length.txt", header=T, sep="", quote = "")

df.list2= list(Tsi_LG1, Tsi_LG2, Tsi_LG3, Tsi_LG4, Tsi_LG5, Tsi_LG6, Tsi_LG7, Tsi_LG8, Tsi_LG9, Tsi_LG10, Tsi_LG11, Tsi_LG12)

table2=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
colnames (table2) [1] <- "LGs"
table2$LG_length <- Tsi_LG_length$length

table2$nb_monomer <- lapply(df.list2, function(x) nrow(x)) #count number of repeat per scaffold

res2=lapply(df.list2, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table2$Repeated_region_size <- lapply(res2, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table2$Counts_median <- lapply(df.list2, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table2$Counts_sum <- lapply(df.list2, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table2$nb_monomer=as.numeric(table2$nb_monomer)
table2$Repeated_region_size=as.numeric(table2$Repeated_region_size)
table2$Counts_sum=as.numeric(table2$Counts_sum)
table2$Counts_median=as.numeric(table2$Counts_median)
table2$Prop_repeated_region= table2$Repeated_region_size/table2$LG_length
table2$LGs <- factor(table2$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
table2$species<-"Tsi"
head(table2)
summary(table2)

ggplot(data=table2, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="blue", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()


Merge1=rbind(table, table2)
Merge1$species <- factor(Merge1$species,levels = c("Tcm", "Tsi"))


ggplot(data=Merge1, aes(x=LGs, y=Prop_repeated_region, fill=species)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()


Tcm$species<-"Tcm"
Tsi$species<-"Tsi"
bindTcmTsi<-rbind(Tcm, Tsi)

sharedTcm<-merge(Tcm, Tsi[,c(1,6)], by="X_remonomsqes", all.x = F, all.y = F)
sharedTsi<-merge(Tsi, Tcm[,c(1,6)], by="X_remonomsqes", all.x = F, all.y = F)
bindsharedTcmTsi<-rbind(sharedTcm, sharedTsi)

ggplot(bindsharedTcmTsi, aes(x=species, y=Counts_s)) + 
  ylim(0,500)+
  geom_violin()


## Tce vs Tms
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_LRv5a")
Tce<-read.table("Tce_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")

Tce_LG1.1<-read.table("Tce_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.2<-read.table("Tce_LRv5a_scf2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.3<-read.table("Tce_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.4<-read.table("Tce_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.5<-read.table("Tce_LRv5a_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.6<-read.table("Tce_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.7<-read.table("Tce_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1.8<-read.table("Tce_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG1<-rbind(Tce_LG1.1, Tce_LG1.2, Tce_LG1.3, Tce_LG1.4, Tce_LG1.5, Tce_LG1.6, Tce_LG1.7, Tce_LG1.8)

Tce_LG2.1<-read.table("Tce_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.2<-read.table("Tce_LRv5a_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.3<-read.table("Tce_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.4<-read.table("Tce_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.5<-read.table("Tce_LRv5a_scf13.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.6<-read.table("Tce_LRv5a_scf14.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.7<-read.table("Tce_LRv5a_scf15.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.8<-read.table("Tce_LRv5a_scf16.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.9<-read.table("Tce_LRv5a_scf17.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.10<-read.table("Tce_LRv5a_scf18.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.11<-read.table("Tce_LRv5a_scf19.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.12<-read.table("Tce_LRv5a_scf20.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.13<-read.table("Tce_LRv5a_scf21.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.14<-read.table("Tce_LRv5a_scf22.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.15<-read.table("Tce_LRv5a_scf23.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.16<-read.table("Tce_LRv5a_scf24.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.17<-read.table("Tce_LRv5a_scf25.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2.18<-read.table("Tce_LRv5a_scf26.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG2<-rbind(Tce_LG2.1, Tce_LG2.2, Tce_LG2.3, Tce_LG2.4, Tce_LG2.5, Tce_LG2.6, Tce_LG2.7, Tce_LG2.8,Tce_LG2.9, Tce_LG2.10, Tce_LG2.11, Tce_LG2.12, Tce_LG2.13, Tce_LG2.14, Tce_LG2.15, Tce_LG2.16, Tce_LG2.17, Tce_LG2.18)

Tce_LG3.1<-read.table("Tce_LRv5a_scf27.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.2<-read.table("Tce_LRv5a_scf28.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.3<-read.table("Tce_LRv5a_scf29.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.4<-read.table("Tce_LRv5a_scf30.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.5<-read.table("Tce_LRv5a_scf31.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.6<-read.table("Tce_LRv5a_scf32.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.7<-read.table("Tce_LRv5a_scf33.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.8<-read.table("Tce_LRv5a_scf34.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.9<-read.table("Tce_LRv5a_scf35.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.10<-read.table("Tce_LRv5a_scf36.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.11<-read.table("Tce_LRv5a_scf37.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.12<-read.table("Tce_LRv5a_scf38.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.13<-read.table("Tce_LRv5a_scf39.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.14<-read.table("Tce_LRv5a_scf40.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3.15<-read.table("Tce_LRv5a_scf41.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG3<-rbind(Tce_LG3.1, Tce_LG3.2, Tce_LG3.3, Tce_LG3.4, Tce_LG3.5, Tce_LG3.6, Tce_LG3.7, Tce_LG3.8,Tce_LG3.9, Tce_LG3.10, Tce_LG3.11, Tce_LG3.12, Tce_LG3.13, Tce_LG3.14, Tce_LG3.15)

Tce_LG4.1<-read.table("Tce_LRv5a_scf42.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4.2<-read.table("Tce_LRv5a_scf43.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4.3<-read.table("Tce_LRv5a_scf44.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4.4<-read.table("Tce_LRv5a_scf45.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4.5<-read.table("Tce_LRv5a_scf47.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4.6<-read.table("Tce_LRv5a_scf48.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4.7<-read.table("Tce_LRv5a_scf49.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG4<-rbind(Tce_LG4.1, Tce_LG4.2, Tce_LG4.3, Tce_LG4.4, Tce_LG4.5, Tce_LG4.6, Tce_LG4.7)

Tce_LG5.1<-read.table("Tce_LRv5a_scf50.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.2<-read.table("Tce_LRv5a_scf51.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.3<-read.table("Tce_LRv5a_scf52.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.4<-read.table("Tce_LRv5a_scf53.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.5<-read.table("Tce_LRv5a_scf54.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.6<-read.table("Tce_LRv5a_scf55.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.7<-read.table("Tce_LRv5a_scf56.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.8<-read.table("Tce_LRv5a_scf57.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.9<-read.table("Tce_LRv5a_scf58.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.10<-read.table("Tce_LRv5a_scf59.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5.11<-read.table("Tce_LRv5a_scf60.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG5<-rbind(Tce_LG5.1, Tce_LG5.2, Tce_LG5.3, Tce_LG5.4, Tce_LG5.5, Tce_LG5.6, Tce_LG5.7, Tce_LG5.8,Tce_LG5.9, Tce_LG5.10, Tce_LG5.11)

Tce_LG6.1<-read.table("Tce_LRv5a_scf61.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.2<-read.table("Tce_LRv5a_scf62.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.3<-read.table("Tce_LRv5a_scf63.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.4<-read.table("Tce_LRv5a_scf64.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.5<-read.table("Tce_LRv5a_scf65.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.6<-read.table("Tce_LRv5a_scf66.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.7<-read.table("Tce_LRv5a_scf67.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.8<-read.table("Tce_LRv5a_scf68.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6.9<-read.table("Tce_LRv5a_scf111.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG6<-rbind(Tce_LG6.1, Tce_LG6.2, Tce_LG6.3, Tce_LG6.4, Tce_LG6.5, Tce_LG6.6, Tce_LG6.7, Tce_LG6.8,Tce_LG6.9)

Tce_LG7.1<-read.table("Tce_LRv5a_scf69.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.2<-read.table("Tce_LRv5a_scf70.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.3<-read.table("Tce_LRv5a_scf71.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.4<-read.table("Tce_LRv5a_scf72.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.5<-read.table("Tce_LRv5a_scf73.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.6<-read.table("Tce_LRv5a_scf74.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.7<-read.table("Tce_LRv5a_scf75.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.8<-read.table("Tce_LRv5a_scf76.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.9<-read.table("Tce_LRv5a_scf77.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.10<-read.table("Tce_LRv5a_scf78.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7.11<-read.table("Tce_LRv5a_scf79.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG7<-rbind(Tce_LG7.1, Tce_LG7.2, Tce_LG7.3, Tce_LG7.4, Tce_LG7.5, Tce_LG7.6, Tce_LG7.7, Tce_LG7.8,Tce_LG7.9, Tce_LG7.10, Tce_LG7.11)

Tce_LG8.1<-read.table("Tce_LRv5a_scf80.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.2<-read.table("Tce_LRv5a_scf81.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.3<-read.table("Tce_LRv5a_scf82.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.4<-read.table("Tce_LRv5a_scf83.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.5<-read.table("Tce_LRv5a_scf84.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.6<-read.table("Tce_LRv5a_scf85.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.7<-read.table("Tce_LRv5a_scf86.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.8<-read.table("Tce_LRv5a_scf87.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.9<-read.table("Tce_LRv5a_scf88.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.10<-read.table("Tce_LRv5a_scf89.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.11<-read.table("Tce_LRv5a_scf90.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.12<-read.table("Tce_LRv5a_scf91.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.13<-read.table("Tce_LRv5a_scf92.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.14<-read.table("Tce_LRv5a_scf93.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.15<-read.table("Tce_LRv5a_scf94.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.16<-read.table("Tce_LRv5a_scf95.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.17<-read.table("Tce_LRv5a_scf96.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.18<-read.table("Tce_LRv5a_scf97.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.19<-read.table("Tce_LRv5a_scf98.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.20<-read.table("Tce_LRv5a_scf100.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8.21<-read.table("Tce_LRv5a_scf101.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG8<-rbind(Tce_LG8.1, Tce_LG8.2, Tce_LG8.3, Tce_LG8.4, Tce_LG8.5, Tce_LG8.6, Tce_LG8.7, Tce_LG8.8,Tce_LG8.9, Tce_LG8.10, Tce_LG8.11, Tce_LG8.12, Tce_LG8.13, Tce_LG8.14, Tce_LG8.15, Tce_LG8.16, Tce_LG8.17, Tce_LG8.18, Tce_LG8.19, Tce_LG8.20, Tce_LG8.21)

Tce_LG9.1<-read.table("Tce_LRv5a_scf102.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.2<-read.table("Tce_LRv5a_scf103.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.3<-read.table("Tce_LRv5a_scf104.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.4<-read.table("Tce_LRv5a_scf105.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.5<-read.table("Tce_LRv5a_scf106.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.6<-read.table("Tce_LRv5a_scf107.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.7<-read.table("Tce_LRv5a_scf108.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.8<-read.table("Tce_LRv5a_scf109.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9.9<-read.table("Tce_LRv5a_scf110.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG9<-rbind(Tce_LG9.1, Tce_LG9.2, Tce_LG9.3, Tce_LG9.4, Tce_LG9.5, Tce_LG9.6, Tce_LG9.7, Tce_LG9.8,Tce_LG9.9)

Tce_LG10.1<-read.table("Tce_LRv5a_scf112.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.2<-read.table("Tce_LRv5a_scf113.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.3<-read.table("Tce_LRv5a_scf114.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.4<-read.table("Tce_LRv5a_scf115.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.5<-read.table("Tce_LRv5a_scf116.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.6<-read.table("Tce_LRv5a_scf117.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.7<-read.table("Tce_LRv5a_scf118.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.8<-read.table("Tce_LRv5a_scf119.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.9<-read.table("Tce_LRv5a_scf120.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.10<-read.table("Tce_LRv5a_scf121.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.11<-read.table("Tce_LRv5a_scf122.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.12<-read.table("Tce_LRv5a_scf123.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10.13<-read.table("Tce_LRv5a_scf124.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG10<-rbind(Tce_LG10.1, Tce_LG10.2, Tce_LG10.3, Tce_LG10.4, Tce_LG10.5, Tce_LG10.6, Tce_LG10.7, Tce_LG10.8,Tce_LG10.9, Tce_LG10.10, Tce_LG10.11, Tce_LG10.12, Tce_LG10.13)

Tce_LG11.1<-read.table("Tce_LRv5a_scf125.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.2<-read.table("Tce_LRv5a_scf126.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.3<-read.table("Tce_LRv5a_scf127.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.4<-read.table("Tce_LRv5a_scf128.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.5<-read.table("Tce_LRv5a_scf129.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.6<-read.table("Tce_LRv5a_scf130.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.7<-read.table("Tce_LRv5a_scf132.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.8<-read.table("Tce_LRv5a_scf133.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11.9<-read.table("Tce_LRv5a_scf134.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG11<-rbind(Tce_LG11.1, Tce_LG11.2, Tce_LG11.3, Tce_LG11.4, Tce_LG11.5, Tce_LG11.6, Tce_LG11.7, Tce_LG11.8,Tce_LG11.9)

Tce_LG12.1<-read.table("Tce_LRv5a_scf135.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.2<-read.table("Tce_LRv5a_scf136.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.3<-read.table("Tce_LRv5a_scf137.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.4<-read.table("Tce_LRv5a_scf138.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.5<-read.table("Tce_LRv5a_scf139.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.6<-read.table("Tce_LRv5a_scf140.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.7<-read.table("Tce_LRv5a_scf141.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.8<-read.table("Tce_LRv5a_scf142.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.9<-read.table("Tce_LRv5a_scf143.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.10<-read.table("Tce_LRv5a_scf144.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12.11<-read.table("Tce_LRv5a_scf145.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG12<-rbind(Tce_LG12.1, Tce_LG12.2, Tce_LG12.3, Tce_LG12.4, Tce_LG12.5, Tce_LG12.6, Tce_LG12.7, Tce_LG12.8,Tce_LG12.9, Tce_LG12.10, Tce_LG12.11)

Tce_LG13.1<-read.table("Tce_LRv5a_scf146.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.2<-read.table("Tce_LRv5a_scf147.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.3<-read.table("Tce_LRv5a_scf148.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.4<-read.table("Tce_LRv5a_scf149.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.5<-read.table("Tce_LRv5a_scf150.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.6<-read.table("Tce_LRv5a_scf151.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.7<-read.table("Tce_LRv5a_scf152.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.8<-read.table("Tce_LRv5a_scf153.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.9<-read.table("Tce_LRv5a_scf154.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.10<-read.table("Tce_LRv5a_scf155.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.11<-read.table("Tce_LRv5a_scf156.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13.12<-read.table("Tce_LRv5a_scf157.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tce_LG13<-rbind(Tce_LG13.1, Tce_LG13.2, Tce_LG13.3, Tce_LG13.4, Tce_LG13.5, Tce_LG13.6, Tce_LG13.7, Tce_LG13.8,Tce_LG13.9, Tce_LG13.10, Tce_LG13.11, Tce_LG13.12)

Tce_LG_length<-read.table("scaffold_length_Tce.txt", header=F, sep="", quote = "")
colnames(Tce_LG_length)<-c("scaffolds", "length")
#Create summary table of LG length, number of distinct repeats and overall repeated region size per scaffold

df.list= list(Tce_LG1, Tce_LG2, Tce_LG3, Tce_LG4, Tce_LG5, Tce_LG6, Tce_LG7, Tce_LG8, Tce_LG9, Tce_LG10, Tce_LG11, Tce_LG12, Tce_LG13)

table=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13"))
colnames (table) [1] <- "LGs"
table$LG_length <- Tce_LG_length$length

table$nb_monomer <- lapply(df.list, function(x) nrow(x)) #count number of repeat per scaffold

res=lapply(df.list, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table$Repeated_region_size <- lapply(res, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table$Counts_median <- lapply(df.list, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table$Counts_sum <- lapply(df.list, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table$nb_monomer=as.numeric(table$nb_monomer)
table$Repeated_region_size=as.numeric(table$Repeated_region_size)
table$Counts_sum=as.numeric(table$Counts_sum)
table$Counts_median=as.numeric(table$Counts_median)
table$Prop_repeated_region= table$Repeated_region_size/table$LG_length
table$LGs <- factor(table$LGs,levels = c("LG3", "LG2", "LG8", "LG13", "LG7", "LG1", "LG5", "LG9", "LG10", "LG6", "LG4", "LG11", "LG12"))
table$species<-"Tce"
head(table)
summary(table)

ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="red", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3.5)+ ylab("Proportion tandem repeats")+
  theme_classic()



setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tms")
Tms<-read.table("Tms_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")

Tms_LG1<-read.table("Tms_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG2<-read.table("Tms_LRv5a_scf2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG3<-read.table("Tms_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG4<-read.table("Tms_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG5<-read.table("Tms_LRv5a_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG6<-read.table("Tms_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG7<-read.table("Tms_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG8<-read.table("Tms_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG9<-read.table("Tms_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG10<-read.table("Tms_LRv5a_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG11<-read.table("Tms_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG12<-read.table("Tms_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG13<-read.table("Tms_LRv5a_scf13.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tms_LG_length<-read.table("scaffold_length_Tms.txt", header=F, sep="", quote = "")
Tms_LG_length<-Tms_LG_length[1:13,]
colnames(Tms_LG_length)<-c("chm", "length")

df.list2= list(Tms_LG1, Tms_LG2, Tms_LG3, Tms_LG4, Tms_LG5, Tms_LG6, Tms_LG7, Tms_LG8, Tms_LG9, Tms_LG10, Tms_LG11, Tms_LG12, Tms_LG13)

table2=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13"))
colnames (table2) [1] <- "LGs"
table2$LG_length <- Tms_LG_length$length

table2$nb_monomer <- lapply(df.list2, function(x) nrow(x)) #count number of repeat per scaffold

res2=lapply(df.list2, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table2$Repeated_region_size <- lapply(res2, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table2$Counts_median <- lapply(df.list2, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table2$Counts_sum <- lapply(df.list2, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table2$nb_monomer=as.numeric(table2$nb_monomer)
table2$Repeated_region_size=as.numeric(table2$Repeated_region_size)
table2$Counts_sum=as.numeric(table2$Counts_sum)
table2$Counts_median=as.numeric(table2$Counts_median)
table2$Prop_repeated_region= table2$Repeated_region_size/table2$LG_length
table2$LGs <- factor(table2$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13"))
table2$species<-"Tms"
head(table2)
summary(table2)

ggplot(data=table2, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="blue", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=2.7)+ ylab("Proportion tandem repeats")+
  theme_classic()


Merge2=rbind(table, table2)
Merge2$species <- factor(Merge2$species,levels = c("Tce", "Tms"))


ggplot(data=Merge2, aes(x=LGs, y=Prop_repeated_region, fill=species)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()




## Tpa vs Tge
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")

Tpa_LG1<-read.table("Tpa_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG2<-read.table("Tpa_LRv5a_scf2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG3<-read.table("Tpa_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG4<-read.table("Tpa_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG5<-read.table("Tpa_LRv5a_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG6<-read.table("Tpa_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG7<-read.table("Tpa_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG8<-read.table("Tpa_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG9<-read.table("Tpa_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG10<-read.table("Tpa_LRv5a_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG11<-read.table("Tpa_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG12<-read.table("Tpa_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG13<-read.table("Tpa_LRv5a_scf13.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG14<-read.table("Tpa_LRv5a_scf14.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tpa_LG_length<-read.table("scaffold_length_Tpa.txt", header=T, sep="", quote = "")

#Create summary table of LG length, number of distinct repeats and overall repeated region size per scaffold

df.list= list(Tpa_LG1, Tpa_LG2, Tpa_LG3, Tpa_LG4, Tpa_LG5, Tpa_LG6, Tpa_LG7, Tpa_LG8, Tpa_LG9, Tpa_LG10, Tpa_LG11, Tpa_LG12, Tpa_LG13, Tpa_LG14)

table=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13", "LG14"))
colnames (table) [1] <- "LGs"
table$LG_length <- Tpa_LG_length$length

table$nb_monomer <- lapply(df.list, function(x) nrow(x)) #count number of repeat per scaffold

res=lapply(df.list, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table$Repeated_region_size <- lapply(res, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table$Counts_median <- lapply(df.list, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table$Counts_sum <- lapply(df.list, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table$nb_monomer=as.numeric(table$nb_monomer)
table$Repeated_region_size=as.numeric(table$Repeated_region_size)
table$Counts_sum=as.numeric(table$Counts_sum)
table$Counts_median=as.numeric(table$Counts_median)
table$Prop_repeated_region= table$Repeated_region_size/table$LG_length
table$LGs <- factor(table$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13", "LG14"))
table$species<-"Tpa"
head(table)
summary(table)

ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="red", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=2.5)+ ylab("Proportion tandem repeats")+
  theme_classic()


setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")

Tge_A_LG1<-read.table("Tge_A_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG2<-read.table("Tge_A_LRv5a_scf2.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG3<-read.table("Tge_A_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG4<-read.table("Tge_A_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG5<-read.table("Tge_A_LRv5a_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG6<-read.table("Tge_A_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG7<-read.table("Tge_A_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG8<-read.table("Tge_A_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG9<-read.table("Tge_A_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG10<-read.table("Tge_A_LRv5a_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG11<-read.table("Tge_A_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG12<-read.table("Tge_A_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG13<-read.table("Tge_A_LRv5a_scf13.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG14<-read.table("Tge_A_LRv5a_scf14.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
Tge_A_LG_length<-read.table("Tge_A_scaffold_length.txt", header=T, sep="", quote = "")

df.list2= list(Tge_A_LG1, Tge_A_LG2, Tge_A_LG3, Tge_A_LG4, Tge_A_LG5, Tge_A_LG6, Tge_A_LG7, Tge_A_LG8, Tge_A_LG9, Tge_A_LG10, Tge_A_LG11, Tge_A_LG12, Tge_A_LG13, Tge_A_LG14)

table2=data.frame(c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13", "LG14"))
colnames (table2) [1] <- "LGs"
table2$LG_length <- Tge_A_LG_length$length

table2$nb_monomer <- lapply(df.list2, function(x) nrow(x)) #count number of repeat per scaffold

res2=lapply(df.list2, function(w) { w$Repeated_region_size <- (w$Counts_s*w$monomer_length); w }) #create new list with Repeated_region_size values (counts*monomer size) for each repeat

table2$Repeated_region_size <- lapply(res2, function(z) sum(z$Repeated_region_size)) #count total repeated region per scaffold
table2$Counts_median <- lapply(df.list2, function(c) median(c$Counts_s)) #count total repeated region per scaffold
table2$Counts_sum <- lapply(df.list2, function(d) sum(d$Counts_s)) #count total repeated region per scaffold
table2$nb_monomer=as.numeric(table2$nb_monomer)
table2$Repeated_region_size=as.numeric(table2$Repeated_region_size)
table2$Counts_sum=as.numeric(table2$Counts_sum)
table2$Counts_median=as.numeric(table2$Counts_median)
table2$Prop_repeated_region= table2$Repeated_region_size/table2$LG_length
table2$LGs <- factor(table2$LGs,levels = c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12", "LG13", "LG14"))
table2$species<-"Tge"
head(table2)
summary(table2)

ggplot(data=table2, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="blue", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()


Merge2=rbind(table, table2)
Merge2$species <- factor(Merge2$species,levels = c("Tpa", "Tge"))


ggplot(data=Merge2, aes(x=LGs, y=Prop_repeated_region, fill=species)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()



####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################



############################
# TR proportion all sps GW #
############################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly")
data<-read.table("TRF_all_species.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data$Repeated_region_size=data$Counts_s*data$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("TRF_Tps.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length
#data_tps<-subset(data_tps, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length
#data_tdi<-subset(data_tdi, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tcm")
data_tcm<-read.table("Tcm_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tcm$Repeated_region_size=data_tcm$Counts_s*data_tcm$monomer_length
#data_tcm<-subset(data_tcm, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length
#data_tsi<-subset(data_tsi, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
data_tce<-read.table("mod_hic_output.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length
#data_tce<-subset(data_tce, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_LRv5a")
data_tce<-read.table("Tce_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tms")
data_tms<-read.table("Tms_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tms$Repeated_region_size=data_tms$Counts_s*data_tms$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length
#data_tpa<-subset(data_tpa, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length
#data_tge<-subset(data_tge, monomer_length<21)

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tbi_LRv4a")
data_tbi<-read.table("Tbi_LRv4a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tbi$Repeated_region_size=data_tbi$Counts_s*data_tbi$monomer_length
#data_tbi<-subset(data_tbi, monomer_length<21)


setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tte_LRv1")
data_Tte<-read.table("Tte_LRv1.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_Tte$Repeated_region_size=data_Tte$Counts_s*data_Tte$monomer_length
#data_Tte<-subset(data_Tte, monomer_length<21)


venn <- list(sex=data_tbi$X_remonomsqes, asex=data_Tte$X_remonomsqes)
ggVennDiagram(venn)

Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1339820256, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1292580112, reproduction = "asex")
Tcm= data.frame(nb_monomer = nrow(data_tcm), monomer_length_median = median(data_tcm$monomer_length), Repeated_region_size= sum(data_tcm$Repeated_region_size), Counts_sum= sum(data_tcm$Counts_s), Counts_mean= mean(data_tcm$Counts_s), Counts_median= median(data_tcm$Counts_s), species = "Tcm", genome_size=1294614351, reproduction = "sex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1324343855, reproduction = "asex")
Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1211418897, reproduction = "sex")#1313940951
Tms= data.frame(nb_monomer = nrow(data_tms), monomer_length_median = median(data_tms$monomer_length), Repeated_region_size= sum(data_tms$Repeated_region_size), Counts_sum= sum(data_tms$Counts_s), Counts_mean= mean(data_tms$Counts_s), Counts_median= median(data_tms$Counts_s), species = "Tms", genome_size=1168089169, reproduction = "asex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1145684399, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1121008509, reproduction = "asex")
Tbi= data.frame(nb_monomer = nrow(data_tbi), monomer_length_median = median(data_tbi$monomer_length), Repeated_region_size= sum(data_tbi$Repeated_region_size), Counts_sum= sum(data_tbi$Counts_s), Counts_mean= mean(data_tbi$Counts_s), Counts_median= median(data_tbi$Counts_s), species = "Tbi", genome_size=1145684399, reproduction = "sex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tcm)
merge_GW=rbind(merge_GW, Tsi)
merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tms)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)
merge_GW=rbind(merge_GW, Tbi)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tcm", "Tsi", "Tce", "Tms", "Tpa", "Tge_A", "Tbi"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=5.5)+
  scale_fill_manual(values=c("blue", "red"))+ ylab("Proportion tandem repeats")+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Repeated_region_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Repeated_region_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=5.5)+
  scale_fill_manual(values=c("blue", "red"))+ ylab("Proportion tandem repeats")+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  scale_fill_manual(values=c("blue", "red"))+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  scale_fill_manual(values=c("blue", "red"))+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  scale_fill_manual(values=c("blue", "red"))+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()



###Proportion between shared monomers
data=data[order(data$Counts_s, decreasing = TRUE),]

#Removing duplicate words in "Scaffold_name"
data$Scaffold_name = as.data.frame(sapply(data$Scaffold_name, function(x) gsub("\"", "", x))) #Removing quotes

reduce_row = function(i) {
  split = strsplit(i, split=",")[[1]]
  paste(unique(split), collapse = ",") 
}

data$Scaffold_name = apply(data, 1, reduce_row) #Removing duplicates

#Counts shared monomers
Tps=data.frame(rowSums(sapply(c("Tps_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tps)<-"Tps"
Tdi=data.frame(rowSums(sapply(c("Tdi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tdi)<-"Tdi"
Tcm=data.frame(rowSums(sapply(c("Tcm_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tcm)<-"Tcm"
Tsi=data.frame(rowSums(sapply(c("Tsi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tsi)<-"Tsi"
Tce=data.frame(rowSums(sapply(c("ScVQC5J_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tce)<-"Tce"
Tpa=data.frame(rowSums(sapply(c("Tpa_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tpa)<-"Tpa"
Tge=data.frame(rowSums(sapply(c("Tge_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tge)<-"Tge"

all_sps<-cbind(data[,c(1,2,3,6)], Tps)
all_sps<-cbind(all_sps,Tdi)
all_sps<-cbind(all_sps,Tcm)
all_sps<-cbind(all_sps,Tsi)
all_sps<-cbind(all_sps,Tce)
all_sps<-cbind(all_sps,Tpa)
all_sps<-cbind(all_sps,Tge)
all_sps$total<-all_sps$Tps+all_sps$Tdi+all_sps$Tcm+all_sps$Tsi+all_sps$Tce+all_sps$Tpa+all_sps$Tge

shared<-subset(all_sps, total=="1")

sharedTps<-subset(data_tps, (X_remonomsqes %in% shared$X_remonomsqes))
sharedTdi<-subset(data_tdi, (X_remonomsqes %in% shared$X_remonomsqes))
sharedTcm<-subset(data_tcm, (X_remonomsqes %in% shared$X_remonomsqes))
sharedTsi<-subset(data_tsi, (X_remonomsqes %in% shared$X_remonomsqes))
sharedTce<-subset(data_tce, (X_remonomsqes %in% shared$X_remonomsqes))
sharedTpa<-subset(data_tpa, (X_remonomsqes %in% shared$X_remonomsqes))
sharedTge<-subset(data_tge, (X_remonomsqes %in% shared$X_remonomsqes))

Tps= data.frame(nb_monomer = nrow(sharedTps), monomer_length_median = median(sharedTps$monomer_length), Repeated_region_size= sum(sharedTps$Repeated_region_size), Counts_sum= sum(sharedTps$Counts_s), Counts_mean= mean(sharedTps$Counts_s), Counts_median= median(sharedTps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(sharedTdi), monomer_length_median = median(sharedTdi$monomer_length), Repeated_region_size= sum(sharedTdi$Repeated_region_size), Counts_sum= sum(sharedTdi$Counts_s), Counts_mean= mean(sharedTdi$Counts_s), Counts_median= median(sharedTdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tcm= data.frame(nb_monomer = nrow(sharedTcm), monomer_length_median = median(sharedTcm$monomer_length), Repeated_region_size= sum(sharedTcm$Repeated_region_size), Counts_sum= sum(sharedTcm$Counts_s), Counts_mean= mean(sharedTcm$Counts_s), Counts_median= median(sharedTcm$Counts_s), species = "Tcm", genome_size=1290000000, reproduction = "sex")
Tsi= data.frame(nb_monomer = nrow(sharedTsi), monomer_length_median = median(sharedTsi$monomer_length), Repeated_region_size= sum(sharedTsi$Repeated_region_size), Counts_sum= sum(sharedTsi$Counts_s), Counts_mean= mean(sharedTsi$Counts_s), Counts_median= median(sharedTsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
Tce= data.frame(nb_monomer = nrow(sharedTce), monomer_length_median = median(sharedTce$monomer_length), Repeated_region_size= sum(sharedTce$Repeated_region_size), Counts_sum= sum(sharedTce$Counts_s), Counts_mean= mean(sharedTce$Counts_s), Counts_median= median(sharedTce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(sharedTpa), monomer_length_median = median(sharedTpa$monomer_length), Repeated_region_size= sum(sharedTpa$Repeated_region_size), Counts_sum= sum(sharedTpa$Counts_s), Counts_mean= mean(sharedTpa$Counts_s), Counts_median= median(sharedTpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(sharedTge), monomer_length_median = median(sharedTge$monomer_length), Repeated_region_size= sum(sharedTge$Repeated_region_size), Counts_sum= sum(sharedTge$Counts_s), Counts_mean= mean(sharedTge$Counts_s), Counts_median= median(sharedTge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")


merge_shared=rbind(Tps, Tdi)
merge_shared=rbind(merge_shared, Tcm)
merge_shared=rbind(merge_shared, Tsi)
merge_shared=rbind(merge_shared, Tce)
merge_shared=rbind(merge_shared, Tpa)
merge_shared=rbind(merge_shared, Tge_A)

merge_shared$species <- factor(merge_shared$species,levels = c("Tps", "Tdi", "Tcm", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_shared$Prop_repeated_region=merge_shared$Repeated_region_size/merge_shared$genome_size
merge_shared$Prop_counts=merge_shared$Counts_sum/merge_shared$genome_size
head(merge_shared)

ggplot(data=merge_shared, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





###Proportion between shared monomers (pair-wise sister species)
#Tps vs Tdi
all_sps$Tps_Tdi<-all_sps$Tps+all_sps$Tdi

sharedTps_Tdi<-subset(all_sps, Tps_Tdi=="2")

sharedTps<-subset(data_tps, (X_remonomsqes %in% sharedTps_Tdi$X_remonomsqes))
sharedTdi<-subset(data_tdi, (X_remonomsqes %in% sharedTps_Tdi$X_remonomsqes))

Tps= data.frame(nb_monomer = nrow(sharedTps), monomer_length_median = median(sharedTps$monomer_length), Repeated_region_size= sum(sharedTps$Repeated_region_size), Counts_sum= sum(sharedTps$Counts_s), Counts_mean= mean(sharedTps$Counts_s), Counts_median= median(sharedTps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(sharedTdi), monomer_length_median = median(sharedTdi$monomer_length), Repeated_region_size= sum(sharedTdi$Repeated_region_size), Counts_sum= sum(sharedTdi$Counts_s), Counts_mean= mean(sharedTdi$Counts_s), Counts_median= median(sharedTdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")

merge_shared=rbind(Tps, Tdi)
merge_shared$species <- factor(merge_shared$species,levels = c("Tps", "Tdi"))
merge_shared$Prop_repeated_region=merge_shared$Repeated_region_size/merge_shared$genome_size
merge_shared$Prop_counts=merge_shared$Counts_sum/merge_shared$genome_size
head(merge_shared)

ggplot(data=merge_shared, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()




#Tcm vs Tsi
all_sps$Tcm_Tsi<-all_sps$Tcm+all_sps$Tsi

sharedTcm_Tsi<-subset(all_sps, Tcm_Tsi=="2")

sharedTcm<-subset(data_tcm, (X_remonomsqes %in% sharedTcm_Tsi$X_remonomsqes))
sharedTsi<-subset(data_tsi, (X_remonomsqes %in% sharedTcm_Tsi$X_remonomsqes))

Tcm= data.frame(nb_monomer = nrow(sharedTcm), monomer_length_median = median(sharedTcm$monomer_length), Repeated_region_size= sum(sharedTcm$Repeated_region_size), Counts_sum= sum(sharedTcm$Counts_s), Counts_mean= mean(sharedTcm$Counts_s), Counts_median= median(sharedTcm$Counts_s), species = "Tcm", genome_size=1290000000, reproduction = "sex")
Tsi= data.frame(nb_monomer = nrow(sharedTsi), monomer_length_median = median(sharedTsi$monomer_length), Repeated_region_size= sum(sharedTsi$Repeated_region_size), Counts_sum= sum(sharedTsi$Counts_s), Counts_mean= mean(sharedTsi$Counts_s), Counts_median= median(sharedTsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")

merge_shared=rbind(Tcm, Tsi)
merge_shared$species <- factor(merge_shared$species,levels = c("Tcm", "Tsi"))
merge_shared$Prop_repeated_region=merge_shared$Repeated_region_size/merge_shared$genome_size
merge_shared$Prop_counts=merge_shared$Counts_sum/merge_shared$genome_size
head(merge_shared)

ggplot(data=merge_shared, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





#Tpa vs Tge
all_sps$Tpa_Tge<-all_sps$Tpa+all_sps$Tge

sharedTpa_Tge<-subset(all_sps, Tpa_Tge=="2")

sharedTpa<-subset(data_tpa, (X_remonomsqes %in% sharedTpa_Tge$X_remonomsqes))
sharedTge<-subset(data_tge, (X_remonomsqes %in% sharedTpa_Tge$X_remonomsqes))

Tpa= data.frame(nb_monomer = nrow(sharedTpa), monomer_length_median = median(sharedTpa$monomer_length), Repeated_region_size= sum(sharedTpa$Repeated_region_size), Counts_sum= sum(sharedTpa$Counts_s), Counts_mean= mean(sharedTpa$Counts_s), Counts_median= median(sharedTpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge= data.frame(nb_monomer = nrow(sharedTge), monomer_length_median = median(sharedTge$monomer_length), Repeated_region_size= sum(sharedTge$Repeated_region_size), Counts_sum= sum(sharedTge$Counts_s), Counts_mean= mean(sharedTge$Counts_s), Counts_median= median(sharedTge$Counts_s), species = "Tge", genome_size=1111346744, reproduction = "asex")

merge_shared=rbind(Tpa, Tge)
merge_shared$species <- factor(merge_shared$species,levels = c("Tpa", "Tge"))
merge_shared$Prop_repeated_region=merge_shared$Repeated_region_size/merge_shared$genome_size
merge_shared$Prop_counts=merge_shared$Counts_sum/merge_shared$genome_size
head(merge_shared)

ggplot(data=merge_shared, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()








####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################



###############################
# TR proportion expressed Tps #
###############################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("TRF_Tps.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

#Removing duplicate words in "Scaffold_name"
data_tps$Scaffold_name = as.data.frame(sapply(data_tps$Scaffold_name, function(x) gsub("\"", "", x))) #Removing quotes

reduce_row = function(i) {
  split = strsplit(i, split=",")[[1]]
  paste(unique(split), collapse =
          ",") 
}

data_tps$Scaffold_name = apply(data_tps, 1, reduce_row) #Removing duplicates
data_tps<-data_tps[order(data_tps$Counts_s, decreasing = TRUE),]

data_tps$Count_LGs <- rowSums(sapply(c("\\bTps_LRv5b_scf1\\b", "\\bTps_LRv5b_scf2\\b", "\\bTps_LRv5b_scf3\\b", "\\bTps_LRv5b_scf4\\b", "\\bTps_LRv5b_scf5\\b", "\\bTps_LRv5b_scf6\\b", "\\bTps_LRv5b_scf7\\b", "\\bTps_LRv5b_scf8\\b", "\\bTps_LRv5b_scf9\\b", "\\bTps_LRv5b_scf10\\b", "\\bTps_LRv5b_scf11\\b", "\\bTps_LRv5b_scf12\\b"),
                                   function(x) grepl(x, data_tps$Scaffold_name)))


setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_transcriptome_long_reads/Tps")
F_G<-read.table("transcripts-2.collapsed.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
M_G<-read.table("transcripts-3.collapsed.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
F_H<-read.table("transcripts-5.collapsed.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
M_H<-read.table("transcripts-6.collapsed.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")

expressed1<-rbind(F_G, M_G, F_H, M_H)
expressed2<-data.frame(unique(expressed1$X_remonomsqes))
colnames(expressed2)<-"X_remonomsqes"
expressed<-subset(data_tps, (X_remonomsqes %in% expressed2$X_remonomsqes))
notexpressed<-subset(data_tps, !(X_remonomsqes %in% expressed2$X_remonomsqes))




###########################################################################
# Sex-asex (Tps vs Tdi) comparison between expressed and non-expressed TR #
###########################################################################

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly")
data<-read.table("TRF_all_species.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data$Repeated_region_size=data$Counts_s*data$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length


###Proportion between shared monomers
data=data[order(data$Counts_s, decreasing = TRUE),]

##Removing duplicate words in "Scaffold_name"
data$Scaffold_name = as.data.frame(sapply(data$Scaffold_name, function(x) gsub("\"", "", x))) #Removing quotes

reduce_row = function(i) {
  split = strsplit(i, split=",")[[1]]
  paste(unique(split), collapse = ",") 
}

data$Scaffold_name = apply(data, 1, reduce_row) #Removing duplicates

##Counts shared monomers
Tps=data.frame(rowSums(sapply(c("Tps_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tps)<-"Tps"
Tdi=data.frame(rowSums(sapply(c("Tdi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tdi)<-"Tdi"

all_sps<-cbind(data[,c(1,2,3,6)], Tps)
all_sps<-cbind(all_sps,Tdi)
all_sps$Tps_Tdi<-all_sps$Tps+all_sps$Tdi

sharedTps_Tdi<-subset(all_sps, Tps_Tdi=="2")

sharedTps<-subset(data_tps, (X_remonomsqes %in% sharedTps_Tdi$X_remonomsqes))
sharedTdi<-subset(data_tdi, (X_remonomsqes %in% sharedTps_Tdi$X_remonomsqes))

expressedTps<-subset(sharedTps, (X_remonomsqes %in% expressed$X_remonomsqes))
expressedTdi<-subset(sharedTdi, (X_remonomsqes %in% expressed$X_remonomsqes))
notexpressedTps<-subset(sharedTps, (X_remonomsqes %in% notexpressed$X_remonomsqes))
notexpressedTdi<-subset(sharedTdi, (X_remonomsqes %in% notexpressed$X_remonomsqes))

##
#Analysis on expressed genes
Tps= data.frame(nb_monomer = nrow(expressedTps), monomer_length_median = median(expressedTps$monomer_length), Repeated_region_size= sum(expressedTps$Repeated_region_size), Counts_sum= sum(expressedTps$Counts_s), Counts_mean= mean(expressedTps$Counts_s), Counts_median= median(expressedTps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(expressedTdi), monomer_length_median = median(expressedTdi$monomer_length), Repeated_region_size= sum(expressedTdi$Repeated_region_size), Counts_sum= sum(expressedTdi$Counts_s), Counts_mean= mean(expressedTdi$Counts_s), Counts_median= median(expressedTdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
#Analysis on non-expressed genes
Tps= data.frame(nb_monomer = nrow(notexpressedTps), monomer_length_median = median(notexpressedTps$monomer_length), Repeated_region_size= sum(notexpressedTps$Repeated_region_size), Counts_sum= sum(notexpressedTps$Counts_s), Counts_mean= mean(notexpressedTps$Counts_s), Counts_median= median(notexpressedTps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(notexpressedTdi), monomer_length_median = median(notexpressedTdi$monomer_length), Repeated_region_size= sum(notexpressedTdi$Repeated_region_size), Counts_sum= sum(notexpressedTdi$Counts_s), Counts_mean= mean(notexpressedTdi$Counts_s), Counts_median= median(notexpressedTdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
##

merge_shared=rbind(Tps, Tdi)
merge_shared$species <- factor(merge_shared$species,levels = c("Tps", "Tdi"))
merge_shared$Prop_repeated_region=merge_shared$Repeated_region_size/merge_shared$genome_size
merge_shared$Prop_counts=merge_shared$Counts_sum/merge_shared$genome_size
head(merge_shared)

ggplot(data=merge_shared, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_shared, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 5)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()











####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################


###############################
# TR proportion all sps X chm #
###############################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf3.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()




###########################################
# TR proportion all sps chm4Tdi - chm4Tge #
###########################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf4.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()



###########################################
# TR proportion all sps chm5Tdi - chm6Tge #
###########################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf5.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()







###########################################
# TR proportion all sps chm6Tdi - chm7Tge #
###########################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf6.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





###########################################
# TR proportion all sps chm7Tdi - chm8Tge #
###########################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf7.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





###########################################
# TR proportion all sps chm8Tdi - chm9Tge #
###########################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf8.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





###########################################
# TR proportion all sps chm9Tdi - chm11Tge #
###########################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf9.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





#############################################
# TR proportion all sps chm10Tdi - chm12Tge #
#############################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf10.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()





#############################################
# TR proportion all sps chm11Tdi - chm13Tge #
#############################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf11.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf13.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf13.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()






#############################################
# TR proportion all sps chm12Tdi - chm14Tge #
#############################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tps")
data_tps<-read.table("Tps_LRv5b_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tps$Repeated_region_size=data_tps$Counts_s*data_tps$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")
data_tdi<-read.table("Tdi_LRv5a_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tdi$Repeated_region_size=data_tdi$Counts_s*data_tdi$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
data_tsi<-read.table("Tsi_LRv5b_scf12.fasta.2.7.7.80.10.100.2000_parse_CMwRC.txt", header=T, sep="\t", quote="")
data_tsi$Repeated_region_size=data_tsi$Counts_s*data_tsi$monomer_length

#setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tce_Nosil/last_genome_assembly")
#data_tce<-read.table("ScVQC5J_9884_HRSCAF11163.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
#data_tce$Repeated_region_size=data_tce$Counts_s*data_tce$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tpa")
data_tpa<-read.table("Tpa_LRv5a_scf14.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tpa$Repeated_region_size=data_tpa$Counts_s*data_tpa$monomer_length

setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tge_A")
data_tge<-read.table("Tge_A_LRv5a_scf14.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
data_tge$Repeated_region_size=data_tge$Counts_s*data_tge$monomer_length


Tps= data.frame(nb_monomer = nrow(data_tps), monomer_length_median = median(data_tps$monomer_length), Repeated_region_size= sum(data_tps$Repeated_region_size), Counts_sum= sum(data_tps$Counts_s), Counts_mean= mean(data_tps$Counts_s), Counts_median= median(data_tps$Counts_s), species = "Tps", genome_size=1340000000, reproduction = "sex")
Tdi= data.frame(nb_monomer = nrow(data_tdi), monomer_length_median = median(data_tdi$monomer_length), Repeated_region_size= sum(data_tdi$Repeated_region_size), Counts_sum= sum(data_tdi$Counts_s), Counts_mean= mean(data_tdi$Counts_s), Counts_median= median(data_tdi$Counts_s), species = "Tdi", genome_size=1290000000, reproduction = "asex")
Tsi= data.frame(nb_monomer = nrow(data_tsi), monomer_length_median = median(data_tsi$monomer_length), Repeated_region_size= sum(data_tsi$Repeated_region_size), Counts_sum= sum(data_tsi$Counts_s), Counts_mean= mean(data_tsi$Counts_s), Counts_median= median(data_tsi$Counts_s), species = "Tsi", genome_size=1320000000, reproduction = "asex")
#Tce= data.frame(nb_monomer = nrow(data_tce), monomer_length_median = median(data_tce$monomer_length), Repeated_region_size= sum(data_tce$Repeated_region_size), Counts_sum= sum(data_tce$Counts_s), Counts_mean= mean(data_tce$Counts_s), Counts_median= median(data_tce$Counts_s), species = "Tce", genome_size=1320000000, reproduction = "sex")
Tpa= data.frame(nb_monomer = nrow(data_tpa), monomer_length_median = median(data_tpa$monomer_length), Repeated_region_size= sum(data_tpa$Repeated_region_size), Counts_sum= sum(data_tpa$Counts_s), Counts_mean= mean(data_tpa$Counts_s), Counts_median= median(data_tpa$Counts_s), species = "Tpa", genome_size=1105020916, reproduction = "sex")
Tge_A= data.frame(nb_monomer = nrow(data_tge), monomer_length_median = median(data_tge$monomer_length), Repeated_region_size= sum(data_tge$Repeated_region_size), Counts_sum= sum(data_tge$Counts_s), Counts_mean= mean(data_tge$Counts_s), Counts_median= median(data_tge$Counts_s), species = "Tge_A", genome_size=1111346744, reproduction = "asex")

merge_GW=rbind(Tps, Tdi)
merge_GW=rbind(merge_GW, Tsi)
#merge_GW=rbind(merge_GW, Tce)
merge_GW=rbind(merge_GW, Tpa)
merge_GW=rbind(merge_GW, Tge_A)

merge_GW$species <- factor(merge_GW$species,levels = c("Tps", "Tdi", "Tsi", "Tce", "Tpa", "Tge_A"))
merge_GW$Prop_repeated_region=merge_GW$Repeated_region_size/merge_GW$genome_size
merge_GW$Prop_counts=merge_GW$Counts_sum/merge_GW$genome_size
head(merge_GW)

ggplot(data=merge_GW, aes(x=species, y=Prop_repeated_region, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=nb_monomer, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(nb_monomer, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_sum, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_sum, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=Counts_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(Counts_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=monomer_length_median, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(monomer_length_median, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()

ggplot(data=merge_GW, aes(x=species, y=genome_size, fill=reproduction)) +
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(genome_size, digits = 3)), vjust=1.6, color="black", position = position_dodge(0.9),size=3.5)+
  theme_minimal()










####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################


##############################
# Shared TRs between species #
##############################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly")

data<-read.table("TRF_all_species.fasta.2.7.7.80.10.100.2000_parse_V2_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
head(data)
data=data[order(data$Counts_s, decreasing = TRUE),]

#Removing duplicate words in "Scaffold_name"
data$Scaffold_name = as.data.frame(sapply(data$Scaffold_name, function(x) gsub("\"", "", x))) #Removing quotes

reduce_row = function(i) {
  split = strsplit(i, split=",")[[1]]
  paste(unique(split), collapse = ",") 
}

data$Scaffold_name = apply(data, 1, reduce_row) #Removing duplicates

#Counts shared monomers
Tps=data.frame(rowSums(sapply(c("Tps_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tps)<-"Tps"
Tdi=data.frame(rowSums(sapply(c("Tdi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tdi)<-"Tdi"
Tcm=data.frame(rowSums(sapply(c("Tcm_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tcm)<-"Tcm"
Tsi=data.frame(rowSums(sapply(c("Tsi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tsi)<-"Tsi"
Tce=data.frame(rowSums(sapply(c("Tce_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tce)<-"Tce"
Tms=data.frame(rowSums(sapply(c("Tms_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tms)<-"Tms"
Tpa=data.frame(rowSums(sapply(c("Tpa_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tpa)<-"Tpa"
Tge=data.frame(rowSums(sapply(c("Tge_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tge)<-"Tge"

all_sps<-cbind(data[,c(1,2,3,6)], Tps)
all_sps<-cbind(all_sps,Tdi)
all_sps<-cbind(all_sps,Tcm)
all_sps<-cbind(all_sps,Tsi)
all_sps<-cbind(all_sps,Tce)
all_sps<-cbind(all_sps,Tms)
all_sps<-cbind(all_sps,Tpa)
all_sps<-cbind(all_sps,Tge)
all_sps$total<-all_sps$Tps+all_sps$Tdi+all_sps$Tcm+all_sps$Tsi+all_sps$Tce+all_sps$Tms+all_sps$Tpa+all_sps$Tge

ggplot(all_sps, aes(x=total)) + 
  geom_bar(color="black", fill="white")+
  xlab("Number of species")+
  ylab("Number of shared repeats")+
  theme_classic()

nrow(all_sps)
nrow(subset(all_sps, total=="1"))


#Proportion shared repeats among species (pairwise comparison)
all_sps$Tps_Tcm<-all_sps$Tps+all_sps$Tcm
Tps_Tcm<-subset(all_sps, Tps_Tcm=="2")

all_sps$Tps_Tce<-all_sps$Tps+all_sps$Tce
Tps_Tce<-subset(all_sps, Tps_Tce=="2")

all_sps$Tps_Tpa<-all_sps$Tps+all_sps$Tpa
Tps_Tpa<-subset(all_sps, Tps_Tpa=="2")

all_sps$Tcm_Tce<-all_sps$Tcm+all_sps$Tce
Tcm_Tce<-subset(all_sps, Tcm_Tce=="2")

all_sps$Tcm_Tpa<-all_sps$Tcm+all_sps$Tpa
Tcm_Tpa<-subset(all_sps, Tcm_Tpa=="2")

all_sps$Tce_Tpa<-all_sps$Tce+all_sps$Tpa
Tce_Tpa<-subset(all_sps, Tce_Tpa=="2")

all_sps$Tdi_Tsi<-all_sps$Tdi+all_sps$Tsi
Tdi_Tsi<-subset(all_sps, Tdi_Tsi=="2")

all_sps$Tdi_Tms<-all_sps$Tdi+all_sps$Tms
Tdi_Tms<-subset(all_sps, Tdi_Tms=="2")

all_sps$Tdi_Tge<-all_sps$Tdi+all_sps$Tge
Tdi_Tge<-subset(all_sps, Tdi_Tge=="2")

all_sps$Tsi_Tms<-all_sps$Tsi+all_sps$Tms
Tsi_Tms<-subset(all_sps, Tsi_Tms=="2")

all_sps$Tsi_Tge<-all_sps$Tsi+all_sps$Tge
Tsi_Tge<-subset(all_sps, Tsi_Tge=="2")

all_sps$Tms_Tge<-all_sps$Tms+all_sps$Tge
Tms_Tge<-subset(all_sps, Tms_Tge=="2")

all_sps$Tps_Tdi<-all_sps$Tps+all_sps$Tdi
Tps_Tdi<-subset(all_sps, Tps_Tdi=="2")

all_sps$Tcm_Tsi<-all_sps$Tcm+all_sps$Tsi
Tcm_Tsi<-subset(all_sps, Tcm_Tsi=="2")

all_sps$Tpa_Tge<-all_sps$Tpa+all_sps$Tge
Tpa_Tge<-subset(all_sps, Tpa_Tge=="2")

all_sps$Tps_Tsi<-all_sps$Tps+all_sps$Tsi
Tps_Tsi<-subset(all_sps, Tps_Tsi=="2")

all_sps$Tps_Tge<-all_sps$Tps+all_sps$Tge
Tps_Tge<-subset(all_sps, Tps_Tge=="2")

all_sps$Tcm_Tdi<-all_sps$Tcm+all_sps$Tdi
Tcm_Tdi<-subset(all_sps, Tcm_Tdi=="2")

all_sps$Tce_Tdi<-all_sps$Tce+all_sps$Tdi
Tce_Tdi<-subset(all_sps, Tce_Tdi=="2")

all_sps$Tpa_Tdi<-all_sps$Tpa+all_sps$Tdi
Tpa_Tdi<-subset(all_sps, Tpa_Tdi=="2")

all_sps$Tcm_Tge<-all_sps$Tcm+all_sps$Tge
Tcm_Tge<-subset(all_sps, Tcm_Tge=="2")

all_sps$Tce_Tsi<-all_sps$Tce+all_sps$Tsi
Tce_Tsi<-subset(all_sps, Tce_Tsi=="2")

all_sps$Tpa_Tsi<-all_sps$Tpa+all_sps$Tsi
Tpa_Tsi<-subset(all_sps, Tpa_Tsi=="2")

all_sps$Tce_Tge<-all_sps$Tce+all_sps$Tge
Tce_Tge<-subset(all_sps, Tce_Tge=="2")

sex<-data.frame(nrow(Tps_Tcm)/nrow(data))
colnames(sex)<-"Prop_shared"
sex[nrow(sex) + 1,] = nrow(Tps_Tce)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tps_Tpa)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tcm_Tce)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tcm_Tpa)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tce_Tpa)/nrow(data)
sex$species<-c("Tps_Tcm", "Tps_Tce", "Tps_Tpa", "Tcm_Tce", "Tcm_Tpa", "Tce_Tpa")
sex$repro<-"sex"

asex<-data.frame(nrow(Tdi_Tsi)/nrow(data))
colnames(asex)<-"Prop_shared"
asex[nrow(asex) + 1,] = nrow(Tdi_Tms)/nrow(data)
asex[nrow(asex) + 1,] = nrow(Tdi_Tge)/nrow(data)
asex[nrow(asex) + 1,] = nrow(Tsi_Tms)/nrow(data)
asex[nrow(asex) + 1,] = nrow(Tsi_Tge)/nrow(data)
asex[nrow(asex) + 1,] = nrow(Tms_Tge)/nrow(data)
asex$species<-c("Tdi_Tsi", "Tdi_Tms", "Tdi_Tge", "Tsi_Tms", "Tsi_Tge", "Tms_Tge")
asex$repro<-"asex"

sex_asex<-rbind(sex,asex)
sex_asex$species <- factor(sex_asex$species,levels = c("Tps_Tcm", "Tdi_Tsi", "Tps_Tce", "Tps_Tpa", "Tdi_Tge", "Tcm_Tce", "Tcm_Tpa", "Tsi_Tge", "Tce_Tpa"))

sister_sp<-data.frame(nrow(Tps_Tdi)/nrow(data))
colnames(sister_sp)<-"Prop_shared"
sister_sp[nrow(sister_sp) + 1,] = nrow(Tcm_Tsi)/nrow(data)
sister_sp[nrow(sister_sp) + 1,] = nrow(Tpa_Tge)/nrow(data)
sister_sp$species<-c("Tps_Tdi", "Tcm_Tsi", "Tpa_Tge")
sister_sp$species <- factor(sister_sp$species,levels = c("Tps_Tdi", "Tcm_Tsi", "Tpa_Tge"))
sister_sp$repro<-"sex_asex"

sp<-data.frame(nrow(Tps_Tsi)/nrow(data))
colnames(sp)<-"Prop_shared"
sp[nrow(sp) + 1,] = nrow(Tps_Tge)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tcm_Tdi)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tce_Tdi)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tpa_Tdi)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tcm_Tge)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tce_Tsi)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tpa_Tsi)/nrow(data)
sp[nrow(sp) + 1,] = nrow(Tce_Tge)/nrow(data)

all<-rbind(sex_asex,sister_sp)
all$species <- factor(all$species,levels = c("Tps_Tdi","Tcm_Tsi","Tpa_Tge","Tps_Tcm", "Tdi_Tsi", "Tps_Tce", "Tps_Tpa", "Tdi_Tge", "Tcm_Tce", "Tcm_Tpa", "Tsi_Tge", "Tce_Tpa"))
all$percent<-all$Prop_shared*100

ggplot(data=all, aes(x=species, y=Prop_shared, fill=repro)) +
  geom_bar(stat="identity")

ggplot(data=sex, aes(x=species, y=Prop_shared)) +
  geom_bar(stat="identity")

ggplot(data=asex, aes(x=species, y=Prop_shared)) +
  geom_bar(stat="identity")

ggplot(data=sex_asex, aes(x=species, y=Prop_shared, fill=repro)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("blue", "red"))

ggplot(data=sister_sp, aes(x=species, y=Prop_shared)) +
  geom_bar(stat="identity")

ggplot(sex_asex, aes(x=repro, y=Prop_shared)) + 
  geom_boxplot()


#Correlation shared motifs vs genetic distance
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly")

corr_table<-read.table("Correlation_GenDist_vs_sharedMotifs.txt", header=T, sep="\t", quote="")
head(corr_table)

ggplot(corr_table, aes(x=genetic.distance, y=prop_shared_motifs, color=reproduction)) + geom_point()

library("RColorBrewer")
ggplot(corr_table, aes(x=(genetic.distance), y=log(prop_shared_motifs), color=clade, label=species)) + geom_point()+geom_text(hjust=-0.1, vjust=0.5)+theme_minimal()+ 
  scale_color_manual(values = brewer.pal(n = 6, name = "RdBu"))+
  geom_line(data = cbind(corr_table, pred = predict(lmm)), aes(y = pred))

glm<-glm(corr_table$prop_shared_motifs~corr_table$genetic.distance+corr_table$reproduction:corr_table$clade)
anova(lm)

library(lme4)

lmm <- lmer(log(corr_table$prop_shared_motifs) ~ corr_table$reproduction * corr_table$genetic.distance + (1|corr_table$clade))
summary(lmm)
plot(lmm)
library(car)
Anova(lmm)
library(MCMCglmm)

####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################


#######################################
# Shared X-linked TRs between species #
#######################################
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly")

data<-read.table("TRF_all_species_chmX.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
head(data)
data=data[order(data$Counts_s, decreasing = TRUE),]

#Removing duplicate words in "Scaffold_name"
data$Scaffold_name = as.data.frame(sapply(data$Scaffold_name, function(x) gsub("\"", "", x))) #Removing quotes

reduce_row = function(i) {
  split = strsplit(i, split=",")[[1]]
  paste(unique(split), collapse = ",") 
}

data$Scaffold_name = apply(data, 1, reduce_row) #Removing duplicates

#Counts shared monomers
Tps=data.frame(rowSums(sapply(c("Tps_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tps)<-"Tps"
Tdi=data.frame(rowSums(sapply(c("Tdi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tdi)<-"Tdi"
Tcm=data.frame(rowSums(sapply(c("Tcm_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tcm)<-"Tcm"
Tsi=data.frame(rowSums(sapply(c("Tsi_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tsi)<-"Tsi"
Tce=data.frame(rowSums(sapply(c("ScVQC5J_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tce)<-"Tce"
Tpa=data.frame(rowSums(sapply(c("Tpa_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tpa)<-"Tpa"
Tge=data.frame(rowSums(sapply(c("Tge_*"),function(x) grepl(x, data$Scaffold_name))))
colnames(Tge)<-"Tge"

all_sps<-cbind(data[,c(1,2,3,6)], Tps)
all_sps<-cbind(all_sps,Tdi)
all_sps<-cbind(all_sps,Tcm)
all_sps<-cbind(all_sps,Tsi)
all_sps<-cbind(all_sps,Tce)
all_sps<-cbind(all_sps,Tpa)
all_sps<-cbind(all_sps,Tge)
all_sps$total<-all_sps$Tps+all_sps$Tdi+all_sps$Tcm+all_sps$Tsi+all_sps$Tce+all_sps$Tpa+all_sps$Tge

ggplot(all_sps, aes(x=total)) + 
  geom_bar()+
  xlab("Number LGs")

nrow(all_sps)
nrow(subset(all_sps, total=="1"))


#Proportion shared repeats among species (pairwise comparison)
all_sps$Tps_Tcm<-all_sps$Tps+all_sps$Tcm
Tps_Tcm<-subset(all_sps, Tps_Tcm=="2")

all_sps$Tps_Tce<-all_sps$Tps+all_sps$Tce
Tps_Tce<-subset(all_sps, Tps_Tce=="2")

all_sps$Tps_Tpa<-all_sps$Tps+all_sps$Tpa
Tps_Tpa<-subset(all_sps, Tps_Tpa=="2")

all_sps$Tcm_Tce<-all_sps$Tcm+all_sps$Tce
Tcm_Tce<-subset(all_sps, Tcm_Tce=="2")

all_sps$Tcm_Tpa<-all_sps$Tcm+all_sps$Tpa
Tcm_Tpa<-subset(all_sps, Tcm_Tpa=="2")

all_sps$Tce_Tpa<-all_sps$Tce+all_sps$Tpa
Tce_Tpa<-subset(all_sps, Tce_Tpa=="2")

all_sps$Tdi_Tsi<-all_sps$Tdi+all_sps$Tsi
Tdi_Tsi<-subset(all_sps, Tdi_Tsi=="2")

all_sps$Tdi_Tge<-all_sps$Tdi+all_sps$Tge
Tdi_Tge<-subset(all_sps, Tdi_Tge=="2")

all_sps$Tsi_Tge<-all_sps$Tsi+all_sps$Tge
Tsi_Tge<-subset(all_sps, Tsi_Tge=="2")

all_sps$Tps_Tdi<-all_sps$Tps+all_sps$Tdi
Tps_Tdi<-subset(all_sps, Tps_Tdi=="2")

all_sps$Tcm_Tsi<-all_sps$Tcm+all_sps$Tsi
Tcm_Tsi<-subset(all_sps, Tcm_Tsi=="2")

all_sps$Tpa_Tge<-all_sps$Tpa+all_sps$Tge
Tpa_Tge<-subset(all_sps, Tpa_Tge=="2")


sex<-data.frame(nrow(Tps_Tcm)/nrow(data))
colnames(sex)<-"Prop_shared"
sex[nrow(sex) + 1,] = nrow(Tps_Tce)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tps_Tpa)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tcm_Tce)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tcm_Tpa)/nrow(data)
sex[nrow(sex) + 1,] = nrow(Tce_Tpa)/nrow(data)
sex$species<-c("Tps_Tcm", "Tps_Tce", "Tps_Tpa", "Tcm_Tce", "Tcm_Tpa", "Tce_Tpa")
sex$repro<-"sex"

asex<-data.frame(nrow(Tdi_Tsi)/nrow(data))
colnames(asex)<-"Prop_shared"
asex[nrow(asex) + 1,] = nrow(Tdi_Tge)/nrow(data)
asex[nrow(asex) + 1,] = nrow(Tsi_Tge)/nrow(data)
asex$species<-c("Tdi_Tsi", "Tdi_Tge", "Tsi_Tge")
asex$repro<-"asex"

sex_asex<-rbind(sex,asex)
sex_asex$species <- factor(sex_asex$species,levels = c("Tps_Tcm", "Tdi_Tsi", "Tps_Tce", "Tps_Tpa", "Tdi_Tge", "Tcm_Tce", "Tcm_Tpa", "Tsi_Tge", "Tce_Tpa"))

sister_sp<-data.frame(nrow(Tps_Tdi)/nrow(data))
colnames(sister_sp)<-"Prop_shared"
sister_sp[nrow(sister_sp) + 1,] = nrow(Tcm_Tsi)/nrow(data)
sister_sp[nrow(sister_sp) + 1,] = nrow(Tpa_Tge)/nrow(data)
sister_sp$species<-c("Tps_Tdi", "Tcm_Tsi", "Tpa_Tge")
sister_sp$species <- factor(sister_sp$species,levels = c("Tps_Tdi", "Tcm_Tsi", "Tpa_Tge"))

ggplot(data=sex, aes(x=species, y=Prop_shared)) +
  geom_bar(stat="identity")

ggplot(data=asex, aes(x=species, y=Prop_shared)) +
  geom_bar(stat="identity")

ggplot(data=sex_asex, aes(x=species, y=Prop_shared, fill=repro)) +
  geom_bar(stat="identity")

ggplot(data=sister_sp, aes(x=species, y=Prop_shared)) +
  geom_bar(stat="identity")

ggplot(sex_asex, aes(x=repro, y=Prop_shared)) + 
  geom_boxplot()




####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################


##############################
# TR distribution among chms #
##############################

####Tdi
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tdi")

Tdi<-read.table("Tdi_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
Tdi$Repeated_region_size<-Tdi$monomer_length*Tdi$Counts_s
Tdi$Tdi_Count_LGs <- rowSums(sapply(c("\\bTdi_LRv5a_scf1\\b", "\\bTdi_LRv5a_scf2\\b", "\\bTdi_LRv5a_scf3\\b", "\\bTdi_LRv5a_scf4\\b", "\\bTdi_LRv5a_scf5\\b", "\\bTdi_LRv5a_scf6\\b", "\\bTdi_LRv5a_scf7\\b", "\\bTdi_LRv5a_scf8\\b", "\\bTdi_LRv5a_scf9\\b", "\\bTdi_LRv5a_scf10\\b", "\\bTdi_LRv5a_scf11\\b", "\\bTdi_LRv5a_scf12\\b"),
                                    function(x) grepl(x, Tdi$Scaffold_name)))

##LG1
Tdi_LG1<-read.table("Tdi_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tdi_LG1$Repeated_region_size<-Tdi_LG1$monomer_length*Tdi_LG1$Counts

ggplot(Tdi_LG1, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,273000000,25000000), limits=c(-1, 273000000))+
  ylab("Copy number")+ ggtitle("LG1_tdi")

ggplot(Tdi_LG1, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,273000000,25000000), limits=c(-1, 273000000))+
  ylab("Repeat region size")+ ggtitle("LG1_tdi")


##LG3
Tdi_LG3<-read.table("Tdi_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tdi_LG3$Repeated_region_size<-Tdi_LG3$monomer_length*Tdi_LG3$Counts

ggplot(Tdi_LG3, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,20000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_tdi")

ggplot(Tdi_LG3, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,25000000), limits=c(-1, 140000000))+
  ylab("Repeat region size")+ ggtitle("LG3_tdi")



####Tcm
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tcm")
Tcm<-read.table("Tcm_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
Tcm$Repeated_region_size<-Tcm$monomer_length*Tcm$Counts_s
Tcm$Tcm_Count_LGs <- rowSums(sapply(c("\\bTcm_LRv5a_scf1\\b", "\\bTcm_LRv5a_scf2\\b", "\\bTcm_LRv5a_scf3\\b", "\\bTcm_LRv5a_scf4\\b", "\\bTcm_LRv5a_scf5\\b", "\\bTcm_LRv5a_scf6\\b", "\\bTcm_LRv5a_scf7\\b", "\\bTcm_LRv5a_scf8\\b", "\\bTcm_LRv5a_scf9\\b", "\\bTcm_LRv5a_scf10\\b", "\\bTcm_LRv5a_scf11\\b", "\\bTcm_LRv5a_scf12\\b"),
                                    function(x) grepl(x, Tcm$Scaffold_name)))

##LG1
Tcm_LG1.1<-read.table("Tcm_LRv5a_scf1.1.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tcm_LG1.2<-read.table("Tcm_LRv5a_scf1.2.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tcm_LG1<-rbind(Tcm_LG1.1, Tcm_LG1.2)
Tcm_LG1$Repeated_region_size<-Tcm_LG1$monomer_length*Tcm_LG1$Counts

ggplot(Tcm_LG1, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,273000000,25000000), limits=c(-1, 273000000))+
  ylab("Copy number")+ ggtitle("LG1_tcm")

ggplot(Tcm_LG1, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,273000000,25000000), limits=c(-1, 273000000))+
  ylab("Copy number")+ ggtitle("LG1_tcm")


##LG3
Tcm_LG3<-read.table("Tcm_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tcm_LG3$Repeated_region_size<-Tcm_LG3$monomer_length*Tcm_LG3$Counts

ggplot(Tcm_LG3, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,20000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_tcm")

ggplot(Tcm_LG3, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,25000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_tcm")






####Tms
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tms")
Tms<-read.table("Tms_LRv5a.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
Tms$Repeated_region_size<-Tms$monomer_length*Tms$Counts_s
Tms$Tms_Count_LGs <- rowSums(sapply(c("\\bTms_LRv5a_scf1\\b", "\\bTms_LRv5a_scf2\\b", "\\bTms_LRv5a_scf3\\b", "\\bTms_LRv5a_scf4\\b", "\\bTms_LRv5a_scf5\\b", "\\bTms_LRv5a_scf6\\b", "\\bTms_LRv5a_scf7\\b", "\\bTms_LRv5a_scf8\\b", "\\bTms_LRv5a_scf9\\b", "\\bTms_LRv5a_scf10\\b", "\\bTms_LRv5a_scf11\\b", "\\bTms_LRv5a_scf12\\b"),
                                    function(x) grepl(x, Tms$Scaffold_name)))

##LG1
Tms_LG1<-read.table("Tms_LRv5a_scf1.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tms_LG1$Repeated_region_size<-Tms_LG1$monomer_length*Tms_LG1$Counts

ggplot(Tms_LG1, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,165000000,25000000), limits=c(-1, 165000000))+
  ylab("Copy number")+ ggtitle("LG1_Tms")

ggplot(Tms_LG1, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,165000000,25000000), limits=c(-1, 165000000))+
  ylab("Copy number")+ ggtitle("LG1_Tms")


##LG3
Tms_LG3<-read.table("Tms_LRv5a_scf3.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tms_LG3$Repeated_region_size<-Tms_LG3$monomer_length*Tms_LG3$Counts

ggplot(Tms_LG3, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,20000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_Tms")

ggplot(Tms_LG3, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,25000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_Tms")



####Tsi
setwd("~/Desktop/Tandem_repeats/TRF-k-seek_processed/TRF_genome_assembly/Tsi")
Tsi<-read.table("Tsi_LRv5b.fasta.2.7.7.80.10.100.2000_parse_CMwRC_minimal_rotation.txt", header=T, sep="\t", quote="")
Tsi$Repeated_region_size<-Tsi$monomer_length*Tsi$Counts_s
Tsi$Tsi_Count_LGs <- rowSums(sapply(c("\\bTsi_LRv5b_scf1\\b", "\\bTsi_LRv5b_scf2\\b", "\\bTsi_LRv5b_scf3\\b", "\\bTsi_LRv5b_scf4\\b", "\\bTsi_LRv5b_scf5\\b", "\\bTsi_LRv5b_scf6\\b", "\\bTsi_LRv5b_scf7\\b", "\\bTsi_LRv5b_scf8\\b", "\\bTsi_LRv5b_scf9\\b", "\\bTsi_LRv5b_scf10\\b", "\\bTsi_LRv5b_scf11\\b", "\\bTsi_LRv5b_scf12\\b"),
                                    function(x) grepl(x, Tsi$Scaffold_name)))

##LG1
Tsi_LG1<-read.table("Tsi_LRv5b_scf1.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tsi_LG1$Repeated_region_size<-Tsi_LG1$monomer_length*Tsi_LG1$Counts

ggplot(Tsi_LG1, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,270000000,25000000), limits=c(-1, 270000000))+
  ylab("Copy number")+ ggtitle("LG1_Tsi")

ggplot(Tsi_LG1, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,270000000,25000000), limits=c(-1, 270000000))+
  ylab("Copy number")+ ggtitle("LG1_Tsi")


##LG3
Tsi_LG3<-read.table("Tsi_LRv5b_scf3.fasta.2.7.7.80.10.100.2000_parse.txt", header=F, sep="\t", quote="", col.names=c("Scaffold_name", "start", "stop", "monomer_length", "Counts", "monomer_seqs"))
Tsi_LG3$Repeated_region_size<-Tsi_LG3$monomer_length*Tsi_LG3$Counts

ggplot(Tsi_LG3, mapping=aes(start, Counts)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,20000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_Tsi")

ggplot(Tsi_LG3, mapping=aes(start, Repeated_region_size)) +
  stat_summary_bin(fun.y = "sum", geom="bar", bins=1000 - 1) +
  scale_x_continuous(name="LG position", breaks=seq(0,140000000,25000000), limits=c(-1, 140000000))+
  ylab("Copy number")+ ggtitle("LG3_Tsi")

