#!/usr/bin/env Rscript

library(stringr)


###Tdi
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/Tdi_LRv5a")
TR_db_GW<-read.delim("Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)



###Tps
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/Tps_LRv5b")
TR_db_GW<-read.delim("Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)


###Tcm
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/Tcm_LRv5a")
TR_db_GW<-read.delim("Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)


###Tce
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/Tce_LRv5a")
TR_db_GW<-read.delim("Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)


###Tpa
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/Tpa_LRv5a")
TR_db_GW<-read.delim("Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)


###Tbi
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/Tbi_LRv4a")
TR_db_GW<-read.delim("Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Tbi_LRv4a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)



###Brossius
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr")
TR_db_GW<-read.delim("Brsri_v3.fasta.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "Brsri_v3.fasta.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "Brsri_v3.fasta.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "Brsri_v3.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "Brsri_v3.fasta.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)



###Locustra migratoria
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr")

TR_db_GW1<-read.delim("GCA_026315105.1_CAU_Lmig_1.0_genomic.fna.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW1$V9, "[;=]", 20))
TR_db_GW1_1<-cbind(TR_db_GW1[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW1_1)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW1_subset<-subset(TR_db_GW1_1, copy_nb >= 5)

TR_db_GW2<-read.delim("GCA_026315105.1_CAU_Lmig_1.0_genomic_chr8.fna.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split2 <- data.frame(str_split_fixed(TR_db_GW2$V9, "[;=]", 20))
TR_db_GW2_1<-cbind(TR_db_GW2[c(1,4,5)], TR_split2[c(6,8,18)])
colnames(TR_db_GW2_1)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW2_subset<-subset(TR_db_GW2_1, copy_nb >= 5)

TR_db_GW3<-read.delim("GCA_026315105.1_CAU_Lmig_1.0_genomic_chr10.fna.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split3 <- data.frame(str_split_fixed(TR_db_GW3$V9, "[;=]", 20))
TR_db_GW3_1<-cbind(TR_db_GW3[c(1,4,5)], TR_split3[c(6,8,18)])
colnames(TR_db_GW3_1)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW3_subset<-subset(TR_db_GW3_1, copy_nb >= 5)

TR_db_GW4<-read.delim("GCA_026315105.1_CAU_Lmig_1.0_genomic_chrX.fna.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split4 <- data.frame(str_split_fixed(TR_db_GW4$V9, "[;=]", 20))
TR_db_GW4_1<-cbind(TR_db_GW4[c(1,4,5)], TR_split4[c(6,8,18)])
colnames(TR_db_GW4_1)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW4_subset<-subset(TR_db_GW4_1, copy_nb >= 5)

TR_db_GW_subset<-rbind(TR_db_GW1_subset, TR_db_GW2_subset, TR_db_GW3_subset, TR_db_GW4_subset)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
write.table(TR_db_GW_subset, file = "GCA_026315105.1_CAU_Lmig_1.0_genomic_chr10.fna.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "GCA_026315105.1_CAU_Lmig_1.0_genomic_chr10.fna.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)



###Ischnura elegans
##Parse TR gff3 file
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr")
TR_db_GW<-read.delim("GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000.gff3", header=F, comment.char="#")
TR_split1 <- data.frame(str_split_fixed(TR_db_GW$V9, "[;=]", 20))

TR_db_GW<-cbind(TR_db_GW[c(1,4,5)], TR_split1[c(6,8,18)])
colnames(TR_db_GW)<-c("chr", "start", "end", "copy_nb", "motif_length","motif")
TR_db_GW_subset<-subset(TR_db_GW, copy_nb >= 5)

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
write.table(TR_db_GW, file = "GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000_parse.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW[c(1,2,3,6)], file = "GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000_parse.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset, file = "GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000_parse_5copies.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(TR_db_GW_subset[c(1,2,3,6)], file = "GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000_parse_5copies.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#These tables are further processed on the cluster (see script_tandem-repeat_annotation.sh)


















