
library( data.table)
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/turnover_orthologs")

##############################
### Exon sequence turnover ###
##############################

CDS<-fread("WillOne2one.CDS.v2.txt")
colnames(CDS)<-c("motif", "chr", "start", "end", "copy_nb", "Tps","Tdi", "Tcm", "Tce", "Tpa", "Tbi", "shared")
CDS$shared<-rowSums(CDS[,c(7,8,9,10,11)])
Tps_Tdi<-subset(CDS, Tdi=="1")
Tps_Tcm<-subset(CDS, Tcm=="1")
Tps_Tce<-subset(CDS, Tce=="1")
Tps_Tpa<-subset(CDS, Tpa=="1")
Tps_Tbi<-subset(CDS, Tbi=="1")

library(venn)
library(ggVennDiagram)
venn <- list(Tps=CDS$motif, Other_Timema=Tps_Tdi$motif)
ggVennDiagram(venn)


################################
### Intron sequence turnover ###
################################

introns<-fread("WillOne2one.intron.v2.txt")
colnames(introns)<-c("motif", "chr", "start", "end", "copy_nb", "Tps","Tdi", "Tcm", "Tce", "Tpa", "Tbi", "shared")
introns$shared<-rowSums(introns[,c(7,8,9,10,11)])
Tps_Tdi<-subset(introns, Tdi=="1")
Tps_Tcm<-subset(introns, Tcm=="1")
Tps_Tce<-subset(introns, Tce=="1")
Tps_Tpa<-subset(introns, Tpa=="1")
Tps_Tbi<-subset(introns, Tbi=="1")

library(venn)
library(ggVennDiagram)
venn <- list(Tps=introns$motif, Other_Timema=Tps_Tbi$motif)
ggVennDiagram(venn)



############################
### TR sequence turnover ###
############################
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/unique_motif_sets")
TR<-read.table("Tps_shared_TR_motifs.txt", header=T, sep="\t", quote="") #In this dataframe, only the shortest TR sequence was kept when TR sequences had the same start or end position in T. poppense TR annotation.
Tps_Tdi<-subset(TR, Tdi=="1")
Tps_Tcm<-subset(TR, Tcm=="1")
Tps_Tce<-subset(TR, Tce=="1")
Tps_Tpa<-subset(TR, Tpa=="1")
Tps_Tbi<-subset(TR, Tbi=="1")

library(venn)
library(ggVennDiagram)
venn <- list(Tps=TR$motif, Other_Timema=Tps_Tdi$motif)
ggVennDiagram(venn)






