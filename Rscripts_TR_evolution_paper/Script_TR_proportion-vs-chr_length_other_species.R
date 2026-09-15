################################################################
###Tandem repeat proportion between chromosomes other species###
################################################################

library(GenomicRanges)

###T. californicum
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tcm<-read.table("Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tcm$array_length<-(Tcm$end-Tcm$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tcm$chr,ranges=IRanges(start=Tcm$start,end=Tcm$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tcm_LRv5a_scf1.1" | seqnames=="Tcm_LRv5a_scf1.2" | seqnames=="Tcm_LRv5a_scf2" | 
                             seqnames=="Tcm_LRv5a_scf3" |seqnames=="Tcm_LRv5a_scf4" | seqnames=="Tcm_LRv5a_scf5.1" | seqnames=="Tcm_LRv5a_scf5.2" | 
                             seqnames=="Tcm_LRv5a_scf6.1" | seqnames=="Tcm_LRv5a_scf6.2" | seqnames=="Tcm_LRv5a_scf7" | seqnames=="Tcm_LRv5a_scf8.1" | 
                             seqnames=="Tcm_LRv5a_scf8.2"| seqnames=="Tcm_LRv5a_scf8.3" | seqnames=="Tcm_LRv5a_scf9" | seqnames=="Tcm_LRv5a_scf10"|
                             seqnames=="Tcm_LRv5a_scf11.1"| seqnames=="Tcm_LRv5a_scf11.2" | seqnames=="Tcm_LRv5a_scf12")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr")

table=data.frame(c("chr1", "chr2", "chrX", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chrX", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))
table$LG_length <- c(23198153+226161097, 177478266, 137281248, 97644469, 3754722+78996087, 4779763+77527622, 79148075, 4715686+3623685+65511858, 72138168, 68730460, 5143100+61968689, 40193497)
table$TR_length <- c(3855880+14535310, 13348953, 6388414, 10441058, 1242988+8511518, 1314250+8690794, 10359313, 1293001+846244+6602734, 8237877, 8769087, 1514396+6829676, 5550062)
table$Prop_repeated_region<-table$TR_length/table$LG_length


ggplot(data=table, aes(x=LG_length, y=Prop_repeated_region)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')


ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()



#Stats

cor.test(table$Prop_repeated_region, table$LG_length, method = "spearman")
cor.test(table$Prop_repeated_region, table$LG_length, method = "pearson")

model<-lm(table$Prop_repeated_region~table$LG_length)
summary(model)





###T. podura
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Tpa<-read.table("Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Tpa$array_length<-(Tpa$end-Tpa$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Tpa$chr,ranges=IRanges(start=Tpa$start,end=Tpa$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Tpa_LRv5a_scf1" | seqnames=="Tpa_LRv5a_scf2" | 
                             seqnames=="Tpa_LRv5a_scf3" |seqnames=="Tpa_LRv5a_scf4" | seqnames=="Tpa_LRv5a_scf5" | 
                             seqnames=="Tpa_LRv5a_scf6" | seqnames=="Tpa_LRv5a_scf7" | seqnames=="Tpa_LRv5a_scf8" | 
                             seqnames=="Tpa_LRv5a_scf9" | seqnames=="Tpa_LRv5a_scf10"|
                             seqnames=="Tpa_LRv5a_scf11"| seqnames=="Tpa_LRv5a_scf12" | seqnames=="Tpa_LRv5a_scf13" | seqnames=="Tpa_LRv5a_scf14")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr")

table=data.frame(c("chrX", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chrX", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"))
table$LG_length <- c(128341569, 128200689, 88008444, 86604445, 77135086, 75690447, 75422829, 71837933, 69097152, 69064364, 68537582, 65766924, 62896878, 38416574)
table$TR_length <- c(4190518, 6723599, 3795640, 6044234, 5524117, 5621056, 5536766, 5337683, 4941944, 3816059, 4455838, 5380789, 4142877, 3972271)
table$Prop_repeated_region<-table$TR_length/table$LG_length


ggplot(data=table, aes(x=LG_length, y=Prop_repeated_region)) +
  geom_point(stat="identity", color="black", fill="white")+
  geom_text(aes(label=LGs), vjust=1.6, size=5.5)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic(base_size = 17)+
  geom_smooth(method='lm')


ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()



#Stats

cor.test(table$Prop_repeated_region, table$LG_length, method = "spearman")
cor.test(table$Prop_repeated_region, table$LG_length, method = "pearson")

model<-lm(table$Prop_repeated_region~table$LG_length)
summary(model)

