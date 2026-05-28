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





###Bacillus
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
Brsri<-read.table("Brsri_v3.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Brsri$array_length<-(Brsri$end-Brsri$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Brsri$chr,ranges=IRanges(start=Brsri$start,end=Brsri$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="Brsri_v3_scf1" | seqnames=="Brsri_v3_scf2" | seqnames=="Brsri_v3_scf3" | 
                             seqnames=="Brsri_v3_scf4_1" |seqnames=="Brsri_v3_scf4_2" | seqnames=="Brsri_v3_scf5" | seqnames=="Brsri_v3_scf6" | 
                             seqnames=="Brsri_v3_scf7" | seqnames=="Brsri_v3_scf8" | seqnames=="Brsri_v3_scf9_1" | seqnames=="Brsri_v3_scf9_2" | 
                             seqnames=="Brsri_v3_scf10"| seqnames=="Brsri_v3_scf11" | seqnames=="Brsri_v3_scf12" | seqnames=="Brsri_v3_scf13"|
                             seqnames=="Brsri_v3_scf14"| seqnames=="Brsri_v3_scf15" | seqnames=="Brsri_v3_scf16" | seqnames=="Brsri_v3_scf17"| seqnames=="Brsri_v3_scf18")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr")

table=data.frame(c("chr1", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18"))
table$LG_length <- c(380301637, 135635846, 132210475, 31528091+53763000, 85288357, 84311336, 72962172, 70899808, 22111569+40511000, 60214734, 59965237, 55548340, 49675683, 48856500, 48679402, 47617215, 40534271, 30374980)
table$TR_length <- c(5992114, 3072112, 2659112, 708488+823283, 1557084, 1960474, 1629936, 1435610, 670245+660204, 1251897, 1707820, 1454847, 1109077, 1209938, 1271349, 1340087, 1084923, 998831)
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




###Ischnura elegans
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
Ieleg<-read.table("GCA_921293095.2_ioIscEleg1.2_genomic.fna.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Ieleg$array_length<-(Ieleg$end-Ieleg$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Ieleg$chr,ranges=IRanges(start=Ieleg$start,end=Ieleg$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="OV121100.1" | seqnames=="OV121101.1" | seqnames=="OV121102.1" | 
                             seqnames=="OV121103.1" |seqnames=="OV121104.1" | seqnames=="OV121105.1" | seqnames=="OV121107.1" | 
                             seqnames=="OV121108.1" | seqnames=="OV121109.1" | seqnames=="OV121110.1" | seqnames=="OV121111.1" | 
                             seqnames=="OV121112.1" | seqnames=="OV121113.1" | seqnames=="OV121106.1")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

table=data.frame(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chrX"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chrX"))
table$LG_length <- c(170575982, 147998192, 139038753, 138073150, 137531555, 126002890, 118302689, 118119564, 115518074, 108619349, 103413653, 94743817, 21320372, 123641648)
table$TR_length <- nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`
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




###Locusta migratoria
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
Lmig<-read.table("GCA_026315105.1_CAU_Lmig_1.0_genomic_chr10.fna.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Lmig$array_length<-(Lmig$end-Lmig$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Lmig$chr,ranges=IRanges(start=Lmig$start,end=Lmig$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="CM048744.1" | seqnames=="CM048745.1" | seqnames=="CM048746.1" | 
                             seqnames=="CM048747.1" |seqnames=="CM048748.1" | seqnames=="CM048749.1" | seqnames=="CM048750.1" | 
                             seqnames=="CM048751.1" | seqnames=="CM048753.1" | seqnames=="CM048755.1")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

table=data.frame(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr10", "chrX"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr10", "chrX"))
table$LG_length <- c(991018619, 901651818, 689659731, 581103692, 485803872, 427868946, 418016564, 377406850, 142112905, 806237231)
table$TR_length <- nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`
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





###Panorpa germanica
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
Pger<-read.table("GCA_963678705.1_ijPanGerm4.1_genomic.fna.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Pger$array_length<-(Pger$end-Pger$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Pger$chr,ranges=IRanges(start=Pger$start,end=Pger$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="OY783204.1" | seqnames=="OY783205.1" | seqnames=="OY783206.1" | 
                             seqnames=="OY783207.1" |seqnames=="OY783208.1" | seqnames=="OY783209.1" | seqnames=="OY783210.1" | 
                             seqnames=="OY783211.1" | seqnames=="OY783212.1" | seqnames=="OY783213.1" | seqnames=="OY783214.1" |
                            seqnames=="OY783215.1" | seqnames=="OY783216.1" | seqnames=="OY783217.1" | seqnames=="OY783218.1" |
                           seqnames=="OY783219.1" | seqnames=="OY783220.1" | seqnames=="OY783221.1" | seqnames=="OY783222.1" |
                           seqnames=="OY783223.1" | seqnames=="OY783203.1")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

table=data.frame(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chrX"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                                         "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chrX"))
table$LG_length <- c(30639171, 29108473, 27763522, 25360444, 24800094, 24531805, 22111551, 20879825, 19178891, 17756717, 17094263, 16625943, 16451125, 16304600,
                     15303954, 14495396, 14118734, 14015732, 13756844, 10115498, 47989005)
table$TR_length <- nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`
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




###Laodelphax striatellus
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
Lstr<-read.table("GCA_017141395.1_ASM1714139v1_genomic.fna.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Lstr$array_length<-(Lstr$end-Lstr$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Lstr$chr,ranges=IRanges(start=Lstr$start,end=Lstr$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="CM029530.1" | seqnames=="CM029531.1" | seqnames=="CM029532.1" | 
                             seqnames=="CM029533.1" |seqnames=="CM029534.1" | seqnames=="CM029535.1" | seqnames=="CM029536.1" | 
                             seqnames=="CM029537.1" | seqnames=="CM029538.1" | seqnames=="CM029539.1" | seqnames=="CM029540.1" |
                             seqnames=="CM029541.1" | seqnames=="CM029542.1" | seqnames=="CM029543.1" | seqnames=="CM029544.1")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/minimal_rotations")

table=data.frame(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                   "chr13", "chr14", "chrX"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                                         "chr13", "chr14", "chrX"))
table$LG_length <- c(50771829, 56807035, 59587101, 40643106, 42735337, 34959101, 34869687, 13928342, 30918797, 29457956, 26969898, 20289471, 17750148, 14510125, 31551787)
table$TR_length <- nonoverlapWind_chr_sum$`nonoverlapWind_chr$width`
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




###Myrmecophilus myrmecophilus
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_other_species_Xchr/parsed_files")
Myr<-read.table("Myrmeco.fasta.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
Myr$array_length<-(Myr$end-Myr$start)+1

#Remove overlaps between repeats
myranges<-GRanges(seqnames=Myr$chr,ranges=IRanges(start=Myr$start,end=Myr$end))
nonoverlapWind<-reduce(myranges)
nonoverlapWind<-as.data.frame(nonoverlapWind)

nonoverlapWind_chr<-subset(nonoverlapWind, seqnames=="scaffold_1" | seqnames=="scaffold_2" | seqnames=="scaffold_3" | 
                             seqnames=="scaffold_4" |seqnames=="scaffold_5" | seqnames=="scaffold_6" | seqnames=="scaffold_7" | 
                             seqnames=="scaffold_8" | seqnames=="scaffold_9" | seqnames=="scaffold_10")

nonoverlapWind_chr_sum<-aggregate(nonoverlapWind_chr$width~nonoverlapWind_chr$seqnames, FUN=sum)


#Calculate proportion

table=data.frame(c("chrX", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9"))
colnames (table) [1] <- "LGs"
table$LGs <- factor(table$LGs,levels = c("chrX", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9"))
table$LG_length <- c(170024704, 98987496, 59505047, 42366615, 41412273, 35887947, 35364768, 29775746, 22055277, 11230474)
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


#Stats

cor.test(table$Prop_repeated_region, table$LG_length, method = "spearman")
cor.test(table$Prop_repeated_region, table$LG_length, method = "pearson")

model<-lm(table$Prop_repeated_region~table$LG_length)
summary(model)

