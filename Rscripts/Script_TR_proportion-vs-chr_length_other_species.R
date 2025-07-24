################################################################
###Tandem repeat proportion between chromosomes other species###
################################################################

###Bacillus
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
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
  geom_text(aes(label=LGs), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()


ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()


###Ischnura elegans
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_annotation_timema/parsed_files")
Ieleg<-read.table("GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000_parse_5copies.txt", header=T, sep="\t", quote="")
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
table$LGs <- factor(table$LGs,levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chrX", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13"))
table$LG_length <- c(170575982, 147998192, 139038753, 138073150, 137531555, 126002890, 118302689, 118119564, 115518074, 108619349, 103413653, 94743817, 21320372, 123641648)
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
  geom_text(aes(label=LGs), vjust=1.6, size=4)+ ylab("Proportion tandem repeats")+ xlab("chromosome size")+
  theme_classic()

ggplot(data=table, aes(x=LGs, y=Prop_repeated_region)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_text(aes(label=round(Prop_repeated_region, digits = 3)), vjust=1.6, size=3)+ ylab("Proportion tandem repeats")+
  theme_classic()




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
