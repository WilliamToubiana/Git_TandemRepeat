#!/bin/bash


#Slurm options:

#SBATCH -p cpu
#SBATCH --time=0-10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB


#Commands start here:

echo "Starting job $SLURM_JOB_NAME with ID $SLURM_JOB_ID".

module load gcc
module load perl
module load picard
module load bedtools2
module load samtools
module load python
module load trimmomatic



##########################################################
## Prepare genome assemblies (example with T. poppense) ##
##########################################################

#cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/genomes

#samtools faidx Tps_LRv5b_mtDNAv350.fasta

#cut -f1,2 Tps_LRv5b_mtDNAv350.fasta.fai > Tps_chm_size_mtDNAv350.txt

#bedtools makewindows -g Tps_chm_size_mtDNAv350.txt -w 250000  > Tps_chm_size_mtDNAv350_w250000.bed


################################################
## Remove reads that failed Casava from files ##
################################################

#### explanation - so reads that failed basic tests at the seq centre were not removed
##### but insteads just flagged with a 'Y' for failed and an 'N' for passed: e.g.

# @WINDU:16:C5UFBANXX:7:1310:11475:72909 1:N:0:GAGTGG ## passed
# @WINDU:16:C5UFBANXX:7:1310:11694:72820 1:Y:0:GAGTGG ## failed

# also seq company breaks reads up into parts of a file so they need cat-ing



### used something like this to get reads that passed - note also need to remove weird -- lines that appear from grepping

cd /scratch/wtoubian/Illumina_reads

#zcat Tps/Tps/ReSeq_Ps14*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind2_R1.fastq 2> Tps-filtered_ind2_R1.err
#zcat Tps/Tps/ReSeq_Ps14*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind2_R2.fastq 2> Tps-filtered_ind2_R2.err
#zcat Tps/Tps/ReSeq_Ps16*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind3_R1.fastq 2> Tps-filtered_ind3_R1.err
#zcat Tps/Tps/ReSeq_Ps16*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind3_R2.fastq 2> Tps-filtered_ind3_R2.err
#zcat Tps/Tps/ReSeq_Ps18*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind4_R1.fastq 2> Tps-filtered_ind4_R1.err #not included
#zcat Tps/Tps/ReSeq_Ps18*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind4_R2.fastq 2> Tps-filtered_ind4_R2.err #not included
#zcat Tps/is_550/*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind1_R1.fastq 2> Tps-filtered_ind1_R1.err
#zcat Tps/is_550/*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tps/Tps-filtered_ind1_R2.fastq 2> Tps-filtered_ind1_R2.err

#zcat Tcm/is_550/*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind1_R1.fastq 2> Tcm-filtered_R1_ind1.err
#zcat Tcm/is_550/*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind1_R2.fastq 2> Tcm-filtered_R2_ind1.err
#zcat Tcm/Tcm/HM217*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind2_R1.fastq 2> Tcm-filtered_R1_ind2.err
#zcat Tcm/Tcm/HM217*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind2_R2.fastq 2> Tcm-filtered_R2_ind2.err
#zcat Tcm/Tcm/HM218*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind3_R1.fastq 2> Tcm-filtered_R1_ind3.err
#zcat Tcm/Tcm/HM218*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind3_R2.fastq 2> Tcm-filtered_R2_ind3.err
#zcat Tcm/Tcm/HM219*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind4_R1.fastq 2> Tcm-filtered_R1_ind4.err
#zcat Tcm/Tcm/HM219*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind4_R2.fastq 2> Tcm-filtered_R2_ind4.err
#zcat Tcm/Tcm/HM220*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind5_R1.fastq 2> Tcm-filtered_R1_ind5.err
#zcat Tcm/Tcm/HM220*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind5_R2.fastq 2> Tcm-filtered_R2_ind5.err
#zcat Tcm/Tcm/HM221*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind6_R1.fastq 2> Tcm-filtered_R1_ind6.err
#zcat Tcm/Tcm/HM221*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tcm/Tcm-filtered_ind6_R2.fastq 2> Tcm-filtered_R2_ind6.err

#zcat Tpa/is_550/*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind1_R1.fastq 2> Tpa-filtered_R1_ind1.err
#zcat Tpa/is_550/*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind1_R2.fastq 2> Tpa-filtered_R2_ind1.err
#zcat Tpa/Tpa/Pa_AB*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind2_R1.fastq 2> Tpa-filtered_R1_ind2.err
#zcat Tpa/Tpa/Pa_AB*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind2_R2.fastq 2> Tpa-filtered_R2_ind2.err
#zcat Tpa/Tpa/PA_CD*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind3_R1.fastq 3> Tpa-filtered_R1_ind3.err
#zcat Tpa/Tpa/PA_CD*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind3_R2.fastq 3> Tpa-filtered_R2_ind3.err
#zcat Tpa/Tpa/PA_E*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind4_R1.fastq 4> Tpa-filtered_R1_ind4.err
#zcat Tpa/Tpa/PA_E*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind4_R2.fastq 4> Tpa-filtered_R2_ind4.err
#zcat Tpa/Tpa/H54*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind5_R1.fastq 5> Tpa-filtered_R1_ind5.err
#zcat Tpa/Tpa/H54*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind5_R2.fastq 5> Tpa-filtered_R2_ind5.err
#zcat Tpa/Tpa/H56*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind6_R1.fastq 6> Tpa-filtered_R1_ind6.err
#zcat Tpa/Tpa/H56*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tpa/Tpa-filtered_ind6_R2.fastq 6> Tpa-filtered_R2_ind6.err

#zcat Tce/is_550/*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind1_R1.fastq 2> Tce-filtered_R1_ind1.err
#zcat Tce/is_550/*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind1_R2.fastq 2> Tce-filtered_R2_ind1.err
#zcat Tce/Tce/CC22B*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind2_R1.fastq 2> Tce-filtered_R1_ind2.err
#zcat Tce/Tce/CC22B*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind2_R2.fastq 2> Tce-filtered_R2_ind2.err
#zcat Tce/Tce/CC22C*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind3_R1.fastq 2> Tce-filtered_R1_ind3.err
#zcat Tce/Tce/CC22C*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind3_R2.fastq 2> Tce-filtered_R2_ind3.err
#zcat Tce/Tce/CC24B*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind4_R1.fastq 2> Tce-filtered_R1_ind4.err
#zcat Tce/Tce/CC24B*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind4_R2.fastq 2> Tce-filtered_R2_ind4.err
#zcat Tce/Tce/CC24C*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind5_R1.fastq 2> Tce-filtered_R1_ind5.err
#zcat Tce/Tce/CC24C*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind5_R2.fastq 2> Tce-filtered_R2_ind5.err
#zcat Tce/Tce/CC25B*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind6_R1.fastq 2> Tce-filtered_R1_ind6.err
#zcat Tce/Tce/CC25B*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tce/Tce-filtered_ind6_R2.fastq 2> Tce-filtered_R2_ind6.err

#zcat Tbi/is_550/*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind1_R1.fastq 2> Tbi-filtered_R1_ind1.err
#zcat Tbi/is_550/*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind1_R2.fastq 2> Tbi-filtered_R2_ind1.err
#zcat Tbi/Tbi/CC86B*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind2_R1.fastq 2> Tbi-filtered_R1_ind2.err
#zcat Tbi/Tbi/CC86B*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind2_R2.fastq 2> Tbi-filtered_R2_ind2.err
#zcat Tbi/Tbi/CC86C*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind3_R1.fastq 2> Tbi-filtered_R1_ind3.err
#zcat Tbi/Tbi/CC86C*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind3_R2.fastq 2> Tbi-filtered_R2_ind3.err
#zcat Tbi/Tbi/CC87B*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind4_R1.fastq 2> Tbi-filtered_R1_ind4.err
#zcat Tbi/Tbi/CC87B*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind4_R2.fastq 2> Tbi-filtered_R2_ind4.err
#zcat Tbi/Tbi/CC87C*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind5_R1.fastq 2> Tbi-filtered_R1_ind5.err
#zcat Tbi/Tbi/CC87C*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind5_R2.fastq 2> Tbi-filtered_R2_ind5.err
#zcat Tbi/Tbi/CC88B*R1*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind6_R1.fastq 2> Tbi-filtered_R1_ind6.err
#zcat Tbi/Tbi/CC88B*R2*.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > Tbi/Tbi-filtered_ind6_R2.fastq 2> Tbi-filtered_R2_ind6.err

##fastq.gz for indG are in /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/Illumina_reads_Will/

##############################################
## read trimming (example with T. poppense) ##
##############################################


#for i in Tps*_R1.fastq ; do
#        foo1=`echo $i`
#                basename=`echo $foo1 | sed 's/_R1.fastq//' | sed 's/.*\///'`
#        infileR1=`echo $foo1`
#        infileR2=`echo $foo1 | sed 's/_R1.fastq/_R2.fastq/'`
#        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq"`
#        outfileR2=`echo "./"$basename"_R2_qtrimmed.fq"`
#        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq"`
#        outfileR2_UP=`echo "./"$basename"_R2_qtrimmed_UNPAIRED.fq"`

#        echo $infileR1
#        echo $infileR2
#        echo $outfileR1
#        echo $outfileR1_UP
#        echo $outfileR2
#        echo $outfileR2_UP

#        trimmomatic PE -threads 16 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
#done






#############################################
## Read mapping (example with T. poppense) ##
#############################################

##change read file names and genome file names accoding to species


cd /scratch/wtoubian
module load minimap2

#minimap2 -ax sr /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/genomes/bwa_index/Tps_LRv5b_mtDNAv350.fasta Tps-filtered_ind1_R1_qtrimmed.fq Tps-filtered_ind1_R2_qtrimmed.fq -t 48 > Tps-filtered_ind1_to_Tps_genome_pe_minimap2.sam
#minimap2 -ax sr /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/genomes/bwa_index/Tps_LRv5b_mtDNAv350.fasta Tps-filtered_ind2_R1_qtrimmed.fq Tps-filtered_ind2_R2_qtrimmed.fq -t 48 > Tps-filtered_ind2_to_Tps_genome_pe_minimap2.sam
#minimap2 -ax sr /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/genomes/bwa_index/Tps_LRv5b_mtDNAv350.fasta Tps-filtered_ind3_R1_qtrimmed.fq Tps-filtered_ind3_R2_qtrimmed.fq -t 48 > Tps-filtered_ind3_to_Tps_genome_pe_minimap2.sam
#minimap2 -ax sr /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/genomes/bwa_index/Tps_LRv5b_mtDNAv350.fasta Tps-filtered_ind4_R1_qtrimmed.fq Tps-filtered_ind4_R2_qtrimmed.fq -t 48 > Tps-filtered_ind4_to_Tps_genome_pe_minimap2.sam






##### remove supp reads
# I filter out Supplemental alignments (secondary alignments are already removed from the bwa output). These are chimeric alignments
# A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
# They have an SA:Z tag. If these are paired reads I also remove the other pair

cd /scratch/wtoubian/

#for i in Tps-filtered*.sam; do
#    outflagstat=`echo $i | sed 's/.sam/_flagstat.txt/'`
#    outbam=`echo $i | sed 's/.sam/.bam/'`
#    badnames=`echo $i | sed 's/.sam/_badnames.txt/'`
#    outbamF=`echo $i | sed 's/.sam/_final.bam/'`
#    outflagstatF=`echo $i | sed 's/.sam/_final_flagstat.txt/'`


#       echo $i
#       echo $outflagstat

#       samtools flagstat $i > $outflagstat     #flagstat
#       samtools view -u $i | samtools sort -o $outbam  #sort bam
#       samtools view $outbam | fgrep SA:Z: | cut -f 1 > $badnames      #remove supp (chimeric) alignments
#       samtools view -h $outbam | fgrep -vf $badnames | samtools view -b > $outbamF    #remove supp (chimeric) alignments
#       samtools flagstat $outbamF > $outflagstatF      #flagstat


#done







        ##remove PCR duplicates

cd /scratch/wtoubian/

#for i in Tps-filtered_ind*_final.bam; do
#    outbam=`echo $i | sed 's/_final.bam/_final_DR.bam/'`
#       flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
#       metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

#       echo $i
#       echo $outbam
#       echo $metric_file
#       echo $flagstat_out_bam

#       picard MarkDuplicates REMOVE_DUPLICATES=true \
#       INPUT=$i \
#    OUTPUT=$outbam \
#    METRICS_FILE=$metric_file

#       samtools flagstat $outbam > $flagstat_out_bam

#done







#######################################################################
## Per-base coverage and sum (example with T. poppense individual 1) ##
#######################################################################

##bed files were generated from the R script "TRF_parsing.R"


#sortBed -faidx genomes/Tps_chm_size_mtDNAv350.txt -i Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse.bed > Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted.bed #sort bed file

#awk '{print $1 "\t" $2 "\t" $3;}' Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted.bed > Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_sorted_noSeq.bed #remove "motif" column


cd /scratch/wtoubian/

#samtools depth -a Tps-filtered_ind1_to_Tps_genome_pe_minimap2_final_DR.bam > Tps-filtered_ind1_minimap2_GW_coverage.txt
#samtools depth -a -b Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies_sorted_noSeq.bed Tps-filtered_ind1_to_Tps_genome_pe_minimap2_final_DR.bam > Tps-filtered_ind1_minimap2_TR_coverage.txt
#awk 'END { print s } { s += $3 }' Tps-filtered_ind1_minimap2_GW_coverage.txt
#awk 'END { print s } { s += $3 }' Tps-filtered_ind1_minimap2_TR_coverage.txt





#########################################################################
## Sum coverage 250kb windows (example with T. poppense individual 1) ##
#########################################################################

#awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_ind1_minimap2_GW_coverage.txt Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_ind1_minimap2_GW_sum_coverage_250kb.txt

#awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_ind1_minimap2_TR_coverage.txt Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_ind1_minimap2_TR_sum_coverage_250kb.txt
