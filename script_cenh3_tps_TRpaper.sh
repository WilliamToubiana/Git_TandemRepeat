#!/bin/bash


#Slurm options:

#SBATCH -p cpu
#SBATCH --time=0-1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB


#Commands start here:

echo "Starting job $SLURM_JOB_NAME with ID $SLURM_JOB_ID".

module load gcc
module load bwa
module load perl
module load bedtools2
module load samtools
module load python

########################
## Prepare chip reads ##
########################

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/

cat *cenh3*R1*.gz > Tps_testes_cenh3_1_R1.fq.gz
cat *cenh3*R3*.gz > Tps_testes_cenh3_1_R3.fq.gz

cat *input*R1*.gz > Tps_testes_input_1_R1.fq.gz
cat *input*R3*.gz > Tps_testes_input_1_R3.fq.gz



################
## Trim reads ##
################

module load trimmomatic

for i in *_R1.fq.gz ; do
        foo1=`echo $i`
                basename=`echo $foo1 | sed 's/_R1.fq.gz*//' | sed 's/.*\///'`
        infileR1=`echo $foo1`
        infileR2=`echo $foo1 | sed 's/_R1.fq.gz/_R3.fq.gz/'`
        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq"`
        outfileR2=`echo "./"$basename"_R3_qtrimmed.fq"`
        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq"`
        outfileR2_UP=`echo "./"$basename"_R3_qtrimmed_UNPAIRED.fq"`

#        echo $infileR1
#        echo $infileR2
#        echo $outfileR1
#        echo $outfileR1_UP
#        echo $outfileR2
#        echo $outfileR2_UP

        trimmomatic PE -threads 16 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
done




###################
## Mapping reads ##
###################

	## bwa index


cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip

bwa index genomes/Tps_LRv5b_mtDNAv350.fasta


	## bwa map
mkdir bwa

bwa mem -t 16 -c 1000000000 -T 30 genomes/Tps_LRv5b_mtDNAv350.fasta Tps_testes_cenh3_1_R1_qtrimmed.fq Tps_testes_cenh3_1_R3_qtrimmed.fq > bwa/chip_cenh3_tps_testes_1_bwa.sam
bwa mem -t 16 -c 1000000000 -T 30 genomes/Tps_LRv5b_mtDNAv350.fasta chip/Tps_testes_input_1_R1_qtrimmed.fq chip/Tps_testes_input_1_R3_qtrimmed.fq > bwa/chip_input_tps_testes_1_bwa.sam


	## flagstat

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip/bwa

samtools flagstat chip_cenh3_tps_testes_1_bwa.sam > chip_cenh3_tps_testes_1_bwa_flagstat.txt
samtools flagstat chip_input_tps_testes_1_bwa.sam > chip_input_tps_testes_1_bwa_flagstat.txt


	## sort bam

samtools view -u chip_cenh3_tps_testes_1_bwa.sam | samtools sort -o chip_cenh3_tps_testes_1_bwa.bam
samtools view -u chip_input_tps_testes_1_bwa.sam | samtools sort -o chip_input_tps_testes_1_bwa.bam

#rm chip_cenh3_tps_testes_1_bwa.sam
#rm chip_input_tps_testes_1_bwa.sam


	## remove supp (chimeric) alignements

samtools view chip_cenh3_tps_testes_1_bwa.bam | fgrep SA:Z: | cut -f 1 > chip_cenh3_tps_testes_1_bwa_badnames.txt
samtools view -h chip_cenh3_tps_testes_1_bwa.bam | fgrep -vf chip_cenh3_tps_testes_1_bwa_badnames.txt | samtools view -b > chip_cenh3_tps_testes_1_bwa_final.bam
#rm chip_cenh3_tps_testes_1_bwa_badnames.txt
samtools flagstat chip_cenh3_tps_testes_1_bwa_final.bam > chip_cenh3_tps_testes_1_bwa_final_flagstat.txt


samtools view chip_input_tps_testes_1_bwa.bam | fgrep SA:Z: | cut -f 1 > chip_input_tps_testes_1_bwa_badnames.txt
samtools view -h chip_input_tps_testes_1_bwa.bam | fgrep -vf chip_input_tps_testes_1_bwa_badnames.txt | samtools view -b > chip_input_tps_testes_1_bwa_final.bam
#rm chip_input_tps_testes_1_bwa_badnames.txt
samtools flagstat chip_input_tps_testes_1_bwa_final.bam > chip_input_tps_testes_1_bwa_final_flagstat.txt


	##remove PCR duplicates
module load picard

for i in  *tps_testes_1*.bam; do
    outbam=`echo $i | sed 's/_bwa_final.bam/_bwa_final_DR.bam/'`
       flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
       metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

#       echo $i
#       echo $outbam
#       echo $metric_file
#       echo $flagstat_out_bam

       picard MarkDuplicates REMOVE_DUPLICATES=true \
       INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file

       samtools flagstat $outbam > $flagstat_out_bam
       mv $flagstat_out_bam /scratch/wtoubian/bwa_timema_genomes/mapping_genomes/BWA_out/flagstat_out_paired

done





#########################################
## Calculate coverage GW (250 kb bins) ##
#########################################


cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/wtoubian/chip
mkdir coverage

bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b bwa/chip_cenh3_tps_testes_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5a_mtDNAv350.fasta.fai -mean > coverage/Tps_cenh3_testes_1_GW_coverage_DR_250kb.txt
bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b bwa/chip_input_tps_testes_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5a_mtDNAv350.fasta.fai -mean > coverage/Tps_input_testes_1_GW_coverage_DR_250kb.txt



