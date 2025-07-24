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



###############################################################
## Tandem repeat annotation on Tdi genome assembly using TRF ##
###############################################################

module load gcc
module load trf
module load r

#run TRF and generate gff3 output file
trf Tps_LRv5b_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h #T. poppense
trf Tdi_LRv5a_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h #T. douglasi
trf Tcm_LRv5a_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h #T. californicum
trf Tce_LRv5a_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h #T. cristinae
trf Tpa_LRv5a_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h #T. podura
trf Tbi_LRv4a_mtDNAv350.fasta 2 7 7 80 10 50 2000 -d -m -h #T. bartmani

trf GCA_026315105.1_CAU_Lmig_1.0_genomic.fna 2 7 7 80 10 50 2000 -d -m -h #L. migratoria #TRF analysis had to be run separately for chr X, 8 and 10 because of time limitations on the cluster
trf GCA_921293095.1_ioIscEleg1.1_genomic.fna 2 7 7 80 10 50 2000 -d -m -h #I. elegans
trf Brsri_v3.fasta 2 7 7 80 10 50 2000 -d -m -h #B. rossius


trf2gff -i Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i Tcm_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i Tce_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i Tpa_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i Tbi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)

trf2gff -i GCA_026315105.1_CAU_Lmig_1.0_genomic.fna.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i GCA_921293095.1_ioIscEleg1.1_genomic.fna.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)
trf2gff -i Brsri_v3.fasta.2.7.7.80.10.50.2000.dat #convert .dat to .gff3 file for parsing (installed from https://github.com/Adamtaranto/TRF2GFF)



#Parse TRF
./TRF_parsing_TR_paper.R #this script was used to parse TRF output Tdi_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000.gff3 file and keep repeat arrays of at least 5 copies
#An example of the files subsequently generated is provided for T. poppense
#Tps_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.txt file was subsequently generated
#Tps_LRv5a_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed file was subsequently generated
