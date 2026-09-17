#!/bin/bash
#SBATCH --account=tschwand_default
#SBATCH --partition cpu
#SBATCH --job-name=new_prop_TR_perwindow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40G
#SBATCH --output=new_prop_TR_perwindow.out
#SBATCH --error=new_prop_TR_perwindow.err
#SBATCH --mail-user=Joao.Souto@unil.ch
#SBATCH --mail-type=ALL
#SBATCH --time=15:00:00

#sbatch new_prop_TR_perwindow.sh

# Tps to Tcm

#conda activate cactus_py3



export PATH=/users/jsouto/software/cactus-bin-v2.8.4/bin:$PATH

source /users/jsouto/software/cactus-bin-v2.8.4/venv-cactus-v2.8.4/bin/activate











IMPDIR=/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/gifts_from_vincent

OUTDIR=/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/WGA

WINDOWDIR=/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/windows

ANNOTATIONDIR=/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/annotation

SP_SOURCE=TPS

SP_TARGET=TCM





## use bedtols subtract, to remove from 250kb everything that is not repeats

if [[ ! -s ${WINDOWDIR}/Tps_Tr.per250000w.sorted.bed  ]]; then

cut -f 1,2,3 ${ANNOTATIONDIR}/Tps_LRv5b_mtDNAv350.fasta.2.7.7.80.10.50.2000_parse_5copies.bed | grep -v "mtDNA" | bedtools sort | bedtools merge | bedtools sort > ${ANNOTATIONDIR}/Tps_Tr.sorted.bed

grep -v "mtDNA" ${WINDOWDIR}/Tps_chm_size_mtDNAv350_w250000.bed | bedtools sort > ${WINDOWDIR}/Tps_w250000.sorted.bed

bedtools subtract -a ${WINDOWDIR}/Tps.fasta.fai.sorted.bed -b ${ANNOTATIONDIR}/Tps_Tr.sorted.bed | bedtools sort > ${WINDOWDIR}/Tps_everythingexcept_Tr.per250000w.sorted.bed

bedtools subtract -a ${WINDOWDIR}/Tps_w250000.sorted.bed -b ${WINDOWDIR}/Tps_everythingexcept_Tr.per250000w.sorted.bed | bedtools sort > ${WINDOWDIR}/Tps_Tr.per250000w.sorted.bed 

fi

## liftOver the TR windows to the target species (TCM)

if [[ ! -s ${OUTDIR}/${SP_SOURCE}_TO_${SP_TARGET}.TR.per250000w.psl  ]]; then

/users/jsouto/software/cactus-bin-v2.8.4/bin/halLiftover --noDupes --outPSL ${IMPDIR}/Timemas_All.hal ${SP_SOURCE} ${WINDOWDIR}/Tps_Tr.per250000w.sorted.bed ${SP_TARGET} ${OUTDIR}/${SP_SOURCE}_TO_${SP_TARGET}.TR.per250000w.psl

fi

awk '{print $10"\t"$12"\t"$13"\t"$1"\t"$2}' ${OUTDIR}/${SP_SOURCE}_TO_${SP_TARGET}.TR.per250000w.psl | sort -k1.14n -k2,2n -k3,3n > ${OUTDIR}/${SP_SOURCE}_TO_${SP_TARGET}.TR.per250000w.sorted.n_aligned.bed


source /users/jsouto/.bashrc

conda activate programs

# compute the proportion of TR aligned to the target species in each window, and the total aligned proportion of the windows

python compute_aligned_proportions.py