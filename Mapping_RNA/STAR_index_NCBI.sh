#!/bin/sh

# Grid Engine options
#$ -N S_index
#$ -cwd
#$ -l h_rt=120:30:00
#$ -l h_vmem=4G
#$ -m baes
#$ -pe sharedmem 12

# If you plan to load any software modules, then you must first initialise the modules framework.
. /etc/profile.d/modules.sh
# Choose the staging environment
export OMP_NUM_THREADS=$NSLOTS

# Then, you must load the modules themselves
module load igmm/apps/STAR/2.7.8a

Species=G_NCBI
ref='/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/WCS/RNA_job/NCBI/GCF_028769735.1_RI_Zleu_2.0_genomic.fna'
anno='/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/WCS/RNA_job/NCBI/GCF_028769735.1_RI_Zleu_2.0_genomic.gtf'

mkdir Star_index_${Species}
STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ./Star_index_${Species} \
--sjdbGTFfile ${anno} \
--genomeFastaFiles ${ref}

###END###
