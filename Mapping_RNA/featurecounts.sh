#!/bin/bash
#$ -N featureCounts
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -l h_rt=48:00:00
#$ -m baes

date
featurecounts=./Install/subread-2.0.6-Linux-x86_64/bin/featureCounts

#In batch
GENOME=./WCS/RNA_job/NCBI/GCF_028769735.1_RI_Zleu_2.0_genomic.fna

anno=./WCS/RNA_job/NCBI/GCF_028769735.1_RI_Zleu_2.0_genomic.gtf

$featurecounts -O -p --countReadPairs -B -d 50 -D 5000 -t exon -g gene_id -T 8 -a ${anno} -G $GENOME \
 -o Output_WCS_RNA_ncbi.txt  \
mapping/*WCS*/*Aligned.sortedByCoord.out.bam  \
  >& Output_WCS_RNA_ncbi.log

#END#
