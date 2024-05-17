#!/bin/sh
# Grid Engine options
#$ -N map
#$ -cwd
#$ -l h_rt=200:00:00
#$ -l h_vmem=8G
#$ -m baes
#$ -pe sharedmem 6

# If you plan to load any software modules, then you must first initialise the modules framework.
. /etc/profile.d/modules.sh
# Choose the staging environment
export OMP_NUM_THREADS=$NSLOTS

# Then, you must load the modules themselves
module load igmm/apps/STAR/2.7.8a
#Usage: qsub -N mapADR STAR_map_multiFQ.sh List_data.GWCS.ADREN

#STAR Index - done
#STAR mapping
mkdir -p mapping

for prefix in `less $1 `;do \
#for prefix in "FAT-47GWCS-70" ;do \  #
#prefix=$1 #Sample_list.txt e.g. K1_1 $i ##  for i in `less ../GWCS_RNA_list.txt`;do sh STAR_map_multiFQ.sh $i ;done

output=./WCS/RNA_job/mapping/${prefix}
mkdir -p $output
ref=./WCS/RNA_job/Star_index_G_NCBI
indir=./WCS/RNA_job/Raw_reads/
cd ${indir}
#ls  ${listR1} ${listR2}
 date
echo "Output: ${output}/${prefix}"

# ${prefix}_mulFQ` should be created before this script
mytext=`cat ${prefix}_mulFQ`
#echo"
STAR --genomeDir ${ref} --runThreadN 6 --readFilesCommand zcat --readFilesIn \
$mytext \
--outFilterType BySJout --outSAMunmapped None --outReadsUnmapped Fastx --outFileNamePrefix ${output}/${prefix}_ncbi. --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 7900000000
#"
done
