#!/bin/bash

# bash hisat2-stringtie.sh -s /data1/STAR_genomeDIR/hg38.p13 -x /data1/hisat2_index/hg38.p13 -g /data1/STAR_genomeDIR/gencode.v33.chr_patch_hapl_scaff.annotation.gtf -d ../test

# bash stringtie_e.sh -s /data1/STAR_genomeDIR/hg38.p13 -x /data1/hisat2_index/hg38.p13 -g ../stringtie/stringtie_merged.gtf -d ../fastq/


# parsing command line options
while getopts "s:x:g:d:" OPTION
do
     case $OPTION in
         s) Path_to_STAR_genomedir=$OPTARG ;;
         x) Path_to_index=$OPTARG ;;
         g) Ref_gtf=$OPTARG ;;
	  d) Path_fq_dir=$OPTARG ;;
         ?)
             echo "incorrect option"
             exit 0
             ;;
     esac
done


mkdir "../Ballgown_stringtie_merged"
STAR_dir="../STAR"
HISAT_dir="../hisat2"
Stringtie_dir="../stringtie"
Bam_dir="../Bam"
BG_dir="../Ballgown_stringtie_merged"

BC=$(ls ${Path_fq_dir}|grep ".fastq.gz$")

for FILE in ${BC[*]}; do
       SAMPLE=${FILE/\.fastq\.gz/}

	echo -e "processing ${FILE} \t ${SAMPLE} " 
	

stringtie -e -B -p 8 -G ${Ref_gtf} \
 -A ${BG_dir}/${SAMPLE}/${SAMPLE}.abund.txt\
 -o ${BG_dir}/${SAMPLE}/${SAMPLE}.gtf \
 ${Bam_dir}/${SAMPLE}.merged.bam


done


