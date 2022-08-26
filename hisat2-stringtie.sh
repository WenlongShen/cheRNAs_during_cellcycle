#!/bin/bash

# bash hisat2-stringtie.sh -s /data1/STAR_genomeDIR/hg38.p13 -x /data1/hisat2_index/hg38.p13 -g /data1/STAR_genomeDIR/gencode.v33.chr_patch_hapl_scaff.annotation.gtf -d ../test

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


mkdir ../STAR
mkdir ../hisat2
mkdir ../stringtie
mkdir ../Bam
STAR_dir="../STAR"
HISAT_dir="../hisat2"
Stringtie_dir="../stringtie"
Bam_dir="../Bam"


BC=$(ls ${Path_fq_dir}|grep ".fastq.gz$")

for FILE in ${BC[*]}; do
       SAMPLE=${FILE/\.fastq\.gz/}

	echo -e "processing ${FILE} \t ${SAMPLE} " 
	
	STAR   --runThreadN 12\
	--genomeDir ${Path_to_STAR_genomedir} \
	--genomeLoad LoadAndKeep \
	--readFilesIn ${Path_fq_dir}/$FILE\
	--readFilesCommand zcat\
       --outReadsUnmapped Fastx\
       --chimSegmentMin 18\
       --chimScoreMin 12\
       --chimOutType WithinBAM\
       --outSAMtype BAM SortedByCoordinate\
       --limitBAMsortRAM 40800000000 \
       --outFileNamePrefix ${STAR_dir}/${SAMPLE}.


	hisat2 -p 8 --dta -x ${Path_to_index} -U ${STAR_dir}/${SAMPLE}.Unmapped.out.mate1  -S ${HISAT_dir}/${SAMPLE}.sam
       samtools sort -@ 16 -o ${HISAT_dir}/${SAMPLE}.hisat2.bam ${HISAT_dir}/${SAMPLE}.sam
      
       samtools merge -h ${HISAT_dir}/${SAMPLE}.hisat2.bam ${Bam_dir}/${SAMPLE}.merged.bam ${STAR_dir}/${SAMPLE}.Aligned.sortedByCoord.out.bam ${HISAT_dir}/${SAMPLE}.hisat2.bam 
       
       rm ${HISAT_dir}/${SAMPLE}.sam
       stringtie -p 8 -G ${Ref_gtf} -o ${Stringtie_dir}/${SAMPLE}.gtf ${Bam_dir}/${SAMPLE}.merged.bam
       echo "${Stringtie_dir}/${SAMPLE}.gtf" >> ${Stringtie_dir}/mergelist.txt
     

done


	
stringtie --merge -p 8 -G ${Ref_gtf} -o ${Stringtie_dir}/stringtie_merged.gtf ${Stringtie_dir}/mergelist.txt
cd ${Stringtie_dir}
gffcompare –r ${Ref_gtf} –G –o merged stringtie_merged.gtf

