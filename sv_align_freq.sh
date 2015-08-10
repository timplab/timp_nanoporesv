#!/bin/bash
##Ok - take in the fastq for alignment from the $1

if [ -z "$1" ]; then
    echo "No arg or file doesn't exist"
    exit 1
else 
    fastq=$1
fi

##And the outdir from $2

##If command line arg empty - or not a dir - die
if [ ! -d "$2" ]; then
    echo "No arg or dir doesn't exist"
    exit 1
else 
    outdir=$2
fi

echo $fastq

echo "$outdir $prefix "

vcfprefix=`basename $fastq`

for bc in {01..12}
#for bc in {01..01}
do 

    prefix=`basename $fastq`_BC${bc}

    if [ -e "${prefix}.fastq.gz" ]
    then

	##Number of reads
	zcat ${fastq}_BC${bc}.fastq.gz | echo $((`wc -l`/4)) >${outdir}/${prefix} >${outdir}/${prefix}.readnum
	
	##BWAMEM
	~/Code/bwa/bwa mem -R "@RG\tID:${fastq}_BC${bc}\tSM:BC${bc}\tPL:ONT" -x ont2d /mithril/Data/NGS/Reference/human/hg19.fa \
	    ${fastq}_BC${bc}.fastq.gz | samtools view -S -b - >${outdir}/${prefix}.bwa.bam
	
	samtools index ${outdir}/${prefix}.bwa.bam

	
	##Error rate
	samtools view -b -F 4 ${outdir}/${prefix}.bwa.bam | samtools calmd -b - /mithril/Data/NGS/Reference/human/hg19.fa >${outdir}/${prefix}.md.bam
	samtools index ${outdir}/${prefix}.md.bam 		
	python ~/Code/timp_nanopore/oxford/bam_extract.py ${outdir}/${prefix}.md.bam


	#samtools view -b -S ${outdir}/${prefix}.bwa.sam | samtools sort - ${outdir}/${prefix}.bwa

	#samtools index ${outdir}/${prefix}.bwa.bam
	
	##From LUMPY readme
	
	# Extract the split-read alignments
	samtools view -h ${outdir}/${prefix}.bwa.bam \
	    | ~/Code/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
	    | samtools view -Sb - \
	    > ${outdir}/${prefix}.bwa.splitters.unsorted.bam
	
	# Sort split alignments
	samtools sort ${outdir}/${prefix}.bwa.splitters.unsorted.bam ${outdir}/${prefix}.bwa.splitters	
	samtools index ${outdir}/${prefix}.bwa.splitters.bam
	
	bamToBed -bed12 -i ${outdir}/${prefix}.bwa.bam >${outdir}/${prefix}.bwa.bed

 	##Calculate Reads in Designed Region
	#samtools flagstat ${outdir}/${prefix}.bwa.bam >${outdir}/${prefix}.filtstats
	coverageBed -a p16_amp.bed -b ${outdir}/${prefix}.bwa.bam >${outdir}/${prefix}.bamreport
	coverageBed -d -a p16_amp.bed -b ${outdir}/${prefix}.bwa.bam >${outdir}/${prefix}.splitbamreport

	##LUMPY
	
	~/Code/lumpy-sv/bin/lumpy \
	    -mw 4 \
	    -tt 0 \
	    -sr id:${prefix}.bwa,bam_file:${outdir}/${prefix}.bwa.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
	    > ${outdir}/${prefix}.bwa.vcf

	~/Code/lumpy-sv/bin/lumpy \
	    -mw 4 \
	    -tt 0 \
	    -sr id:${prefix}.bwa,bam_file:${outdir}/${prefix}.sorted.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
	    > ${outdir}/${prefix}.nosplit.vcf



    fi
done
    
    




