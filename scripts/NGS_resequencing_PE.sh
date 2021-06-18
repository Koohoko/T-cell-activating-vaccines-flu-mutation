#!/bin/bash

# Maireid 2020-12-06
mkdir ../results
mkdir ../FastQC
mkdir ../results/MEGAHIT

cd ../data

~/softwares/samtools-1.11/bin/samtools faidx reference.fasta
~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index reference.fasta
~/softwares/bwa/bwa index reference.fasta
java -jar ~/softwares/picard.jar CreateSequenceDictionary -R reference.fasta -O reference.dict
~/softwares/bioawk/bioawk -c fastx '{print $name"\t1\t"length($seq)"\t"$name}' reference.fasta > reference.bed


for fwdread in `ls | grep '_1.fastq.gz'`
do
	sample=$(echo $fwdread | cut -d"_" -f 1)
	rwsread=$sample"_2.fastq.gz"
	bamfile=$sample"_sorted.bam"
	# mapping
	~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 46 reference.fasta $fwdread $rwsread | ~/softwares/samtools-1.11/bin/samtools view -b --threads 46 - | ~/softwares/samtools-1.11/bin/samtools sort -o aln.sorted.bam --threads 46 -m 4G - 
	mv aln.sorted.bam $sample"_sorted.bam"
	~/softwares/samtools-1.11/bin/samtools index $sample"_sorted.bam"

	# de novo assembly	
	# ~/softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -t 46 -1 $fwdread -2 $rwsread -o ../results/MEGAHIT/$sample
	# ## readcount
	bedtools makewindows -w 650 -b reference.bed > ref_slide_window.bed
	cat ref_slide_window.bed | awk '{print $1":"$2"-"$3}' | parallel -j 47 -k "~/softwares/samtools-1.11/bin/samtools mpileup -aa -d 1000000000 -f reference.fasta --region {} $sample"_sorted.bam" " | ~/softwares/mpileup2readcounts/build/mpileup2readcounts > ../FastQC/$sample"_bam.readcount.mpileup.txt"	

	#FAST-QC
	# ~/FastQC/fastqc $sample"_marked_duplicates.bam" -t 48 -o ../FastQC/

	## freebayes
	~/softwares/freebayes/scripts/freebayes-parallel <(~/softwares/freebayes/scripts/fasta_generate_regions.py reference.fasta.fai 650) 47 -f reference.fasta -F 0.01 -C 1 -p 1 -q 30 -K --min-coverage 5 $bamfile > ../results/"freebayes_"$sample".vcf"
			
	## VarDict
	bedtools makewindows -w 650 -b reference.bed > ref_slide_window.bed
	cat ref_slide_window.bed | awk '{print $1":"$2"-"$3}' | parallel -j 47 -k "~/softwares/VarDictJava/build/install/VarDict/bin/VarDict -G reference.fasta -f 0.01 -N samplename -th 48 -b $bamfile -R {} " | ~/softwares/VarDictJava/build/install/VarDict/bin/teststrandbias.R | ~/softwares/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl -N samplename -E -f 0.01 > ../results/"vardict_"$sample".vcf"

	## lofreq
	/home/hggu/softwares/lofreq/src/lofreq/lofreq call-parallel --pp-threads 47 -f reference.fasta -o ../results/"lofreq_"$sample".vcf" $bamfile
	
done

