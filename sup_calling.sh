#!/bin/sh

#This  bascript pipeline for snp calling and  msmsc2 estimate of effective population size

#change directory
cd /gondor/dieunel/Austreococcus_tauri/snp_calling


#Download
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR072/SRR4026808/SRR4026808_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR072/SRR4026808/SRR4026808_2.fastq.gz

#unzip
gunzip SRR4026808_2.fastq.gz
gunzip SRR4026808_1.fastq.gz

#filter
fastp -i SRR4026808_1.fastq -I SRR4026808_2.fastq -o SRR4026808_1.fq -O SRR4026808_2.fq

#repair
repair.sh in=SRR4026808_1.fq in2=SRR4026808_2.fq out=SRR4026808_f1.fq out2=SRR4026808_f2.fq

#mapping with bbmap
bbmap.sh -Xmx121g t=16 in1=SRR4026808_f1.fq in2=SRR4026808_f2.fq fast=f ref=GCF_000214015.3_version_140606_genomic.fna out=SRR4026808.sam

## To get the summary statistics about the alignment
samtools flagstat SRR4026808.sam > SRR4026808_sam_summary.txt

#convert sam to bam
samtools view -bS SRR4026808.sam -@ 16 > SRR4026808.bam

# sort the bam file
samtools sort SRR4026808.bam -o SRR4026808.sorted.bam -@ 16

#Creating the fasta sequence dictionary file ( this scriot should generate a "*.dict" file
#java -jar /gondor/dieunel/picard.jar CreateSequenceDictionary R= GCF_000214015.3_version_140606_genomic.fna O= GCF_000214015.3_version_140606_genomic.dict


#Creating the fasta index file (this script will generate a  ("*.fai" file)

#Mark and remove duplicates
java -XX:ParallelGCThreads=16 -jar /gondor/dieunel/picard.jar MarkDuplicates INPUT=SRR4026808.sorted.bam OUTPUT=SRR4026808_Nodup.sorted.bam METRICS_FILE=metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

# Replace all read groups in the INPUT (sorted.bam)file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
java -jar  /gondor/dieunel/picard.jar AddOrReplaceReadGroups I=SRR4026808_Nodup.sorted.bam O=SRR4026808_NoDup_addrplced.bam RGID=4 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=20

# inedxing the output bam file
samtools index SRR4026808_NoDup_addrplced.bam

#Define intervals
java -XX:ParallelGCThreads=16 -jar /gondor/dieunel/GenomeAnalysisTK.jar -T RealignerTargetCreator -R GCF_000214015.3_version_140606_genomic.fna -I SRR4026808_NoDup_addrplced.bam -o Intervals.list

#Realign against indels
java -XX:ParallelGCThreads=16 -jar /gondor/dieunel/GenomeAnalysisTK.jar -T IndelRealigner -R GCF_000214015.3_version_140606_genomic.fna -I SRR4026808_NoDup_addrplced.bam -targetIntervals Intervals.list -o SRR4026808_NoDupRealigned.bam


samtools mpileup -C50 -Q 20  -uf GCF_000214015.3_version_140606_genomic.fna SRR4026808_NoDupRealigned.bam | bcftools call -c > SRR4026808.vcf
vcftools --vcf SRR4026808.vcf --recode --recode-INFO-all --out SRR4026808 --minQ 20.00 --minDP 5
bcftools stats  SRR4026808.recode.vcf >  SRR4026808_FilteredVCFStats.txt





#For effective population size inference, we need snps in each chromosome
grep '>' GCF_000214015.3_version_140606_genomic.fna  #printing the chromosome name (header)

#Create input snps for msmc2
./generate_multihetsep.py chr1.vcf.gz > chr1.msmc.input.txt #Repeated for all the chromosomes
#Effective population size history
./msmc2 -p 1*2+15*1+1*2 -o final_11 *input.txt #All chromosomes that were generated in previous step

#Plotting the msmc2 output in R
f10=read.table('/Users/mzillur/thesis/msmc/final_11.final.txt',header = T)
plot(f10$left_time_boundary/mu*g, (1/f10$lambda)/(2*mu),log='x', type = 'n', xlab='Years ago',ylab = 'Effective population size')
lines(f10$left_time_boundary/mu*g,(1/f10$lambda)/(2*mu),type='s',col='red')
