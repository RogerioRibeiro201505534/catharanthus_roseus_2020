#!/bin/bash 



#Run from /mnt/Disk1/rogerio

#File tree 
#--00_raw_data
#	|
#	---idioblastoma
#	|	|
#	|	---idio|meso|leaves|totalptt
#	|		|
#	|		---sample_1|sample_2|sample3
#	|				|	
#	|				R1_cut_paired
#	|				R2_cut_paired
#	|
#	---InOut
#		|
#		---Ext|Int
#			|
#			---F1|F4
#					|
#					---P1|P2|P3
#						|
#						R1_cut_paired
#						R2_cut_paired
#
#
#--catharanthus_genome_v2
#	|
#	cro_v2.gene_models.gff3 
#	cro_v2_asm.fasta
#


##Note: Idioblast, mesophyll, totalptt and leaf samples were trimmed with the following comand (prior to this dissertation)
#java -jar $trimomaticPath/trimmomatic-0.39.jar PE -threads 10 -phred33 140817Can73-16_S1_L001_R1_001.fastq.gz 140817Can73-16_S1_L001_R2_001.fastq.gz R1_cut_paired.fastq.gz  R1_cut_unpaired.fastq.gz  R2_cut_paired.fastq.gz  R2_cut_unpaired.fastq.gz ILLUMINACLP:TruSeq3-PE.fa:2:30:10 ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50

##Note2: Interior and Exterior leaves were trimmed with the following comand (prior to this dissertation)
#java -jar $trimomaticPath/trimmomatic-0.39.jar PE -threads 10 -phred33 P1F1Int_R1_001.fastq.gz P1F1Int_R2_001.fastq.gz R1_cut_paired.fastq.gz  R1_cut_unpaired.fastq.gz  R2_cut_paired.fastq.gz  R2_cut_unpaired.fastq.gz ILLUMINACLP:TruSeq3-PE.fa:2:30:10 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50


###Step 0) Build an index out of catharanthus genome

gffread -E catharanthus_genome_v2/cro_v2.gene_models.gff3 -T > catharanthus_genome_v2/cro_v2.gene_models.gtf
hisat2_extract_exons.py catharanthus_genome_v2/cro_v2.gene_models.gtf > catharanthus_genome_v2/catharanthus_exons.txt
hisat2_extract_splice_sites.py catharanthus_genome_v2/cro_v2.gene_models.gtf > catharanthus_genome_v2/catharanthus_ss.txt

mkdir -p mapping_to_genome
mkdir -p mapping_to_genome/index 

echo 'generating index for hisat2'

hisat2-build --seed 1 -p 8 -f --exon catharanthus_genome_v2/catharanthus_exons.txt --ss catharanthus_genome_v2/catharanthus_ss.txt catharanthus_genome_v2/cro_v2_asm.fasta mapping_to_genome/index/cro_index


###Step 1) Map read with hisat2 

mkdir -p mapping_to_genome/01_hisat
mkdir -p mapping_to_genome/01_hisat/summary_files

echo 'mapping with hisat2' 

##Map idioblastos##

mkdir -p mapping_to_genome/01_hisat/idio

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/idio/sample_1/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/idio/sample_1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/idio_1.txt -S idio_1.sam
samtools view -Su idio_1.sam | samtools sort -o mapping_to_genome/01_hisat/idio/idio_1_sorted.bam
rm idio_1.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/idio/sample_2/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/idio/sample_2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/idio_2.txt -S idio_2.sam
samtools view -Su idio_2.sam | samtools sort -o mapping_to_genome/01_hisat/idio/idio_2_sorted.bam
rm idio_2.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/idio/sample_3/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/idio/sample_3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/idio_3.txt -S idio_3.sam
samtools view -Su idio_3.sam | samtools sort -o mapping_to_genome/01_hisat/idio/idio_3_sorted.bam
rm idio_3.sam

##Map leaves## 

mkdir -p mapping_to_genome/01_hisat/leaves

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/leaves/sample_1/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/leaves/sample_1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/leaves_1.txt -S leaves_1.sam
samtools view -Su leaves_1.sam | samtools sort -o mapping_to_genome/01_hisat/leaves/leaves_1_sorted.bam
rm leaves_1.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/leaves/sample_2/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/leaves/sample_2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/leaves_2.txt -S leaves_2.sam
samtools view -Su leaves_2.sam | samtools sort -o mapping_to_genome/01_hisat/leaves/leaves_2_sorted.bam
rm leaves_2.sam

##Map meso##

mkdir -p mapping_to_genome/01_hisat/meso

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/meso/sample_1/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/meso/sample_1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/meso_1.txt -S meso_1.sam
samtools view -Su meso_1.sam | samtools sort -o mapping_to_genome/01_hisat/meso/meso_1_sorted.bam
rm meso_1.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/meso/sample_2/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/meso/sample_2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/meso_2.txt -S meso_2.sam
samtools view -Su meso_2.sam | samtools sort -o mapping_to_genome/01_hisat/meso/meso_2_sorted.bam
rm meso_2.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/meso/sample_3/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/meso/sample_3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/meso_3.txt -S meso_3.sam
samtools view -Su meso_3.sam | samtools sort -o mapping_to_genome/01_hisat/meso/meso_3_sorted.bam
rm meso_3.sam

##Map totalptt##
mkdir -p mapping_to_genome/01_hisat/totalptt

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/totalptt/sample_1/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/totalptt/sample_1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/totalptt_1.txt -S totalptt_1.sam
samtools view -Su totalptt_1.sam | samtools sort -o mapping_to_genome/01_hisat/totalptt/totalptt_1_sorted.bam
rm totalptt_1.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/totalptt/sample_2/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/totalptt/sample_2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/totalptt_2.txt -S totalptt_2.sam
samtools view -Su totalptt_2.sam | samtools sort -o mapping_to_genome/01_hisat/totalptt/totalptt_2_sorted.bam
rm totalptt_2.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/idioblastome/totalptt/sample_3/R1_cut_paired.fastq.gz -2 00_raw_data/idioblastome/totalptt/sample_3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/totalptt_3.txt -S totalptt_3.sam
samtools view -Su totalptt_3.sam | samtools sort -o mapping_to_genome/01_hisat/totalptt/totalptt_3_sorted.bam
rm totalptt_3.sam

##Map F1_ext## 

mkdir -p mapping_to_genome/01_hisat/F1_EXT

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/EXT/F1/P1/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/EXT/F1/P1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F1P1_EXT.txt -S F1P1_EXT.sam
samtools view -Su F1P1_EXT.sam | samtools sort -o mapping_to_genome/01_hisat/F1_EXT/F1P1_EXT_sorted.bam
rm F1P1_EXT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/EXT/F1/P2/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/EXT/F1/P2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F1P2_EXT.txt -S F1P2_EXT.sam
samtools view -Su F1P2_EXT.sam | samtools sort -o mapping_to_genome/01_hisat/F1_EXT/F1P2_EXT_sorted.bam
rm F1P2_EXT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/EXT/F1/P3/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/EXT/F1/P3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F1P3_EXT.txt -S F1P3_EXT.sam
samtools view -Su F1P3_EXT.sam | samtools sort -o mapping_to_genome/01_hisat/F1_EXT/F1P3_EXT_sorted.bam
rm F1P3_EXT.sam


##Map F4_ext##

mkdir -p mapping_to_genome/01_hisat/F4_EXT

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/EXT/F4/P1/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/EXT/F4/P1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F4P1_EXT.txt -S F4P1_EXT.sam
samtools view -Su F4P1_EXT.sam | samtools sort -o mapping_to_genome/01_hisat/F4_EXT/F4P1_EXT_sorted.bam
rm F4P1_EXT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/EXT/F4/P2/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/EXT/F4/P2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F4P2_EXT.txt -S F4P2_EXT.sam
samtools view -Su F4P2_EXT.sam | samtools sort -o mapping_to_genome/01_hisat/F4_EXT/F4P2_EXT_sorted.bam
rm F4P2_EXT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/EXT/F4/P3/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/EXT/F4/P3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F4P3_EXT.txt -S F4P3_EXT.sam
samtools view -Su F4P3_EXT.sam | samtools sort -o mapping_to_genome/01_hisat/F4_EXT/F4P3_EXT_sorted.bam
rm F4P3_EXT.sam


##Map F1_INT## 

mkdir -p mapping_to_genome/01_hisat/F1_INT

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/INT/F1/P1/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/INT/F1/P1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F1P1_INT.txt -S F1P1_INT.sam
samtools view -Su F1P1_INT.sam | samtools sort -o mapping_to_genome/01_hisat/F1_INT/F1P1_INT_sorted.bam
rm F1P1_INT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/INT/F1/P2/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/INT/F1/P2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F1P2_INT.txt -S F1P2_INT.sam
samtools view -Su F1P2_INT.sam | samtools sort -o mapping_to_genome/01_hisat/F1_INT/F1P2_INT_sorted.bam
rm F1P2_INT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/INT/F1/P3/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/INT/F1/P3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F1P3_INT.txt -S F1P3_INT.sam
samtools view -Su F1P3_INT.sam | samtools sort -o mapping_to_genome/01_hisat/F1_INT/F1P3_INT_sorted.bam
rm F1P3_INT.sam


##Map F4_INT##

mkdir -p mapping_to_genome/01_hisat/F4_INT

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/INT/F4/P1/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/INT/F4/P1/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F4P1_INT.txt -S F4P1_INT.sam
samtools view -Su F4P1_INT.sam | samtools sort -o mapping_to_genome/01_hisat/F4_INT/F4P1_INT_sorted.bam
rm F4P1_INT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/INT/F4/P2/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/INT/F4/P2/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F4P2_INT.txt -S F4P2_INT.sam
samtools view -Su F4P2_INT.sam | samtools sort -o mapping_to_genome/01_hisat/F4_INT/F4P2_INT_sorted.bam
rm F4P2_INT.sam

hisat2 -q --dta -p 4 -x mapping_to_genome/index/cro_index -1 00_raw_data/InOut/INT/F4/P3/R1_cut_paired.fastq.gz -2 00_raw_data/InOut/INT/F4/P3/R2_cut_paired.fastq.gz --summary-file mapping_to_genome/01_hisat/summary_files/F4P3_INT.txt -S F4P3_INT.sam
samtools view -Su F4P3_INT.sam | samtools sort -o mapping_to_genome/01_hisat/F4_INT/F4P3_INT_sorted.bam
rm F4P3_INT.sam


###Step 2)use countFeatures before stringTie. This will give us an idea of the difference of reads used before and after assembly

mkdir -p mapping_to_genome/02_FeatureCounts_before_assembly

echo 'running feature counts with catharanthus genome gene models as reference'

##Idioblastome experience 
featureCounts mapping_to_genome/01_hisat/idio/idio_1_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/idio_1_counts.txt
featureCounts mapping_to_genome/01_hisat/idio/idio_2_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/idio_2_counts.txt
featureCounts mapping_to_genome/01_hisat/idio/idio_3_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/idio_3_counts.txt

featureCounts mapping_to_genome/01_hisat/meso/meso_1_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/meso_1_counts.txt
featureCounts mapping_to_genome/01_hisat/meso/meso_2_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/meso_2_counts.txt
featureCounts mapping_to_genome/01_hisat/meso/meso_3_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/meso_3_counts.txt

featureCounts mapping_to_genome/01_hisat/leaves/leaves_1_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/leaves_1_counts.txt
featureCounts mapping_to_genome/01_hisat/leaves/leaves_2_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/leaves_2_counts.txt

featureCounts mapping_to_genome/01_hisat/totalptt/totalptt_1_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/totalptt_1_counts.txt
featureCounts mapping_to_genome/01_hisat/totalptt/totalptt_2_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/totalptt_2_counts.txt
featureCounts mapping_to_genome/01_hisat/totalptt/totalptt_3_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/totalptt_3_counts.txt

#InOut experience
featureCounts mapping_to_genome/01_hisat/F1_EXT/F1P1_EXT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F1P1_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_EXT/F1P2_EXT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F1P2_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_EXT/F1P3_EXT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F1P3_EXT_counts.txt

featureCounts mapping_to_genome/01_hisat/F4_EXT/F4P1_EXT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F4P1_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_EXT/F4P2_EXT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F4P2_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_EXT/F4P3_EXT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F4P3_EXT_counts.txt

featureCounts mapping_to_genome/01_hisat/F1_INT/F1P1_INT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F1P1_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_INT/F1P2_INT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F1P2_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_INT/F1P3_INT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F1P3_INT_counts.txt

featureCounts mapping_to_genome/01_hisat/F4_INT/F4P1_INT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F4P1_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_INT/F4P2_INT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F4P2_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_INT/F4P3_INT_sorted.bam -a catharanthus_genome_v2/cro_v2.gene_models.gtf  -g transcript_id -p -o mapping_to_genome/02_FeatureCounts_before_assembly/F4P3_INT_counts.txt


###Step 3)
#3.1) Use stringtie for transcriptome assembly

mkdir -p mapping_to_genome/03_stringtie
mkdir -p mapping_to_genome/03_stringtie/single_assembly

echo 'Building individual assemblies with stringtie'

#Idioblastoma experience
stringtie mapping_to_genome/01_hisat/idio/idio_1_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/idio_1.gtf -p 4 -l idio_1
stringtie mapping_to_genome/01_hisat/idio/idio_2_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/idio_2.gtf -p 4 -l idio_2
stringtie mapping_to_genome/01_hisat/idio/idio_3_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/idio_3.gtf -p 4 -l idio_3

stringtie mapping_to_genome/01_hisat/leaves/leaves_1_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/leaves_1.gtf -p 4 -l leaves_1
stringtie mapping_to_genome/01_hisat/leaves/leaves_2_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/leaves_2.gtf -p 4 -l leaves_1

stringtie mapping_to_genome/01_hisat/meso/meso_1_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/meso_1.gtf -p 4 -l meso_1
stringtie mapping_to_genome/01_hisat/meso/meso_2_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/meso_2.gtf -p 4 -l meso_2
stringtie mapping_to_genome/01_hisat/meso/meso_3_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/meso_3.gtf -p 4 -l meso_3

stringtie mapping_to_genome/01_hisat/totalptt/totalptt_1_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/totalptt_1.gtf -p 4 -l totalptt_1
stringtie mapping_to_genome/01_hisat/totalptt/totalptt_2_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/totalptt_2.gtf -p 4 -l totalptt_2
stringtie mapping_to_genome/01_hisat/totalptt/totalptt_3_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/totalptt_3.gtf -p 4 -l totalptt_3


#Inout experience
stringtie mapping_to_genome/01_hisat/F1_EXT/F1P1_EXT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F1P1E.gtf -p 4 -l F1P1E
stringtie mapping_to_genome/01_hisat/F1_EXT/F1P2_EXT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F1P2E.gtf -p 4 -l F1P2E
stringtie mapping_to_genome/01_hisat/F1_EXT/F1P3_EXT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F1P3E.gtf -p 4 -l F1P3E

stringtie mapping_to_genome/01_hisat/F4_EXT/F4P1_EXT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F4P1E.gtf -p 4 -l F4P1E
stringtie mapping_to_genome/01_hisat/F4_EXT/F4P2_EXT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F4P2E.gtf -p 4 -l F4P2E
stringtie mapping_to_genome/01_hisat/F4_EXT/F4P3_EXT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F4P3E.gtf -p 4 -l F4P3E

stringtie mapping_to_genome/01_hisat/F1_INT/F1P1_INT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F1P1I.gtf -p 4 -l F1P1I
stringtie mapping_to_genome/01_hisat/F1_INT/F1P2_INT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F1P2I.gtf -p 4 -l F1P2I
stringtie mapping_to_genome/01_hisat/F1_INT/F1P3_INT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F1P3I.gtf -p 4 -l F1P3I

stringtie mapping_to_genome/01_hisat/F4_INT/F4P1_INT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F4P1I.gtf -p 4 -l F4P1I
stringtie mapping_to_genome/01_hisat/F4_INT/F4P2_INT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F4P2I.gtf -p 4 -l F4P2I
stringtie mapping_to_genome/01_hisat/F4_INT/F4P3_INT_sorted.bam -o mapping_to_genome/03_stringtie/single_assembly/F4P3I.gtf -p 4 -l F4P3I


#3.2) Use stringtie merge to merge assemblies
echo 'merging stringTie assembly'

#Assembly with all stringtie files
mkdir -p mapping_to_genome/03_stringtie/merged_transcriptome/
stringtie --merge mapping_to_genome/03_stringtie/single_assembly/* -G catharanthus_genome_v2/cro_v2.gene_models.gff3 -o mapping_to_genome/03_stringtie/merged_transcriptome/merged_transcriptome.gtf -l CATHA

#3.3) Use gffread to retrieve transcripts before filtering
gffread -F -w mapping_to_genome/03_stringtie/merged_transcriptome/merged_transcriptome.fasta -g catharanthus_genome_v2/cro_v2_asm.fasta mapping_to_genome/03_stringtie/merged_transcriptome/merged_transcriptome.gtf


#Step 4) obtain statistics about the merged file and an annotation file which will be used downstream 
echo 'obtaining statistics for Stringtie assembly'

mkdir -p mapping_to_genome/04_gffcompare

gffcompare mapping_to_genome/03_stringtie/merged_transcriptome/merged_transcriptome.gtf -o mapping_to_genome/04_gffcompare/merged_transcriptome -r catharanthus_genome_v2/cro_v2.gene_models.gff3 -s catharanthus_genome_v2/cro_v2_asm.fasta -V

mv mapping_to_genome/03_stringtie/merged_transcriptome/*map mapping_to_genome/04_gffcompare

#4.2) #run the filter python script in the gff compare 
py filter_gff.py

#4.3) run gffread to retrieve the transcripts after filtering 
gffread -F -w mapping_to_genome/04_gffcompare/merged_transcriptome_filtered.fasta -g catharanthus_genome_v2/cro_v2_asm.fasta mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf


#Step 5)Count features after 

echo 'running feature counts with stringTie merged assembly gene models as reference'

mkdir -p mapping_to_genome/05_featureCounts

featureCounts mapping_to_genome/01_hisat/idio/idio_1_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/idio_1_counts.txt
featureCounts mapping_to_genome/01_hisat/idio/idio_2_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/idio_2_counts.txt
featureCounts mapping_to_genome/01_hisat/idio/idio_3_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/idio_3_counts.txt

featureCounts mapping_to_genome/01_hisat/meso/meso_1_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/meso_1_counts.txt
featureCounts mapping_to_genome/01_hisat/meso/meso_2_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/meso_2_counts.txt
featureCounts mapping_to_genome/01_hisat/meso/meso_3_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/meso_3_counts.txt

featureCounts mapping_to_genome/01_hisat/leaves/leaves_1_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/leaves_1_counts.txt
featureCounts mapping_to_genome/01_hisat/leaves/leaves_2_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/leaves_2_counts.txt

featureCounts mapping_to_genome/01_hisat/totalptt/totalptt_1_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/totalptt_1_counts.txt
featureCounts mapping_to_genome/01_hisat/totalptt/totalptt_2_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/totalptt_2_counts.txt
featureCounts mapping_to_genome/01_hisat/totalptt/totalptt_3_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/totalptt_3_counts.txt

#InOut experience
featureCounts mapping_to_genome/01_hisat/F1_EXT/F1P1_EXT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F1P1_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_EXT/F1P2_EXT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F1P2_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_EXT/F1P3_EXT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F1P3_EXT_counts.txt

featureCounts mapping_to_genome/01_hisat/F4_EXT/F4P1_EXT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F4P1_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_EXT/F4P2_EXT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F4P2_EXT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_EXT/F4P3_EXT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F4P3_EXT_counts.txt

featureCounts mapping_to_genome/01_hisat/F1_INT/F1P1_INT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F1P1_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_INT/F1P2_INT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F1P2_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F1_INT/F1P3_INT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F1P3_INT_counts.txt

featureCounts mapping_to_genome/01_hisat/F4_INT/F4P1_INT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F4P1_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_INT/F4P2_INT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F4P2_INT_counts.txt
featureCounts mapping_to_genome/01_hisat/F4_INT/F4P3_INT_sorted.bam -a mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf  -g gene_name -p -o mapping_to_genome/05_featureCounts/F4P3_INT_counts.txt

echo 'input these files into the deSeq2 script after cleaning up the outputs (counts)' 


#mapping with kallisto 

#6.1) build index
kallisto index -i mapping_to_genome/index/cro_kal_index mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf

#target de Bruijn graph has 525804 contigs and contains 55484603 k-mers

#6.2) Map with kallisto quant

mkdir -p mapping_to_genome/07_kallisto

##idioblastoma experience

#Map idioblastos
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/idio_1 -t 8 00_raw_data/idioblastome/idio/sample_1/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/idio/sample_1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/idio_2 -t 8 00_raw_data/idioblastome/idio/sample_2/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/idio/sample_2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/idio_3 -t 8 00_raw_data/idioblastome/idio/sample_3/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/idio/sample_3/R2_cut_paired.fastq.gz


#map leaves 
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/leaves_1 -t 8 00_raw_data/idioblastome/leaves/sample_1/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/leaves/sample_1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/leaves_2 -t 8 00_raw_data/idioblastome/leaves/sample_2/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/leaves/sample_2/R2_cut_paired.fastq.gz


#map mesophyll 

kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/meso_1 -t 8 00_raw_data/idioblastome/meso/sample_1/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/meso/sample_1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/meso_2 -t 8 00_raw_data/idioblastome/meso/sample_2/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/meso/sample_2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/meso_3 -t 8 00_raw_data/idioblastome/meso/sample_3/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/meso/sample_3/R2_cut_paired.fastq.gz


#Map protoplast 

kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/totalptt_1 -t 8 00_raw_data/idioblastome/totalptt/sample_1/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/totalptt/sample_1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/totalptt_2 -t 8 00_raw_data/idioblastome/totalptt/sample_2/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/totalptt/sample_2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/totalptt_3 -t 8 00_raw_data/idioblastome/totalptt/sample_3/R1_cut_paired.fastq.gz 00_raw_data/idioblastome/totalptt/sample_3/R2_cut_paired.fastq.gz


##InOut

#Map F1_Ext
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F1P1_Ext -t 8 00_raw_data/InOut/EXT/F1/P1/R1_cut_paired.fastq.gz 00_raw_data/InOut/EXT/F1/P1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F1P2_Ext -t 8 00_raw_data/InOut/EXT/F1/P2/R1_cut_paired.fastq.gz 00_raw_data/InOut/EXT/F1/P2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F1P3_Ext -t 8 00_raw_data/InOut/EXT/F1/P3/R1_cut_paired.fastq.gz 00_raw_data/InOut/EXT/F1/P3/R2_cut_paired.fastq.gz


#Map F4_Ext
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F4P1_Ext -t 8 00_raw_data/InOut/EXT/F4/P1/R1_cut_paired.fastq.gz 00_raw_data/InOut/EXT/F4/P1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F4P2_Ext -t 8 00_raw_data/InOut/EXT/F4/P2/R1_cut_paired.fastq.gz 00_raw_data/InOut/EXT/F4/P2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F4P3_Ext -t 8 00_raw_data/InOut/EXT/F4/P3/R1_cut_paired.fastq.gz 00_raw_data/InOut/EXT/F4/P3/R2_cut_paired.fastq.gz


#Map F1_Int
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F1P1_Int -t 8 00_raw_data/InOut/INT/F1/P1/R1_cut_paired.fastq.gz 00_raw_data/InOut/INT/F1/P1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F1P2_Int -t 8 00_raw_data/InOut/INT/F1/P2/R1_cut_paired.fastq.gz 00_raw_data/InOut/INT/F1/P2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F1P3_Int -t 8 00_raw_data/InOut/INT/F1/P3/R1_cut_paired.fastq.gz 00_raw_data/InOut/INT/F1/P3/R2_cut_paired.fastq.gz


#Map F4_Int
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F4P1_Int -t 8 00_raw_data/InOut/INT/F4/P1/R1_cut_paired.fastq.gz 00_raw_data/InOut/INT/F4/P1/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F4P2_Int -t 8 00_raw_data/InOut/INT/F4/P2/R1_cut_paired.fastq.gz 00_raw_data/InOut/INT/F4/P2/R2_cut_paired.fastq.gz
kallisto quant -i mapping_to_genome/index/cro_kal_index -o mapping_to_genome/07_kallisto/F4P3_Int -t 8 00_raw_data/InOut/INT/F4/P3/R1_cut_paired.fastq.gz 00_raw_data/InOut/INT/F4/P3/R2_cut_paired.fastq.gz

#step 7) run DESeq2 for both counts

