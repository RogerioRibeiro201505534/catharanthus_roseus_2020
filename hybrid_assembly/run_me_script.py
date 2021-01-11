#!bin/bash 

#Dont run
#mv barcode01.fastq.gz P1F2Int.fastq.gz
#mv barcode02.fastq.gz P2F2Int.fastq.gz
#mv barcode03.fastq.gz P3F2Int.fastq.gz
#mv barcode04.fastq.gz P1F3Int.fastq.gz
#mv barcode05.fastq.gz P2F3Int.fastq.gz
#mv barcode06.fastq.gz P3F3Int.fastq.gz
#mv barcode07.fastq.gz P1F2Ext.fastq.gz
#mv barcode08.fastq.gz P2F2Ext.fastq.gz
#mv barcode09.fastq.gz P3F2Ext.fastq.gz
#mv barcode10.fastq.gz P1F3Ext.fastq.gz
#mv barcode11.fastq.gz P2F3Ext.fastq.gz
#mv barcode12.fastq.gz P3F3Ext.fastq.gz

#gzip -d * 

##############

#change this accordingly 
trimomaticPath='/home/up201505534/software/Trimmomatic-0.39'
nano_folder='00_data' 


#Step one, pre-process read data
#Orient reads, trim small reads and output fastqc files 
#runfastqc


mkdir -p 01_process_reads
mkdir -p 01_process_reads/fasqc_raw/

fastqc $nano_folder/P1F2Int.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P2F2Int.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P3F2Int.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P1F3Int.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P2F3Int.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P3F3Int.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P1F2Ext.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P2F2Ext.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P3F2Ext.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P1F3Ext.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P2F3Ext.fastq -o 01_process_reads/fasqc_raw/ -f fastq
fastqc $nano_folder/P3F3Ext.fastq -o 01_process_reads/fasqc_raw/ -f fastq


#Orient reads
mkdir -p 01_process_reads/pychopper/
mkdir -p 01_process_reads/pychopper/reports

cdna_classifier.py $nano_folder/P1F2Int.fastq 01_process_reads/pychopper/P1F2Int.oriented.fastq -r 01_process_reads/pychopper/reports/P1F2Int_report.pdf -w 01_process_reads/pychopper/P1F2Int_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P1F2Int_stat.tsv -u 01_process_reads/pychopper/P1F2Int_unclass_reads.fastq
cdna_classifier.py $nano_folder/P2F2Int.fastq 01_process_reads/pychopper/P2F2Int.oriented.fastq -r 01_process_reads/pychopper/reports/P2F2Int_report.pdf -w 01_process_reads/pychopper/P2F2Int_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P2F2Int_stat.tsv -u 01_process_reads/pychopper/P2F2Int_unclass_reads.fastq
cdna_classifier.py $nano_folder/P3F2Int.fastq 01_process_reads/pychopper/P3F2Int.oriented.fastq -r 01_process_reads/pychopper/reports/P3F2Int_report.pdf -w 01_process_reads/pychopper/P3F2Int_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P3F2Int_stat.tsv -u 01_process_reads/pychopper/P3F2Int_unclass_reads.fastq
cdna_classifier.py $nano_folder/P1F3Int.fastq 01_process_reads/pychopper/P1F3Int.oriented.fastq -r 01_process_reads/pychopper/reports/P1F3Int_report.pdf -w 01_process_reads/pychopper/P1F3Int_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P1F3Int_stat.tsv -u 01_process_reads/pychopper/P1F3Int_unclass_reads.fastq
cdna_classifier.py $nano_folder/P2F3Int.fastq 01_process_reads/pychopper/P2F3Int.oriented.fastq -r 01_process_reads/pychopper/reports/P2F3Int_report.pdf -w 01_process_reads/pychopper/P2F3Int_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P2F3Int_stat.tsv -u 01_process_reads/pychopper/P2F3Int_unclass_reads.fastq
cdna_classifier.py $nano_folder/P3F3Int.fastq 01_process_reads/pychopper/P3F3Int.oriented.fastq -r 01_process_reads/pychopper/reports/P3F3Int_report.pdf -w 01_process_reads/pychopper/P3F3Int_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P3F3Int_stat.tsv -u 01_process_reads/pychopper/P3F3Int_unclass_reads.fastq
cdna_classifier.py $nano_folder/P1F2Ext.fastq 01_process_reads/pychopper/P1F2Ext.oriented.fastq -r 01_process_reads/pychopper/reports/P1F2Ext_report.pdf -w 01_process_reads/pychopper/P1F2Ext_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P1F2Ext_stat.tsv -u 01_process_reads/pychopper/P1F2Ext_unclass_reads.fastq
cdna_classifier.py $nano_folder/P2F2Ext.fastq 01_process_reads/pychopper/P2F2Ext.oriented.fastq -r 01_process_reads/pychopper/reports/P2F2Ext_report.pdf -w 01_process_reads/pychopper/P2F2Ext_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P2F2Ext_stat.tsv -u 01_process_reads/pychopper/P2F2Ext_unclass_reads.fastq
cdna_classifier.py $nano_folder/P3F2Ext.fastq 01_process_reads/pychopper/P3F2Ext.oriented.fastq -r 01_process_reads/pychopper/reports/P3F2Ext_report.pdf -w 01_process_reads/pychopper/P3F2Ext_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P3F2Ext_stat.tsv -u 01_process_reads/pychopper/P3F2Ext_unclass_reads.fastq
cdna_classifier.py $nano_folder/P1F3Ext.fastq 01_process_reads/pychopper/P1F3Ext.oriented.fastq -r 01_process_reads/pychopper/reports/P1F3Ext_report.pdf -w 01_process_reads/pychopper/P1F3Ext_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P1F3Ext_stat.tsv -u 01_process_reads/pychopper/P1F3Ext_unclass_reads.fastq
cdna_classifier.py $nano_folder/P2F3Ext.fastq 01_process_reads/pychopper/P2F3Ext.oriented.fastq -r 01_process_reads/pychopper/reports/P2F3Ext_report.pdf -w 01_process_reads/pychopper/P2F3Ext_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P2F3Ext_stat.tsv -u 01_process_reads/pychopper/P2F3Ext_unclass_reads.fastq
cdna_classifier.py $nano_folder/P3F3Ext.fastq 01_process_reads/pychopper/P3F3Ext.oriented.fastq -r 01_process_reads/pychopper/reports/P3F3Ext_report.pdf -w 01_process_reads/pychopper/P3F3Ext_rescued_reads.fastq -S 01_process_reads/pychopper/reports/P3F3Ext_stat.tsv -u 01_process_reads/pychopper/P3F3Ext_unclass_reads.fastq

#concatenate the rescue reads with the oriented reads 
#rescue reads are also oriented aparently
mkdir -p 01_process_reads/pychopper_oriented_all

cat 01_process_reads/pychopper/P1F2Int.oriented.fastq 01_process_reads/pychopper/P1F2Int_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P1F2Int.all.oriented.fastq
cat 01_process_reads/pychopper/P2F2Int.oriented.fastq 01_process_reads/pychopper/P2F2Int_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P2F2Int.all.oriented.fastq
cat 01_process_reads/pychopper/P3F2Int.oriented.fastq 01_process_reads/pychopper/P3F2Int_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P3F2Int.all.oriented.fastq
cat 01_process_reads/pychopper/P1F3Int.oriented.fastq 01_process_reads/pychopper/P1F3Int_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P1F3Int.all.oriented.fastq
cat 01_process_reads/pychopper/P2F3Int.oriented.fastq 01_process_reads/pychopper/P2F3Int_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P2F3Int.all.oriented.fastq
cat 01_process_reads/pychopper/P3F3Int.oriented.fastq 01_process_reads/pychopper/P3F3Int_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P3F3Int.all.oriented.fastq
cat 01_process_reads/pychopper/P1F2Ext.oriented.fastq 01_process_reads/pychopper/P1F2Ext_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P1F2Ext.all.oriented.fastq
cat 01_process_reads/pychopper/P2F2Ext.oriented.fastq 01_process_reads/pychopper/P2F2Ext_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P2F2Ext.all.oriented.fastq
cat 01_process_reads/pychopper/P3F2Ext.oriented.fastq 01_process_reads/pychopper/P3F2Ext_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P3F2Ext.all.oriented.fastq
cat 01_process_reads/pychopper/P1F3Ext.oriented.fastq 01_process_reads/pychopper/P1F3Ext_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P1F3Ext.all.oriented.fastq
cat 01_process_reads/pychopper/P2F3Ext.oriented.fastq 01_process_reads/pychopper/P2F3Ext_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P2F3Ext.all.oriented.fastq
cat 01_process_reads/pychopper/P3F3Ext.oriented.fastq 01_process_reads/pychopper/P3F3Ext_rescued_reads.fastq > 01_process_reads/pychopper_oriented_all/P3F3Ext.all.oriented.fastq


##There still are small reads < 50 bps in this dataset
#run trimmomatic 
mkdir -p 01_process_reads/trimmed_annotation/


java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P1F2Int.all.oriented.fastq 01_process_reads/trimmed_annotation/P1F2Int.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P2F2Int.all.oriented.fastq 01_process_reads/trimmed_annotation/P2F2Int.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P3F2Int.all.oriented.fastq 01_process_reads/trimmed_annotation/P3F2Int.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P1F3Int.all.oriented.fastq 01_process_reads/trimmed_annotation/P1F3Int.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P2F3Int.all.oriented.fastq 01_process_reads/trimmed_annotation/P2F3Int.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P3F3Int.all.oriented.fastq 01_process_reads/trimmed_annotation/P3F3Int.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P1F2Ext.all.oriented.fastq 01_process_reads/trimmed_annotation/P1F2Ext.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P2F2Ext.all.oriented.fastq 01_process_reads/trimmed_annotation/P2F2Ext.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P3F2Ext.all.oriented.fastq 01_process_reads/trimmed_annotation/P3F2Ext.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P1F3Ext.all.oriented.fastq 01_process_reads/trimmed_annotation/P1F3Ext.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P2F3Ext.all.oriented.fastq 01_process_reads/trimmed_annotation/P2F3Ext.oriented.filtered.annot.fastq MINLEN:50
java -jar $trimomaticPath/trimmomatic-0.39.jar SE -threads 8 -phred33 01_process_reads/pychopper_oriented_all/P3F3Ext.all.oriented.fastq 01_process_reads/trimmed_annotation/P3F3Ext.oriented.filtered.annot.fastq MINLEN:50


#Run fastqc on these results
mkdir -p 01_process_reads/trimmed_annotation_fastqc

fastqc 01_process_reads/trimmed_annotation/P1F2Int.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P2F2Int.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P3F2Int.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P1F3Int.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P2F3Int.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P3F3Int.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P1F2Ext.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P2F2Ext.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P3F2Ext.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P1F3Ext.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P2F3Ext.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq
fastqc 01_process_reads/trimmed_annotation/P3F3Ext.oriented.filtered.annot.fastq -o 01_process_reads/trimmed_annotation_fastqc/ -f fastq



###For ILLUMINA reads----
##SRpipeline computed on the fac computers, needs stringtie2 
mkdir -p 02_mapping/ 

#idio
create_rna_sr.py -1 00_data/idio_1/Nextera_R1_cut_paired.fastq -2 00_data/idio_1/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/idio_1 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/idio_2/Nextera_R1_cut_paired.fastq -2 00_data/idio_2/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/idio_2 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/idio_3/Nextera_R1_cut_paired.fastq -2 00_data/idio_3/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/idio_3 --frag-len 315 --frag-std 145

#meso 
create_rna_sr.py -1 00_data/meso_1/Nextera_R1_cut_paired.fastq -2 00_data/meso_1/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/meso_1 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/meso_2/Nextera_R1_cut_paired.fastq -2 00_data/meso_2/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/meso_2 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/meso_3/Nextera_R1_cut_paired.fastq -2 00_data/meso_3/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/meso_3 --frag-len 315 --frag-std 145

#totalptt
create_rna_sr.py -1 00_data/totalptt_1/Nextera_R1_cut_paired.fastq -2 00_data/totalptt_1/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/totalptt_1 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/totalptt_2/Nextera_R1_cut_paired.fastq -2 00_data/totalptt_2/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/totalptt_2 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/totalptt_3/Nextera_R1_cut_paired.fastq -2 00_data/totalptt_3/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/totalptt_3 --frag-len 315 --frag-std 145

#leaves
create_rna_sr.py -1 00_data/leaves_1/Nextera_R1_cut_paired.fastq -2 00_data/leaves_1/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/leaves_1 --frag-len 315 --frag-std 145
create_rna_sr.py -1 00_data/leaves_2/Nextera_R1_cut_paired.fastq -2 00_data/leaves_2/Nextera_R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 36 -o 02_mapping/leaves_2 --frag-len 315 --frag-std 145


#F1Int 
create_rna_sr.py -1 00_data/P1F1Int/R1_cut_paired.fastq -2 00_data/P1F1Int/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P1F1Int --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P2F1Int/R1_cut_paired.fastq -2 00_data/P2F1Int/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P2F1Int --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P3F1Int/R1_cut_paired.fastq -2 00_data/P3F1Int/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P3F1Int --frag-len 300 --frag-std 120

#F1Ext 
create_rna_sr.py -1 00_data/P1F1Ext/R1_cut_paired.fastq -2 00_data/P1F1Ext/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P1F1Ext --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P2F1Ext/R1_cut_paired.fastq -2 00_data/P2F1Ext/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P2F1Ext --frag-len 300 --frag-std 120 
create_rna_sr.py -1 00_data/P3F1Ext/R1_cut_paired.fastq -2 00_data/P3F1Ext/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P3F1Ext --frag-len 300 --frag-std 120


#F4Int 
create_rna_sr.py -1 00_data/P1F4Int/R1_cut_paired.fastq -2 00_data/P1F4Int/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P1F4Int --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P2F4Int/R1_cut_paired.fastq -2 00_data/P2F4Int/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P2F4Int --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P3F4Int/R1_cut_paired.fastq -2 00_data/P3F4Int/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P3F4Int --frag-len 300 --frag-std 120

#F4Ext 
create_rna_sr.py -1 00_data/P1F4Ext/R1_cut_paired.fastq -2 00_data/P1F4Ext/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P1F4Ext --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P2F4Ext/R1_cut_paired.fastq -2 00_data/P2F4Ext/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P2F4Ext --frag-len 300 --frag-std 120
create_rna_sr.py -1 00_data/P3F4Ext/R1_cut_paired.fastq -2 00_data/P3F4Ext/R2_cut_paired.fastq -H indexes/cro_index -G cro_index_gmap -p 4 -o 02_mapping/P3F4Ext --frag-len 300 --frag-std 120


#merge mapping files use stringtie
samtools merge 02_mapping/Illumina_all_samples_merged.bam 02_mapping/idio_1/sr_merge.bam 02_mapping/idio_2/sr_merge.bam 02_mapping/idio_3/sr_merge.bam 02_mapping/leaves_1/sr_merge.bam 02_mapping/leaves_2/sr_merge.bam 02_mapping/meso_1/sr_merge.bam 02_mapping/meso_2/sr_merge.bam 02_mapping/meso_3/sr_merge.bam 02_mapping/totalptt_1/sr_merge.bam 02_mapping/totalptt_2/sr_merge.bam 02_mapping/totalptt_3/sr_merge.bam 02_mapping/P1F1Ext/sr_merge.bam 02_mapping/P2F1Ext/sr_merge.bam 02_mapping/P3F1Ext/sr_merge.bam 02_mapping/P1F4Ext/sr_merge.bam 02_mapping/P2F4Ext/sr_merge.bam 02_mapping/P3F4Ext/sr_merge.bam 02_mapping/P1F1Int/sr_merge.bam P2F1Int/sr_merge.bam 02_mapping/P3F1Int/sr_merge.bam 02_mapping/P1F4Int/sr_merge.bam P2F4Int/sr_merge.bam 02_mapping/P3F4Int/sr_merge.bam -@ 20
samtools sort 02_mapping/Illumina_all_samples_merged.bam > 02_mapping/Illumina_all_samples_sorted.bam

#run stringtie 
stringtie 02_mapping/Illumina_all_samples_sorted.bam  -j 50 --conservative -o 03_stringtie/Illumina_all_samples.gtf -p 20 -l Illumina_Catha -G cro_v2.gene_models.gff3



####FOR NANOPORE READS
#cocatenate all samples in one file
cat 01_process_reads/trimmed_annotation/* > 01_process_reads/trimmed_annotation/all_reads.oriented.filtered.fastq

#build minimap index and align (-uf option for splicing sites only on the same strand) 
minimap2 -I 1000G -d catharanthus_genome_v2/cro_v2_asm.mmi catharanthus_genome_v2/cro_v2_asm.fasta

minimap2 -t 8 -ax splice -uf catharanthus_genome_v2/cro_v2_asm.mmi 01_process_reads/trimmed_annotation/all_reads.oriented.filtered.fastq > all_reads.aligned.sam
samtools view all_reads.aligned.sam -q 40 -F 2304 -Sb | samtools sort -@ 8 -o 02_mapping/all_reads.mapped.sorted.bam

#Run stringtie with the nanoporeSample

#Run stringtie in the assembly of the nanopore reads
stringtie --rf -l CATHA_NANO -L -v -p 8 -c 1.5 02_mapping/all_reads.mapped.sorted.bam -o 03_stringtie/nanopore_transcriptome.gtf -G catharanthus_genome_v2/cro_v2.gene_models.gff3



#####MERGING ASSEMBLIES
#Run --merge with the illumina assemblies to get the merged_transcriptome
stringtie --merge 03_stringtie/Illumina_all_samples.gtf 03_stringtie/nanopore_transcriptome.gtf -l CATHA_NANO -G catharanthus_genome_v2/cro_v2.gene_models.gff3 -o 03_stringtie/hybrid_transcriptome.gtf

#Run gff compare and filter based on the class codes
mkdir -p 04_annotation/hybrid

gffcompare -o 04_annotation/hybrid_transcriptome.gtf -r catharanthus_genome_v2/cro_v2.gene_models.gff3 -V  03_stringtie/hybrid_transcriptome.gtf 

python3 scripts/filter_gff.py -in 04_annotation/hybrid/hybrid_transcriptome.annotated.gtf -out 04_annotation/hybrid


#run gffread to obtain the transcriptome
gffread -g catharanthus_genome_v2/cro_v2_asm.fasta -w 04_annotation/hybrid_transcriptome.fasta 04_annotation/hybrid/hybrid_transcriptome.filtered.gtf


####assemble illumina and nanopore transcriptome to test with busco

##Illumina with reference 

stringtie --merge 03_stringtie/Illumina_all_samples.gtf -l CATHA_Illumina -G catharanthus_genome_v2/cro_v2.gene_models.gff3 -o 03_stringtie/Illumina_merged_with_reference.gtf

Run gff compare and filter based on the class codes
mkdir -p 04_annotation/Illumina

gffcompare -o 04_annotation/illumina_transcriptome.gtf -r catharanthus_genome_v2/cro_v2.gene_models.gff3 -V  03_stringtie/Illumina_merged_with_reference.gtf

python3 scripts/filter_gff.py -in 04_annotation/gffcompare/illumina_transcriptome.gtf -out 04_annotation/illumina


run gffread to obtain the transcriptome
gffread -g catharanthus_genome_v2/cro_v2_asm.fasta -w 04_annotation/illumina_transcriptome.fasta 04_annotation/illumina/illumina_transcriptome.filtered.gtf

##Nanopore transcriptome
stringtie --merge 03_stringtie/nanopore_transcriptome.gtf -l CATHA_NANO -G catharanthus_genome_v2/cro_v2.gene_models.gff3 -o 03_stringtie/nanopore_with_referece.gtf

#Run gff compare and filter based on the class codes
mkdir -p 04_annotation/nanopore

gffcompare -o 04_annotation/nanopore/nanopore_transcriptome.gtf -r catharanthus_genome_v2/cro_v2.gene_models.gff3 -V  03_stringtie/nanopore_with_referece.gtf

python3 scripts/filter_gff.py -in 04_annotation/nanopore/nanopore_transcriptome.gtf -out 04_annotation/nanopore


#run gffread to obtain the transcriptome
gffread -g catharanthus_genome_v2/cro_v2_asm.fasta -w 04_annotation/nanopore_transcriptome.fasta 04_annotation/nanopore/nanopore_transcriptome.filtered.gtf



#####Quantification##### 


#align with minimap2 to a transcriptome
mkdir -p 05_quantification/ 
mkdir -p 05_quantification/alignments 

minimap2 -t 8 -I 1000G -d 05_quantification/transcriptome_index.onlyStrand.mmi 04_annotation/hybrid_transcriptome.fasta

minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P1F2Int.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P1F2Int.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P2F2Int.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P2F2Int.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P3F2Int.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P3F2Int.oriented.aligned.bam

minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P1F3Int.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P1F3Int.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P2F3Int.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P2F3Int.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P3F3Int.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P3F3Int.oriented.aligned.bam


minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P1F2Ext.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P1F2Ext.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P2F2Ext.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P2F2Ext.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P3F2Ext.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P3F2Ext.oriented.aligned.bam

minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P1F3Ext.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P1F3Ext.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P2F3Ext.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P2F3Ext.oriented.aligned.bam
minimap2 -t 8 -ax map-ont -p 0.99 -N 100 05_quantification/transcriptome_index.onlyStrand.mmi 01_process_reads/trimmed_annotation/P3F3Ext.oriented.filtered.annot.fastq | samtools view -Sb > 05_quantification/alignments/P3F3Ext.oriented.aligned.bam


#quantify using salmon

mkdir -p 05_quantification/salmonQuant


salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P1F2Int.oriented.aligned.bam -o 05_quantification/salmonQuant/P1F2Int
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P2F2Int.oriented.aligned.bam -o 05_quantification/salmonQuant/P2F2Int
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P3F2Int.oriented.aligned.bam -o 05_quantification/salmonQuant/P3F2Int

salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P1F2Ext.oriented.aligned.bam -o 05_quantification/salmonQuant/P1F2Ext
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P2F2Ext.oriented.aligned.bam -o 05_quantification/salmonQuant/P2F2Ext
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P3F2Ext.oriented.aligned.bam -o 05_quantification/salmonQuant/P3F2Ext

salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P1F3Int.oriented.aligned.bam -o 05_quantification/salmonQuant/P1F3Int
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P2F3Int.oriented.aligned.bam -o 05_quantification/salmonQuant/P2F3Int
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P3F3Int.oriented.aligned.bam -o 05_quantification/salmonQuant/P3F3Int

salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P1F3Ext.oriented.aligned.bam -o 05_quantification/salmonQuant/P1F3Ext
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P2F3Ext.oriented.aligned.bam -o 05_quantification/salmonQuant/P2F3Ext
salmon quant --noErrorModel  -p 8 -t 04_annotation/hybrid_transcriptome.fasta -l SF -a 05_quantification/alignments/P3F3Ext.oriented.aligned.bam -o 05_quantification/salmonQuant/P3F3Ext


