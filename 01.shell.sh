#! /bin/bash

###############################################################################################################################
# work dir: /home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq
###############################################################################################################################
fastp="/home/sll/miniconda3/bin/fastp"
hisat2="/home/sll/miniconda3/bin/hisat2"
samtools="/home/sll/miniconda3/bin/samtools"
featureCounts="/home/sll/miniconda3/bin/featureCounts"
genomefa="/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic"
genomegtf="/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic.gtf"
###############################################################################################################################

# build index for hisat2
hisat2-build -p 4 GCF_910594005.1_mCerEla1.1_genomic.fna  GCF_910594005.1_mCerEla1.1_genomic > hisat2.log

ls ML*/ML* | while read id;
do echo $id
sample=${id%%/*}

# QC use fastp
$fastp -i ${sample}/${sample}_1.clean.fq.gz -I ${sample}/${sample}_2.clean.fq.gz -g -q 15 -n 5 -l 150 -u 50 -o ${sample}_1.filter.fq.gz -O ${sample}_2.filter.fq.gz -h ${sample}.fastp.html

# reads mapping  ----hisat2
$hisat2 -p 8 -x $genomefa -1 ${sample}_1.filter.fq.gz -2 ${sample}_2.filter.fq.gz -S hismap/sam/${sample}.hismap.sam  

# sam to bam and sorted and index for *sort.bam file
$samtools sort -@ 8 -O bam -o ./hismap/${sample}_sort.bam ./hismap/sam/${sample}.hismap.sam
$samtools index ./hismap/ ${sample}_sort.bam ./hismap/${sample}_sort.bam.index

done

# featurecount count the reads number
$featureCounts -p -t exon -g gene_id -M -T 8 -a $genomegtf -o all.featurecounts.txt ./hismap.sam/*_sort.bam

# multiqc 
$multiqc all.featurecounts.txt.summary -o  all.counts.summary

# get counts to R 
awk -F '\t' '{print $1,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' all.featurecounts.txt > all_fcount.matrix.txt
