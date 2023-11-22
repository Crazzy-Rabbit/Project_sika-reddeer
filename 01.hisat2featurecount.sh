#! /bin/bash
###############################################################################################################################
# /home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq
# ./hismap/sam
# ./hismap
###############################################################################################################################
# build index for hisat2

hisat2-build -p 4 GCF_910594005.1_mCerEla1.1_genomic.fna  GCF_910594005.1_mCerEla1.1_genomic > hisat2.log

###############################################################################################################################
# filter low quality reads and cut adapter use fastp
fastp="/home/sll/miniconda3/bin/fastp"

ls ML*/ML* | cut -d/ -f1 | uniq | while read id; do
    sample=${id}
    $fastp -i ${sample}/${sample}_1.clean.fq.gz -I ${sample}/${sample}_2.clean.fq.gz \
    -g -q 15 -n 5 -l 150 -u 50 \
    -o ${sample}_1.filter.fq.gz -O ${sample}_2.filter.fq.gz -h ${sample}.fastp.html
done

###############################################################################################################################
# align use hisat2
hisat2="/home/sll/miniconda3/bin/hisat2"
genomefa="/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic"
MAX_PROCESSES=4
RUNNING_PROCESSES=0
mkdir -p ./hismap/sam/

ls ML*/ML* | cut -d/ -f1 | uniq | while read id; do
    sample=${id}
    # if not exist sam that after map, then perform hisat2
    if [ ! -f ./hismap/sam/${sample}.hismap.sam ]; then
        $hisat2 -p 8 -x $genomefa -1 ${sample}_1.filter.fq.gz -2 ${sample}_2.filter.fq.gz -S ./hismap/sam/${sample}.hismap.sam &
        RUNNING_PROCESSES=$((RUNNING_PROCESSES+1))
    fi
    # if the maximum number of processes is reached, wait for them to finish
    if [ $RUNNING_PROCESSES -eq $MAX_PROCESSES ]; then
        wait
        RUNNING_PROCESSES=0
    fi
done

###############################################################################################################################
# sort and reformat use samtools
samtools="/home/sll/miniconda3/bin/samtools"
MAX_PROCESSES=4
RUNNING_PROCESSES=0

ls ML*/ML* | cut -d/ -f1 | uniq | while read id; do
    sample=${id}
    # if exist sam that after map, then perform samtools
    if [ -f ./hismap/sam/${sample}.hismap.sam ]; then
        $samtools sort -@ 8 -O bam -o ./hismap/${sample}_sort.bam ./hismap/sam/${sample}.hismap.sam &
        pid=$!
        RUNNING_PROCESSES=$((RUNNING_PROCESSES+1))
        if not kill -0 $pid 2>/dev/null; then
            $samtools index ./hismap/${sample}_sort.bam ./hismap/${sample}_sort.bam.index &
            RUNNING_PROCESSES=$((RUNNING_PROCESSES+1))
        fi
    fi
    # if the maximum number of processes is reached, wait for them to finish
    if [ $RUNNING_PROCESSES -eq $MAX_PROCESSES ]; then
        wait
        RUNNING_PROCESSES=0
    fi
done

###############################################################################################################################
# featurecount count the reads number
featureCounts="/home/sll/miniconda3/bin/featureCounts"
multiqc="/home/sll/miniconda3/bin/multiqc"
genomegtf="/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic.gtf"

$featureCounts -p -t exon -g gene_id -M -T 8 -a $genomegtf -o all.featurecounts.txt ./hismap.sam/*_sort.bam
# multiqc 
$multiqc all.featurecounts.txt.summary -o  all.counts.summary
# get counts to R
awk -F '\t' '{print $1,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' all.featurecounts.txt > all_fcount.matrix.txt
