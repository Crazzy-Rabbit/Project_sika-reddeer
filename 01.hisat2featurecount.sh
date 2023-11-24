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
sam_dir="/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap/sam"

ls ML*/ML* | cut -d/ -f1 | uniq | while read id; do
    sample=${id}
    # if not exist sam that after map, then perform hisat2
    if [ ! -f ./hismap/sam/${sample}.hismap.sam ]; then
        $hisat2 -p 8 -x $genomefa -1 ${sample}_1.filter.fq.gz -2 ${sample}_2.filter.fq.gz -S ${sam_dir}/${sample}.hismap.sam &
        RUNNING_PROCESSES=$((RUNNING_PROCESSES+1))
    fi
    if [ $RUNNING_PROCESSES -eq $MAX_PROCESSES ]; then
        wait
        RUNNING_PROCESSES=0
    fi
done

###############################################################################################################################
# sort and reformat use samtools
samtools="/home/sll/miniconda3/bin/samtools"

# set dir that save sam and bam file
sam_dir="/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap/sam"
bam_dir="/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap"
MAX_PROCESSES=4
RUNNING_PROCESSES=0

ls ML*/ML* | cut -d/ -f1 | uniq | while read id; do
    sample=${id}
    if [ -f ${sam_dir}/${sample}.hismap.sam ]; then
        $samtools sort -@ 8 -O bam -o ${bam_dir}/${sample}_sort.bam ${sam_dir}/${sample}.hismap.sam &
        pid=$!
        RUNNING_PROCESSES=$((RUNNING_PROCESSES+1))
        if ! kill -0 $pid 2>/dev/null; then
            $samtools index ${bam_dir}/${sample}_sort.bam ${bam_dir}/${sample}_sort.bam.index &
        else
            RUNNING_PROCESSES=$((RUNNING_PROCESSES+1))
        fi
    fi
    if [ $RUNNING_PROCESSES -eq $MAX_PROCESSES ]; then
        wait
        RUNNING_PROCESSES=0
    fi
done    

#####################################################

##################################
### python methon for muti process
### hisat2 align 
### samtools sort and index
### max process = 4
##################################
#! usr/bin/python
import os
import subprocess
# muti process
from concurrent.futures import ThreadPoolExecutor

# set dir that save sam and bam file
sam_dir = "/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap/sam"
bam_dir = "/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap"
hisat2 = "/home/sll/miniconda3/bin/hisat2"
samtools = "/home/sll/miniconda3/bin/samtools"
featureCounts = "/home/sll/miniconda3/bin/featureCounts"
multiqc = "/home/sll/miniconda3/bin/multiqc"
genomefa = "/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic"
genomegtf = "/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic.gtf"

def process_sample(sample):
    fq_dir = f"/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/{sample}"
    fq1_file = os.path.join(fq_dir, f"{sample}_1.filter.fq.gz")
    fq2_file = os.path.join(fq_dir, f"{sample}_2.filter.fq.gz")
    sam_file = os.path.join(sam_dir, f"{sample}.hismap.sam")
    bam_file = os.path.join(bam_dir, f"{sample}_sort.bam")
    index_file = os.path.join(bam_dir, f"{sample}_sort.bam.index")

    # if there are fq1_file and fq2_file
    if os.path.isfile(fq1_file) and os.path.isfile(fq2_file):
        # align use hisat2
        subprocess.run([hisat2, "-p", "8", "-x", genomefa, "-1", fqt1_file, "-2", fqt2_file, "-S", sam_file])
        # sort sam file
        subprocess.run([samtools, "sort", "-@", "8", "-O", "bam", "-o", bam_file, sam_file])
        # index bam file
        subprocess.run([samtools, "index", bam_file, index_file]) 
    else:
        print("Please provide the filter.fq.gz file")

# samples = [id for id in os.listdir(sam_dir) if id.startswith("ML")]
# samples = subprocess.check_output("ls ML*/ML* | cut -d/ -f1 | uniq", shell=True).decode().splitlines()
samples = [id.split(".")[0] for id in os.listdir(sam_dir) if id.startswith("ML")]

os.system("mkdir -p ./hismap/sam/")

# set max processes and run
MAX_PROCESSES = 4
with ThreadPoolExecutor(max_workers=MAX_PROCESSES) as executor:
    executor.map(process_sample, samples) 

# run featurecount and multiQC after all file index
os.system("mkdir /home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap/counts")
os.system(f"$featureCounts -p -t exon -g gene_id -M -T 8 -a $genomegtf -o ./hismap/counts/all.featurecounts.txt *_sort.bam")
os.system(f"$multiqc ./hismap/counts/all.featurecounts.txt.summary -o  ./hismap/counts/all.counts.summary")
###############################################################################################################################
cd /home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap
# featurecount count the reads number
featureCounts="/home/sll/miniconda3/bin/featureCounts"
multiqc="/home/sll/miniconda3/bin/multiqc"
genomegtf="/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic.gtf"

$featureCounts -p -t exon -g gene_id -M -T 8 -a $genomegtf -o ./counts/all.featurecounts.txt *_sort.bam
# multiqc 
mkdir counts && cd counts
$multiqc all.featurecounts.txt.summary -o  all.counts.summary
# get counts to R
awk -F '\t' '{printf $1; for(i=6;i<=NF;i++) printf FS $i; printf "\n"}' OFS='\t' all.featurecounts.txt > all_fcount.matrix.txt
