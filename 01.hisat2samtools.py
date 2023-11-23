mkdir -p ./hismap/sam/ # remove this code when use the python script
#! usr/bin/python
##################################
### python method for muti process
### hisat2 align 
### samtools sort and index
### max process = 4
##################################
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

# set dir that save sam and bam file
sam_dir = "/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap/sam"
bam_dir = "/home/sll/5t_wgs_20230814_bam/20231118-deer-rna-seq/hismap"
hisat2 = "/home/sll/miniconda3/bin/hisat2"
samtools = "/home/sll/miniconda3/bin/samtools"
genomefa = "/home/ysq/20221108-deer-depth/20231007-deer-sift-data/reddeer-ref-mCerEla1.1/GCF_910594005.1_mCerEla1.1_genomic"

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

# set max processes
MAX_PROCESSES = 4

with ThreadPoolExecutor(max_workers=MAX_PROCESSES) as executor:
    executor.map(process_sample, samples)
