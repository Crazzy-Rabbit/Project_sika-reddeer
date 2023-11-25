# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# if (!requireNamespace("DESeq2", quietly = TRUE))
# BiocManager::install("DESeq2")
# install.packages("tidyverse")
# install.packages("data.table")

#################################
# calculate FPKM
#################################
setwd("D:/百度网盘同步文件/BaiduSyncdisk/21-24吉大硕士/2021-2024实验及文章/2021-2024其他实验/2023.11-梅花鹿马鹿转录组/hisat2-featurecount-deseq2")
rm(list=ls()) 
options(stringsAsFactors = F)  

library(data.table) 

a1 <- fread('all_fcount.matrix.txt', header = T, data.table = F)
# reads count for gene
counts <- a1[,3:ncol(a1)] 
rownames(counts) <- a1$Geneid
# gene length
geneid_efflen <- subset(a1, select = c("Geneid","Length"))        
colnames(geneid_efflen) <- c("geneid","efflen")   
geneid_efflen_fc <- geneid_efflen
# order same as a1
efflen <- geneid_efflen[match(rownames(counts),                               
                              geneid_efflen$geneid),"efflen"] 

# FPKM/RPKM (Fragments/Reads Per Kilobase Million)
counts2FPKM <- function(count=count, efflength=efflen){    
  PMSC_counts <- sum(count)/1e6   # per million scaling factor    depth scale  
  FPM <- count/PMSC_counts        # Reads/Fragments Per Million   length scale  
  FPM/(efflength/1000)                                       
}

# apply counts2FPKM to each sample
FPKM <- as.data.frame(apply(counts,2,counts2FPKM))
FPKM <- FPKM[rowSums(FPKM)>1,]

#################################
# DEG use DEseq2
#################################
library(DESeq2)

rm(list=ls()) 
rowcount <- read.table("all_fcount.matrix.txt", header=TRUE, row.names = 1)
count <- rowcount[,2:ncol(rowcount)] 
##  filter: the sum count of gene < 10 
count.filter<-count[rowSums(count)>=1&apply(count,1,function(x){all(x>=10)}),]
write.table(count.filter,"read.count.filter.txt", quote = F, sep = "\t")

## establish dds matrix
condition <- factor(c("HM", "MH", "MH", "HM", "MH", "MH", rep("HL",5), "HM","HM","HL", rep("ML",4)))
coldata <- data.frame(row.names=colnames(count.filter), condition)
dds <- DESeqDataSetFromMatrix(countData=count.filter, colData=coldata, design=~condition)

## PCA use rld
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup=c("condition"))

## DEseq2
deseq <- DESeq(dds)
resdata<- results(deseq,contrast = c("condition","HL","HM"))
# show Benjamini-Hochberg correct adjust p < 0.05
table(resdata$padj<0.05) 
resdata <- as.data.frame(resdata)
res_padj <- resdata[order(resdata$padj), ]
# define down and up gene
res_padj <- res_padj %>% 
  mutate(change = case_when(padj < 0.05 & log2FoldChange > 1 ~ "UP",
                            padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
                            TRUE ~ "NOT"))
write.table(res_padj,"HL-HM_diffexpr_padj_results.txt",quote = F,sep = '\t')

# vocano plot
library(ggplot2)
library(patchwork)
degdata1 <- read.table("ML-HL_diffexpr_padj_results.txt",header = T,sep = '\t',row.names = 1)
degdata1$label <- c(rownames(degdata1)[1:10] ,rep(NA,(nrow(degdata1)-10)))
degdata2 <- read.table("HL-MH_diffexpr_padj_results.txt",header = T,sep = '\t',row.names = 1)
degdata2$label <- c(rownames(degdata2)[1:10] ,rep(NA,(nrow(degdata2)-10)))
degdata3 <- read.table("HL-HM_diffexpr_padj_results.txt",header = T,sep = '\t',row.names = 1)
degdata3$label <- c(rownames(degdata3)[1:10] ,rep(NA,(nrow(degdata3)-10)))

p1 = ggplot(degdata1, aes(x=log2FoldChange, y=-log10(padj), color=change))+
        geom_point(alpha=0.8, size=3)+
        geom_text(aes(label=c(label),color = "red"), 
                  size = 2, vjust = 1, hjust = 1)+ 
        labs(x="log2FC", y="-log10(FDR pvalue)")+
        scale_color_manual(values=c('#a121f0','#bebebe',"red",'#ffad21'))+
        theme_bw(base_size = 15)+
        theme(legend.position="none")+
        geom_vline(xintercept = 1,lty="dashed")+
        geom_vline(xintercept = -1,lty="dashed")+
        ggtitle("ML vs HL") 
p2 = ggplot(degdata2, aes(x=log2FoldChange, y=-log10(padj), color=change))+
        geom_point(alpha=0.8, size=3)+
        geom_text(aes(label=c(label),color = "red"), 
                  size = 2, vjust = 1, hjust = 1)+ 
        labs(x="log2FC", y="-log10(FDR pvalue)")+
        scale_color_manual(values=c('#a121f0','#bebebe',"red",'#ffad21'))+
        theme_bw(base_size = 15)+
        theme(legend.position="none")+
        geom_vline(xintercept = 1,lty="dashed")+
        geom_vline(xintercept = -1,lty="dashed")+
        ggtitle("HL vs MH") 

p3 = ggplot(degdata3, aes(x=log2FoldChange, y=-log10(padj), color=change))+
        geom_point(alpha=0.8, size=3)+
        geom_text(aes(label=c(label),color = "red"), 
                  size = 2, vjust = 1, hjust = 1)+ 
        labs(x="log2FC", y="-log10(FDR pvalue)")+
        scale_color_manual(values=c('#a121f0','#bebebe',"red",'#ffad21'))+
        theme_bw(base_size = 15)+
        theme(legend.position="none")+
        geom_vline(xintercept = 1,lty="dashed")+
        geom_vline(xintercept = -1,lty="dashed")+
        ggtitle("HL vs HM") 
  
pdf(file = "DESeq2.pdf",
    width = 15, height = 6)
p1 + p2 + p3 +
  plot_layout(guides = "collect")
dev.off()

# expression heatmap profile
library(pheatmap)
read_count <- read.table("read.count.filter.txt",header = T,sep = '\t',row.names = 1)
choose_gene=head(rownames(degdata1),100)  
choose_matrix=read_count[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))

png(filename = "DEG_pheatmap.png", width = 600, height = 1000)
pheatmap(choose_matrix,show_rownames=FALSE,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(20))
dev.off()






