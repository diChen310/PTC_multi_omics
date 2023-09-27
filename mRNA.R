
#################mRNA analysis, 2022-08-15###################
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggsci)
library(DESeq2)

####Read sample type, tumor or normal

sample.type = read.csv(file='./data/rnaseq/sampleInfo.csv',row.names = 1,check.names = F)

samples.tp = rownames(sample.type[sample.type$Type == 'P',])
samples.nt = rownames(sample.type[sample.type$Type == 'N',])

####Read gene count####
mrna.count = read.delim(file='data/rnaseq/mrna_count.txt',check.names = F)
# b = as.double(mrna.count[,'19-711'])
# a = as.double(mrna.count[,'19-711t'])
# c = as.double(mrna.count[,'146'])
# 
# length(a[a==0])#15887,5573
# length(b[b==0])#7427,4117
# length(c[c!=0])#7225 
mrna.count = mrna.count[,colnames(mrna.count) != '19-711t'] ####Remove poor data

many0.rows = apply(mrna.count, 1, function(a){length(a[a==0])})
many0.cols = apply(mrna.count, 2, function(a){length(a[a==0])})
many0.cols/nrow(mrna.count)
mrna.count = mrna.count[many0.rows/ncol(mrna.count)<0.6,]

dds <- DESeqDataSetFromMatrix(countData = mrna.count,
                              colData = sample.type[colnames(mrna.count),],
                              design= ~ Type)
dds <- DESeq(dds)
res.deseq <- results(dds)
#resLFC <- lfcShrink(dds, coef=2, type="apeglm")
mrna.count.n = counts(dds,normalized = T)
saveRDS(mrna.count.n,file='./variables/deseqNorm_count.Rds')
mrna.count.n = readRDS(file='./variables/deseqNorm_count.Rds')
# DESeq()
# counts(,normalized = TRUE)

res.c = data.frame(id = rownames(res.deseq),gene = rownames(res.deseq),
                   p.value = res.deseq$pvalue,FDR=res.deseq$padj,
                   Log2FC=res.deseq$log2FoldChange)
res = res.c[!is.na(res.c$p.value),]
res$significance = 'Not'
res[res$p.value<0.01 & res$Log2FC > 1,'significance']='Up'
res[res$p.value<0.01 & res$Log2FC < -1,'significance']='Down'
top25 <- res %>% filter(significance != 'Not') %>% group_by(significance) %>% top_n(-25, p.value)%>%  top_n(35, Log2FC)
rownames(res.tcga)=res.tcga$gene
top25$tcga.fc = res.tcga[as.character(top25$gene),'FC']

res.sig = filter(res,significance != 'Not')
res.sig$tcga.fc = res.tcga[as.character(res.sig$gene),'FC']
res.sig = filter(res.sig,!is.na(tcga.fc))
res.up.genes = unique(as.character(res[res$p.value<0.01 & res$Log2FC>0,'gene']))
res.down.genes = unique(as.character(res[res$p.value<0.01 & res$Log2FC<0,'gene']))


ggplot(res,aes(x=Log2FC,y=-log10(p.value),color=significance))+geom_point()+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+geom_hline(yintercept = 2,lty=3)+
  geom_text_repel(aes(x=Log2FC,y=-log10(p.value),label = gene,color=significance),top25,size = 2)+
  theme_cowplot()+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )
ggplot(res.sig,aes(x=Log2FC,y=tcga.fc,color=significance))+geom_point()+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))
#deseq2 count vs tcga.fc.tiff

#dds = DESeq2::counts(dds,normalized = T)
pcr.int = prcomp(t(log2(dds+1)),center = T,scale. = T)

pcr.cluster = data.frame(cluster = ifelse(colnames(dds) %in% samples.tp,'T','N'),
                         row.names  = colnames(dds),
                         PC1 = pcr.int$x[,1],
                         PC2 = pcr.int$x[,2])
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=cluster),size = 1.3)+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()

#dds = dds[,colnames(dds) != '146'] 
dds.r <- DESeqDataSetFromMatrix(countData = mrna.count[,c(samples.tp.r,samples.nt.r)],
                                colData = sample.type[c(samples.tp.r,samples.nt.r),],
                                design= ~ Type)
dds.r <- DESeq(dds.r)
res.r <- results(dds.r)
resLFC <- lfcShrink(dds.r, coef=2, type="apeglm")
res.mrna = data.frame(id = rownames(res.r),gene = rownames(res.r),
                      p.value = res.r$pvalue,FDR=res.r$padj,
                      Log2FC=res.r$log2FoldChange)
saveRDS(res.mrna,file = './variables/res.mrna.paired.Rds')





