setwd('H:/project/multiOmics')

res.mrna = readRDS(file = './variables/res.mrna.paired.Rds')
res.protein = readRDS(file = './variables/res.protein.paired.Rds')

res.protein.diffTrends = res.protein[res.protein$mRNA.fc * res.protein$log2FC < 0,]
res.protein.diffTrends = res.protein.diffTrends[!is.na(res.protein.diffTrends$gene),]
rownames(res.mrna)=res.mrna$id
res.protein.diffTrends$mRNA.fc.p = res.mrna[res.protein.diffTrends$gene,'p.value']
res.protein.diffTrends.sig = res.protein.diffTrends[res.protein.diffTrends$p.value<0.01 & res.protein.diffTrends$mRNA.fc.p<0.01,]
res.protein.diffTrends.sig$E3 = ifelse(as.character(res.protein.diffTrends.sig$gene) %in% e3.genes,'E3','Others')
####Spearman correlations####
mrna.count = read.delim(file='data/rnaseq/mrna_count.txt',check.names = F)
mrna.count = mrna.count[,colnames(mrna.count) != '19-711t']
sample.type = read.csv(file='H:/project/multiOmics/data/rnaseq/sampleInfo.csv',row.names = 1,check.names = F)
dds <- DESeqDataSetFromMatrix(countData = mrna.count,
colData = sample.type[colnames(mrna.count),],
design= ~ Type)
dds <- DESeq(dds)

data.mrna = counts(dds,normalized = T)
data.mrna = log2(data.mrna+1)

data.pro = read.csv(file='./data/proteomics/MWY-22-157D_武汉协和医院74例人组织蛋白组lablefree检测服务/MWY-22-157D/4.Quantitation/Expressed_annotation.csv',
                    row.names=1,check.names=F)
data.pro.num.0 = apply(data.pro[,4:77], 1, function(a){length(a[a==0])})

data.pro = data.pro[data.pro.num.0/74<0.3,]#nrow:3147
pair.info = read.csv(file='H:/project/multiOmics/data/paired sample info.csv')

samples.tp.p = as.character(pair.info$Tumor)
samples.nt.p = as.character(pair.info$Normal)

pro.exp = as.matrix(data.pro[,c(samples.tp.p,samples.nt.p)])
rownames(pro.exp)=data.pro$Gene
pro.exp = log2(pro.exp+1)

####sample-wise correlations####

both.items = intersect(rownames(pro.exp),rownames(data.mrna))
both.cols.t = intersect(samples.tp.p,colnames(data.mrna))
both.cols.n = intersect(samples.nt.p,colnames(data.mrna))

corr.t = numeric(length = length(both.cols.t))
for(i in 1:length(both.cols.t)){
  sample.i = both.cols.t[i]
  sample.i.mrna = as.double(data.mrna[both.items,sample.i])
  sample.i.pro = as.double(pro.exp[both.items,sample.i])
  corr.t[i]=cor(sample.i.mrna,sample.i.pro,method = 'spearman')
}
boxplot(corr.t)


corr.n = numeric(length = length(both.cols.n))
for(i in 1:length(both.cols.n)){
  sample.i = both.cols.n[i]
  sample.i.mrna = as.double(data.mrna[both.items,sample.i])
  sample.i.pro = as.double(pro.exp[both.items,sample.i])
  corr.n[i]=cor(sample.i.mrna,sample.i.pro,method = 'spearman')
}
boxplot(corr.t,corr.n)

boxplot(c(corr.t,corr.n),col="#4DBBD5FF")
sample.corrs = c(corr.t,corr.n)
median(sample.corrs)#0.2937053
hist(sample.corrs,breaks = 30,col = "#4DBBD5FF",freq = F)
lines(density(na.omit(sample.corrs)),lwd=2,col="black")
sample.corrs.df = data.frame(row.names = both.cols.t,corr=corr.t)
clin.d.in$sample.corr = sample.corrs.df[clin.d.in$ID,'corr']
library(ggpubr)
library(ggplot2)
ggplot(clin.d.in,aes(x=RecurRisk,y=sample.corr,fill=RecurRisk))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))

ggplot(clin.d.in,aes(x=N.stage,y=sample.corr,fill=N.stage))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))


####Gene-wise correlations#####

corr.genes = cor(t(data.mrna[both.items,c(both.cols.n,both.cols.t)]),
                 t(pro.exp[both.items,c(both.cols.n,both.cols.t)]),method = 'spearman')

corr.genes.unique = diag(corr.genes)
boxplot(corr.genes.unique)
median(corr.genes.unique,na.rm = T)#0.08561368
hist(corr.genes.unique,breaks = 30,col = "#E64B35FF",freq = F)
lines(density(na.omit(corr.genes.unique)),lwd=2,col="black")

corr.genes.p = unlist(lapply(both.items, function(a){
  cor.test(as.double(data.mrna[a,c(both.cols.n,both.cols.t)]),
           as.double(pro.exp[a,c(both.cols.n,both.cols.t)]))$p.value
}))

corr.genes.df = data.frame(row.names = both.items,corr=corr.genes.unique,p.value=corr.genes.p)
corr.genes.df = corr.genes.df[!is.na(corr.genes.df$corr),]
corr.genes.df$FDR = p.adjust(corr.genes.df$p.value)

corr.genes.n = cor(t(data.mrna[both.items,both.cols.n]),
                 t(pro.exp[both.items,both.cols.n]),method = 'spearman')

corr.genes.unique.n = diag(corr.genes.n)
boxplot(corr.genes.unique.n)
hist(corr.genes.unique.n,breaks = 30,col = "#4DBBD5FF",freq = F)
lines(density(na.omit(corr.genes.unique.n)),lwd=2,col="black")

####pathway analysis####
load('E:/project/metaboliteProteinInteraction/variables/KEGG_path2Gene.RData')
path2cate <- read.csv(file='E:/data/KEGG/KEGG hsa pathway list.csv',row.names = 1)
rownames(path2cate) = path2cate$Name

path.corr.sum = numeric(length = length(unique(path2Gene$Pathway)))

for(i in 1:length(path.corr.sum)){
  path = unique(path2Gene$Pathway)[i]
  genes.i = path2Gene[path2Gene$Pathway == path,'Item']
  genes.i.in = intersect(genes.i,both.items)
  path.corr.sum[i] = median(corr.genes.df[genes.i.in,'corr'],na.rm = T)
}
path.corr.sum = data.frame(pathway = unique(path2Gene$Pathway),score = path.corr.sum)

path2gene.in = path2Gene[path2Gene$Item %in% both.items,]
path2gene.in$corr = corr.genes.df[as.character(path2gene.in$Item),'corr']

sig.genes.input = as.character(rownames(corr.genes.df[corr.genes.df$p.value<0.05,]))
library(clusterProfiler)
path.enrich.res = enricher(sig.genes.input,pvalueCutoff = 1.1,TERM2GENE = path2Gene[,c(1,2)],
                           minGSSize = 3)@result

corr.genes.df = corr.genes.df[order(corr.genes.df$corr,decreasing = T),]

gsea.input = corr.genes.df$corr
names(gsea.input)=rownames(corr.genes.df)
disease.pathways = as.character(path2cate[path2cate$cate == '6. Human Diseases','Name'])
path.gsea = GSEA(gsea.input,pvalueCutoff = 1.1,TERM2GENE = path2Gene[path2Gene$Pathway %in% disease.pathways == F,c(1,2)],
                 minGSSize = 3,maxGSSize =400,seed = 1357)
path.gsea.res = path.gsea@result

#gseaplot(path.gsea,geneSetID = c(path.gsea.res$ID[1:10]))
path.gsea.res.sig = path.gsea.res[path.gsea.res$pvalue<0.05,]
top = path.gsea.res.sig$ID
path.gsea.2 = GSEA(gsea.input,pvalueCutoff = 1.1,TERM2GENE = path2Gene[path2Gene$Pathway %in% top,c(1,2)],
                 minGSSize = 3,maxGSSize =400,seed = 1357)

enrichplot::ridgeplot(path.gsea.2,showCategory = length(top))+theme(axis.text.y = element_text(size=8),
                                       axis.text.x = element_text(size=8))
write.csv(path.gsea.res,file = './variables/GSEA gene-wise correlations.csv')
write.csv(corr.genes.df,file = './variables/gene-wise correlations.csv')

####Find E3 or DUB####
E3.dbi = read.delim(file='E:/project/ubq prediction/R space/R space/e3net ubp pairs.txt')
dub.dbi = read.csv(file='E:/project/USP analysis/R space/merged DUB-SUB pairs.csv')
E3.genes = unique(as.character(E3.dbi$e3.genes))
dub.genes = unique(as.character(dub.dbi$dub))
res.protein$DUB = ifelse(res.protein$gene %in% dub.genes,'DUB','Others')
res.protein.ubi.sig = res.protein[res.protein$p.value<0.05 & (res.protein$E3 == 'E3'|res.protein$DUB == 'DUB'),]
res.protein.ubi.sig = res.protein.ubi.sig[order(res.protein.ubi.sig$E3),]
sig.items = as.character(res.protein.ubi.sig$gene)

pro.exp.sig.items = pro.exp[sig.items,]

library(reshape2)

pro.exp.sig.items.data = melt(pro.exp.sig.items,measure.vars = colnames(pro.exp))
pro.exp.sig.items.data$Type = ifelse(pro.exp.sig.items.data$Var2 %in% samples.tp.p,'Tumor','Normal')


ggplot(pro.exp.sig.items.data,aes(x=Type,y=value,fill=Type))+geom_boxplot(outlier.size = 1)+
  facet_wrap(~Var1,ncol = 3)+
  stat_compare_means(comparisons = list(c(1,2)),size=3)+
  theme_cowplot()+ylim(c(0,5.5))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1))
####Find PPI network####

ppis = read.delim(file='E:/data/PPI/BIOGRID/BIOGRID-MV-Physical-2023/BIOGRID-MV-Physical-4.4.224.tab2.txt')

diff.items = as.character(res.protein.diffTrends.sig$gene)

ppis.re = ppis[ppis$Official.Symbol.Interactor.A %in% diff.items | ppis$Official.Symbol.Interactor.B %in% diff.items,]

ppis.re = ppis.re[ppis.re$Official.Symbol.Interactor.A != ppis.re$Official.Symbol.Interactor.B,]
ppis.re$two = paste(ppis.re$Official.Symbol.Interactor.A,ppis.re$Organism.Interactor.B)
ppis.re = ppis.re[!duplicated(ppis.re$two),]
ppis.re = ppis.re[ppis.re$Throughput == 'Low Throughput',]
library(igraph)
library(ggraph)
library(tidygraph)
sub.g = graph_from_data_frame(ppis.re[,c('Official.Symbol.Interactor.A','Official.Symbol.Interactor.B')])
sub.g = as_tbl_graph(sub.g,directed = F)
sub.g.degree = degree(sub.g)
degree.res.df = data.frame(gene = names(sub.g.degree),degree = as.numeric(sub.g.degree))
degree.res.df$inDiff = ifelse(degree.res.df$gene %in% diff.items,'DiffTrend','No')

ppis.res.PCNA = ppis.re[ppis.re$Official.Symbol.Interactor.B == 'PCNA' & ppis.re$Modification != '-',]
sub.g %>% 
  activate(nodes) %>%
  ggraph(layout = "gem") + 
  geom_edge_link(aes(colour = weight),width=3,edge_alpha=0.6)+
  scale_edge_colour_gradient2(low = 'blue',high = 'red')+
  geom_node_point(color='cyan',size=10) +
  geom_node_text(aes(label = name),size=8, repel = TRUE)+    theme_graph()
