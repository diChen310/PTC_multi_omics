####
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggsci)
library(MOVICS) #R4.0.4 only
library(survival)
library(survminer)
####Differential across RecurRisk
runDiff = function(omics,cols.1,cols.2){
  p.values <- lapply(1:nrow(omics), function(i){
    wilcox.test(as.double(omics[i,cols.1]),as.double(omics[i,cols.2]))$p.value
  })
  
  FCs <- lapply(1:nrow(omics), function(i){
    mean(as.double(omics[i,cols.1]))-mean(as.double(omics[i,cols.2]))
  })
  
  res <- data.frame(id = rownames(omics),
                    p.value = as.double(unlist(p.values)),
                    FDR = p.adjust(as.double(unlist(p.values)),method = 'fdr'),
                    FC = as.double(unlist(FCs)))
  return(res)
}

findDiffItems = function(omics,clin.both){
  
  clin.both = clin.both[rownames(clin.both) %in% colnames(omics),]
  high.cols = clin.both[clin.both$RecurRisk=='High','ID']
  m.cols = clin.both[clin.both$RecurRisk=='Median','ID']
  l.cols = clin.both[clin.both$RecurRisk=='Low','ID']
  
  res.high = runDiff(omics,high.cols,c(m.cols,l.cols))
  res.m = runDiff(omics,m.cols,c(high.cols,l.cols))
  res.l = runDiff(omics,l.cols,c(m.cols,high.cols))
  
  res.high$Risk = 'H'
  res.m$Risk='M'
  res.l$Risk='L'
  
  res = rbind(res.high,res.m,res.l)
  
  return(res)
  
}

runGSEA.i = function(path2Gene,res){
  res = arrange(res,desc(FC),p.value)
  FC.gsea = res$FC
  names(FC.gsea)=as.character(res$id)
  
  diff.pathway = GSEA(FC.gsea,minGSSize = 2,maxGSSize = 500,pvalueCutoff = 1.2,
                      TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                      seed=1234)
  #minGSSize=5 for omics other than metabolomics
  return(diff.pathway@result)
}

findDiffPathway = function(diff.res,path2Gene){
  res.high = diff.res[diff.res$Risk == 'H',]
  res.m = diff.res[diff.res$Risk == 'M',]
  res.l = diff.res[diff.res$Risk == 'L',]
  
  path.high =runGSEA.i(path2Gene,res.high)
  path.m = runGSEA.i(path2Gene,res.m)
  path.l = runGSEA.i(path2Gene,res.l)
  path.high$Risk = 'H'
  path.m$Risk='M'
  path.l$Risk='L'
  
  res = rbind(path.high,path.m,path.l)
  
  return(res)
  
}
getKEGGIDAndLogFC = function(res.r.nona){
  kegg.ids = c()
  log.fcs = c()
  froms = unique(res.r.nona$From)
  froms.new = c()
  risks = c('H','M','L')
  risks.new = c()
  for(fr in froms){
    for(ri in risks){
 
  kegg2logfc = filter(res.r.nona,From == fr & Risk == ri & kegg %in% c('','-') == F) %>% group_by(kegg) %>% mutate(mLog2FC = ifelse(mean(FC)>0,max(FC),min(FC)))
  
  for( i in 1:nrow(kegg2logfc)){
    kegg.id = kegg2logfc$kegg[i]
    logfc = kegg2logfc$mLog2FC[i]
    if(grepl(',',kegg.id)){
      ids = strsplit(kegg.id,',')[[1]]
      for(id in ids){
        kegg.ids = append(kegg.ids,id)
        log.fcs = append(log.fcs,logfc)
        froms.new = append(froms.new,fr)
        risks.new = append(risks.new,ri)
      }
    }else{
      kegg.ids = append(kegg.ids,kegg.id)
      log.fcs = append(log.fcs,logfc)
      froms.new = append(froms.new,fr)
      risks.new = append(risks.new,ri)
    }
   
  }
  
  }}
  
  rr = data.frame(id = kegg.ids,FC = log.fcs,From=froms.new,Risk =risks.new)
  return(rr)
}
data.meta = readRDS(file = './variables/normalizedMetabolism.Rds')##The names were changed into the same as mrna-count
data.mrna = readRDS(file = './variables/deseqNorm_count.Rds')
data.mrna = log2(data.mrna+1)
sample.type = read.csv(file='H:/project/multiOmics/data/rnaseq/sampleInfo.csv',row.names = 1,check.names = F)
samples.tp = rownames(sample.type[sample.type$Type == 'P',])

clin.input = read.csv(file='./data/clin.info.csv')
rownames(clin.input)=clin.input[,1]

both.samples = intersect(colnames(data.meta),colnames(data.mrna))
both.samples = intersect(both.samples,samples.tp)####92 samples

clin.both = clin.input[both.samples,]

gene.aov.p.values = unlist(lapply(1:nrow(data.mrna), function(a){
  clin.both$gene = data.mrna[a,both.samples]
  res = aov(gene ~ RecurRisk,data=clin.both)
  return(summary(res)[[1]][1,5])
}))
meta.aov.p.values = unlist(lapply(1:nrow(data.meta), function(a){
  clin.both$gene = data.meta[a,both.samples]
  res = aov(gene ~ RecurRisk,data=clin.both)
  return(summary(res)[[1]][1,5])
}))
#####Proteomics####
data.pro = read.csv(file='./data/proteomics/MWY-22-157D_武汉协和医院74例人组织蛋白组lablefree检测服务/MWY-22-157D/4.Quantitation/Expressed_annotation.csv',
                    row.names=1,check.names=F)

data.pro2 = read.csv(file='./data/proteomics/MWY-22-160D_武汉协和医院74例人组织蛋白磷酸化修饰lablefree检测服务/MWY-22-160D/4.Quantitation/Expressed_annotation.csv',
                     row.names=1,check.names=F)


data.pro.num.0 = apply(data.pro[,4:77], 1, function(a){length(a[a==0])})

data.pro = data.pro[data.pro.num.0/74<0.3,]#nrow:3147

data.pro2.num.0= apply(data.pro2[,4:77], 1, function(a){length(a[a==0])})

data.pro2 = data.pro2[data.pro2.num.0/74<0.3,]#nrow:652

pair.info = read.csv(file='./data/paired sample info.csv')

samples.tp.p = as.character(pair.info$Tumor)
samples.nt.p = as.character(pair.info$Normal)

pro.exp = as.matrix(data.pro[,c(samples.tp.p,samples.nt.p)])
rownames(pro.exp)=data.pro$Gene
pro.exp = log2(pro.exp+1)

pro2.exp = as.matrix(data.pro2[,c(samples.tp.p,samples.nt.p)])
rownames(pro2.exp) = data.pro2$Gene
pro2.exp = log2(pro2.exp+1)
# gene.res.aov = data.frame(gene = rownames(data.mrna),
#                           p.value=gene.aov.p.values)
# meta.res.aov = data.frame(gene = rownames(data.meta),
#                           p.value=meta.aov.p.values)

gene.res.risk = findDiffItems(data.mrna,clin.input)
meta.res.risk = findDiffItems(data.meta,clin.input)
pro.res.risk = findDiffItems(pro.exp,clin.input)
ph.res.risk = findDiffItems(pro2.exp,clin.input)
gene.res.risk$From = 'mRNA'
meta.res.risk$From = 'metabolite'
pro.res.risk$From = 'protein'
ph.res.risk$From = 'phosph'

gene.res.risk$name = gene.res.risk$id
meta.res.risk$name = meta.res.risk$id
pro.res.risk$name = paste0('pro',pro.res.risk$id)
ph.res.risk$name = paste0('ph',ph.res.risk$id)


res.risk = rbind(gene.res.risk,meta.res.risk,pro.res.risk,ph.res.risk)

risk.top = rbind(filter(res.risk,p.value<0.05 & FC >0) %>% group_by(Risk,From) %>% top_n(5,abs(FC)),
             filter(res.risk,p.value<0.05 & FC <0) %>% group_by(Risk,From) %>% top_n(5,abs(FC)))
risk.top$id = as.character(risk.top$id)
risk.top$id[risk.top$id ==  "TG(16:0_18:1_22:5);TG(18:1_18:2_20:3)"] = 'TG(16:0_18:1_22:5)' #short name
names.sig = arrange(risk.top,Risk,0-FC,p.value)$id
risk.top$id = factor(risk.top$id,levels = names.sig[!duplicated(names.sig)])
risk.top$Risk = factor(risk.top$Risk,levels = c('H','M','L'))

##Figure 3a
ggplot(risk.top,aes(x=id,y=Risk,colour=FC,size = -log10(p.value)))+geom_point()+
  facet_wrap(~From,nrow = 4,scales = 'free')+
  theme_bw()+
  theme(axis.text = element_text(size=8,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_colour_gradientn(colours = hcl.colors(10,palette = 'Blue-Red 3'))+
  scale_size_continuous(range = c(1,5))

#Supplementary figure Tumor to normal diff barplot
t2n.diff.pro = res.protein[paste(sig.pros2,'Pro'),c('gene','p.value','Log2FC')]
t2n.diff.pro$gene = factor(t2n.diff.pro$gene,levels = rev(sig.pros2))
ggplot(t2n.diff.pro,aes(fill=Log2FC,y=gene,x=-log10(p.value)))+geom_bar(stat = 'identity')+
  theme_bw()+
  theme(axis.text = element_text(size=8,colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(10,palette = 'Blue-Red 3'))

t2n.diff.phos = res.phosph[paste(sig.phos2,'Phosph'),c('gene','p.value','Log2FC')]
t2n.diff.phos$gene = factor(t2n.diff.phos$gene,levels = rev(sig.phos2))
ggplot(t2n.diff.phos,aes(fill=Log2FC,y=gene,x=-log10(p.value)))+geom_bar(stat = 'identity')+
  theme_bw()+
  theme(axis.text = element_text(size=8,colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(10,palette = 'Blue-Red 3'))
  
t2n.diff = res.mrna[sig.genes2,c('gene','p.value','Log2FC')]
t2n.diff$gene = factor(t2n.diff$gene,levels = rev(sig.genes2))
ggplot(t2n.diff,aes(fill=Log2FC,y=gene,x=-log10(p.value)))+geom_bar(stat = 'identity')+
  theme_bw()+
  theme(axis.text = element_text(size=8,colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(10,palette = 'Blue-Red 3'))

t2n.diff.m = res.meta[top.meta,c('name','p.value','Log2FC')]
t2n.diff.m$gene = factor(t2n.diff.m$name,levels = rev(top.meta))
ggplot(t2n.diff.m,aes(fill=Log2FC,y=gene,x=-log10(p.value)))+geom_bar(stat = 'identity')+
  theme_bw()+
  theme(axis.text = element_text(size=8,colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(10,palette = 'Blue-Red 3'))

#####RecurRisk relevant pathway ####
# MSIGDB.FILE = 'H:/project/multiomics/MSigDB_Hallmark.gmt'
# 
# path2Gene = read.gmt(MSIGDB.FILE)
load('./variables/KEGG_path2Gene.RData')
load(file = './variables/KEGG_path2Meta.RData')

path2cate <- read.csv(file='./data/KEGG/KEGG hsa pathway list.csv',row.names = 1)
rownames(path2cate) = path2cate$Name

gene.res.risk.path = findDiffPathway(gene.res.risk,path2Gene)
pro.res.risk.path = findDiffPathway(pro.res.risk,path2Gene)
ph.res.risk.path = findDiffPathway(ph.res.risk,path2Gene)

meta.res.risk$kegg = meta.info[as.character(meta.res.risk$id),'KEGG']
meta.res.risk.ch = getKEGGIDAndLogFC(meta.res.risk)
meta.res.risk.ch$p.value = 'N'
meta.res.risk.path = findDiffPathway(meta.res.risk.ch,path2Meta)



gene.res.risk.path$From = 'mRNA'
pro.res.risk.path$From = 'protein'
ph.res.risk.path$From = 'phosph'
meta.res.risk.path$From = 'metabolite'
res.risk.path = rbind(gene.res.risk.path,pro.res.risk.path,ph.res.risk.path,meta.res.risk.path)
res.risk.path.sig = res.risk.path[res.risk.path$pvalue<0.01,]
res.risk.path.sig.up = filter(res.risk.path.sig,abs(NES)>0) %>% group_by(Risk,From) %>% top_n(10,abs(NES))
res.risk.path.sig.up = arrange(res.risk.path.sig.up,From,Risk,NES)
up.paths = unique(res.risk.path.sig.up$ID)

up.paths.plot = res.risk.path[res.risk.path$ID %in% up.paths,]

up.paths.plot$ID = factor(up.paths.plot$ID ,levels = up.paths)
up.paths.plot$Risk = factor(up.paths.plot$Risk,levels = c('H','M','L'))
up.paths.plot$From = factor(up.paths.plot$From,levels = c('metabolite','mRNA','protein','phosph'))
up.paths.plot$cate = path2cate[as.character(up.paths.plot$ID),'cate']
up.paths.plot=up.paths.plot[up.paths.plot$cate!='6. Human Diseases',]

####Figure 3e and Supplementary figure
ggplot(up.paths.plot,aes(x=Risk,y=ID,fill=NES ))+geom_tile(colour='grey')+
  geom_text(aes(label=ifelse(pvalue<0.01 ,'*','')),size=3)+
  facet_grid(cate~From,scales = 'free',space = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,angle = 45),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(10,palette = 'Blue-Red 3'))



###############molecular network#################
library(reshape2)

top.meta = unique(risk.top[risk.top$From == 'metabolite','id'])
top.meta = as.character(top.meta$id)
top.meta[top.meta ==  'TG(16:0_18:1_22:5)'] ="TG(16:0_18:1_22:5);TG(18:1_18:2_20:3)" 

sig.genes2 = unique(risk.top[risk.top$From == 'mRNA','id'])
sig.genes2 = as.character(as.character(sig.genes2$id))
#sig.metas2 = unique(meta.res.risk[meta.res.risk$p.value<0.05,'id'])
sig.pros2 = unique(risk.top[risk.top$From == 'protein','id'])
sig.pros2 = as.character(as.character(sig.pros2$id))

sig.phos2 = unique(risk.top[risk.top$From == 'phosph','id'])
sig.phos2 = as.character(as.character(sig.phos2$id))
all.samples.t = as.character(pair.info.all$Tumor)
all.samples.n = as.character(pair.info.all$Normal)

high.cols = intersect(all.samples.t,clin.input$ID[clin.input$RecurRisk == 'High'])
m.cols = intersect(all.samples.t,clin.both$ID[clin.both$RecurRisk == 'Median'])
l.cols = intersect(all.samples.t,clin.both$ID[clin.both$RecurRisk == 'Low'])

X1 = data.mrna[sig.genes,c(high.cols,m.cols,l.cols,all.samples.n)]
X2 = data.meta[sig.metas,c(high.cols,m.cols,l.cols,all.samples.n)]

X3 = pro.exp[sig.pros2,colnames(X1)]
X4=pro2.exp[sig.phos2,colnames(X1)]
rownames(X3) = paste(rownames(X3),'Pro')
rownames(X4) = paste(rownames(X4),'Phosph')
ht.input = rbind(data.mrna[sig.genes2,colnames(X1)],
data.meta[top.meta,colnames(X1)],
X3,
X4)
ht.input = t(scale(t(ht.input)))
row.type = data.frame(Item = rownames(ht.input),
Type = c(rep('mRNA',length(sig.genes2)),
rep('metabolite',length(top.meta)),
rep('protein',length(sig.pros2)),
rep('phosphprotein',length(sig.phos2))),row.names = rownames(ht.input))
corr.mat =cor(t(ht.input),method = 'spearman')
corr.meta.mat = melt(corr.mat,measure.vars = colnames(corr.mat))
corr.meta.mat$type1 = row.type[corr.meta.mat$Var1,'Type']
corr.meta.mat$type2 = row.type[corr.meta.mat$Var2,'Type']

corr.meta.high = corr.meta.mat[abs(corr.meta.mat$value)>0.65,]
corr.meta.high = corr.meta.high[corr.meta.high$type1 != corr.meta.high$type2,]
library(igraph)
library(ggraph)
library(tidygraph)

g.net = graph_from_data_frame(corr.meta.high[,c(1,2)])
sub.g = as_tbl_graph(g.net,directed = F)
####Figure 3b, f
sub.g = activate(sub.g,nodes) %>%
mutate(From = ifelse( names(V(sub.g)) %in% sig.genes,'aGene','Metabolite')) %>%
mutate(Type = 
         row.type[names(V(sub.g)),'Type']) %>%
activate(edges) %>%
mutate(edge_weights = ifelse(corr.meta.high$value<0,'-','+'))
ggraph(sub.g,layout = "graphopt") +
geom_edge_link(aes(edge_colour = edge_weights),width=0.5,edge_alpha=0.6)+
scale_edge_colour_manual(values = c('blue','red'))+
geom_node_point(aes(colour = Type,shape=From),size=4) +
geom_node_text(aes(label = name),size=2, repel = TRUE)+
scale_color_manual(values =pal_npg('nrc')(4))

####Figure 3c,d,g,h
FFA=as.double(ht.input['FFA 24:2',])
ALOX5_pro=as.double(ht.input['ALOX5 Pro',])
coef = cor.test(FFA,ALOX5_pro,method = 'spearman')$estimate
p = cor.test(FFA,ALOX5_pro,method = 'spearman')$p.value
qplot(FFA,ALOX5_pro,geom=c('point','smooth'),method='lm',se=F)+ggtitle("Correlations")+theme_bw()+
annotate("text", label = paste0("R = ",signif(coef,2), "\nP-value = ",signif(p,2)), x = median(FFA)+1, y = min(ALOX5_pro)+0.3, size = 4)



