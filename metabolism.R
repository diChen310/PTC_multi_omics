###Metabolism
options(stringsAsFactors = F)
mergeDuplicated = function(meta.exp,metas){
  dup.metas = metas[duplicated(metas)]
  dup.inx = c(1:length(metas))[metas %in% dup.metas]
  not.dup.inx = c(1:length(metas))[metas %in% dup.metas == F]
  not.dup.metas = metas[not.dup.inx]
  not.dup.mat = meta.exp[not.dup.inx,]
  rownames(not.dup.mat)=not.dup.metas
  
  new.mat = matrix(,nrow = length(dup.metas),ncol = ncol(meta.exp),
                   dimnames = list(dup.metas,
                                   colnames(meta.exp)))
  for(meta in dup.metas){
    dup.inx = c(1:length(metas))[metas == meta]
    new.mat[meta,]= apply(meta.exp[dup.inx,], 2, mean)
  }
  
  new.mat = rbind(new.mat,not.dup.mat)
  
  return(new.mat)
}


mergeSummary = function(meta.exp,short.name){
  metas = short.name[rownames(meta.exp),'short']
  dup.metas = unique(metas[duplicated(metas)])
  dup.inx = c(1:length(metas))[metas %in% dup.metas]
  not.dup.inx = c(1:length(metas))[metas %in% dup.metas == F]
  not.dup.metas = metas[not.dup.inx]
  not.dup.mat = meta.exp[not.dup.inx,]
  rownames(not.dup.mat)=not.dup.metas
  
  new.mat = matrix(,nrow = length(dup.metas),ncol = ncol(meta.exp),
                   dimnames = list(dup.metas,
                                   colnames(meta.exp)))
  for(meta in dup.metas){
    dup.inx = c(1:length(metas))[metas == meta]
    new.mat[meta,]= apply(meta.exp[dup.inx,], 2, sum)
  }
  
  new.mat = rbind(new.mat,not.dup.mat)
  
  return(new.mat)
}

library(randomForest)
itemImportance<-function(omics,cols,both.res,clu.s){
  
  X=t(omics[,cols])
  
  set.seed(71)
  rf <- randomForest(as.matrix(X),as.factor(both.res[,clu.s]), importance=TRUE,
                     proximity=TRUE)
  varImpPlot(rf,n.var=50)
  ## Look at variable importance:
  imp.s<-round(importance(rf), 2)
  return(imp.s)
}

getKEGGIDAndLogFC = function(res.r.nona){
  kegg.ids = c()
  log.fcs = c()
  
  kegg2logfc = filter(res.r.nona,kegg %in% c('','-') == F) %>% group_by(kegg) %>% mutate(mLog2FC = mean(Log2FC))
  
  for( i in 1:nrow(kegg2logfc)){
    kegg.id = kegg2logfc$kegg[i]
    logfc = kegg2logfc$mLog2FC[i]
    if(grepl(',',kegg.id)){
      ids = strsplit(kegg.id,',')[[1]]
      for(id in ids){
        kegg.ids = append(kegg.ids,id)
        log.fcs = append(log.fcs,logfc)
      }
    }else{
      kegg.ids = append(kegg.ids,kegg.id)
      log.fcs = append(log.fcs,logfc)
    }
  }
  
  rr = data.frame(kegg = kegg.ids,Log2FC = log.fcs)
  return(rr)
}


getKEGGIDAndName = function(res.r.nona){
  kegg.ids = c()
  names = c()
  

  for( i in 1:nrow(res.r.nona)){
    kegg.id = res.r.nona$kegg[i]
    name = res.r.nona$name[i]
    if(grepl(',',kegg.id)){
      ids = strsplit(kegg.id,',')[[1]]
      for(id in ids){
        kegg.ids = append(kegg.ids,id)
        names = append(names,name)
      }
    }else{
      kegg.ids = append(kegg.ids,kegg.id)
      names = append(names,name)
    }
  }
  
  rr = data.frame(kegg = kegg.ids,name = names)
  return(rr)
}

setwd('H:/project/multiOmics')

library(DMwR)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(ggsci)
library(dplyr)
data.meta = read.csv(file='./data/metabolism/Fudan_metabolites_send.csv',check.names = F)
short.name = read.csv(file='./data/metabolism/short name.csv',check.names = F)
short.name = short.name[!duplicated(short.name$`Compound name`),]
rownames(short.name)=short.name$`Compound name`

sample.info = read.csv(file='./data/metabolism/sample info.csv',check.names = F)

samples.tp.m = as.character(sample.info[sample.info$class=='T',2])
samples.nt.m = as.character(sample.info[sample.info$class!='T',2])


meta.exp = data.meta[,c(-1,-2,-3,-4)]
meta.exp = apply(meta.exp, 2, as.double)
meta.exp = matrix(data = as.double(meta.exp),nrow = nrow(meta.exp),
                  dimnames = list(data.meta$`Compound name`,
                                  c(colnames(data.meta)[c(-1,-2,-3,-4)])))
metas = data.meta$`Compound name`
rownames(meta.exp)=data.meta$`Compound name`

new.mat = mergeDuplicated(meta.exp,metas)
meta.info = data.meta[,c(1,2,3,4)]
meta.info = meta.info[!duplicated(meta.info$`Compound name`),]
meta.info$KEGG[meta.info$KEGG == '-']=meta.info$`Compound name`[meta.info$KEGG == '-']
meta.info$KEGG[meta.info$KEGG == '']=meta.info$`Compound name`[meta.info$KEGG == '']
rownames(meta.info)=meta.info$`Compound name`

manyNAs(meta.exp,0.1)
cc = colnames(meta.exp)[manyNAs(t(meta.exp),0.1)] ###"19-407" "19-709" "3449"  !!!!!!!!!!!!!!!!!

#notna.exp = meta.exp[,-manyNAs(t(meta.exp),0.0001)] ####27 samples without any NA

new.mat.nona = knnImputation(new.mat,scale = F)

new.mat.nona.s = mergeSummary(new.mat.nona,short.name )



################paired, fill the NA values################################
p.values.r <- lapply(1:nrow(new.mat.nona), function(i){
  wilcox.test(as.double(new.mat.nona[i,samples.tp.r]),as.double(new.mat.nona[i,samples.nt.r]),paired = T)$p.value
})

FCs.r <- lapply(1:nrow(new.mat.nona), function(i){
  mean(as.double(new.mat.nona[i,samples.tp.r]))/mean(as.double(new.mat.nona[i,samples.nt.r]))
})

res.r.nona <- data.frame(name = rownames(new.mat.nona),
                    kegg = meta.info[rownames(new.mat.nona),'KEGG'],
                    p.value = as.double(unlist(p.values.r)),
                    FDR = p.adjust(as.double(unlist(p.values.r)),method = 'fdr'),
                    Log2FC = log2(as.double(unlist(FCs.r))))

#res.r.nona = res.r.nona[!is.na(res.r$p.value),]
res.r.nona$significance = 'Not'
res.r.nona[res.r.nona$p.value<0.05 & res.r.nona$Log2FC > 1,'significance']='Up'
res.r.nona[res.r.nona$p.value<0.05 & res.r.nona$Log2FC < -1,'significance']='Down'
top25.r.nona <- arrange(res.r.nona,p.value,desc(Log2FC),kegg) %>% filter(significance != 'Not' & !duplicated(kegg)) %>% group_by(significance) %>% top_n(-10, p.value) %>% top_n(10, Log2FC)



ggplot(res.r.nona,aes(x=Log2FC,y=-log10(p.value),color=significance))+geom_point()+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+geom_hline(yintercept = -log10(0.05),lty=3)+
  geom_text_repel(aes(x=Log2FC,y=-log10(p.value),label = name,color=significance),top25.r.nona,size = 2)+
  theme_cowplot()+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )

rownames(res.r.nona)=res.r.nona$name
res.r = res.r[!duplicated(res.r$name),]
rownames(res.r)=res.r$name

res.r.nona$Log2FC.na.ignore = res.r[rownames(res.r.nona),'Log2FC']
plot(res.r.nona$Log2FC.na.ignore,res.r.nona$Log2FC)

table(res.r.nona$significance)
sig.metas = rownames(res.r.nona[res.r.nona$p.value<0.01,])


saveRDS(res.r.nona,file = './variables/all37Paired_diffMetaRes.Rds')

###########Pathway enrichment ############
load(file = './variables/KEGG_path2Meta.RData')

changeKegg2logfc = getKEGGIDAndLogFC(res.r.nona)
changeKegg2logfc = changeKegg2logfc[!duplicated(changeKegg2logfc$kegg),]
changeKegg2logfc = changeKegg2logfc[order(changeKegg2logfc$Log2FC,decreasing = T),]
FCs = changeKegg2logfc$Log2FC
names(FCs)=changeKegg2logfc$kegg
metabolism.pathways = path2cate[path2cate$cate =='1. Metabolism','Name']
diff.pathway.meta = GSEA(FCs,minGSSize = 3,maxGSSize = 200,pvalueCutoff = 1.2,pAdjustMethod = 'fdr',
                         TERM2GENE = path2Meta[path2Meta$Pathway %in% metabolism.pathways,c(1,2)],seed=1234)
gsea.res.m = diff.pathway.meta@result
#gsea.res.m = gsea.res.m[gsea.res.m$pvalue<0.05,]

gsea.top = gsea.res.m[gsea.res.m$pvalue<0.1,]
ggplot(gsea.top,aes(y=ID,x= NES,fill = -log10(pvalue)))+geom_bar(stat = 'identity',width = 0.6)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))


