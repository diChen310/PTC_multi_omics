library(DMwR)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(ggsci)
library(dplyr)

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
keggid2name = getKEGGIDAndName(res.r.nona)



#####Merge the differential analysis results from different omics#### 

###The differential results were from mRNA.R metabolism.R proteom.R
res.meta = readRDS(file = './variables/all37Paired_diffMetaRes.Rds')
res.mrna = readRDS(file = './variables/res.mrna.paired.Rds')
res.protein = readRDS(file = './variables/res.protein.paired.Rds')
res.phosph = readRDS(file = './variables/res.phosph.paired.Rds')

res = res.protein[!is.na(res.protein$p.value),]
res$significance = 'Not'
res[res$p.value<0.01 & res$log2FC > 1,'significance']='Up'
res[res$p.value<0.01 & res$log2FC < -1,'significance']='Down'
#top25 <- res%>% group_by(significance) %>% filter(significance != 'Not' & !duplicated(kegg))  %>% top_n(-10, p.value)%>%  top_n(10, Log2FC)
top25 <- res%>% group_by(significance) %>% filter(significance != 'Not')  %>%  top_n(10, abs(log2FC))

###Supplementary Figure S2
ggplot(res,aes(x=log2FC,y=-log10(p.value),color=significance))+geom_point()+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+geom_hline(yintercept = 2,lty=3)+
  geom_text_repel(aes(x=log2FC,y=-log10(p.value),label = gene,color=significance),top25,size = 2)+
  theme_cowplot()+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )

###############sig diff molecules pathway enrichment analysis,2022-12-09####


###################

sig.mrna = res.mrna[res.mrna$p.value<0.01,'gene']
sig.metas = res.meta[res.meta$p.value<0.01,'name']
sig.pros = res.protein[res.protein$p.value<0.01,'gene']
sig.phops = res.phosph[res.phosph$p.value<0.01,'gene']

sig.metas.ids = meta.info[sig.metas,'KEGG']
sig.metas.ids = unique(sig.metas.ids)
sig.metas.ids = unlist(lapply(sig.metas.ids, function(a){
  if(grepl(',',a)){return(strsplit(a,',')[[1]][1])}else{return(a)}
}))
sig.meta.paths.k = enricher(sig.metas.ids,minGSSize = 2,TERM2GENE = path2Meta[,c(1,2)])@result
length(sig.pros)
sig.meta.paths.k$cate = path2cate[sig.meta.paths.k$Description,'cate']
sig.meta.paths.sm = enricher(sig.metas.ids,minGSSize = 2,TERM2GENE = path2meta.smpdb.c)@result
length(sig.metas.ids)

disease.pathways.sm = as.character(path2meta.smpdb[path2meta.smpdb$Pathway.Subject=='Disease',
                                                   'name.short'])
sig.meta.paths.sm[sig.meta.paths.sm$ID %in% disease.pathways.sm,]
path.file <- './data/pathways/MsigDBC2-curated.all.v6.2.symbols.gmt'
path2gene.all = read.gmt(gmtfile = path.file)
path2gene.all = path2gene.all[!grepl('_UP',path2gene.all$ont),]
path2gene.all = path2gene.all[!grepl('_DN',path2gene.all$ont),]

sig.mrna.paths = enricher(sig.mrna,minGSSize = 3,TERM2GENE = path2gene.all)@result
sig.pro.paths = enricher(sig.pros,minGSSize = 3,TERM2GENE = path2gene.all)@result
sig.phops.paths = enricher(sig.phops,minGSSize = 3,TERM2GENE = path2gene.all)@result

sig.mrna.path.names = sig.mrna.paths[sig.mrna.paths$pvalue<0.05,'ID']
sig.pro.path.names = sig.pro.paths[sig.pro.paths$pvalue<0.05,'ID']
sig.phops.path.names = sig.phops.paths[sig.phops.paths$pvalue<0.05,'ID']

intersect(intersect(sig.mrna.path.names,sig.pro.path.names),sig.phops.path.names)


sig.mrna.paths.k = enricher(sig.mrna,minGSSize = 3,TERM2GENE = path2Gene[,c(1,2)])@result
sig.pro.paths.k = enricher(sig.pros,minGSSize = 3,TERM2GENE = path2Gene[,c(1,2)])@result
sig.phops.paths.k = enricher(sig.phops,minGSSize = 3,TERM2GENE = path2Gene[,c(1,2)])@result
sig.mrna.paths.k$cate = path2cate[as.character(sig.mrna.paths.k$Description),'cate']
sig.pro.paths.k$cate = path2cate[as.character(sig.pro.paths.k$Description),'cate']
sig.phops.paths.k$cate = path2cate[as.character(sig.phops.paths.k$Description),'cate']


####Figure 2a
pie(table(as.character(sig.mrna.paths.k$cate[sig.mrna.paths.k$cate != '6. Human Diseases' & sig.mrna.paths.k $pvalue < 0.05])),
    col = pal_npg()(6),main = 'Transcriptome')
pie(table(as.character(sig.pro.paths.k$cate[sig.pro.paths.k$cate != '6. Human Diseases'& sig.pro.paths.k $pvalue < 0.05])),
    col = pal_npg()(6),main = 'Proteome')
pie(table(as.character(sig.phops.paths.k$cate[sig.phops.paths.k$cate != '6. Human Diseases'& sig.phops.paths.k $pvalue < 0.05])),
    col = pal_npg()(6),main = 'Phosphorylated proteomics')

pie(table(as.character(sig.meta.paths.k$cate[sig.meta.paths.k$cate != '6. Human Diseases'& sig.meta.paths.k $pvalue < 0.05])),
    col = pal_npg()(6),main = 'Metabolism')

sig.mrna.paths.sm = enricher(sig.mrna,minGSSize = 3,TERM2GENE = path2gene.smpdb.c[,c('name.short','Gene.Name')])@result
sig.pro.paths.sm = enricher(sig.pros,minGSSize = 3,TERM2GENE = path2gene.smpdb.c[,c('name.short','Gene.Name')])@result
sig.phops.paths.sm = enricher(sig.phops,minGSSize = 3,TERM2GENE = path2gene.smpdb.c[,c('name.short','Gene.Name')])@result





sig.mrna.path.names.k = sig.mrna.paths.k[sig.mrna.paths.k$pvalue<0.05,'ID']
sig.pro.path.names.k = sig.pro.paths.k[sig.pro.paths.k$pvalue<0.05,'ID']
sig.phops.path.names.k = sig.phops.paths.k[sig.phops.paths.k$pvalue<0.05,'ID']
sig.meta.path.names.k = sig.meta.paths.k[round(sig.meta.paths.k$pvalue,2)<=0.05,'ID']

intersect(intersect(sig.mrna.path.names.k,sig.pro.path.names.k),sig.phops.path.names.k)

intersect(intersect(intersect(sig.mrna.path.names.k,sig.pro.path.names.k),sig.phops.path.names.k)
,sig.meta.path.names.k)

sig.mrna.paths.k$From = 'mRNA'
sig.pro.paths.k$From = 'protein'
sig.phops.paths.k$From = 'phosphProtein'
sig.meta.paths.k$From = 'metabolite'


sig.paths = rbind(sig.mrna.paths.k,sig.pro.paths.k,
                  sig.phops.paths.k,sig.meta.paths.k)
write.csv(sig.paths,file = './variables/sig molecule enriched pathways.csv')
sig.paths = read.csv(file = './variables/sig molecule enriched pathways.csv')
disease.pathways = as.character(path2cate[path2cate$cate == '6. Human Diseases','Name'])

pathways.input = unique(c(sig.mrna.path.names.k,
                          sig.meta.path.names.k,
                          sig.pro.path.names.k,
                          sig.phops.path.names.k))
path.res = rbind(sig.mrna.paths.k[intersect(rownames(sig.mrna.paths.k),pathways.input),],
                 sig.pro.paths.k[intersect(rownames(sig.pro.paths.k),pathways.input),],
                 sig.phops.paths.k[intersect(rownames(sig.phops.paths.k),pathways.input),],
                 sig.meta.paths.k[intersect(rownames(sig.meta.paths.k),pathways.input),])
dup.pathways = plyr::count(path.res[path.res$pvalue<0.05,]$ID)
dup.pathways = dup.pathways[dup.pathways$freq>=3,'x']

ggplot(path.res[path.res$ID %in% dup.pathways & path.res$ID %in% disease.pathways == F,],aes(x=ID,y=-log10(pvalue),fill=From,color=From))+
  geom_bar(stat = 'identity',width = 0.1)+
  geom_point(aes(x=ID,y=-log10(pvalue),color=From,size=Count))+
  geom_hline(yintercept = -log10(0.05),lty=3)+
  facet_grid(~From,scales = 'free_x')+coord_flip()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = pal_npg('nrc')(5))+
  scale_color_manual(values = pal_npg('nrc')(5))+
  scale_size_continuous(range = c(1,3))
#Figure 2d
ggplot(path.res[path.res$ID %in% as.character(sig.meta.paths.k$ID[1:10]),],aes(x=ID,y=-log10(pvalue),fill=From,color=From))+
  geom_bar(stat = 'identity',width = 0.1)+
  geom_point(aes(x=ID,y=-log10(pvalue),color=From,size=Count))+
  geom_hline(yintercept = -log10(0.05),lty=3)+
  facet_grid(~From,scales = 'free_x')+coord_flip()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = pal_npg('nrc')(5))+
  scale_color_manual(values = pal_npg('nrc')(5))+
  scale_size_continuous(range = c(1,3))

#Figure 2b
library(VennDiagram)
library(grid)
venn.p <- venn.diagram(x=list(Transcriptome = sig.mrna.path.names.k[sig.mrna.path.names.k %in% disease.pathways == F],
                              Proteome = sig.pro.path.names.k[sig.pro.path.names.k %in% disease.pathways == F],
                              PhosphProteome = sig.phops.path.names.k[sig.phops.path.names.k %in% disease.pathways == F]),
                       filename=NULL,
                       fill = pal_npg('nrc')(3))
grid.draw(venn.p)



#####Multi-omics visulization#####Figure 2e

library(pathview)
library(org.Hs.eg.db)

genes.all = union(union(as.character(res.mrna$gene),as.character(res.protein$gene)),as.character(res.phosph$gene))
genes.info = select(org.Hs.eg.db,keys = genes.all,keytype = 'SYMBOL',columns = c('SYMBOL','ENTREZID'))
genes.info = genes.info[!duplicated(genes.info$SYMBOL),]
genes.info$ENTREZID[is.na(genes.info$ENTREZID)] = genes.info$SYMBOL[is.na(genes.info$ENTREZID)]
rownames(genes.info)=genes.info$SYMBOL
res.mrna$id = genes.info[as.character(res.mrna$gene),'ENTREZID']
mrna.fc = res.mrna$Log2FC
names(mrna.fc)=res.mrna$id

res.protein$id = genes.info[as.character(res.protein$gene),'ENTREZID']
pro.fc = res.protein$log2FC
names(pro.fc)=res.protein$id

res.phosph$id = genes.info[as.character(res.phosph$gene),'ENTREZID']
phosphpro.fc = res.phosph$log2FC
names(phosphpro.fc)=res.phosph$id


res.meta = readRDS(file = './variables/all37Paired_diffMetaRes.Rds')
res.meta = res.meta[order(abs(res.meta$Log2FC),decreasing = T),]

res.meta$KEGG=meta.info[rownames(res.meta),'KEGG2']
res.meta.un = res.meta[!duplicated(res.meta$KEGG),]
rownames(res.meta.un)=res.meta.un$KEGG
fc.meta = res.meta.un$Log2FC
names(fc.meta)=rownames(res.meta.un)

rownames(res.mrna)=res.mrna$id
rownames(res.protein)=res.protein$id
rownames(res.phosph)=res.phosph$id

three.m = data.frame(mrna = res.mrna[unique(genes.info$ENTREZID),'Log2FC'],
                     protein = res.protein[unique(genes.info$ENTREZID),'log2FC'],
                     phosph =res.phosph[unique(genes.info$ENTREZID),'log2FC'] ,
                     row.names = unique(genes.info$ENTREZID))
three.m = as.matrix(three.m)
three.m[is.na(three.m)]=0

three.m.2 = data.frame(meta1 = fc.meta,
                       meta2 = fc.meta,
                       meta3 = fc.meta,
                       row.names = names(fc.meta))

int.path.res = path.res[path.res$Description == 'Glycolysis / Gluconeogenesis',]
int.path.items.g = unlist(lapply(int.path.res[1:3,]$geneID, function(a){
  strsplit(a,'/')[[1]]
}))
int.path.items.g = unique(int.path.items.g)
int.path.items.m = strsplit(int.path.res[4,'geneID'],'/')[[1]]
pv.out <- pathview(gene.data=three.m,cpd.data  = as.matrix(three.m.2), pathway.id =
                     "00010", species = "hsa", out.suffix = "Glycolysis",
                   same.layer = F)
