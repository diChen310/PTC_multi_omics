
library(survminer)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(MOVICS)
######################Clustering#################
#gene.res.aov, meta.res.aov from recurRiskAnalysis.R
genes.sig = gene.res.aov[order(gene.res.aov$p.value),'gene'][1:500]
metas.sig = meta.res.aov[meta.res.aov$p.value<0.05,'gene']

mo.data = list(omics1 = data.mrna[genes.sig,both.samples],
               omics2 = data.meta[metas.sig,both.samples])
set.seed(112334)
moic.res.try = getMOIC(data  = mo.data,
                        methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster",
                                           "ConsensusClustering", "CIMLR", "MoCluster",'iClusterBayes','IntNMF'),
                        N.clust     = 3,
                        type        = c("gaussian", "gaussian"))

cmoic <- getConsensusMOIC(moic.res.list = moic.res.try,
                          fig.name      = "CONSENSUS HEATMAP",
                          distance      = "euclidean")
ss.mat = cmoic$similarity.matrix
ss.dist = cmoic$clust.dend

hc.cs.4 = cutree(ss.dist,k=4)
cmoic$clust.res$clust = hc.cs.4
clin.both = clin.both[both.samples,]

# pheatmap::pheatmap(ss.mat[ss.dist$order,ss.dist$order],cluster_cols = F,cluster_rows = F,
#                    col=hcl.colors(100),
#                    show_rownames = F,
#                    show_colnames = F,
#                    annotation_col = data.frame(Cluster =paste0('CS', hc.cs.4[ss.dist$order]),
#                                                RecurRisk = clin.both[ss.dist$order,'RecurRisk'],
#                                                row.names = rownames(ss.mat)[ss.dist$order]),
#                    annotation_colors = list(RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red')))


#clin.both$Cluster = cmoic$clust.res[,'clust']#
clin.both$Cluster = hc.cs.4
#clin.both$Cluster = moic.res$ConsensusClustering$clust.res[,'clust']
survdiff(Surv(days_RFS, Recur) ~ Cluster, data = clin.both,rho = 0)
#clin.both = clin.both[!(clin.both$Cluster == '2' & clin.both$RecurRisk == 'Low'),]

fisher.test(table(clin.both$Cluster,clin.both$RecurRisk))   
clin.display$Cluster = clin.both[rownames(clin.display),'Cluster']



##############Cluster survival analysis################## Figure 5c


surv = compSurv(moic.res  = cmoic,
                surv.info        = clin.both,
                convt.time       = "m", # convert day unit to month
                surv.median.line = "h", # draw horizontal line at median survival
                fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

####Cluster markers identification####
runDEA(dea.method = "deseq2",
       expr       = mrna.count[,both.samples],
       moic.res   = cmoic,
       prefix     = "Thyroid")
data.meta.retain = retainCorrItem(data.meta[,both.samples],th=0.7)
runDEA(dea.method = "limma",
       expr       = data.meta[data.meta.retain,both.samples], 
       moic.res   = cmoic,
       prefix     = "Thyroid metabolism mo")


meta.info$KEGG[meta.info$KEGG == '-'] = meta.info$`Compound name`[meta.info$KEGG == '-']
meta.info$KEGG[meta.info$KEGG == ''] = meta.info$`Compound name`[meta.info$KEGG == '']
meta.info$KEGG2 = unlist(lapply(meta.info$KEGG,function(a){
  if(grepl(',',a)){return(strsplit(a,',')[[1]][1])}else{return(a)}
}))

data.meta.un.rows.ids = meta.info[rownames(data.meta.un),'KEGG2']
rownames(data.meta.un)[1:18]
data.meta.un.rows.ids[1:18] = c('C01194','C02862','C00195','Hex1Cer','ChE',
                                'C00165','C04233','LPC-O','C00157','C00350',
                                'C05973','C02737','C00550','C00319','C00422',
                                'C00344','C00162','OAHFA')
data.meta.un.rows.ids[duplicated(data.meta.un.rows.ids)]
rownames(data.meta.un)=data.meta.un.rows.ids

runDEA(dea.method = "limma",
       expr       = data.meta.un[,both.samples], 
       moic.res   = cmoic,
       prefix     = "Thyroid metabolism kegg ids2")


# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic,
                       dea.method    = "deseq2", # name of DEA method
                       prefix        = "Thyroid", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.01, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.01, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 10, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = data.mrna[,both.samples], # use normalized expression as heatmap input
                       
                       show_rownames = T, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

marker.down <- runMarker(moic.res      = cmoic,
                         dea.method    = "deseq2", # name of DEA method
                         prefix        = "Thyroid", # MUST be the same of argument in runDEA()
                         dat.path      = getwd(), # path of DEA files
                         res.path      = getwd(), # path to save marker files
                         p.cutoff      = 0.01, # p cutoff to identify significant DEGs
                         p.adj.cutoff  = 0.01, # padj cutoff to identify significant DEGs
                         dirct         = "down", # direction of dysregulation in expression
                         n.marker      = 10, # number of biomarkers for each subtype
                         doplot        = TRUE, # generate diagonal heatmap
                         norm.expr     = data.mrna[,both.samples], # use normalized expression as heatmap input
                         
                         show_rownames = T, # show no rownames (biomarker name)
                         fig.name      = "UPREGULATED BIOMARKER HEATMAP")
marker.up.meta <- runMarker(moic.res      = cmoic,
                            dea.method    = "limma", # name of DEA method
                            prefix        = "Thyroid metabolism mo", # MUST be the same of argument in runDEA()
                            dat.path      = getwd(), # path of DEA files
                            res.path      = getwd(), # path to save marker files
                            p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                            p.adj.cutoff  = 1.1, # padj cutoff to identify significant DEGs
                            dirct         = "up", # direction of dysregulation in expression
                            n.marker      = 10, # number of biomarkers for each subtype
                            doplot        = TRUE, # generate diagonal heatmap
                            norm.expr     = data.meta[data.meta.retain,both.samples], # use normalized expression as heatmap input
                            
                            show_rownames = T, # show no rownames (biomarker name)
                            fig.name      = "UPREGULATED BIOMARKER HEATMAP_metabolism2")
marker.d.meta <- runMarker(moic.res      = cmoic,
                           dea.method    = "limma", # name of DEA method
                           prefix        = "Thyroid metabolism mo", # MUST be the same of argument in runDEA()
                           dat.path      = getwd(), # path of DEA files
                           res.path      = getwd(), # path to save marker files
                           p.cutoff      = 0.1, # p cutoff to identify significant DEGs
                           p.adj.cutoff  = 1, # padj cutoff to identify significant DEGs
                           dirct         = "down", # direction of dysregulation in expression
                           n.marker      = 10, # number of biomarkers for each subtype
                           doplot        = TRUE, # generate diagonal heatmap
                           norm.expr     = data.meta[data.meta.retain,both.samples], # use normalized expression as heatmap input
                           
                           show_rownames = T, # show no rownames (biomarker name)
                           fig.name      = "UPREGULATED BIOMARKER HEATMAP_metabolism")
##############mrna markers for each cluster############ 
cosmic.hallmark.genes = read.delim(file='E:/data/COSMIC/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv')
ch.genes.list = as.character(unique(cosmic.hallmark.genes$GENE_NAME))
mrna.markers = read.delim(file='./consensusMOIC_deseq2_upregulated_marker_templates.txt')
res.cs1 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS1_vs_Others.txt',row.names = 1)
res.cs2 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS2_vs_Others.txt',row.names = 1)
res.cs3 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS3_vs_Others.txt',row.names = 1)
res.cs4 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS4_vs_Others.txt',row.names = 1)

res.cs1=res.cs1[as.character(mrna.markers$probe),]
res.cs1$gene = rownames(res.cs1)
res.cs2=res.cs2[as.character(mrna.markers$probe),]
res.cs2$gene = rownames(res.cs2)
res.cs3=res.cs3[as.character(mrna.markers$probe),]
res.cs3$gene = rownames(res.cs3)
res.cs4=res.cs4[as.character(mrna.markers$probe),]
res.cs4$gene = rownames(res.cs4)

res.cs1$Cluster = 'CS1'
res.cs2$Cluster = 'CS2'
res.cs3$Cluster = 'CS3'
res.cs4$Cluster = 'CS4'

res.cs.all = rbind(res.cs1,res.cs2,res.cs3,res.cs4)
res.cs.all$gene = factor(res.cs.all$gene,levels = as.character(mrna.markers$probe))
####Figure 6a
ggplot(res.cs.all,aes(x=gene,y=Cluster,color=log2fc,size = -log10(pvalue)))+geom_point()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_colour_gradientn(colours = hcl.colors(30,palette = 'Blue-Red 3'))+
  scale_size_continuous(range = c(1,5))

meta.markers = read.delim(file='./consensusMOIC_limma_upregulated_marker_templates.txt')
meta.markers2 = read.delim(file='./consensusMOIC_limma_downregulated_marker_templates.txt')
meta.markers = rbind(meta.markers,meta.markers2)
meta.markers = meta.markers[!duplicated(meta.markers$probe),]
res.cs1 = read.delim(file='./consensusMOIC_Thyroid metabolism mo_limma_test_result.CS1_vs_Others.txt',row.names = 1)
res.cs2 = read.delim(file='./consensusMOIC_Thyroid metabolism mo_limma_test_result.CS2_vs_Others.txt',row.names = 1)
res.cs3 = read.delim(file='./consensusMOIC_Thyroid metabolism mo_limma_test_result.CS3_vs_Others.txt',row.names = 1)
res.cs4 = read.delim(file='./consensusMOIC_Thyroid metabolism mo_limma_test_result.CS4_vs_Others.txt',row.names = 1)

res.cs1=res.cs1[as.character(meta.markers$probe),]
res.cs1$gene = rownames(res.cs1)
res.cs2=res.cs2[as.character(meta.markers$probe),]
res.cs2$gene = rownames(res.cs2)
res.cs3=res.cs3[as.character(meta.markers$probe),]
res.cs3$gene = rownames(res.cs3)
res.cs4=res.cs4[as.character(meta.markers$probe),]
res.cs4$gene = rownames(res.cs4)

res.cs1$Cluster = 'CS1'
res.cs2$Cluster = 'CS2'
res.cs3$Cluster = 'CS3'
res.cs4$Cluster = 'CS4'

res.cs.all = rbind(res.cs1,res.cs2,res.cs3,res.cs4)
res.cs.all$gene = factor(res.cs.all$gene,levels = as.character(meta.markers$probe))

####Figure 6b
ggplot(res.cs.all,aes(x=gene,y=Cluster,color=log2fc,size = -log10(pvalue)))+geom_point()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_colour_gradientn(colours = hcl.colors(30,palette = 'Blue-Red 3'))+
  scale_size_continuous(range = c(1,5))




# MUST locate ABSOLUTE path of msigdb file
#MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
MSIGDB.FILE = 'H:/project/multiomics/MSigDB_Hallmark.gmt'
# run GSEA to identify up-regulated GO pathways using results from edgeR
hallmark.up <- runGSEA(moic.res     = cmoic,
                       dea.method   = "deseq2", # name of DEA method
                       prefix       = "Thyroid", # MUST be the same of argument in runDEA()
                       dat.path     = getwd(), # path of DEA files
                       res.path     = getwd(), # path to save GSEA files
                       msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                       norm.expr    = data.mrna[,both.samples], # use normalized expression to calculate enrichment score
                       dirct        = "up", # direction of dysregulation in pathway
                       n.path = 10,
                       p.cutoff     = 0.05, # p cutoff to identify significant pathways
                       p.adj.cutoff = 0.5, # padj cutoff to identify significant pathways
                       gsva.method  = "gsva", # method to calculate single sample enrichment score
                       norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                       fig.name     = "UPREGULATED Msigdb_HallMark PATHWAY HEATMAP")

hallmark.d <- runGSEA(moic.res     = cmoic,
                      dea.method   = "deseq2", # name of DEA method
                      prefix       = "Thyroid", # MUST be the same of argument in runDEA()
                      dat.path     = getwd(), # path of DEA files
                      res.path     = getwd(), # path to save GSEA files
                      msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                      norm.expr    = data.mrna[,both.samples], # use normalized expression to calculate enrichment score
                      dirct        = "down", # direction of dysregulation in pathway
                      n.path = 15,
                      p.cutoff     = 0.05, # p cutoff to identify significant pathways
                      p.adj.cutoff = 0.5, # padj cutoff to identify significant pathways
                      gsva.method  = "gsva", # method to calculate single sample enrichment score
                      norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                      fig.name     = "DOWN-REGULATED Msigdb_HallMark PATHWAY HEATMAP")
res.cs2 = res.cs2[order(res.cs2$log2fc,decreasing = T),]
FC = res.cs2$log2fc
names(FC) = rownames(res.cs2)
gsea.hallmark.cs2 = GSEA(geneList = FC,TERM2GENE = path2Gene.Hall)
gseaplot(gsea.hallmark.cs2,'EPITHELIAL_MESENCHYMAL_TRANSITION')

metakegg.FILE <- 'H:/project/multiOmics/metabolism_kegg2symbol.gmt'
#setwd('H:/project/multiomics/temp/')
# run GSEA to identify up-regulated GO pathways using results from edgeR
meta.gsea.up <- runGSEA(moic.res     = cmoic,
                        dea.method   = "limma", # name of DEA method
                        prefix       = "Thyroid metabolism kegg ids2", # MUST be the same of argument in runDEA()
                        dat.path     = getwd(), # path of DEA files
                        res.path     = getwd(), # path to save GSEA files
                        msigdb.path  = metakegg.FILE, # MUST be the ABSOLUTE path of msigdb file
                        norm.expr    = data.meta.un[,both.samples], # use normalized expression to calculate enrichment score
                        dirct        = 'up', # direction of dysregulation in pathway
                        p.cutoff     = 0.05, # p cutoff to identify significant pathways
                        p.adj.cutoff = 0.9, # padj cutoff to identify significant pathways
                        gsva.method  = "gsva", # method to calculate single sample enrichment score
                        norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                        fig.name     = " Metabolism Pathway Heatmap")
meta.gsea.d <- runGSEA(moic.res     = cmoic,
                       dea.method   = "limma", # name of DEA method
                       prefix       = "Thyroid metabolism kegg ids2", # MUST be the same of argument in runDEA()
                       dat.path     = getwd(), # path of DEA files
                       res.path     = getwd(), # path to save GSEA files
                       msigdb.path  = metakegg.FILE, # MUST be the ABSOLUTE path of msigdb file
                       norm.expr    = data.meta.un[,both.samples], # use normalized expression to calculate enrichment score
                       dirct        = 'down', # direction of dysregulation in pathway
                       p.cutoff     = 0.05, # p cutoff to identify significant pathways
                       p.adj.cutoff = 0.9, # padj cutoff to identify significant pathways
                       gsva.method  = "gsva", # method to calculate single sample enrichment score
                       norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                       fig.name     = " Metabolism Pathway Heatmap")


print(gsea.up$gsea.list$CS1[1:6,3:6])
# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = moic.res$IntNMF,
                   dea.method   = "deseq2",
                   prefix       = "Thyroid",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = data.mrna[,both.samples],
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.5,
                   gsva.method  = "gsva", # switch to ssgsea
                   norm.method  = "mean", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 



##########getMoHeatmap with annRow_mrna/meta.marker#############
# data normalization for heatmap
clin.both$futime = clin.both$days_RFS*30
clin.both$fustat = clin.both$Recur
####Figure 5c
surv = compSurv(moic.res  = cmoic,
                surv.info        = clin.both,
                convt.time       = "m", # convert day unit to month
                surv.median.line = "h", # draw horizontal line at median survival
                fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
clin.both.high = clin.both[clin.both$Cluster != '1' &
                             clin.both$RecurRisk == 'High',]

####Figure 5d,e
km.by.cluster4 = survfit(Surv(days_RFS, Recur) ~ Cluster=='3', data = clin.both.high)
survdiff(Surv(days_RFS, Recur) ~ Cluster=='3', data = clin.both.high,rho = 0)

ggsurvplot(km.by.cluster4,
           pval = T, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() ,# Change ggplot2 theme
           palette = pal_npg('nrc')(10)
)

plotdata <- getStdiz(data       = mo.data,
                     halfwidth  = c(2,2), # no truncation for mutation
                     centerFlag = c(T,T), # no center for mutation
                     scaleFlag  = c(T,T)) # no scale for mutation
mRNA.col   <- brewer.pal(6,'RdBu')
meta.col   <- brewer.pal(6,'PiYG')
col.list   <- list(mRNA.col, meta.col)
CCDC6_RET.samples = clin.d[clin.d$CCDC6_RET==1,'Tumor_Sample_Barcode']
NCOA4_RET.samples= clin.d[clin.d$NCOA4_RET==1,'Tumor_Sample_Barcode']
NCOA4_RET.samples = gsub('_','-',gsub('S','',NCOA4_RET.samples))
CCDC6_RET.samples = gsub('_','-',gsub('S','',CCDC6_RET.samples))
MUC16.mut.samples = as.character(subsetMaf(var.annovar.maf,genes = c('MUC16'))@data$Tumor_Sample_Barcode)
MUC16.mut.samples = gsub('_','-',gsub('S','',MUC16.mut.samples))
BRAF.mut.samples = as.character(subsetMaf(var.annovar.maf,genes = c('BRAF'))@data$Tumor_Sample_Barcode)
BRAF.mut.samples = gsub('_','-',gsub('S','',BRAF.mut.samples))
BRAF.mut.samples=gsub('--','-',BRAF.mut.samples)

clin.display$CCDC6_RET = ifelse(rownames(clin.display) %in% CCDC6_RET.samples, 'TRUE', 'FALSE')
clin.display$NCOA4_RET = ifelse(rownames(clin.display) %in% NCOA4_RET.samples, 'TRUE', 'FALSE')
clin.display$BRAF = ifelse(rownames(clin.display) %in% BRAF.mut.samples, 'TRUE', 'FALSE')
clin.display$MUC16 = ifelse(rownames(clin.display) %in% MUC16.mut.samples, 'TRUE', 'FALSE')
####Figure 5a
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","Metabolism"),
             is.binary     = c(F,F), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.normalized","Metabolism.normalized"),
             clust.res     = cmoic$clust.res, # cluster results
             clust.dend    = cmoic$clust.dend, 
             annRow = list(mrna = sample(genes.sig,20),
                           meta = c(union(marker.up.meta$templates$probe,marker.d.meta$templates$probe))),
             show.rownames = c(F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             color         = col.list,
             annCol        = clin.display[,c('RecurRisk','Recur','TERT_promoter','BRAF','MUC16','CCDC6_RET','NCOA4_RET')], # no annotation for samples
             annColors     = list(RecurRisk = c(
               "Median" = 'cyan',
               "Low" = 'purple',
               "High" = 'red'
             ),
             Recur=c("TRUE"='black',"FALSE"='lightgrey'),
             TERT_promoter = c('TERT C228T' = 'black', "No" = 'lightgrey'),
             CCDC6_RET = c("TRUE" = 'black', "FALSE" = 'lightgrey'),
             NCOA4_RET = c("TRUE" = 'black', "FALSE" = 'lightgrey'),
             BRAF = c("TRUE" = 'black', "FALSE" = 'lightgrey'),
             MUC16 = c("TRUE" = 'black', "FALSE" = 'lightgrey')), # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP")

#####Compare TDS BRAF RAS scores###########
setwd('H:/project/multiOmics')
clin.display = read.csv(file = './variables/clin.display_TDSBSRS.csv',row.names = 1)
clin.display$Cluster = paste0('CS',cmoic$clust.res[rownames(clin.display),'clust'])
c.data = clin.display[,c('Age','TumorSize','LNM.No','Cluster','TDS','RAS.Score','BRAF.Score')]
library(reshape2)
c.data$Cluster2 = ifelse(c.data$Cluster == 'CS2','CS2','Others')

c.data = melt(c.data,measure.vars = c('Age','TumorSize','LNM.No','TDS','RAS.Score','BRAF.Score'))

ggplot(c.data,aes(x=Cluster,y=value,fill=Cluster))+
  geom_boxplot()+
  facet_wrap(~variable,scales = 'free')+
  theme_cowplot()+scale_fill_manual(values = pal_npg('nrc')(4))

ggplot(c.data,aes(x=Cluster,y=value,fill=Cluster))+geom_boxplot()+facet_grid(~variable,scales = 'free',space = 'free_y')+
  stat_compare_means(comparisons = list(c('CS1','CS2'),c('CS1','CS3'),c('CS2','CS3'),
                                        c('CS1','CS4'),c('CS4','CS3'),c('CS2','CS4')),
                     method = 't.test')+theme_bw()+scale_fill_manual(values = pal_npg('nrc')(4))

ggplot(c.data.plot,aes(x=Cluster2,y=value,fill=Cluster2))+
  geom_boxplot(scale='width')+facet_grid(~variable,scales = 'free')+
  stat_compare_means(method = 't.test')+theme_bw()+
  scale_fill_manual(values = pal_npg('nrc')(6)[c(1,6)])



####Compare clinical and mutation ####
check = function(a,b){
  cc = table(a,b)
  #总共有N件产品，其中M件次品，现在从中抽取n件做检查，抽到k件次品的概率分布服从超几何分布
  k1=cc[2,1]#k
  n1=sum(cc[2,])#n
  m1=sum(cc[,1])#M
  n=sum(cc)#N
  # q = the number of white balls drawn from the urn (without replacement)
  # q对应到抽样问题，为k
  # 
  # m = the number of white balls in the urn
  # m对应到抽样问题，为M
  # 
  # n = the number of black balls in the urn
  # n对应到抽样问题，为N-M
  # 
  # k = the number of balls drawn from the urn (sample size)
  # k对应到抽样问题，为n
  # ————————————————phyper(q-1, m, n, k, lower.tail=F)
  p1=phyper(k1, m1, n-m1, n1, lower.tail=F)
  k2=cc[2,2]#k
  n2=sum(cc[2,])#n
  m2=sum(cc[,2])#M
  
  k3=cc[2,3]#k
  n3=sum(cc[2,])#n
  m3=sum(cc[,3])#M
  p3=phyper(k3, m3, n-m3, n3, lower.tail=F)
  
  
  p2=phyper(k2, m2, n-m2, n2, lower.tail=F)
  
  k4=cc[2,4]#k
  n4=sum(cc[2,])#n
  m4=sum(cc[,4])#M
  p4=phyper(k4, m4, n-m4, n4, lower.tail=F)
  
  return(c(p1,p2,p3,p4))
}

check2 = function(a,b){
  cc = table(a,b)
  #总共有N件产品，其中M件次品，现在从中抽取n件做检查，抽到k件次品的概率分布服从超几何分布
  k1=cc[2,1]#k
  n1=sum(cc[2,])#n
  m1=sum(cc[,1])#M
  n=sum(cc)#N
  # q = the number of white balls drawn from the urn (without replacement)
  # q对应到抽样问题，为k
  # 
  # m = the number of white balls in the urn
  # m对应到抽样问题，为M
  # 
  # n = the number of black balls in the urn
  # n对应到抽样问题，为N-M
  # 
  # k = the number of balls drawn from the urn (sample size)
  # k对应到抽样问题，为n
  # ————————————————phyper(q-1, m, n, k, lower.tail=F)
  p1=dhyper(k1, m1, n-m1, n1)
  k2=cc[2,2]#k
  n2=sum(cc[2,])#n
  m2=sum(cc[,2])#M
  
  k3=cc[2,3]#k
  n3=sum(cc[2,])#n
  m3=sum(cc[,3])#M
  p3=dhyper(k3, m3, n-m3, n3)
  
  
  p2=dhyper(k2, m2, n-m2, n2)
  
  k4=cc[2,4]#k
  n4=sum(cc[2,])#n
  m4=sum(cc[,4])#M
  p4=dhyper(k4, m4, n-m4, n4)
  
  return(c(p1,p2,p3,p4))
}

for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter','CCDC6_RET','NCOA4_RET',
                 'BRAF','MUC16')){
  print(c(feature,chisq.test(clin.display$Cluster,clin.display[,feature])$p.value))
  
}

check(clin.display$TERT_promoter,clin.display$Cluster)

#[1] 0.352678124 0.714572810 0.842829199 0.003487506
check(clin.display$BRAF,clin.display$Cluster)
#0.1233816302 0.7131736013 0.9915548449 0.0004093522
check(clin.display$Recur,clin.display$Cluster)
check(clin.display$CCDC6_RET,clin.display$Cluster)

check(clin.display$NCOA4_RET,clin.display$Cluster)

check(clin.display$MUC16,clin.display$Cluster)
chisq.test(table(clin.display$TERT_promoter,clin.display$Cluster=='CS4'))
chisq.test(table(clin.display$BRAF,clin.display$Cluster=='CS4'))

check2(clin.display$BRAF,clin.display$Cluster)


check2(clin.display$TERT_promoter,clin.display$Cluster)

for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter','CCDC6_RET','NCOA4_RET',
                 'BRAF','MUC16')){
  print(c(feature,chisq.test(clin.display$Cluster %in% c('CS2'),clin.display[,feature])$p.value))
  
}

clin.display$RecurRisk = factor(clin.display$RecurRisk,levels = c('High','Median','Low'))
plot.data = melt(clin.display,measure.vars = c('RecurRisk','exNInvasion','N.stage','T.stage','BRAF','TERT_promoter','MUC16'))
plot.data$value = factor(plot.data$value,levels = c('High','Median','Low','TRUE','FALSE',
                                                    "N0","N1a","N1b","T1" ,"T2","T3" ,"T4",'TERT C228T','No'))
plot.data$variable = factor(plot.data$variable,levels = c('RecurRisk','T.stage','N.stage','BRAF','exNInvasion','TERT_promoter','MUC16'))

####Figure 5b
ggplot(plot.data,aes(x=Cluster,fill=value))+geom_bar(width = 0.75,position = 'fill')+
  facet_wrap(~variable,scales = 'free',nrow = 2)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c('red','cyan','purple','black','grey',pal_npg()(10)[c(-1,-2)],'lightgrey'))
#"Median"='cyan',"Low"='purple',"High"='red'



####Pathway check####'
marker.up <- runMarker(moic.res      = cmoic,
                       dea.method    = "deseq2", # name of DEA method
                       prefix        = "Thyroid", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.01, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 50, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = data.mrna[,both.samples], # use normalized expression as heatmap input
                       
                       show_rownames = T, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
####Diff mRNA####
res.cs1 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS1_vs_Others.txt',row.names = 1)
res.cs2 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS2_vs_Others.txt',row.names = 1)
res.cs3 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS3_vs_Others.txt',row.names = 1)
res.cs4 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS4_vs_Others.txt',row.names = 1)
mrna.markers = read.delim(file='./consensusMOIC_deseq2_upregulated_marker_templates.txt')
cs2.path = rbind(read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS2_unique_downexpr_pathway.txt',row.names = 1),
                 read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS2_unique_upexpr_pathway.txt',row.names = 1))


res.cs2.up = mrna.markers[mrna.markers$class == 'CS2','probe']

cs2.path.up = enricher(res.cs2.up,TERM2GENE = path2Gene[,c(1,2)])@result


res.cs2.d = rownames(res.cs2[res.cs2$log2fc<0 & res.cs2$padj<0.01,])

cs2.path.d = enricher(res.cs2.d,TERM2GENE = path2Gene[,c(1,2)])@result


res.cs3.up = rownames(res.cs3[res.cs3$log2fc>2 & res.cs3$padj<0.01,])

cs3.path.up = enricher(res.cs3.up,TERM2GENE = path2Gene[,c(1,2)])@result


res.cs2.d = rownames(res.cs2[res.cs2$log2fc<0 & res.cs2$padj<0.01,])

cs2.path.d = enricher(res.cs2.d,TERM2GENE = path2Gene[,c(1,2)])@result


####GSEA kegg pathway####
res.cs1 = res.cs1[order(res.cs1$log2fc,decreasing = T),]
res.cs1.fc = res.cs1$log2fc
names(res.cs1.fc)= rownames(res.cs1)
cs1.path.gsea = GSEA(res.cs1.fc,
                     minGSSize = 2,
                     TERM2GENE = path2Gene[,c(1,2)],
                     pvalueCutoff = 1.5)
cs1.path.gsea.res = cs1.path.gsea@result


res.cs2 = res.cs2[order(res.cs2$log2fc,decreasing = T),]
res.cs2.fc = res.cs2$log2fc
names(res.cs2.fc)= rownames(res.cs2)
cs2.path.gsea = GSEA(res.cs2.fc,
                     minGSSize = 2,
                     TERM2GENE = path2Gene[,c(1,2)],
                     pvalueCutoff = 1.5)
cs2.path.gsea.res = cs2.path.gsea@result

res.cs3 = res.cs3[order(res.cs3$log2fc,decreasing = T),]
res.cs3.fc = res.cs3$log2fc
names(res.cs3.fc)= rownames(res.cs3)
cs3.path.gsea = GSEA(res.cs3.fc,
                     minGSSize = 2,
                     TERM2GENE = path2Gene[,c(1,2)],
                     pvalueCutoff = 1.5)
cs3.path.gsea.res = cs3.path.gsea@result

res.cs4 = res.cs4[order(res.cs4$log2fc,decreasing = T),]
res.cs4.fc = res.cs4$log2fc
names(res.cs4.fc)= rownames(res.cs4)
cs4.path.gsea = GSEA(res.cs4.fc,
                     minGSSize = 2,
                     TERM2GENE = path2Gene[,c(1,2)],
                     pvalueCutoff = 1.5)
cs4.path.gsea.res = cs4.path.gsea@result



cs1.path.gsea.res$From = 'CS1'
cs2.path.gsea.res$From = 'CS2'
cs3.path.gsea.res$From = 'CS3'
cs4.path.gsea.res$From = 'CS4'


kegg.gsea.css = rbind(cs1.path.gsea.res,
                      cs2.path.gsea.res,
                      cs3.path.gsea.res,
                      cs4.path.gsea.res)
saveRDS(kegg.gsea.css,file = './variables/cs1tocs4 kegg gsea.Rds')
kegg.gsea.css = readRDS(file = './variables/cs1tocs4 kegg gsea.Rds')

gsea.paths = as.character(kegg.gsea.css[kegg.gsea.css$pvalue<0.05,'ID'])
kegg.gsea.css = kegg.gsea.css[kegg.gsea.css$ID %in% gsea.paths,]
kegg.gsea.css$cate = path2cate[kegg.gsea.css$ID,'cate']
kegg.gsea.css$subcate = path2cate[kegg.gsea.css$ID,'subcate']

kegg.gsea.css = kegg.gsea.css[kegg.gsea.css$cate != '6. Human Diseases',]

ggplot(kegg.gsea.css,
       aes(y=ID,x=From))+
  geom_raster(aes(fill=NES))+
  facet_grid(rows = vars(cate),scales = 'free',space = 'free')+
  geom_text(aes(label= ifelse(pvalue<0.05,'*','')))+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(100,'Blue-Red 3'))


####Figure 6h
ggplot(kegg.gsea.css[kegg.gsea.css$cate == '1. Metabolism',],
       aes(y=ID,x=From))+
  geom_raster(aes(fill=NES))+
  geom_text(aes(label= ifelse(pvalue<0.05,'*','')))+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(100,'Blue-Red 3'))

####Figure 6i
ggplot(kegg.gsea.css[kegg.gsea.css$cate == '5. Organismal Systems',],
       aes(y=ID,x=From))+
  facet_grid(rows = vars(subcate),scales = 'free',space = 'free')+
  geom_raster(aes(fill=NES))+
  geom_text(aes(label= ifelse(pvalue<0.05,'*','')))+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(100,'Blue-Red 3'))

ggplot(kegg.gsea.css[kegg.gsea.css$cate %in% c('1. Metabolism','5. Organismal Systems') == F ,],
       aes(y=ID,x=From))+
  geom_raster(aes(fill=NES))+
  facet_grid(rows = vars(cate),scales = 'free',space = 'free')+
  geom_text(aes(label= ifelse(pvalue<0.05,'*','')))+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_gradientn(colours = hcl.colors(100,'Blue-Red 3'))



