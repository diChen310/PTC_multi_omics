####Find the biomarkers for prediction of the models
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggsci)
library(survival)
library(randomForest)
itemImportance<-function(omics,cols,both.res,clu.s){
  
  X=t(omics[,cols])
  
  set.seed(71)
  rf <- randomForest(as.matrix(X),as.factor(both.res[,clu.s]), importance=TRUE,
                     proximity=TRUE)
  varImpPlot(rf,n.var=50)
  ## Look at variable importance:
  imp.s<-round(importance(rf), 2)
  
  #write.csv(round(importance(rf), 2),file='E:/project/omics integration/proteomics/gene importance score.csv')
  return(imp.s)
}

imp.genes <-
  itemImportance(t(scale(t(mo.data$omics1))), rownames(clin.both), clin.both, 'Cluster')
imp.pp.genes <-
  row.names(imp.genes[rev(order(as.double(imp.genes[, 'MeanDecreaseGini']))), ])[1:30]
imp.metas <-
  itemImportance(t(scale(t(mo.data$omics2))), rownames(clin.both), clin.both, 'Cluster')
imp.pp.metas <-
  c( row.names(imp.metas[rev(order(as.double(imp.metas[, 'MeanDecreaseGini']))), ])[1:30]
  )

######machine learning methods: problem: dataset-imbalance

library(e1071)
library(caret)
omics.merge = rbind(t(scale(t(mo.data$omics1[,both.samples]))),
                    t(scale(t(mo.data$omics2[,both.samples]))))
sel.mols = c(imp.pp.genes,imp.pp.metas)
input.data = t( omics.merge[sel.mols,])
input.data=as.data.frame(input.data)
input.data$Cluster = as.factor(paste0('CS',clin.both[rownames(input.data),'Cluster']))
trainIndex <- createDataPartition(input.data$Cluster, p = .6, 
                                  list = FALSE, 
                                  times = 1)

train= input.data[trainIndex,]
test= input.data[-trainIndex,]
control = trainControl(method = 'cv',number = 10,classProbs = TRUE)
#t

rf.train <- train(x=train[,colnames(train) != 'Cluster'],y=train$Cluster, method = 'rf',
                   trControl = control)

pred.p = predict(rf.train,test,type = 'prob')
pred.p.t = predict(rf.train,train,type = 'prob')

pred.y = predict(rf.train,test,type = 'raw')

# Check accuracy:
table(pred.y, test$Cluster)

roc.res.cs1 = roc(response = ifelse(test$Cluster == 'CS1','CS1','Others'),
                    predictor = pred.p[,1],
                    levels = c('CS1','Others'))
roc.res.cs2 = roc(response = ifelse(test$Cluster == 'CS2','CS2','Others'),
                  predictor = pred.p[,2],
                  levels = c('CS2','Others'))
roc.res.cs3 = roc(response = ifelse(test$Cluster == 'CS3','CS3','Others'),
                  predictor = pred.p[,3],
                  levels = c('CS3','Others'))
roc.res.cs4 = roc(response = ifelse(test$Cluster == 'CS4','CS4','Others'),
                  predictor = pred.p[,4],
                  levels = c('CS4','Others'))


####Figure 7c
plot(roc.res.cs1,type = "S",col = pal_npg()(4)[1])
plot(roc.res.cs2,add=T,col = pal_npg()(4)[2])
plot(roc.res.cs3,add=T,col = pal_npg()(4)[3])
plot(roc.res.cs4,add=T,col = pal_npg()(4)[4])
text(0.4,0.2,labels = paste('AUC(CS1):',signif(roc.res.cs1$auc,2)),col=pal_npg()(4)[1])
text(0.4,0.3,labels = paste('AUC(CS2):',signif(roc.res.cs2$auc,2)),col=pal_npg()(4)[2])
text(0.4,0.4,labels = paste('AUC(CS3):',signif(roc.res.cs3$auc,2)),col=pal_npg()(4)[3])
text(0.4,0.5,labels = paste('AUC(CS4):',signif(roc.res.cs4$auc,2)),col=pal_npg()(4)[4])

saveRDS(rf.train,file = './variables/rf_train_twoOmics_top30.Rds')
####Check on the probability for prognosis####
pred.all = predict(rf.train,input.data[,sel.mols],type = 'prob')
saveRDS(pred.all,file = './variables/subType_prob-all.Rds')

clin.both$prob = pred.all[,3]
coxph(Surv(days_RFS, Recur) ~ prob, data = clin.both)

####Figure 7d
clin.both$probCS2_l = ifelse(pred.all[,2]>0.25,'L','S')
km.by.cluster = survfit(Surv(days_RFS, Recur) ~ probCS2_l, data = clin.both)
survdiff(Surv(days_RFS, Recur) ~ probCS2_l, data = clin.both,rho = 0)

library(survminer)
library(RColorBrewer)
ggsurvplot(km.by.cluster,
           pval = T, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() ,# Change ggplot2 theme
           palette = pal_npg('nrc')(10)
)

#####Find biomarkers or targets based on databases####

res.cs1 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS1_vs_Others.txt',row.names = 1)
res.cs2 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS2_vs_Others.txt',row.names = 1)
res.cs3 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS3_vs_Others.txt',row.names = 1)
res.cs4 = read.delim(file='./consensusMOIC_Thyroid_deseq2_test_result.CS4_vs_Others.txt',row.names = 1)
mrna.markers = read.delim(file='./consensusMOIC_deseq2_upregulated_marker_templates.txt')

res.cs1$gene = rownames(res.cs1)
res.cs2$gene = rownames(res.cs2)
res.cs3$gene = rownames(res.cs3)
res.cs4$gene = rownames(res.cs4)

res.cs1$Cluster = 'C1'
res.cs2$Cluster = 'C2'
res.cs3$Cluster = 'C3'
res.cs4$Cluster = 'C4'

res.cs.all = rbind(res.cs1,res.cs2,res.cs3,res.cs4)

res.cs.all.sig = res.cs.all[res.cs.all$padj<0.001 & res.cs.all$log2fc>2,]
sig.genes = unique(res.cs.all.sig$gene)
length(sig.genes)
library(maftools)
diffCluster.dgi = drugInteractions(genes = sig.genes, drugs = TRUE)
dg.cs1 = drugInteractions(genes = mrna.markers[mrna.markers$class == 'CS1','probe'])
dg.cs2 = drugInteractions(genes = mrna.markers[mrna.markers$class == 'CS2','probe'])
dg.cs3 = drugInteractions(genes = mrna.markers[mrna.markers$class == 'CS3','probe'])
dg.cs4 = drugInteractions(genes = mrna.markers[mrna.markers$class == 'CS4','probe'])

dg.cs1.genes = unique(dg.cs1$Gene)
dg.cs2.genes = unique(dg.cs2$Gene)
dg.cs3.genes = unique(dg.cs3$Gene)
dg.cs4.genes = unique(dg.cs4$Gene)

diffCluster.dgi.cs1 = drugInteractions(genes = dg.cs1.genes, drugs = TRUE)
diffCluster.dgi.cs2 = drugInteractions(genes = dg.cs2.genes, drugs = TRUE)
diffCluster.dgi.cs3 = drugInteractions(genes = dg.cs3.genes, drugs = TRUE)
diffCluster.dgi.cs4 = drugInteractions(genes = dg.cs4.genes, drugs = TRUE)

diffCluster.dgi.cs1 = diffCluster.dgi.cs1[diffCluster.dgi.cs1$interaction_types !='agonist',]
diffCluster.dgi.cs2 = diffCluster.dgi.cs2[diffCluster.dgi.cs2$interaction_types !='agonist',]
diffCluster.dgi.cs3 = diffCluster.dgi.cs3[diffCluster.dgi.cs3$interaction_types !='agonist',]
diffCluster.dgi.cs4 = diffCluster.dgi.cs4[diffCluster.dgi.cs4$interaction_types !='agonist',]

dg.cs1.genes = unique(diffCluster.dgi.cs1$Gene)
dg.cs2.genes = unique(diffCluster.dgi.cs2$Gene)
dg.cs3.genes = unique(diffCluster.dgi.cs3$Gene)
dg.cs4.genes = unique(diffCluster.dgi.cs4$Gene)

dg.genes = data.frame(Gene = c(dg.cs1.genes,dg.cs2.genes,dg.cs3.genes,dg.cs4.genes),
                      From = c(rep('CS1',length(dg.cs1.genes)),
                               rep('CS2',length(dg.cs2.genes)),
                               rep('CS3',length(dg.cs3.genes)),
                               rep('CS4',length(dg.cs4.genes))))
rownames(dg.genes)=dg.genes$Gene
res.mrna = readRDS(file = './variables/res.mrna.paired.Rds')
rownames(res.mrna)=res.mrna$id
# mat.dg = t(scale(t(data.mrna[sig.genes,rownames(clin.both)])))
# plsda.res = plsda(t(mat.dg),as.factor(clin.both$Cluster),ncomp = 2)
# plotIndiv(plsda.res,ind.names = F,ellipse=F,col.per.group =pal_npg()(4),style = 'graphics')
# plotLoadings(plsda.res, comp = 1, contrib = "max", ndisplay = 50)
# 

####Figure 7a
pheatmap::pheatmap(t(scale(t(data.mrna[dg.genes$Gene,both.samples[ss.dist$order]]))),cluster_cols = F,cluster_rows = F,
                   col=hcl.colors(100,'Blue-Red 3'),border_color = NA,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_col = data.frame(Cluster =paste0('CS', hc.cs.4[ss.dist$order]),
                                               row.names = rownames(ss.mat)[ss.dist$order]),
                   annotation_row = data.frame(From = dg.genes$From,row.names = rownames(dg.genes),
                                               CS1=res.cs1[dg.genes$Gene,'log2fc'],
                                               CS2=res.cs2[dg.genes$Gene,'log2fc'],
                                               CS3=res.cs3[dg.genes$Gene,'log2fc'],
                                               CS4=res.cs4[dg.genes$Gene,'log2fc']),
                   annotation_colors = list(Cluster=c('CS1'=pal_npg()(4)[1],
                                                      'CS2'=pal_npg()(4)[2],
                                                      'CS3'=pal_npg()(4)[3],
                                                      'CS4'=pal_npg()(4)[4]),
                                            From = c('CS1'=pal_npg()(4)[1],
                                                     'CS2'=pal_npg()(4)[2],
                                                     'CS3'=pal_npg()(4)[3],
                                                     'CS4'=pal_npg()(4)[4]),
                                            CS1=hcl.colors(50,'Blue-Red 2'),
                                            CS2=hcl.colors(50,'Blue-Red 2'),
                                            CS3=hcl.colors(50,'Blue-Red 2'),
                                            CS4=hcl.colors(50,'Blue-Red 2')),
                   gaps_row = c(6,14,32),
                   gaps_col = c(26,41,58),
                   fontsize_row = 8
                   )
res.cs.all.sig = res.cs.all[res.cs.all$gene %in% dg.genes$Gene,]
res.cs.all.sig$gene = factor(res.cs.all.sig$gene,levels = rev(dg.genes$Gene))

ggplot(res.cs.all.sig,aes(y=gene,x=Cluster,color=log2fc,size = -log10(pvalue)))+geom_point()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_colour_gradientn(colours = hcl.colors(30,palette = 'Blue-Red 3'))+
  scale_size_continuous(range = c(1,5))
             

####Correlations between mrna and metabolism#### Figure 7b
sig.genes = as.character(dg.genes$Gene)
meta.mrna.corrs = cor(t(data.mrna[sig.genes,both.samples]),
                      t(data.meta[,both.samples]),
                      method = 'pearson')
library(reshape2)
meta.mrna.corrs = melt(meta.mrna.corrs,measure.vars = colnames(meta.mrna.corrs))
fivenum(meta.mrna.corrs$value)
meta.mrna.corrs.high = meta.mrna.corrs[abs(meta.mrna.corrs$value)>quantile(abs(meta.mrna.corrs$value),1-50/nrow(meta.mrna.corrs)),]
library(igraph)
library(ggraph)
library(tidygraph)
#library(ggnewscale)
g.net = graph_from_data_frame(meta.mrna.corrs.high[,c(1,2)],directed = F)
sub.g = as_tbl_graph(g.net,directed = F)
E(sub.g)$weight = meta.mrna.corrs.high$value
subtype = function(a,input){
  
  if(a %in% sig.genes){
    
    type.m = 'mRNA'
  }else{
    type.m = input
  }
}
sub.g=
  activate(sub.g,nodes) %>% 
  mutate(From = ifelse( names(V(sub.g)) %in% sig.genes,'aGene','Metabolite')) %>%
  mutate(Type = unlist(lapply(names(V(sub.g)), function(i){
    subtype(i,'Metabolite')
  }))) %>% 
  activate(edges) %>% 
  mutate(edge_weights = ifelse(meta.mrna.corrs.high$value<0,'-','+'))

  ggraph(sub.g,layout = "graphopt") + 
  geom_edge_link(aes(edge_colour = edge_weights),width=0.5,edge_alpha=0.6)+
  scale_edge_colour_manual(values = c('blue','red'))+
  geom_node_point(aes(colour = Type,shape=From),size=4) +
  geom_node_text(aes(label = name),size=2, repel = TRUE)+
  scale_color_manual(values =pal_npg('nrc')(4))+
  theme_graph()

  
meta.i = c('N8-Acetylspermidine','Glutamic acid','Taurine')

data.plot = data.meta[meta.i,c(cs1.samples,cs2.samples,cs3.samples)]
data.plot = melt(data.plot,measure.vars = colnames(data.plot))
data.plot$cluster = paste0('CS',clin.display[data.plot$Var2,'Cluster'])

ggplot(data.plot,aes(x=cluster,y=value,fill=cluster))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('CS1','CS2'),
                                        c('CS1','CS3'),
                                        c('CS2','CS3')),method = 't.test')+
  facet_wrap(~Var1)+theme_bw()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = pal_npg()(3))
  

