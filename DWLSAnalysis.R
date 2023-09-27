library(DWLS)#############################################Codes from Ubuntu###############################
scdata = readRDS(file='./scdata_matrix.Rds')
labels = readRDS(file='./scdata_labels.Rds')
labels = as.character(labels)
labels[labels == 'T/NK']='TorNK'
unique(labels)
dataBulk = readRDS(file='./deseqNorm_count.Rds')
scdata = as.matrix(scdata)
sig <- buildSignatureMatrixMAST(scdata,id=as.character(labels),path='./DWLS_result')

cellComp.matrix = matrix(,nrow = length(unique(labels)),ncol=ncol(dataBulk),
                         dimnames = list(unique(labels),colnames(dataBulk)))
for(i in 1:ncol(dataBulk)){
  i.dataBulk <- dataBulk[,i]
  trimmed <- trimData(sig, i.dataBulk)
  S <- trimmed$sig
  B <- trimmed$bulk
  res <- solveDampenedWLS(S, B)
  cellComp.matrix[,i]=res
}
heatmap(cellComp.matrix)
###################################End of codes in Ubuntu########################
setwd('H:/project/multiOmics')
cellComp.matrix = readRDS(file='./data/rnaseq/cellCompMatrix.Rds')

cols = rownames(ss.mat)[ss.dist$order]

pheatmap::pheatmap( t(scale(t(cellComp.matrix[, cols]))),cluster_cols = F,cluster_rows = T,
                   border_color = NA,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_col = data.frame(Cluster =paste0('CS', hc.cs.4[ss.dist$order]),
                                               RecurRisk = clin.both[ss.dist$order,'RecurRisk'],
                                               row.names = rownames(ss.mat)[ss.dist$order]),
                   annotation_colors = list(RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red')))

c.data.plot = cellComp.matrix[, cols]

c.data.plot = data.frame(t(c.data.plot))
c.data.plot$Immune = apply(c.data.plot[,c('B','TorNK','Myeloid')],1,sum)
c.data.plot$Cluster = paste('C',clin.both[rownames(c.data.plot),'Cluster'])
c.data.plot$Cluster2 = ifelse(c.data.plot$Cluster == 'C 1','C 1','Others')

c.data.plot = melt(c.data.plot,measure.vars = c(rownames(cellComp.matrix),'Immune'))

ggplot(c.data.plot,aes(x=Cluster2,y=value,fill=Cluster2))+
  geom_boxplot(scale='width')+facet_grid(~variable,scales = 'free')+
  stat_compare_means(method = 't.test')+theme_bw()+
  scale_fill_manual(values = pal_npg('nrc')(6)[c(1,6)])

ggplot(c.data.plot,aes(x=Cluster,y=value,fill=Cluster))+
  geom_boxplot(scale='width')+facet_grid(~variable,scales = 'free')+
  stat_compare_means(method = 'anova')+theme_bw()+
  scale_fill_manual(values = pal_npg('nrc')(4))

ggplot(c.data.plot,aes(x=Cluster,y=value,fill=Cluster))+geom_violin(scale='width')+facet_grid(~variable,scales = 'free')+
  stat_compare_means(method = 't.test',comparisons = list(c('C 1','C 2'),
                                                          c('C 2','C 3'),
                                                          c('C 1','C 3'),
                                                          c('C 1','C 4'),
                                                          c('C 2','C 4'),
                                                          c('C 4','C 3')))+
  theme_bw()+
  scale_fill_manual(values = pal_npg('nrc')(4))



cellComp.matrix = readRDS(file='./data/rnaseq/cellCompMatrix2.Rds')
pheatmap::pheatmap( t(scale(t(cellComp.matrix[, cols]))),cluster_cols = F,cluster_rows = T,
                    border_color = NA,
                    show_rownames = T,
                    show_colnames = F,
                    annotation_col = data.frame(Cluster =paste0('CS', hc.cs.4[ss.dist$order]),
                                                RecurRisk = clin.both[ss.dist$order,'RecurRisk'],
                                                row.names = rownames(ss.mat)[ss.dist$order]),
                    annotation_colors = list(RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red')))

cell.sig.pvalues = unlist(lapply(1:nrow(cellComp.matrix), function(i){
  data.i = data.frame(cell=cellComp.matrix[i,cols],
                      Cluster=paste0('CS', hc.cs.4[ss.dist$order]))
  res = aov(cell ~ Cluster,data=data.i)
  return(summary(res)[[1]][1,5])
}))

cell.sig.pvalues.cs1 = unlist(lapply(1:nrow(cellComp.matrix), function(i){
  data.i = data.frame(cell=cellComp.matrix[i,cols],
                      Cluster=ifelse(hc.cs.4[ss.dist$order]=='1','CS1','Others'))
  res = kruskal.test(cell ~ Cluster,data=data.i)
  return(res$p.value)
}))
cell.sig.pvalues.cs2 = unlist(lapply(1:nrow(cellComp.matrix), function(i){
  data.i = data.frame(cell=cellComp.matrix[i,cols],
                      Cluster=ifelse(hc.cs.4[ss.dist$order]=='2','CS2','Others'))
  res = kruskal.test(cell ~ Cluster,data=data.i)
  return(res$p.value)
}))

cell.sig.pvalues.cs3 = unlist(lapply(1:nrow(cellComp.matrix), function(i){
  data.i = data.frame(cell=cellComp.matrix[i,cols],
                      Cluster=ifelse(hc.cs.4[ss.dist$order]=='3','CS3','Others'))
  res = kruskal.test(cell ~ Cluster,data=data.i)
  return(res$p.value)
}))

cell.sig.pvalues.cs4 = unlist(lapply(1:nrow(cellComp.matrix), function(i){
  data.i = data.frame(cell=cellComp.matrix[i,cols],
                      Cluster=ifelse(hc.cs.4[ss.dist$order]=='4','CS4','Others'))
  res = kruskal.test(cell ~ Cluster,data=data.i)
  return(res$p.value)
}))
cell.sig.res = data.frame(cell = rownames(cellComp.matrix),p.value = cell.sig.pvalues,
                          p1 = cell.sig.pvalues.cs1,
                          p2 = cell.sig.pvalues.cs2,
                          p3 = cell.sig.pvalues.cs3,
                          p4= cell.sig.pvalues.cs4)
sig.cells = cell.sig.res[cell.sig.res$p.value<0.05 |
                           cell.sig.res$p1<0.05 |
                           cell.sig.res$p2<0.05 |
                           cell.sig.res$p3<0.05 |
                           cell.sig.res$p4<0.05,'cell']

c.data.plot = cellComp.matrix[sig.cells, cols]

c.data.plot = data.frame(t(c.data.plot))
c.data.plot$Cluster = paste0('CS', hc.cs.4[ss.dist$order])
c.data.plot$Cluster2 = ifelse(c.data.plot$Cluster == 'CS1','CS1','Others')

c.data.plot = melt(c.data.plot,measure.vars = sig.cells)

c.data.plot$variable = factor(c.data.plot$variable,levels = c('c03','c04','c05','c06','c07','c08',
                                                              'CD4_c6_FOXP3','arteries','veins',
                                                              'tip_cell'))
ggplot(c.data.plot,aes(x=Cluster,y=value,fill=Cluster))+
  geom_boxplot(scale='width')+facet_wrap(~variable,scales = 'free',nrow = 2)+
  theme_bw()+theme(strip.background = element_blank(),axis.text = element_text(size=8) )+
  scale_fill_manual(values = pal_npg('nrc')(4))
ggplot(c.data.plot,aes(x=Cluster,y=value,fill=Cluster))+
  geom_boxplot(scale='width')+
  facet_wrap(~variable,scales = 'free',nrow = 2)+
  stat_compare_means(method = 't.test',comparisons = list(c('CS1','CS2'),
                                                          c('CS2','CS3'),
                                                          c('CS1','CS3'),
                                                          c('CS1','CS4'),
                                                          c('CS2','CS4'),
                                                          c('CS4','CS3')),
                     label = 'p.signif',
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                          symbols = c("****", "***", "**", "*", "")))+
  theme_bw()+
  scale_fill_manual(values = pal_npg('nrc')(4))
pheatmap( t(scale(t(cellComp.matrix[sig.cells, cols]))),cluster_cols = F,cluster_rows = T,
          color = hcl.colors(100),
                    border_color = NA,
                    show_rownames = T,
                    show_colnames = F,
                    annotation_col = data.frame(Cluster =paste0('CS', hc.cs.4[ss.dist$order]),
                                                RecurRisk = clin.both[ss.dist$order,'RecurRisk'],
                                                row.names = rownames(ss.mat)[ss.dist$order]),
                    annotation_colors = list(RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red')))
