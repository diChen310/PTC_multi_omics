library(mixOmics)
library(ggsci)
FSbyMAD = function(omics,topN){
  mad.res = apply(omics, 1, mad)
  omics = omics[order(mad.res,decreasing = T),]
  return(omics[1:topN,])
}
change2HML = function(merge.risk.res,tag){
  
  
  res.h = merge.risk.res[merge.risk.res$Risk == 'H',]
  rownames(res.h)=res.h$name
  
  res.m = merge.risk.res[merge.risk.res$Risk == 'M',]
  rownames(res.m)=res.m$name
  
  res.l = merge.risk.res[merge.risk.res$Risk == 'L',]
  rownames(res.l)=res.l$name
  
  res.new = res.h
  
  res.new$H=res.new[,tag]
  res.new$M = res.m[rownames(res.h),tag]
  res.new$L = res.l[rownames(res.h),tag]
  
  return(res.new[,c('id','name','H','M','L')])
  
}
#res.mrna,res.meta,res.protein, res.phosph from mergeDiffResAmongOmics.R
sig.genes = res.mrna[res.mrna$p.value<0.01,'gene']
sig.metas = res.meta[res.meta$p.value<0.01,'name']
sig.pros = as.character(res.protein[res.protein$p.value<0.01,'gene'])
sig.phos = as.character(res.phosph[res.phosph$p.value<0.01,'gene'])

#gene.res.risk, meta.res.risk,pro.res.risk, ph.res.risk from  recurRishAnalysis.R
sig.genes2 = unique(gene.res.risk[gene.res.risk$p.value<0.01,'id'])
sig.metas2 = unique(meta.res.risk[meta.res.risk$p.value<0.01,'id'])
sig.pros2 = unique(pro.res.risk[pro.res.risk$p.value<0.01,'id'])
sig.phos2 = unique(ph.res.risk[ph.res.risk$p.value<0.01,'id'])

sig.genes = union(sig.genes,sig.genes2)#5907
sig.metas = union(sig.metas,sig.metas2)#376
sig.pros = union(sig.pros,sig.pros2)#1987
sig.phos = union(sig.phos,sig.phos2)#462


both.samples = intersect(intersect(colnames(data.meta),colnames(data.mrna)),colnames(pro.exp))

pair.info.all = pair.info[pair.info$Tumor %in% both.samples & pair.info$Normal %in% both.samples,]

all.samples.t = as.character(pair.info.all$Tumor)
all.samples.n = as.character(pair.info.all$Normal)

high.cols = intersect(all.samples.t,clin.input$ID[clin.input$RecurRisk == 'High'])
m.cols = intersect(all.samples.t,clin.both$ID[clin.both$RecurRisk == 'Median'])
l.cols = intersect(all.samples.t,clin.both$ID[clin.both$RecurRisk == 'Low'])

X1 = data.mrna[sig.genes,c(high.cols,m.cols,l.cols,all.samples.n)]
X2 = data.meta[sig.metas,c(high.cols,m.cols,l.cols,all.samples.n)]
X3 = pro.exp[sig.pros,c(high.cols,m.cols,l.cols,all.samples.n)]
X4 = pro2.exp[sig.phos,c(high.cols,m.cols,l.cols,all.samples.n)]

X1 = FSbyMAD(X1, topN = 500)
X3 = FSbyMAD(X3, topN=500)

rownames(X3) = paste(rownames(X3),'Pro')
rownames(X4) = paste(rownames(X4),'Phosph')
dim(X1)
dim(X2)
dim(X3)
dim(X4)




X <- list(mrna = t(X1), metabolite = t(X2), protein = t(X3),phosph = t(X4))
Y <- c(rep('H',length(high.cols)),rep('M',length(m.cols)),rep('L',length(l.cols)),
       rep('N',length(all.samples.n)))
Y = factor(Y,levels = c('H','M','L','N'))

design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

list.keepX = list(metabolite = rep(20,3),
                  mrna = rep(20,3),
                  protein= rep(20,3),
                  phosph = rep(20,3))


# set the optimised DIABLO model
final.diablo.model = block.splsda(X, Y , ncomp = 3, 
                                  keepX = list.keepX, design = design)

# the features selected to form the first component
selectVar(final.diablo.model, block = 'mrna', comp = 1)$mrna$name 


plotIndiv(final.diablo.model, ind.names = FALSE,style = 'graphics',legend = T,
          col = pal_npg()(4))
####Figure 4a
plotDiablo(final.diablo.model, ncomp = 1,col=pal_npg()(4))
circosPlot(final.diablo.model, cutoff=0.8)
# extended example:
cimDiablo(final.diablo.model, 
          color.blocks = pal_npg('nrc')(4), 
          comp = c(1,2), margin=c(8,20), 
          legend.position = "right",
          transpose = T)
plotLoadings(final.diablo.model, comp = 2, contrib = "max")
plotLoadings(final.diablo.model, comp = 1, contrib = "max")
plotLoadings(final.diablo.model, comp = 3, contrib = "max")
#plotLoadings(final.diablo.model, comp = 4, contrib = "max")


#####Re-analyze the heatmap####

top.mrna = c(selectVar(final.diablo.model, block = 'mrna', comp = 1)$mrna$name,
             selectVar(final.diablo.model, block = 'mrna', comp = 2)$mrna$name)
top.mrna = unique(top.mrna)

top.meta=  c(selectVar(final.diablo.model, block = 'metabolite', comp = 1)$metabolite$name,
             selectVar(final.diablo.model, block = 'metabolite', comp = 2)$metabolite$name)
top.meta = unique(top.meta)

top.pro=  c(selectVar(final.diablo.model, block = 'protein', comp = 1)$protein$name,
            selectVar(final.diablo.model, block = 'protein', comp = 2)$protein$name)
top.pro = unique(top.pro)

top.pro2=  c(selectVar(final.diablo.model, block = 'phosph', comp = 1)$phosph$name,
             selectVar(final.diablo.model, block = 'phosph', comp = 2)$phosph$name)
top.pro2 = unique(top.pro2)


ht.input = rbind(X1[top.mrna,],
                 X2[top.meta,],
                 X3[top.pro,],
                 X4[top.pro2,])
row.type = data.frame(Item = rownames(ht.input),
                      Type = c(rep('mRNA',length(top.mrna)),
                               rep('metabolite',length(top.meta)),
                               rep('protein',length(top.pro)),
                               rep('phosphprotein',length(top.pro2))),row.names = rownames(ht.input))
ht.input = t(scale(t(ht.input)))
library(pheatmap)
library(circlize)
library(reshape2)
corr.mat =cor(t(ht.input),method = 'pearson')

hc = hclust(d =dist(ht.input))

hc.cs.4 = cutree(hc,k=4)
hc.cs.4[hc$order]

pheatmap(ht.input[hc$order,c(high.cols,m.cols,l.cols,all.samples.n)],cluster_rows = F,cluster_cols = F,
         color = colorRamp2(c(-2,0,2),c('blue','white','red'))(seq(-4,4,0.02)),show_rownames = T,
         show_colnames = F,
         annotation_row = data.frame(Type = row.type$Type[hc$order],
                                     cluster = paste0('C',hc.cs.4[hc$order]),
                                     row.names = rownames(ht.input[hc$order,])),
         annotation_col = data.frame(From = Y,
                                     row.names = c(high.cols,m.cols,l.cols,all.samples.n)),
         annotation_colors = list(From = c('H'='#E64B35FF','M'='#4DBBD5FF','L'="#00A087FF",'N'="#3C5488FF"),
                                  cluster=c('C1'="#00468BFF",'C2'="#ED0000FF",'C3'="#42B540FF",'C4'="#0099B4FF"),
                                  Type = c('metabolite'="#00468BFF",'mRNA'="#ED0000FF",'phosphprotein'="#42B540FF",'protein'="#0099B4FF")),
         fontsize_row = 6,
         fontsize_col = 6)

pro.res.risk$name = paste(pro.res.risk$id,'Pro')
ph.res.risk$name = paste(ph.res.risk$id,'Phosph')
rownames(res.protein) = paste(res.protein$gene,'Pro')
rownames(res.phosph) = paste(res.phosph$gene,'Phosph')
rownames(res.mrna)=res.mrna$gene
colnames(res.protein)=c("gene","p.value", "FDR"  ,"Log2FC" ,      "significance" ,"mRNA.fc"     )
colnames(res.phosph)=c("gene","p.value", "FDR"  ,"Log2FC" ,      "significance" ,"pro.fc"     )

merge.t2n.res = rbind(res.mrna[,c('p.value','Log2FC')],
                      res.meta[,c('p.value','Log2FC')],
                      res.protein[,c('p.value','Log2FC')],
                      res.phosph[,c('p.value','Log2FC')]
                      )

merge.risk.res = rbind(meta.res.risk[,colnames(pro.res.risk)],
                       gene.res.risk[,colnames(pro.res.risk)],
                       pro.res.risk[,colnames(pro.res.risk)],
                       ph.res.risk[,colnames(pro.res.risk)])
merge.fc = change2HML(merge.risk.res ,'FC')
merge.p = change2HML(merge.risk.res,'p.value')
merge.fc$N = merge.t2n.res[rownames(merge.fc),'Log2FC']
merge.p$N = merge.t2n.res[rownames(merge.p),'p.value']

rows.input = rownames(ht.input[hc$order,])
mat.fc = merge.fc[rows.input,c('H','M','L','N')]
mat.p = merge.p[rows.input,c('H','M','L','N')]
mat.fc[mat.p>0.05]=NA


pheatmap(mat.fc,cluster_rows = F,cluster_cols = F,
         color = hcl.colors(100,palette = 'Blue-Red 3'),
         gaps_col = c(1,2,3))

####Figure 4b
pheatmap(ht.input[hc$order,c(high.cols,m.cols,l.cols,all.samples.n)],cluster_rows = F,cluster_cols = F,
         color = colorRamp2(c(-2,0,2),c('blue','white','red'))(seq(-4,4,0.02)),show_rownames = F,
         show_colnames = F,
         annotation_row = data.frame(Type = row.type$Type[hc$order],
                                           cluster = paste0('C',hc.cs.4[hc$order]),
                                           row.names = rownames(ht.input[hc$order,])),
         annotation_col = data.frame(From = Y,
                                     row.names = c(high.cols,m.cols,l.cols,all.samples.n)),
         annotation_colors = list(From = c('H'='#E64B35FF','M'='#4DBBD5FF','L'="#00A087FF",'N'="#3C5488FF"),
                                  cluster=c('C1'="#00468BFF",'C2'="#ED0000FF",'C3'="#42B540FF",'C4'="#0099B4FF"),
                                  Type = c('metabolite'="#00468BFF",'mRNA'="#ED0000FF",'phosphprotein'="#42B540FF",'protein'="#0099B4FF")
                                 ),
         fontsize_row = 6,
         fontsize_col = 6)
pheatmap(mat.fc,cluster_rows = F,cluster_cols = F,
         color = hcl.colors(100,palette = 'Blue-Red 3'),show_rownames = F,
         annotation_col = data.frame(From = colnames(mat.fc),
                                     row.names = colnames(mat.fc)),
         gaps_col = c(1,2,3),na_col = 'white')

#####The molecule module network#### Figure 4C
cluster.items = data.frame(Type = row.type$Type[hc$order],
                           cluster = paste0('C',hc.cs.4[hc$order]),
                           row.names = rownames(ht.input[hc$order,]))

cluster.i = 'C3' ####Try different clusters

items = rownames(cluster.items[cluster.items$cluster == cluster.i,])

corr.items = cor(t(ht.input[items,]),method = 'spearman')
corr.items = melt(corr.items,measure.vars = colnames(corr.items))
corr.items = corr.items[corr.items$Var1 != corr.items$Var2,]
corr.items.high = corr.items[corr.items$value > quantile(corr.items$value,0.3),]
corr.items.high = group_by(corr.items.high,Var1) %>% top_n(5,value)                        

library(igraph)
library(ggraph)
library(tidygraph)
#library(ggnewscale)
g.net = graph_from_data_frame(corr.items.high[,c(1,2)])
sub.g = as_tbl_graph(g.net,directed = F)

sub.g %>%
  activate(nodes) %>%
  mutate(group=cluster.items[names(V(sub.g)),'Type']) %>%
  mutate(name2=unlist(lapply(names(V(sub.g)), function(a){
    a=gsub(' Pro','',a)
    a=gsub(' Phosph','',a)
    return(a)
  }))) %>%
  activate(edges) %>%
  mutate(weight=corr.items.high$value) %>%
  ggraph(layout = "graphopt") + 
  geom_edge_link(aes(colour = weight),width=0.5,edge_alpha=0.6)+
  scale_edge_colour_gradient2(low = 'blue',high = 'red')+
  geom_node_point(aes(colour = group),size=2) +
  geom_node_text(aes(label = name2),size=2, repel = TRUE)+
  scale_color_manual(values =c('metabolite'="#00468BFF",'mRNA'="#ED0000FF",'phosphprotein'="#42B540FF",'protein'="#0099B4FF") )+
  theme_graph()
