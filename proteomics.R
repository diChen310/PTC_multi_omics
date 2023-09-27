#############Proteomics analysis, 2020-08-18

data.pro = read.csv(file='./data/proteomics/MWY-22-157D_武汉协和医院74例人组织蛋白组lablefree检测服务/MWY-22-157D/4.Quantitation/Expressed_annotation.csv',
row.names=1,check.names=F)
#nrow:5612
data.pro2 = read.csv(file='./data/proteomics/MWY-22-160D_武汉协和医院74例人组织蛋白磷酸化修饰lablefree检测服务/MWY-22-160D/4.Quantitation/Expressed_annotation.csv',
                    row.names=1,check.names=F)


data.pro.num.0 = apply(data.pro[,4:77], 1, function(a){length(a[a==0])})

data.pro = data.pro[data.pro.num.0/74<0.3,]#nrow:3147

data.pro2.num.0= apply(data.pro2[,4:77], 1, function(a){length(a[a==0])})

data.pro2 = data.pro2[data.pro2.num.0/74<0.3,]#nrow:652

pair.info = read.csv(file='H:/project/multiOmics/data/paired sample info.csv')

samples.tp.p = as.character(pair.info$Tumor)
samples.nt.p = as.character(pair.info$Normal)

pro.exp = as.matrix(data.pro[,c(samples.tp.p,samples.nt.p)])
rownames(pro.exp)=data.pro$Gene
pro.exp = log2(pro.exp+1)


pro2.exp = as.matrix(data.pro2[,c(samples.tp.p,samples.nt.p)])
rownames(pro2.exp) = data.pro2$Gene
pro2.exp = log2(pro2.exp+1)

###############Proteomics###############
####Differential analysis,  paired####

p.values.p <- lapply(1:nrow(pro.exp), function(i){
  wilcox.test(as.double(pro.exp[i,samples.tp.p]),as.double(pro.exp[i,samples.nt.p]),paired = T)$p.value
})

log2FCs.p <- lapply(1:nrow(pro.exp), function(i){
  mean(as.double(pro.exp[i,samples.tp.p]))-mean(as.double(pro.exp[i,samples.nt.p]))
})

res.p <- data.frame(gene = rownames(pro.exp),
                  p.value = as.double(unlist(p.values.p)),
                  FDR = p.adjust(as.double(unlist(p.values.p)),method = 'fdr'),
                  log2FC = as.double(unlist(log2FCs.p)))
res.p = res.p[!is.na(res.p$p.value),]
res.p$significance = 'Not'
res.p[res.p$p.value<0.01 & res.p$log2FC > 1,'significance']='Up'
res.p[res.p$p.value<0.01 & res.p$log2FC < -1,'significance']='Down'
top25 <- res.p %>% filter(significance != 'Not') %>% group_by(significance) %>% top_n(-10, p.value) %>% top_n(10, abs(log2FC))


res.g = res[order(res$p.value),]
res.g = res.g[!duplicated(res.g$gene),]
rownames(res.g)=res.g$gene
#rownames(res)=res$gene
top25$mRNA.fc = res.g[as.character(top25$gene),'Log2FC']

res.sig.p = filter(res.p,significance != 'Not')
res.sig.p$mRNA.fc = res.g[as.character(res.sig.p$gene),'Log2FC']
res.sig.p = filter(res.sig.p,!is.na(mRNA.fc))


ggplot(res.p,aes(x=log2FC,y=-log10(p.value),color=significance))+geom_point()+
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
res.p$mRNA.fc = res.g[as.character(res.p$gene),'Log2FC']

ggplot(res.p,aes(x=log2FC,y=mRNA.fc,color=significance))+geom_point()+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+
  geom_text_repel(aes(x=log2FC,y=mRNA.fc,label = gene,color=significance),top25,size = 2)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))

res.p.sig = res.p[res.p$significance != 'Not',]

cor.test(res.p.sig$log2FC,res.p.sig$mRNA.fc)
#Pearson's product-moment correlation
# 
# data:  res.p.sig$log2FC and res.p.sig$mRNA.fc
# t = 8.7558, df = 279, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.3672127 0.5513363
# sample estimates:
# cor 
# 0.4642759 
cor.test(res.p$log2FC,res.p$mRNA.fc)

####Proteomics
# Pearson's product-moment correlation
# 
# data:  res.p.sig$log2FC and res.p.sig$mRNA.fc
# t = 9.5501, df = 281, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.4015931 0.5782352
# sample estimates:
# cor 
# 0.4950118 

####Phosph
# Pearson's product-moment correlation
# 
# data:  res.p.sig$log2FC and res.p.sig$mRNA.fc
# t = 5.4157, df = 129, p-value = 2.882e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.2795082 0.5605429
# sample estimates:
# cor 
# 0.4303992 

write.csv(res.p,file = './variables/diffProteinExpPaired.csv')


#####PCA

t.pro.exp = t(pro.exp)
t.pro.exp[is.na(t.pro.exp)]=0
t.pro.exp = t.pro.exp[,apply(t.pro.exp, 2, sd)!=0]
colnames(t.pro.exp)=1:ncol(t.pro.exp)
t.pro.exp[is.na(t.pro.exp)]=0


pcr.int = prcomp(t.pro.exp)

pcr.cluster = data.frame(cluster = ifelse(rownames(t.pro.exp) %in% samples.tp,'T','N'),
                         row.names  = rownames(t.pro.exp),
                         PC1 = pcr.int$x[,1],
                         PC2 = pcr.int$x[,2])
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=cluster),size = 1.3)+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()


####Pathway enrichment

res.p = res.p[order(res.p$p.value),]
res.p = res.p[!duplicated(res.p$gene),]
res.p = res.p[order(res.p$log2FC,decreasing = T),]
FC.gsea.p = res.p$log2FC
names(FC.gsea.p)=as.character(res.p$gene)
path2Gene.p = path2Gene[path2Gene$Item %in% names(FC.gsea.p),]
diff.pathway.pro = GSEA(FC.gsea.p,minGSSize = 5,maxGSSize = 200,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene.p[,c(1,2)],pAdjustMethod = 'fdr')
gsea.res.p = diff.pathway.pro@result
gsea.res.p = gsea.res.p[gsea.res.p$pvalue<0.05,]
gsea.res.p$cate = path2cate[as.character(gsea.res.p$Description),'cate']
gsea.res.p = filter(gsea.res.p,cate != '6. Human Diseases')
gsea.top = gsea.res.p %>% group_by(cate) %>% top_n(-15, pvalue)

ggplot(gsea.top,aes(y=ID,x= NES,fill = cate))+geom_bar(stat = 'identity',width = 0.6)+
  theme_cowplot()+facet_grid(vars(cate),space = 'free',scales = 'free')+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = pal_npg('nrc')(5))


######
#https://github.com/saezlab/COSMOS_MSB/tree/main/scripts


###############################################Phosph proteomics######################################################


###############Proteomics###############
####Differential analysis,  paired####

p.values.ph <- lapply(1:nrow(pro2.exp), function(i){
  wilcox.test(as.double(pro2.exp[i,samples.tp.p]),as.double(pro2.exp[i,samples.nt.p]),paired = T)$p.value
})

log2FCs.ph <- lapply(1:nrow(pro2.exp), function(i){
  mean(as.double(pro2.exp[i,samples.tp.p]))-mean(as.double(pro2.exp[i,samples.nt.p]))
})

res.p2 <- data.frame(gene = rownames(pro2.exp),
                    p.value = as.double(unlist(p.values.ph)),
                    FDR = p.adjust(as.double(unlist(p.values.ph)),method = 'fdr'),
                    log2FC = as.double(unlist(log2FCs.ph)))
res.p2 = res.p2[!is.na(res.p2$p.value),]
res.p2$significance = 'Not'
res.p2[res.p2$p.value<0.01 & res.p2$log2FC > 1,'significance']='Up'
res.p2[res.p2$p.value<0.01 & res.p2$log2FC < -1,'significance']='Down'
top25.2 <- res.p2 %>% filter(significance != 'Not') %>% group_by(significance) %>% top_n(25, abs(log2FC))


rownames(res.p)=res.p$gene
#rownames(res)=res$gene
top25.2$pro.fc = res.p[as.character(top25.2$gene),'log2FC']

res.sig.p2 = filter(res.p2,significance != 'Not')
res.sig.p2$pro.fc = res.p[as.character(res.sig.p2$gene),'log2FC']
res.sig.p2 = filter(res.sig.p2,!is.na(pro.fc))


ggplot(res.p2,aes(x=log2FC,y=-log10(p.value),color=significance))+geom_point()+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+geom_hline(yintercept = 2,lty=3)+
  geom_text_repel(aes(x=log2FC,y=-log10(p.value),label = gene,color=significance),top25.2,size = 2)+
  theme_cowplot()+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )
res.p2$pro.fc = res.p[as.character(res.p2$gene),'log2FC']

ggplot(res.p2,aes(x=log2FC,y=pro.fc,color=significance))+geom_point()+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 1,lty=3)+geom_vline(xintercept = -1,lty=3)+
  geom_text_repel(aes(x=log2FC,y=pro.fc,label = gene,color=significance),top25.2,size = 2)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8) )+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))

res.p2.sig = res.p2[res.p2$significance != 'Not',]

cor.test(res.p2.sig$log2FC,res.p2.sig$pro.fc)
# Pearson's product-moment correlation
# 
# data:  res.p2.sig$log2FC and res.p2.sig$pro.fc
# t = 16.919, df = 86, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.8176441 0.9177929
# sample estimates:
# cor 
# 0.8769153 


write.csv(res.p2,file = './variables/diffProteinphosphorylateExpPaired.csv')


#####PCA

t.pro.exp = t(pro.exp2)
t.pro.exp[is.na(t.pro.exp)]=0
t.pro.exp = t.pro.exp[,apply(t.pro.exp, 2, sd)!=0]
colnames(t.pro.exp)=1:ncol(t.pro.exp)
t.pro.exp[is.na(t.pro.exp)]=0


pcr.int = prcomp(t.pro.exp)

pcr.cluster = data.frame(cluster = ifelse(rownames(t.pro.exp) %in% samples.tp,'T','N'),
                         row.names  = rownames(t.pro.exp),
                         PC1 = pcr.int$x[,1],
                         PC2 = pcr.int$x[,2])
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=cluster),size = 1.3)+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()


##########################check correlations

both.genes = intersect(rownames(pro.exp),rownames(pro2.exp))
print(c(nrow(pro.exp),nrow(pro2.exp),length(both.genes)))#3147  652  606

corrs.p2ph = cor(t(pro.exp),t(pro2.exp))

cis.corrs = diag(corrs.p2ph[both.genes,both.genes])
not.p.genes = rownames(pro2.exp)[rownames(pro2.exp) %in% both.genes == F]
not.ph.genes = rownames(pro.exp)[rownames(pro.exp) %in% both.genes == F]
trans.corrs.1 = corrs.p2ph[both.genes,both.genes]
diag(trans.corrs.1)=NA

hist(trans.corrs.1)

trans.corrs.2 = corrs.p2ph[not.ph.genes,both.genes]
trans.corrs.3 = corrs.p2ph[both.genes,not.p.genes]
trans.corrs.4 = corrs.p2ph[not.ph.genes,not.p.genes]

trans.corrs = c(as.double(trans.corrs.1),
                as.double(trans.corrs.2),
                as.double(trans.corrs.3),
                as.double(trans.corrs.4))
diff.c = wilcox.test(trans.corrs,cis.corrs)
boxplot(trans.corrs,cis.corrs,col=c('blue','red'),main=paste('p-value=',signif(diff.c$p.value,4)),
        names=c('same','different'))

#################################################Clustering based on the two omics#####################
clin.p = clin.input[intersect(rownames(clin.input),colnames(pro.exp)),]

clin.display = clin.p[,c('Age','Gender','RecurRisk','Recur','TumorSize',
                            'exTInvasion','LNM','LNM.No','exNInvasion','LNM.3cm','T.stage',
                            'N.stage','M.stage','TNM.stage','TERT_promoter')]
clin.display$CCDC6_RET = ifelse(rownames(clin.display) %in% CCDC6_RET.samples,'TRUE','FALSE')
clin.display$NCOA4_RET = ifelse(rownames(clin.display) %in% NOCA4_RET.samples,'TRUE','FALSE')

clin.display = arrange(clin.display,TERT_promoter,T.stage,N.stage,M.stage,TNM.stage,
                       Age,Gender,RecurRisk,Recur,TumorSize)
clin.display$LNM = as.character(clin.display$LNM)
clin.display$exTInvasion = as.character(clin.display$exTInvasion)
clin.display$exNInvasion = as.character(clin.display$exNInvasion)
clin.display$Recur = as.character(clin.display$Recur)
clin.display$TERT_promoter[clin.display$TERT_promoter==""]='No'
clin.display$LNM.3cm = as.character(clin.display$LNM.3cm)
cols = rownames(clin.display)
ha.pp<-HeatmapAnnotation(df=clin.display[cols,],annotation_name_gp = gpar(fontsize = 6),
                         show_legend = TRUE,
                         show_annotation_name = TRUE,
                         simple_anno_size = unit(0.2, "cm"),
                         col = list(Cluster2=c('P'='red',"B"='blue',"M"='gray'),
                                    Age = col_fun,Cluster=cols.t8,
                                    Gender=c("Male"='cyan',"Female"='pink'),
                                    TNM.stage=c("I"="#00A087FF",'II'="#3C5488FF","III"="#8491B4FF","IV"="#DC0000FF"),
                                    T.stage=c("T1"="#00A087FF",'T2'="#3C5488FF","T3"="#8491B4FF","T4"="#DC0000FF"),
                                    M.stage=c('M0'="#00A087FF","M1"="#3C5488FF"),
                                    N.stage=c("N0"="#00A087FF",'N1a'="#3C5488FF","N1b"="#8491B4FF"),
                                    exNInvasion=c("TRUE"='black',"FALSE"='white'),
                                    LNM.No=col_fun3,
                                    LNM.3cm = c("TRUE"='black',"FALSE"='white'),
                                    LNM=c("TRUE"='black',"FALSE"='white'),
                                    exTInvasion=c("TRUE"='black',"FALSE"='white'),
                                    TumorSize=col_fun2,
                                    Recur=c("TRUE"='black',"FALSE"='white'),
                                    RecurRisk=c("Median"='cyan',"Low"='white',"High"='red'),
                                    TERT_promoter=c('TERT C228T'='black',"No"='white'),
                                    CCDC6_RET=c("TRUE"='black',"FALSE"='white'),
                                    NCOA4_RET=c("TRUE"='black',"FALSE"='white'))
)


ha = Heatmap(rbind(t(scale(t(data1[,cols]))),t(scale(t(data2[,cols])))),row_split = 4,
        column_split = 3,
        top_annotation = ha.pp,
        show_row_dend = F,
        show_column_dend = T,
        show_row_names=F,
        show_column_names=F,
        cluster_rows = T,
        cluster_columns = T,
        row_names_gp = gpar(fontsize = 6),
        width=unit(8,'cm'),
        height=unit(9,'cm'),
        clustering_distance_columns = c('pearson'))
ha.clusters= column_order(ha)

results.cluster = unlist(lapply(1:length(cols), function(a){
  for(i in 1:length(ha.clusters)){
    if(a %in% ha.clusters[[i]]){
      return(i)
    }
  }
}))

clin.p[cols,]$Cluster=results.cluster
#clin.p$status = clin.p$Recur

km.by.cluster = survfit(Surv(days_RFS, status) ~ Cluster, data = clin.p)
survdiff(Surv(days_RFS, status) ~ Cluster, data = clin.p,rho = 0)
for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter')){
  print(c(feature,chisq.test(clin.p$Cluster,clin.p[,feature])$p.value))
}

# [1] "Gender"            "0.550288600710807"
# [1] "RecurRisk"          "0.0157919670792115"
# [1] "Recur"             "0.327629128581164"
# [1] "exTInvasion"       "0.185445132805489"
# [1] "LNM"               "0.509421675314032"
# [1] "exNInvasion"         "0.00939907191651143"
# [1] "LNM.3cm"          "0.41351979798438"
# [1] "T.stage"           "0.972314378476886"
# [1] "N.stage"           "0.314432682203122"
# [1] "M.stage"           "0.658253944633061"
# [1] "TNM.stage"         "0.644983079123938"
# [1] "TERT_promoter"    "0.90730934288795"
c.data = clin.p[,c('Age','TumorSize','LNM.No','Cluster')]
c.data$Cluster = factor(c.data$Cluster)
library(reshape2)
c.data=melt(c.data,measure.vars = c('Age','TumorSize','LNM.No'))
ggplot(c.data,aes(x=Cluster,y=value,fill=Cluster))+geom_boxplot()+facet_grid(~variable,scales = 'free')+
  stat_compare_means()+theme_bw()+scale_fill_manual(values = pal_npg('nrc')(3))


# data1 = FSbyMAD(pro.exp[,cols], cut.type="topk",value=700)
# data2 = FSbyMAD(pro2.exp[,cols], cut.type="topk",value=150)
# TCA = list(proExp = t(scale(t(data1))),
#            phExp = t(scale(t(data2))))
# set.seed(123)
# result=ExecuteSNF(datasets=TCA, clusterNum =4,K=25)
# clin.p[cols,]$Cluster=result$group
# clin.p$status = clin.p$Recur
# 
# km.by.cluster = survfit(Surv(days_RFS, status) ~ Cluster, data = clin.p)
# survdiff(Surv(days_RFS, status) ~ Cluster, data = clin.p,rho = 0)
# ggsurvplot(km.by.cluster,
#            pval = T, conf.int = FALSE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw() ,# Change ggplot2 theme
#            palette = pal_npg('nrc')(10)
# )
# 
# clin.display$Cluster = paste('M',clin.p[cols,'Cluster'])
# cols.t8 = pal_npg('nrc')(length(unique(clin.display$Cluster)))
# names(cols.t8) = unique(clin.display$Cluster)
# clin.display = arrange(clin.display,Cluster,TERT_promoter,T.stage,N.stage,M.stage,TNM.stage,
#                        Age,Gender,RecurRisk,Recur,TumorSize)
# 
# 
# cols = rownames(clin.display)
# ha.pp<-HeatmapAnnotation(df=clin.display[cols,],annotation_name_gp = gpar(fontsize = 6),
#                          show_legend = TRUE,
#                          show_annotation_name = TRUE,
#                          simple_anno_size = unit(0.2, "cm"),
#                          col = list(Cluster2=c('P'='red',"B"='blue',"M"='gray'),
#                                     Age = col_fun,Cluster=cols.t8,
#                                     Gender=c("Male"='cyan',"Female"='pink'),
#                                     TNM.stage=c("I"="#00A087FF",'II'="#3C5488FF","III"="#8491B4FF","IV"="#DC0000FF"),
#                                     T.stage=c("T1"="#00A087FF",'T2'="#3C5488FF","T3"="#8491B4FF","T4"="#DC0000FF"),
#                                     M.stage=c('M0'="#00A087FF","M1"="#3C5488FF"),
#                                     N.stage=c("N0"="#00A087FF",'N1a'="#3C5488FF","N1b"="#8491B4FF"),
#                                     exNInvasion=c("TRUE"='black',"FALSE"='white'),
#                                     LNM.No=col_fun3,
#                                     LNM.3cm = c("TRUE"='black',"FALSE"='white'),
#                                     LNM=c("TRUE"='black',"FALSE"='white'),
#                                     exTInvasion=c("TRUE"='black',"FALSE"='white'),
#                                     TumorSize=col_fun2,
#                                     Recur=c("TRUE"='black',"FALSE"='white'),
#                                     RecurRisk=c("Median"='cyan',"Low"='white',"High"='red'),
#                                     TERT_promoter=c('TERT C228T'='black',"No"='white'),
#                                     CCDC6_RET=c("TRUE"='black',"FALSE"='white'),
#                                     NCOA4_RET=c("TRUE"='black',"FALSE"='white'))
# )
# 
# imp.pros<-itemImportance(t(scale(t(data1))),rownames(clin.p),clin.p,'Cluster')
# imp.pp.pros<-row.names(imp.pros[rev(order(as.double(imp.pros[,'MeanDecreaseGini']))),])[1:40]
# rownames(data2)=paste('p',rownames(data2))
# imp.phs<-itemImportance(t(scale(t(data2))),rownames(clin.p),clin.p,'Cluster')
# imp.pp.phs<-row.names(imp.phs[rev(order(as.double(imp.phs[,'MeanDecreaseGini']))),])[1:40]
# 
# Heatmap(rbind(t(scale(t(data1[imp.pp.pros,cols]))),t(scale(t(data2[imp.pp.phs,cols])))),row_split = 4,
#         top_annotation = ha.pp,
#         show_row_dend = F,
#         show_column_dend = F,
#         show_row_names=T,
#         show_column_names=F,
#         cluster_rows = T,
#         cluster_columns = F,
#         row_names_gp = gpar(fontsize = 6),
#         width=unit(10,'cm'),
#         height=unit(13,'cm'))

save.image(file = 'mRNA&metabolism&Pro&Phosph.RData')
