library(maftools)
library(dplyr)
setwd('H:/project/multiOmics')


####Read clinical information
clin.d = read.csv(file = './data/clin.info.csv',row.names = 1)
clin.d$Tumor_Sample_Barcode = rownames(clin.d)
clin.d$Status = ifelse(clin.d$Recur,1,0)
####Read gene fusions recognized by starFusion and Arriba, and merged with clinical information
sfc = read.csv(file='./variables/fusion_sfc.csv')
sfc$Sample = unlist(lapply(sfc$Sample,function(a){
  a=gsub('-','_',a)
  if(substr(a,1,1) %in% as.character(1:9)){
    a=paste0('S',a)
  }
  return(a)
}))
fusion.tumor = table(sfc$FusionName)
fusion.tumor.freq = sfc %>% filter(Gene1A != Gene1B) %>% group_by(Fusion_Type,FusionName) %>% dplyr::count()
fusion.tumor.freq.top = fusion.tumor.freq %>% group_by(Fusion_Type) %>% top_n(10,n) %>% arrange(n) %>% filter(n>1)
fusion.tumor.freq.top$FusionName = factor(fusion.tumor.freq.top$FusionName,
                                          levels = fusion.tumor.freq.top$FusionName)

top.fusion.mat = matrix(0,nrow = nrow(clin.d),ncol = nrow(fusion.tumor.freq.top),
                        dimnames = list(clin.d$Tumor_Sample_Barcode,fusion.tumor.freq.top$FusionName))

for(i in 1:ncol(top.fusion.mat)){
  item = colnames(top.fusion.mat)[i]
  item.samples = unique(sfc[sfc$FusionName == item,'Sample'])
  item.samples = intersect(item.samples,rownames(top.fusion.mat))
  top.fusion.mat[item.samples,i]=1
}
colnames(top.fusion.mat)=gsub('--','_',colnames(top.fusion.mat))
top.fusion.mat = top.fusion.mat[,-1]

rownames(clin.d)=clin.d$Tumor_Sample_Barcode
clin.d = cbind(clin.d,as.data.frame(top.fusion.mat))
# clin.cccc = read.csv(file = './variables/clin.display_TDSBSRS.csv',row.names = 1)
# 
# rownames(clin.cccc)=gsub('-','_',rownames(clin.cccc))
# rownames(clin.cccc)[substr(rownames(clin.cccc),1,1) %in% c('1','2','3','4')] = paste0('S',rownames(clin.cccc)[substr(rownames(clin.cccc),1,1) %in% c('1','2','3','4')])
# 
# rownames(clin.cccc)[rownames(clin.cccc) %in% rownames(clin.d) == F] = 'S20__11'
# clin.d$Cluster = clin.cccc[rownames(clin.d),'Cluster']

var.annovar.maf <- annovarToMaf(annovar = "I:/project/fudanOmics/results/annovar.varscan_all_v2.vcf", #annovar.merge_all_v2_strick.vcf
                             refBuild = 'hg38',
                             tsbCol = 'Tumor_Sample_Barcode', 
                             table = 'refGene',
                             MAFobj = T,
                             sampleAnno = clin.d)

cosmic.genes = read.csv(file='./data/cancer_gene_census.csv')
cosmic.genes.list = cosmic.genes$Gene.Symbol
sig.genes = cosmic.genes.list[cosmic.genes.list %in% var.annovar.maf@gene.summary$Hugo_Symbol]
ptc = subsetMaf(var.annovar.maf,genes=sig.genes)
oncostrip(ptc,top=30,
          clinicalFeatures = c('RecurRisk'
          ),
          removeNonMutated=F,
          anno_height = 0.3,
          annoBorderCol='White',
          showTumorSampleBarcodes=T,
          sortByAnnotation=T
)



library(circlize)
col_fun = colorRamp2(c(0, 70), c("white", "red"))
col_fun2 = colorRamp2(c(0, 11), c("white", "red"))
col_fun3 = colorRamp2(c(0, 39), c("white", "red"))
order.list = read.table(file='./variables/mutation orders.txt')
order.list = as.character(order.list$V1)
####Figure 1a
oncostrip(var.annovar.maf,genes=top.30.cosmic,
          clinicalFeatures = c('TERT_promoter'    
                              ),
          removeNonMutated=F,fontSize=0.6,
          anno_height = 0.3,
          annoBorderCol='White',
          showTumorSampleBarcodes=F,
          sampleOrder = order.list,
          annotationColor=list(Age = col_fun,
                                     Gender=c("Male"='cyan',"Female"='pink'),
                                     TNM.stage=c("I"="#00A087FF",'II'="#3C5488FF","III"="#8491B4FF","IV"="#DC0000FF"),
                                     T.stage=c("T1"="#00A087FF",'T2'="#3C5488FF","T3"="#8491B4FF","T4"="#DC0000FF"),
                                     M.stage=c('M0'="#00A087FF","M1"="#3C5488FF"),
                                     N.stage=c("N0"="#00A087FF",'N1a'="#3C5488FF","N1b"="#8491B4FF"),
                                     exNInvasion=c("TRUE"='black',"FALSE"='lightgrey'),
                                     LNM.No=col_fun3,
                                     LNM.3cm = c("TRUE"='black',"FALSE"='lightgrey'),
                                     LNM=c("TRUE"='black',"FALSE"='lightgrey'),
                                     exTInvasion=c("TRUE"='black',"FALSE"='lightgrey'),
                                     TumorSize=col_fun2,
                                     Recur=c("TRUE"='black',"FALSE"='lightgrey'),
                                     RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red'),
                                     TERT_promoter=c('TERT C228T'='blue',"No"='white'),
                                     TLK2_FAM157A=c('1'='#8491B4FF','0'='lightgrey'),
                               CCDC6_RET=c('1'='#8491B4FF','0'='lightgrey'),
                               ZNF33B_NCOA4=c('1'='#8491B4FF','0'='lightgrey'),
                               NTRK3_ETV6=c('1'='#8491B4FF','0'='lightgrey'),
                               TPR_NTRK1=c('1'='#8491B4FF','0'='lightgrey'),
                               ETV6_NTRK3=c('1'='#8491B4FF','0'='lightgrey'),
                               FBXO25_SEPTIN14=c('1'='#8491B4FF','0'='lightgrey'),
                               NCOA4_RET=c('1'='#8491B4FF','0'='lightgrey')
                               ))
oncostrip(var.annovar.maf,genes=top.30.cosmic,
          clinicalFeatures = c('RecurRisk','Gender',
                               "Recur", "exTInvasion", "LNM"   , "exNInvasion" ,
                               "LNM.3cm" ,"T.stage" ,"N.stage" ,"M.stage", "TNM.stage"      
          ),
          removeNonMutated=F,fontSize=0.6,
          anno_height = 4,
          annoBorderCol='White',
          showTumorSampleBarcodes=F,
          sampleOrder = order.list,
          annotationColor=list(Age = col_fun,
                               Gender=c("Male"='cyan',"Female"='pink'),
                               TNM.stage=c("I"="#00A087FF",'II'="#3C5488FF","III"="#8491B4FF","IV"="#DC0000FF"),
                               T.stage=c("T1"="#00A087FF",'T2'="#3C5488FF","T3"="#8491B4FF","T4"="#DC0000FF"),
                               M.stage=c('M0'="#00A087FF","M1"="#3C5488FF"),
                               N.stage=c("N0"="#00A087FF",'N1a'="#3C5488FF","N1b"="#8491B4FF"),
                               exNInvasion=c("TRUE"='black',"FALSE"='lightgrey'),
                               LNM.No=col_fun3,
                               LNM.3cm = c("TRUE"='black',"FALSE"='lightgrey'),
                               LNM=c("TRUE"='black',"FALSE"='lightgrey'),
                               exTInvasion=c("TRUE"='black',"FALSE"='lightgrey'),
                               TumorSize=col_fun2,
                               Recur=c("TRUE"='black',"FALSE"='lightgrey'),
                               RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red'),
                               TERT_promoter=c('TERT C228T'='blue',"No"='white'),
                               TLK2_FAM157A=c('1'='#8491B4FF','0'='lightgrey'),
                               CCDC6_RET=c('1'='#8491B4FF','0'='lightgrey'),
                               ZNF33B_NCOA4=c('1'='#8491B4FF','0'='lightgrey'),
                               NTRK3_ETV6=c('1'='#8491B4FF','0'='lightgrey'),
                               TPR_NTRK1=c('1'='#8491B4FF','0'='lightgrey'),
                               ETV6_NTRK3=c('1'='#8491B4FF','0'='lightgrey'),
                               FBXO25_SEPTIN14=c('1'='#8491B4FF','0'='lightgrey'),
                               NCOA4_RET=c('1'='#8491B4FF','0'='lightgrey')
          ))
oncostrip(var.annovar.maf,genes=top.30.cosmic,
          clinicalFeatures = c("CCDC6_RET","NCOA4_RET" ,"FBXO25_SEPTIN14"  , 
                               "TLK2_FAM157A" ,        "ZNF33B_NCOA4"   ,      "NTRK3_ETV6" ,         
                               "TPR_NTRK1" ,           "ETV6_NTRK3"          
          ),
          removeNonMutated=F,fontSize=0.6,
          anno_height = 2.8,
          annoBorderCol='White',
          showTumorSampleBarcodes=F,
          sampleOrder = order.list,
          annotationColor=list(Age = col_fun,
                               Gender=c("Male"='cyan',"Female"='pink'),
                               TNM.stage=c("I"="#00A087FF",'II'="#3C5488FF","III"="#8491B4FF","IV"="#DC0000FF"),
                               T.stage=c("T1"="#00A087FF",'T2'="#3C5488FF","T3"="#8491B4FF","T4"="#DC0000FF"),
                               M.stage=c('M0'="#00A087FF","M1"="#3C5488FF"),
                               N.stage=c("N0"="#00A087FF",'N1a'="#3C5488FF","N1b"="#8491B4FF"),
                               exNInvasion=c("TRUE"='black',"FALSE"='lightgrey'),
                               LNM.No=col_fun3,
                               LNM.3cm = c("TRUE"='black',"FALSE"='lightgrey'),
                               LNM=c("TRUE"='black',"FALSE"='lightgrey'),
                               exTInvasion=c("TRUE"='black',"FALSE"='lightgrey'),
                               TumorSize=col_fun2,
                               Recur=c("TRUE"='black',"FALSE"='lightgrey'),
                               RecurRisk=c("Median"='cyan',"Low"='purple',"High"='red'),
                               TERT_promoter=c('TERT C228T'='blue',"No"='white'),
                               TLK2_FAM157A=c('1'='#8491B4FF','0'='lightgrey'),
                               CCDC6_RET=c('1'='#8491B4FF','0'='lightgrey'),
                               ZNF33B_NCOA4=c('1'='#8491B4FF','0'='lightgrey'),
                               NTRK3_ETV6=c('1'='#8491B4FF','0'='lightgrey'),
                               TPR_NTRK1=c('1'='#8491B4FF','0'='lightgrey'),
                               ETV6_NTRK3=c('1'='#8491B4FF','0'='lightgrey'),
                               FBXO25_SEPTIN14=c('1'='#8491B4FF','0'='lightgrey'),
                               NCOA4_RET=c('1'='#8491B4FF','0'='lightgrey')
          ))

prog_geneset = survGroup(maf = ptc, top = 40, geneSetSize = 1, time = "days_RFS", Status = "Status", verbose = FALSE)

mafSurvival(maf = var.annovar.maf, genes = c('MUC16'), time = 'days_RFS', Status = 'Status', isTCGA = F)
prog_geneset = survGroup(maf = ptc, top = 40, geneSetSize = 2, time = "days_RFS", Status = "Status", verbose = FALSE)

mafSurvGroup(maf = var.annovar.maf, geneSet = c("BRAF", "MUC16"), time = "days_RFS", Status = "Status")

OncogenicPathways(maf=var.annovar.maf)
PlotOncogenicPathways(maf=var.annovar.maf,pathways = "RTK-RAS")



####Compare to TCGA-THCA####

load('./variables/THCA.mut.maf.RData')

coBarplot(m1 = ptc, m2 = THCA.mut,
          m1Name = 'Our', m2Name = 'TCGA')
####Clinical enrichment####
for(item in c('Gender',"RecurRisk" ,"Recur","exTInvasion" ,    "LNM",        
              "exNInvasion" ,"LNM.3cm","T.stage" ,"N.stage",             
              "M.stage", "TNM.stage","TERT_promoter")){
  ce = clinicalEnrichment(maf = var.annovar.maf, clinicalFeature = item)
  print(ce)
}

####mutation enrichment in one risk type#### Figure 1b
clin.d.short = clin.d[rownames(clin.d) %in% var.annovar.maf@variants.per.sample$Tumor_Sample_Barcode ,]
risk.p.mat = matrix(nrow = 3,ncol = 30,dimnames = list(names(table(clin.d.short$RecurRisk)),top.30.cosmic))
for(gene in top.30.cosmic){
  mut.samples =  as.character(subsetMaf(var.annovar.maf,genes =gene)@data$Tumor_Sample_Barcode)
  clin.d.short$geneMut = rownames(clin.d.short) %in% mut.samples
  table(clin.d.short$geneMut,clin.d.short$RecurRisk)
  # print(c(gene,chisq.test(clin.d.short$geneMut,clin.d.short$RecurRisk)$p.value))
  # print(c(gene,chisq.test(clin.d.short$geneMut,clin.d.short$RecurRisk %in% c('High','Median'))$p.value))
  
  risk.p.mat[,gene]=check(clin.d.short$geneMut,clin.d.short$RecurRisk)
}
risk.p.mat2 = matrix(nrow = 3,ncol = 9,dimnames = list(names(table(clin.d.short$RecurRisk)),c('TERT_promoter',"CCDC6_RET","NCOA4_RET" ,"FBXO25_SEPTIN14"  , 
                                                                                              "TLK2_FAM157A" ,        "ZNF33B_NCOA4"   ,      "NTRK3_ETV6" ,         
                                                                                              "TPR_NTRK1" ,           "ETV6_NTRK3")))
for(feature in c('TERT_promoter',"CCDC6_RET","NCOA4_RET" ,"FBXO25_SEPTIN14"  , 
                 "TLK2_FAM157A" ,        "ZNF33B_NCOA4"   ,      "NTRK3_ETV6" ,         
                 "TPR_NTRK1" ,           "ETV6_NTRK3")){
  risk.p.mat2[,feature]=check(clin.d[,feature],clin.d$RecurRisk)
}
risk.p.mat = cbind(risk.p.mat,risk.p.mat2)

risk.p.mat =0- log10(risk.p.mat)
risk.p.mat[is.infinite(risk.p.mat)]=max(risk.p.mat[!is.infinite(risk.p.mat)])+1
pheatmap(risk.p.mat)

risk.p.mat = melt(risk.p.mat,measure.vars = colnames(risk.p.mat))
risk.p.mat.high = risk.p.mat[risk.p.mat$value> -log10(0.05),]
ggplot(risk.p.mat,aes(x=Var2,y=Var1,fill=value,size=value))+geom_raster(width=0.6,height=0.6)+
  geom_text(aes(label=ifelse(value> -log10(0.05),'*','')),size=3)+
  theme_bw()+
  theme(axis.text = element_text(size=7),
        axis.text.x = element_text(angle = 30, hjust=1),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))+
  scale_fill_gradient2(low = 'blue',high = 'red')


#####MUC16 relevant clin#### Figure 1c
clin.d = clin.d[clin.d$Tumor_Sample_Barcode %in% var.annovar.maf@data$Tumor_Sample_Barcode,]
MUC16.mut.samples = unique(as.character(subsetMaf(var.annovar.maf,genes = c('MUC16'))@data$Tumor_Sample_Barcode))
clin.d$MUC16 = clin.d$Tumor_Sample_Barcode %in% MUC16.mut.samples
for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter')){
  print(c(feature,check(clin.d$MUC16,clin.d[,feature])))
}
plot.data = melt(clin.d,measure.vars = c('RecurRisk','Recur','LNM.3cm',
                                               'T.stage','N.stage','M.stage'))
ggplot(plot.data,aes(x=value,fill=MUC16))+geom_bar(width = 0.75,position = position_fill(reverse = TRUE))+
  facet_grid(~variable,space='free',scales = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c(
    'grey','red'
  ))
c.data = clin.d[,c('Age','MUC16')]
library(reshape2)
c.data = melt(c.data,measure.vars = c('Age'))
ggplot(c.data,aes(x=MUC16,y=value,fill=MUC16))+geom_boxplot()+facet_grid(~variable,scales = 'free')+
  stat_compare_means()+theme_bw()+scale_fill_manual(values = pal_npg('nrc')(3))+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8))

####TERT_promoter relevant #### Figure 1d

for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage')){
 
  print(c(feature,check(clin.d$TERT_promoter,clin.d[,feature])))
}
plot.data = melt(clin.d,measure.vars = c('RecurRisk','Recur','LNM.3cm',
                                         'exTInvasion','LNM','exNInvasion'))

plot.data$TERT_promoter = factor(plot.data$TERT_promoter,levels = c('','TERT C228T'))
ggplot(plot.data,aes(x=value,fill=TERT_promoter))+geom_bar(width = 0.75,position = position_fill(reverse = TRUE))+
  facet_grid(~variable,space='free',scales = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
 scale_fill_manual(values = c('grey','blue'))

c.data = clin.d[,c('Age','TumorSize','LNM.No','TERT_promoter')]
library(reshape2)
c.data = melt(c.data,measure.vars = c('Age','TumorSize','LNM.No'))
ggplot(c.data,aes(x=TERT_promoter,y=value,fill=TERT_promoter))+geom_boxplot()+facet_grid(~variable,scales = 'free')+
  stat_compare_means()+theme_bw()+scale_fill_manual(values = pal_npg('nrc')(3))

####CCDC^_RET relevant####
clin.d$CCDC6_RET=factor(clin.d$CCDC6_RET,levels = c(0,1))
for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage')){
  print(c(feature,check(clin.d$CCDC6_RET,clin.d[,feature])))
}
plot.data = melt(clin.d,measure.vars = c('LNM.3cm','M.stage'))

ggplot(plot.data,aes(x=value,fill=CCDC6_RET))+geom_bar(width = 0.75,position = position_fill(reverse = TRUE))+
  facet_grid(~variable,space='free',scales = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c('grey','#8491B4FF'))

####BRAF clin ####

BRAF.mut.samples = as.character(subsetMaf(var.annovar.maf,genes = c('BRAF'))@data$Tumor_Sample_Barcode)
clin.d$BRAF = ifelse(clin.d$Tumor_Sample_Barcode %in% BRAF.mut.samples,'Mut',"WT")
for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter')){
  print(c(feature,fisher.test(clin.d$BRAF,clin.d[,feature])$p.value))
}
plot.data = melt(clin.d,measure.vars = c('TERT_promoter'))

ggplot(plot.data,aes(x=value,fill=BRAF))+geom_bar(width = 0.75,position = position_fill(reverse = TRUE))+
  facet_grid(~variable,space='free',scales = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))+scale_fill_manual(values = pal_npg('nrc')(3))

####BRAF and TERT
tert.samples = clin.d[clin.d$TERT_promoter!='','Tumor_Sample_Barcode']
tert.braf.samples = union(tert.samples,BRAF.mut.samples)
clin.d$TERT_BRAF=ifelse(clin.d$Tumor_Sample_Barcode %in% BRAF.mut.samples,'Both',"No")
for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter')){
  print(c(feature,fisher.test(clin.d$TERT_BRAF,clin.d[,feature])$p.value))
}
survdiff(Surv(days_RFS, Recur) ~ TERT_promoter, data = clin.d,rho = 0)
####CNV in NRAS RET NCOA4, et al.####


#item.samples = sfc[sfc$Gene1A == item | sfc$Gene1B == item,'Sample']
item.samples = as.character(subsetMaf(ptc.plus.gistic,genes = c('AK2'))@data$Tumor_Sample_Barcode)
item.samples = item.samples[item.samples %in% c('S19_613') == F]
clin.d$CNV = ifelse(clin.d$Tumor_Sample_Barcode %in% item.samples,'del','No')
clin.d.short = clin.d[clin.d$Tumor_Sample_Barcode %in% mut.all.samples,]
for(feature in c('Gender','RecurRisk','Recur','exTInvasion','LNM','exNInvasion','LNM.3cm','T.stage',
                 'N.stage','M.stage','TNM.stage','TERT_promoter')){
  print(c(feature,fisher.test(clin.d.short$CNV,clin.d.short[,feature])$p.value))
}
####

plot.data = melt(clin.d.short,measure.vars = c('exTInvasion','exNInvasion','M.stage',
                                         'T.stage','N.stage'))
plot.data$CNV = factor(plot.data$CNV,levels = c('del','No'))
ggplot(plot.data,aes(x=value,fill=CNV))+geom_bar(width = 0.75,position = position_fill(reverse = TRUE))+
  facet_grid(~variable,space='free',scales = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c(
    'Not'='grey','Up'='red','Down'='blue'
  ))+scale_fill_manual(values = c('blue','grey'))

c.data = clin.d.short[,c('Age','TumorSize','LNM.No','CNV')]
library(reshape2)
c.data = melt(c.data,measure.vars = c('Age','TumorSize','LNM.No'))
ggplot(c.data,aes(x=CNV,y=value,fill=CNV))+geom_boxplot()+facet_grid(~variable,scales = 'free')+
  stat_compare_means()+theme_bw()+scale_fill_manual(values =c('blue','grey'))
####Compare S19_126




####################Find mutation relevant molecules by multi-omics
mut.all.samples = unique(as.character(var.annovar.maf@data$Tumor_Sample_Barcode))#97 samples
findRelItems = function(gene,omics.mat,mut.all.samples){
  colnames(omics.mat)= unlist(lapply(colnames(omics.mat),function(a){
    a=gsub('-','_',a)
    if(substr(a,1,1) %in% as.character(1:9)){
      a=paste0('S',a)
    }
    if(a=='S20_11'){
      a='S20__11'
    }
    return(a)
  }
                                     ))
  cols = intersect(mut.all.samples,colnames(omics.mat)) 
  if(gene == 'TERT_promoter'){
    mut.samples =as.character( clin.d[clin.d$TERT_promoter != '','Tumor_Sample_Barcode'])
  }else{
    mut.samples = as.character(subsetMaf(var.annovar.maf,genes = gene)@data$Tumor_Sample_Barcode)
    
  }
  
  cols.mut = intersect(cols,mut.samples)
  cols.wt = cols[cols %in% mut.samples == F]
  
  p.values <- lapply(1:nrow(omics.mat), function(i){
    wilcox.test(as.double(omics.mat[i,cols.mut]),as.double(omics.mat[i,cols.wt]),paired = F)$p.value
  })
  
  log2FCs <- lapply(1:nrow(omics.mat), function(i){
    mean(as.double(omics.mat[i,cols.mut]))-mean(as.double(omics.mat[i,cols.wt]))
  })
  
  res.diff <- data.frame(gene = rownames(omics.mat),
                      p.value = as.double(unlist(p.values)),
                      FDR = p.adjust(as.double(unlist(p.values)),method = 'fdr'),
                      Log2FC = as.double(unlist(log2FCs)))
   
  return(res=list(cols.mut,cols.wt,res.diff))
  
}

seleItemFromPathSig = function(pathway,gsea.braf,BRAF.gene.rel,BRAF.pro.rel,BRAF.ph.rel){
  #pathway = 'Oxidative phosphorylation'
  path.res = gsea.braf[gsea.braf$ID == pathway,]
  genes.sel = c()
  BRAF.gene.rel = arrange(BRAF.gene.rel,p.value)
  BRAF.pro.rel = arrange(BRAF.pro.rel,p.value)
  BRAF.ph.rel = arrange(BRAF.ph.rel,p.value)
  
  for(i in 1:nrow(path.res)){
    genes.1 = strsplit(path.res$core_enrichment[i],'/')[[1]]
    if(path.res$From[i] == 'mrna'){
      genes.1 = BRAF.gene.rel[BRAF.gene.rel$gene %in% genes.1,]
    }else if(path.res$From[i] == 'protein'){
      genes.1 = BRAF.pro.rel[BRAF.pro.rel$gene %in% genes.1,]
    }else{
      genes.1 = BRAF.ph.rel[BRAF.ph.rel$gene %in% genes.1,]
    }
    
    genes.1.sel = genes.1[genes.1$p.value<0.05,'gene']
    if(length(genes.1.sel)>10){
      genes.1.sel = genes.1.sel[1:10]
    }
    if(length(genes.1.sel)<5){
      genes.1.sel = strsplit(path.res$core_enrichment[i],'/')[[1]][1:5]
    }
    genes.sel = union(genes.sel,genes.1.sel)
  }
  return(genes.sel)
}

data.meta = readRDS(file = './variables/normalizedMetabolism.Rds')##The names were changed into the same as mrna-count
data.mrna = readRDS(file = './variables/deseqNorm_count.Rds')
data.mrna = log2(data.mrna+1)
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

BRAF.pro.rel.res = findRelItems('BRAF',pro.exp,mut.all.samples)
BRAF.pro.rel = arrange(BRAF.pro.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = BRAF.pro.rel$Log2FC
names(FC.gsea)=as.character(BRAF.pro.rel$gene)

diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.pro = diff.pathway@result
BRAF.pro.sig = BRAF.pro.rel[BRAF.pro.rel$p.value<0.05,'gene']
BRAF.pro.sig.paths.res = enricher(BRAF.pro.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                  TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result

BRAF.ph.rel = findRelItems('BRAF',pro2.exp,mut.all.samples)
BRAF.ph.rel = arrange(BRAF.ph.rel[[3]],desc(Log2FC),p.value)
FC.gsea = BRAF.ph.rel$Log2FC
names(FC.gsea)=as.character(BRAF.ph.rel$gene)

diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.ph = diff.pathway@result
BRAF.ph.sig = BRAF.ph.rel[BRAF.ph.rel$p.value<0.05,'gene']
BRAF.ph.sig.paths.res = enricher(BRAF.ph.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                  TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result


BRAF.gene.rel = findRelItems('BRAF',data.mrna,mut.all.samples)
BRAF.gene.rel = arrange(BRAF.gene.rel[[3]],desc(Log2FC),p.value)
FC.gsea = BRAF.gene.rel$Log2FC
names(FC.gsea)=as.character(BRAF.gene.rel$gene)
diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.gene = diff.pathway@result

BRAF.gene.sig = BRAF.gene.rel[BRAF.gene.rel$p.value<0.05,'gene']
BRAF.gene.sig.paths.res = enricher(BRAF.gene.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                  TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result

BRAF.meta.rel = findRelItems('BRAF',data.meta,mut.all.samples)
BRAF.meta.rel = arrange(BRAF.meta.rel[[3]],desc(Log2FC),p.value)
BRAF.meta.sig = BRAF.meta.rel[BRAF.meta.rel$p.value<0.05,'gene']
BRAF.meta.sig.keggids = meta.info[BRAF.meta.sig,'KEGG2']
BRAF.meta.sig.paths.res = enricher(BRAF.meta.sig.keggids,minGSSize = 2,maxGSSize = 500,pvalueCutoff = 1.2,
                                   TERM2GENE = path2Meta[,c(1,2)],pAdjustMethod = 'fdr')@result

gsea.res.gene$From = 'mrna'
gsea.res.pro$From = 'protein'
gsea.res.ph$From = 'ph'

gsea.braf = rbind(gsea.res.gene,gsea.res.pro,gsea.res.ph)
disease.pathways = path2cate[path2cate$cate=='6. Human Diseases','Name']
gsea.braf = gsea.braf[gsea.braf$ID %in% disease.pathways == F & gsea.braf$ID != 'Metabolic pathways',]
gsea.braf = gsea.braf[gsea.braf$pvalue<0.01,]
gsea.braf.paths.dup = unique(gsea.braf[duplicated(gsea.braf$ID),'ID'])
gsea.dup.braf = gsea.braf[gsea.braf$ID %in% gsea.braf.paths.dup,]

MUC16.pro.rel.res = findRelItems('MUC16',pro.exp,mut.all.samples)
MUC16.pro.rel = arrange(MUC16.pro.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = MUC16.pro.rel$Log2FC
names(FC.gsea)=as.character(MUC16.pro.rel$gene)

diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gseaplot(diff.pathway,'Oxidative phosphorylation')
gsea.res.pro.muc16 = diff.pathway@result
MUC16.pro.sig = MUC16.pro.rel[MUC16.pro.rel$p.value<0.05,'gene']
MUC16.pro.sig.paths.res = enricher(MUC16.pro.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                 TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result


MUC16.ph.rel.res = findRelItems('MUC16',pro2.exp,mut.all.samples)
MUC16.ph.rel = arrange(MUC16.ph.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = MUC16.ph.rel$Log2FC
names(FC.gsea)=as.character(MUC16.ph.rel$gene)

diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.ph.muc16 = diff.pathway@result
MUC16.ph.sig = MUC16.ph.rel[MUC16.ph.rel$p.value<0.05,'gene']
MUC16.ph.sig.paths.res = enricher(MUC16.ph.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                   TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result


MUC16.gene.rel.res = findRelItems('MUC16',data.mrna,mut.all.samples)
MUC16.gene.rel = arrange(MUC16.gene.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = MUC16.gene.rel$Log2FC
names(FC.gsea)=as.character(MUC16.gene.rel$gene)

diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.gene.muc16 = diff.pathway@result
MUC16.gene.sig = MUC16.gene.rel[MUC16.gene.rel$p.value<0.05,'gene']
MUC16.gene.sig.paths.res = enricher(MUC16.gene.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                   TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result


MUC16.meta.rel.res = findRelItems('MUC16',data.meta,mut.all.samples)
MUC16.meta.rel = arrange(MUC16.meta.rel.res[[3]],desc(Log2FC),p.value)
MUC16.meta.sig = MUC16.meta.rel[MUC16.meta.rel$p.value<0.05,'gene']
MUC16.meta.sig.keggids = unique(meta.info[MUC16.meta.sig,'KEGG2'])
MUC16.meta.sig.paths.res = enricher(MUC16.meta.sig.keggids,minGSSize = 2,maxGSSize = 500,pvalueCutoff = 1.2,
                                   TERM2GENE = path2Meta[,c(1,2)],pAdjustMethod = 'fdr')@result

gsea.res.gene.muc16$From = 'mrna'
gsea.res.pro.muc16$From = 'protein'
gsea.res.ph.muc16$From = 'ph'

gsea.muc16 = rbind(gsea.res.gene.muc16,gsea.res.pro.muc16,gsea.res.ph.muc16)

gsea.muc16 = gsea.muc16[gsea.muc16$ID %in% disease.pathways == F & gsea.muc16$ID != 'Metabolic pathways',]
gsea.muc16 = gsea.muc16[gsea.muc16$pvalue<0.01,]
gsea.muc16.paths.dup = unique(gsea.muc16[duplicated(gsea.muc16$ID),'ID'])
gsea.dup.muc16 = gsea.muc16[gsea.muc16$ID %in% gsea.muc16.paths.dup,]


TERT.pro.rel.res = findRelItems('TERT_promoter',pro.exp,mut.all.samples)
TERT.pro.rel = arrange(TERT.pro.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = TERT.pro.rel$Log2FC
names(FC.gsea)=as.character(TERT.pro.rel$gene)
diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.pro.tert = diff.pathway@result
TERT.pro.sig = TERT.pro.rel[TERT.pro.rel$p.value<0.05,'gene']
TERT.pro.sig.paths.res = enricher(TERT.pro.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                   TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result

TERT.ph.rel.res = findRelItems('TERT_promoter',pro2.exp,mut.all.samples)
TERT.ph.rel = arrange(TERT.ph.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = TERT.ph.rel$Log2FC
names(FC.gsea)=as.character(TERT.ph.rel$gene)
diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
gsea.res.ph.tert = diff.pathway@result
TERT.ph.sig = TERT.ph.rel[TERT.ph.rel$p.value<0.05,'gene']
TERT.ph.sig.paths.res = enricher(TERT.ph.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                  TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result

TERT.gene.rel.res = findRelItems('TERT_promoter',data.mrna,mut.all.samples)
TERT.gene.rel = arrange(TERT.gene.rel.res[[3]],desc(Log2FC),p.value)
FC.gsea = TERT.gene.rel$Log2FC
names(FC.gsea)=as.character(TERT.gene.rel$gene)
diff.pathway = GSEA(FC.gsea,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                    TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr',
                    seed=1234)
TERT.gene.sig = TERT.gene.rel[TERT.gene.rel$p.value<0.05,'gene']
TERT.gene.sig.paths.res = enricher(TERT.gene.sig,minGSSize = 5,maxGSSize = 500,pvalueCutoff = 1.2,
                                  TERM2GENE = path2Gene[,c(1,2)],pAdjustMethod = 'fdr')@result

TERT.meta.rel.res = findRelItems('TERT_promoter',data.meta,mut.all.samples)
TERT.meta.rel = arrange(TERT.meta.rel.res[[3]],desc(Log2FC),p.value)
TERT.meta.sig =TERT.meta.rel[TERT.meta.rel$p.value<0.05,'gene']
TERT.meta.sig.keggids = unique(meta.info[TERT.meta.sig,'KEGG2'])
TERT.meta.sig.paths.res = enricher(TERT.meta.sig.keggids,minGSSize = 2,maxGSSize = 500,pvalueCutoff = 1.2,
                                    TERM2GENE = path2Meta[,c(1,2)],pAdjustMethod = 'fdr')@result

gsea.res.gene.tert = diff.pathway@result
gsea.res.gene.tert$From = 'mrna'
gsea.res.pro.tert$From = 'protein'
gsea.res.ph.tert$From = 'ph'

gsea.tert = rbind(gsea.res.gene.tert,gsea.res.pro.tert,gsea.res.ph.tert)

gsea.tert = gsea.tert[gsea.tert$ID %in% disease.pathways == F & gsea.tert$ID != 'Metabolic pathways',]
gsea.tert = gsea.tert[gsea.tert$pvalue<0.01,]
gsea.tert.paths.dup = unique(gsea.tert[duplicated(gsea.tert$ID),'ID'])
gsea.dup.tert = gsea.tert[gsea.tert$ID %in% gsea.tert.paths.dup,]


#####Plot pathway res

BRAF.pro.sig.paths.res$From = 'protein'
BRAF.ph.sig.paths.res$From = 'phospho'
BRAF.gene.sig.paths.res$From = 'mRNA'
BRAF.meta.sig.paths.res$From = 'metabolite'
disease.pathways = as.character(path2cate[path2cate$cate == '6. Human Diseases','Name'])
BRAF.meta.sig.paths.res = BRAF.meta.sig.paths.res[BRAF.meta.sig.paths.res$ID %in% metabolism.pathways,]

BRAF.path.enrich = rbind(BRAF.pro.sig.paths.res,BRAF.ph.sig.paths.res,BRAF.gene.sig.paths.res,BRAF.meta.sig.paths.res)
BRAF.path.enrich.top = arrange(BRAF.path.enrich,From,pvalue) %>% filter(ID %in% metabolism.pathways & pvalue <0.05 & Count >1) %>% group_by (From) %>% top_n(-10,pvalue)
BRAF.path.enrich.top$ID = factor(BRAF.path.enrich.top$ID,levels = rev(as.character(unique(BRAF.path.enrich.top$ID))))
BRAF.path.enrich.top$From = factor(BRAF.path.enrich.top$From,levels = c('metabolite','mRNA','protein','phospho'))



MUC16.pro.sig.paths.res$From = 'protein'
MUC16.ph.sig.paths.res$From = 'phospho'
MUC16.gene.sig.paths.res$From = 'mRNA'
MUC16.meta.sig.paths.res$From = 'metabolite'
MUC16.meta.sig.paths.res = MUC16.meta.sig.paths.res[MUC16.meta.sig.paths.res$ID %in% metabolism.pathways,]
MUC16.path.enrich = rbind(MUC16.pro.sig.paths.res,MUC16.ph.sig.paths.res,MUC16.gene.sig.paths.res,MUC16.meta.sig.paths.res)
MUC16.path.enrich.top = arrange(MUC16.path.enrich,From,pvalue) %>% filter(ID %in% metabolism.pathways & pvalue <0.05 & Count >1) %>% group_by (From) %>% top_n(-10,pvalue)
MUC16.path.enrich.top$ID = factor(MUC16.path.enrich.top$ID,levels = rev(as.character(unique(MUC16.path.enrich.top$ID))))
MUC16.path.enrich.top$From = factor(MUC16.path.enrich.top$From,levels = c('metabolite','mRNA','protein','phospho'))


TERT.pro.sig.paths.res$From = 'protein'
TERT.ph.sig.paths.res$From = 'phospho'
TERT.gene.sig.paths.res$From = 'mRNA'
TERT.meta.sig.paths.res$From = 'metabolite'
TERT.meta.sig.paths.res = TERT.meta.sig.paths.res[TERT.meta.sig.paths.res$ID %in% metabolism.pathways,]
TERT.path.enrich = rbind(TERT.pro.sig.paths.res,TERT.ph.sig.paths.res,TERT.gene.sig.paths.res,TERT.meta.sig.paths.res)
TERT.path.enrich.top = arrange(TERT.path.enrich,From,pvalue) %>% filter(ID %in% metabolism.pathways & pvalue <0.05 & Count >1) %>% group_by (From) %>% top_n(-10,pvalue)
TERT.path.enrich.top$ID = factor(TERT.path.enrich.top$ID,levels = rev(as.character(unique(TERT.path.enrich.top$ID))))
TERT.path.enrich.top$From = factor(TERT.path.enrich.top$From,levels = c('metabolite','mRNA','protein','phospho'))
BRAF.path.enrich.top$Mutation = 'BRAF'
MUC16.path.enrich.top$Mutation = 'MUC16'
TERT.path.enrich.top$Mutation = 'TERT promoter'
mut.path.m = rbind(BRAF.path.enrich.top,MUC16.path.enrich.top,TERT.path.enrich.top)

ggplot(mut.path.m,aes(y=ID,x=From,color=Count,size=-log10(pvalue)))+
  geom_point()+scale_color_gradientn(colours = hcl.colors(100,"Blues 3"))+
  facet_wrap(~Mutation,ncol = 3,scales = 'free')+scale_size_continuous(range = c(1,3))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(size=6),
        strip.background = element_blank(),
        strip.text = element_text(size=8),
        legend.direction = 'vertical',
        legend.position = 'right',
        legend.key.size = unit(8, "pt"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

BRAF.path.enrich$Mutation = 'BRAF'
MUC16.path.enrich$Mutation = 'MUC16'
TERT.path.enrich$Mutation = 'TERT'

path.m = rbind(BRAF.path.enrich,MUC16.path.enrich,TERT.path.enrich)
path.m = path.m[path.m$pvalue<0.05,]
write.csv(path.m,file = './variables/three mutation relevant pathways.csv')
####Print the top five items

findTopItems = function(BRAF.gene.rel,BRAF.pro.rel,BRAF.ph.rel){
  
  BRAF.gene.rel = arrange(BRAF.gene.rel,p.value)
  BRAF.pro.rel = arrange(BRAF.pro.rel,p.value)
  BRAF.ph.rel = arrange(BRAF.ph.rel,p.value)
  
  top.items = unique(c(BRAF.gene.rel$gene[1:10],
                    BRAF.pro.rel$gene[1:10],
                    BRAF.ph.rel$gene[1:10]))
  
  BRAF.gene.rel$From = 'mrna'
  BRAF.pro.rel$From = 'protein'
  BRAF.ph.rel$From = 'phosphprotein'
  rownames(BRAF.gene.rel)=BRAF.gene.rel$gene
  rownames(BRAF.pro.rel)=BRAF.pro.rel$gene
  rownames(BRAF.ph.rel)=BRAF.ph.rel$gene
  res.m = rbind(BRAF.gene.rel[intersect(rownames(BRAF.gene.rel),top.items),],
                BRAF.pro.rel[intersect(rownames(BRAF.pro.rel),top.items),],
                BRAF.ph.rel[intersect(rownames(BRAF.ph.rel),top.items),])
  
  
  res.m$gene = factor(res.m$gene,levels=top.items)
  return(res.m)
}

####BRAF\MUC16\TERT relevant gene protein or phosph items#### Supplementary  Figure 
res.m = findTopItems(BRAF.gene.rel,BRAF.pro.rel,BRAF.ph.rel)
ggplot(res.m,aes(x=gene,y=From,color=Log2FC,size=-log10(p.value)+0.1))+
  geom_point()+scale_color_gradient2(low = 'blue',high = 'red')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=90),
        legend.direction= 'horizontal',
        legend.position = 'bottom')
res.m2 = findTopItems(MUC16.gene.rel,MUC16.pro.rel,MUC16.ph.rel)
ggplot(res.m2,aes(x=gene,y=From,color=Log2FC,size=-log10(p.value)+0.1))+
  geom_point()+scale_color_gradient2(low = 'blue',high = 'red')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=90),
        legend.direction= 'horizontal',
        legend.position = 'bottom')

res.m3= findTopItems(TERT.gene.rel,TERT.pro.rel,TERT.ph.rel)
ggplot(res.m3,aes(x=gene,y=From,color=Log2FC,size=-log10(p.value)+0.1))+
  geom_point()+scale_color_gradient2(low = 'blue',high = 'red')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=90),
        legend.direction= 'horizontal',
        legend.position = 'bottom')


res.m$Mutation = 'BRAF'
res.m2$Mutation = 'MUC16'
res.m3$Mutation = 'TERT promoter'
res.all = rbind(res.m,res.m2,res.m3)
res.all$Log2FC[res.all$Log2FC>2]=2
res.all$Log2FC[res.all$Log2FC< -2]= -2

ggplot(res.all,aes(x=gene,y=From,color=Log2FC,size=-log10(p.value)+0.1))+
  geom_point()+scale_color_gradientn(colours = hcl.colors(100,"Blue-Red 2"))+
  facet_wrap(~Mutation,ncol = 1,scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=90),
        axis.text = element_text(size=8),
        strip.background = element_blank(),
        legend.direction= 'horizontal',
        legend.position = 'bottom',
        legend.key.size = unit(8, "pt"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))



#################metabolism relevances#######
BRAF.meta.rel$Mut = 'BRAF'
MUC16.meta.rel$Mut = 'MUC16'
TERT.meta.rel$Mut='TERT_promoter'

meta.rel.muts = rbind(BRAF.meta.rel,MUC16.meta.rel,TERT.meta.rel)
meta.rel.muts.top = arrange(meta.rel.muts,Mut,p.value) %>% group_by (Mut) %>% top_n(10,abs(Log2FC))
metas = unique(meta.rel.muts.top[meta.rel.muts.top$p.value<0.05,'gene'])

res.m4 = meta.rel.muts[meta.rel.muts$gene %in% as.character(metas$gene),]
res.m4$gene = factor(res.m4$gene,levels = rev(as.character(metas$gene)))
# ggplot(res.m4,aes(y=gene,fill=Log2FC,x=-log10(p.value)))+
#   facet_grid(~Mut)+
#   geom_bar(stat = 'identity',position = 'dodge',width = 0.6)+
#   scale_fill_gradient2(low = 'blue',high = 'red')+
#   geom_vline(xintercept = 2,lty=3,color='red')+
#   theme_cowplot()+
#   theme(axis.text = element_text(size=8),
#         strip.background = element_blank())

ggplot(res.m4,aes(y=gene,x=Mut,color=Log2FC,size=-log10(p.value)+0.1))+
  geom_point()+scale_color_gradientn(colours = hcl.colors(100,"Blue-Red 2"))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(size=8),
        strip.background = element_blank(),
        legend.direction = 'vertical',
        legend.position = 'right',
        legend.key.size = unit(8, "pt"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))


prepare_violin = function(mols.res,N=10){
  
  mols.res$type = 'NS'
  
  mols.res[mols.res$pvalue<0.05 & mols.res$log2FC > 0,'type']='Up'
  mols.res[mols.res$pvalue<0.05 & mols.res$log2FC < 0, 'type']= 'Down'
  
  mols.res$DisplayName = ''
  
  mols.res[order(mols.res$log2FC)[1:N],'DisplayName']=mols.res[order(mols.res$log2FC)[1:N],'name']
  mols.res[order(mols.res$log2FC,decreasing = T)[1:N],'DisplayName']=mols.res[order(mols.res$log2FC,decreasing = T)[1:N],'name']
  
  mols.res[mols.res$type == 'NS','DisplayName']=''
  return(mols.res)
  
}
colnames(MUC16.meta.rel)=c('name','pvalue','FDR','log2FC','Mut')
ggplot(prepare_violin(MUC16.meta.rel),aes(x=log2FC,y=-log10(pvalue),color = type))+
  geom_point(size=2)+
  geom_text_repel(aes(label=DisplayName),size=2)+
  geom_hline(yintercept =-log10(0.05),color='red',lty=3)+
  theme_bw()+
  scale_color_manual(values = c('blue','grey','red'))

