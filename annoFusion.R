
library(annoFuse)
library(annoFuseData)
cases = read.delim(file = './RNASeq_cases.txt')

sf.merge = data.frame()
for(case in cases){
  
  fileStarFusion = read_starfusion_calls(paste0('./cleandata/',case,'SFOut/star-fusion.fusion_predictions.abridged.coding_effect.tsv'))
  fileStarFusion$Tumor_ID = case
  formattedStarFusion <- fusion_standardization(fileStarFusion,
                                                caller = "STARFUSION",
                                                tumorID = strsplit(case,'_')[[1]][1]
  )
  sf.merge = rbind(sf.merge,formattedStarFusion)
  print(case)
  
}
geneListReferenceDataTab <- read.delim(
  system.file("extdata", "genelistreference.txt", package = "annoFuseData"),
  stringsAsFactors = FALSE
)
fusionReferenceDataTab <- read.delim(
  system.file("extdata", "fusionreference.txt", package = "annoFuseData"),
  stringsAsFactors = FALSE
)
bioMartDataPfam <-
  readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData"))
kinaseid <- unique(bioMartDataPfam$pfam_id[grep("kinase", bioMartDataPfam$NAME)])
sf.merge <- annotate_fusion_calls(sf.merge,
                                  geneListReferenceDataTab = geneListReferenceDataTab,
                                  fusionReferenceDataTab = fusionReferenceDataTab,
                                  checkReciprocal = TRUE)
fusion_driver_df <- fusion_driver(sf.merge,
                                  annotated = T,
                                  geneListReferenceDataTab = geneListReferenceDataTab,
                                  fusionReferenceDataTab = fusionReferenceDataTab,
                                  checkDomainStatus = FALSE,
                                  domainsToCheck = kinaseid
)

plot_recurrent_fusions(fusion_driver_df,
                       groupby = "Fusion_Type",
                       countID = "Sample",
                       plotn=30
)
plot_recurrent_genes(
  fusion_driver_df,
  groupby='Fusion_Type',
  plotn = 30,
  countID='Sample',
  palette_rec = NULL,
  base_size = 20
)
plot_summary(sf.merge, groupby = "Fusion_Type")
#reportFuse(out_annofuse = sf.merge)
saveRDS(sf.merge,file='./02.fusion/merge.Rds')


##############################20221020, arriba#######################################

sf.merge2 = data.frame()
for(case in cases){
  
  fileStarFusion = read_starfusion_calls(paste0('./02.fusion/',case,'_fusion.tsv'))
  fileStarFusion$Tumor_ID = case
  formattedStarFusion <- fusion_standardization(fileStarFusion,
                                                caller = "ARRIBA",
                                                tumorID = strsplit(case,'_')[[1]][1]
  )
  sf.merge2 = rbind(sf.merge2,formattedStarFusion)
  print(case)
  
}


sf.merge2 <- annotate_fusion_calls(sf.merge2,
                                  geneListReferenceDataTab = geneListReferenceDataTab,
                                  fusionReferenceDataTab = fusionReferenceDataTab,
                                  checkReciprocal = TRUE)
fusion_driver_df2 <- fusion_driver(sf.merge2,
                                  annotated = T,
                                  geneListReferenceDataTab = geneListReferenceDataTab,
                                  fusionReferenceDataTab = fusionReferenceDataTab,
                                  checkDomainStatus = FALSE,
                                  domainsToCheck = kinaseid
)

plot_recurrent_fusions(sf.merge2,
                       groupby = "Fusion_Type",
                       countID = "Sample",
                       plotn=30
)
plot_recurrent_genes(
  fusion_driver_df2,
  groupby='Fusion_Type',
  plotn = 30,
  countID='Sample',
  palette_rec = NULL,
  base_size = 20
)
plot_summary(sf.merge, groupby = "Fusion_Type")
#reportFuse(out_annofuse = sf.merge)
saveRDS(sf.merge,file='./02.fusion/merge.Rds')
saveRDS(sf.merge2,file='./02.fusion/merge2.Rds')


##############combine two methods#################
standardFusioncalls = rbind(sf.merge,sf.merge2)
fusionQCFiltered <- fusion_filtering_QC(
  standardFusioncalls = standardFusioncalls,
  readingFrameFilter = "in-frame|frameshift|other",
  artifactFilter = "GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap|ConjoinG",
  junctionReadCountFilter = 3,
  spanningFragCountFilter = 0,
  readthroughFilter = TRUE
)


fusion_driver <- fusion_driver(fusionQCFiltered,
                                   annotated = T,
                                   geneListReferenceDataTab = geneListReferenceDataTab,
                                   fusionReferenceDataTab = fusionReferenceDataTab,
                                   checkDomainStatus = FALSE,
                                   domainsToCheck = kinaseid
)
sfc <- as.data.frame(
  fusionQCFiltered[which(fusionQCFiltered$Fusion_Type == "in-frame" & fusionQCFiltered$BreakpointLocation == "Genic"), ]
)
plot_recurrent_fusions(sfc,
                       groupby = "Fusion_Type",
                       countID = "Sample",
                       plotn=30
)+theme_publication(base_size = 8)
plot_recurrent_genes(
  sfc,
  groupby='Fusion_Type',
  plotn = 30,
  countID='Sample',
  palette_rec = NULL,
  base_size = 8
)
plot_summary(sf.merge, groupby = "Fusion_Type")

write.csv(sfc,file = './variables/fusion_sfc.csv')
