# PTC_multi_omics
Multi omics integration analysis of papillary thyroid cancer

annoFusion.R: annotation of the fusions identified based on arriba and starFusion.
mutation.R: somatic mutations and gene fusions
mRNA.R: basic analysis about transcriptomics data
metabolism.R: basic analysis about metabolomics data
proteomics.R: basic analysis about proteomics and phosphorylation-proteomics data
mergeDiffResAmongOmics.R: find the differences between tumor and normal samples based on the transcriptomics.R, metabolomics and proteomics data.
recurRiskAnalysis.R: find the molecules and pathways correlated with PTC recurrence risks.
diablo_omics_correlations.R: supervised multi-omics correlation analysis
subtypesClustering.R:find the new subtypes based on clustering of omics data
findBiomarkerOrTargets.R: find the potential targets of the subtyes
differentBetweenTrans&Prot.R: compute the sample-wise and gene wise mRNA-protein correlations
DWLSAnalysis.R: integrate scRNA-seq with bulk RNA-seq
