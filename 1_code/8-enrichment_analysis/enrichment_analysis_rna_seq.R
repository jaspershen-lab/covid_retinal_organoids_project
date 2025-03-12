library(r4projects)
library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)  
library(ggplot2)
library(stringr)
library(enrichplot) 

setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

# read data
load("3_data_analysis/2-data_cleaning/1-transcriptome/transcriptome_data.RData")


dir.create("3_data_analysis/8-enrichment_analysis/rna_seq", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/8-enrichment_analysis/rna_seq")


# 读取cluster基因
read_cluster_genes <- function(cluster_num) {
  project_wd <- get_project_wd()
  file_path <- file.path(project_wd,
                         "3_data_analysis/7-heatmap_analysis/rna_seq/clustering_results",  
                         paste0("cluster_", cluster_num, "_genes.csv"))
  
  if(!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  genes_df <- read.csv(file_path)
  return(genes_df)
}



############### Enrichment Analysis for Each Cluster ###############
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# Create function to perform enrichment analysis
perform_enrichment_analysis <- function(cluster_num) {
  # Read cluster genes
  cluster_genes <- read_cluster_genes(cluster_num)
  
  # Get ENTREZ IDs (remove NA)
  entrez_ids <- na.omit(cluster_genes$ENTREZID)
  
  # Create directory for results
  dir.create(paste0("cluster", cluster_num), showWarnings = FALSE)
  
  # 1. GO Enrichment Analysis
  ego_bp <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
  
  ego_mf <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
  
  ego_cc <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
  
  # 2. KEGG Pathway Analysis
  ekegg <- enrichKEGG(gene = entrez_ids,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
  
  # 3. Reactome Pathway Analysis
  ereactome <- enrichPathway(gene = entrez_ids,
                             organism = "human",
                             pvalueCutoff = 0.05)
  
  # Save results
  ## Save tables
  if(nrow(as.data.frame(ego_bp)) > 0) {
    write.xlsx(as.data.frame(ego_bp), 
               paste0("cluster", cluster_num, "/GO_BP_enrichment.xlsx"))
  }
  
  if(nrow(as.data.frame(ego_mf)) > 0) {
    write.xlsx(as.data.frame(ego_mf), 
               paste0("cluster", cluster_num, "/GO_MF_enrichment.xlsx"))
  }
  
  if(nrow(as.data.frame(ego_cc)) > 0) {
    write.xlsx(as.data.frame(ego_cc), 
               paste0("cluster", cluster_num, "/GO_CC_enrichment.xlsx"))
  }
  
  if(nrow(as.data.frame(ekegg)) > 0) {
    write.xlsx(as.data.frame(ekegg), 
               paste0("cluster", cluster_num, "/KEGG_enrichment.xlsx"))
  }
  
  if(nrow(as.data.frame(ereactome)) > 0) {
    write.xlsx(as.data.frame(ereactome), 
               paste0("cluster", cluster_num, "/Reactome_enrichment.xlsx"))
  }
  
  ## Save plots
  # GO plots
  if(nrow(as.data.frame(ego_bp)) > 0) {
    pdf(paste0("cluster", cluster_num, "/GO_BP_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(ego_bp, showCategory = 20))
    dev.off()
  }
  
  if(nrow(as.data.frame(ego_mf)) > 0) {
    pdf(paste0("cluster", cluster_num, "/GO_MF_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(ego_mf, showCategory = 20))
    dev.off()
  }
  
  if(nrow(as.data.frame(ego_cc)) > 0) {
    pdf(paste0("cluster", cluster_num, "/GO_CC_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(ego_cc, showCategory = 20))
    dev.off()
  }
  
  # KEGG plot
  if(nrow(as.data.frame(ekegg)) > 0) {
    pdf(paste0("cluster", cluster_num, "/KEGG_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(ekegg, showCategory = 20))
    dev.off()
  }
  
  # Reactome plot
  if(nrow(as.data.frame(ereactome)) > 0) {
    pdf(paste0("cluster", cluster_num, "/Reactome_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(ereactome, showCategory = 20))
    dev.off()
  }
  
  return(list(
    GO_BP = ego_bp,
    GO_MF = ego_mf,
    GO_CC = ego_cc,
    KEGG = ekegg,
    Reactome = ereactome
  ))
}

# Perform analysis for each cluster
for(i in 1:6) {
  cat(paste("\nProcessing Cluster", i, "...\n"))
  results <- perform_enrichment_analysis(i)
  cat("Done!\n")
}





