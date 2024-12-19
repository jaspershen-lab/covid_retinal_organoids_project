library(r4projects)
library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)  
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

# read data
load("3_data_analysis/2-data_cleaning/1-transcriptome/transcriptome_data.RData")


dir.create("3_data_analysis/8-enrichment_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/8-enrichment_analysis")


# 读取cluster基因
read_cluster_genes <- function(cluster_num) {
  project_wd <- get_project_wd()
  file_path <- file.path(project_wd,
                         "3_data_analysis/7-heatmap_analysis", 
                         paste0("cluster_", cluster_num), 
                         paste0("cluster", cluster_num, ".xlsx"))
  
  if(!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  genes_df <- readxl::read_excel(file_path)
  return(genes_df)
}

# 创建基因ID映射
symbol_to_ensembl <- as.data.frame(transcriptome_data@variable_info) %>%
  dplyr::select(SYMBOL, ENSEMBL) %>%
  na.omit() %>%
  distinct()

# 结果目录
results_dir <- getwd()

for(cluster_num in 1:6) {
  cat(sprintf("\nProcessing cluster %d...\n", cluster_num))
  
  # 读取并映射基因
  genes_df <- read_cluster_genes(cluster_num)
  mapped_genes <- genes_df %>%
    left_join(symbol_to_ensembl, by = c("variable_id" = "SYMBOL"))
  
  cat(sprintf("Found %d genes with ENSEMBL IDs\n", sum(!is.na(mapped_genes$ENSEMBL))))
  
  # 创建cluster结果目录
  cluster_result_dir <- file.path(results_dir, paste0("cluster_", cluster_num))
  dir.create(cluster_result_dir, showWarnings = FALSE)
  
  # GO enrichment analysis
  ego <- enrichGO(gene = unique(mapped_genes$ENSEMBL),
                  universe = unique(transcriptome_data@variable_info$ENSEMBL),
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL",
                  ont = "ALL",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.1,
                  qvalueCutoff = 0.1)
  
  if(!is.null(ego)) {
    ego_df <- as.data.frame(ego)
    if(nrow(ego_df) > 0) {
      write.xlsx(ego_df,
                 file = file.path(cluster_result_dir, "GO_enrichment.xlsx"),
                 rowNames = FALSE)
      cat(sprintf("GO analysis: Found %d enriched terms\n", nrow(ego_df)))
    } else {
      cat("GO analysis: No significant enrichment found\n")
    }
  }
  
  # KEGG enrichment analysis
  # 先转换ENSEMBL ID到ENTREZ ID
  entrez_ids <- bitr(mapped_genes$ENSEMBL,
                     fromType = "ENSEMBL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  
  ekegg <- enrichKEGG(gene = unique(entrez_ids$ENTREZID),
                      organism = 'hsa',  # 人类的KEGG代码
                      pvalueCutoff = 0.1,
                      qvalueCutoff = 0.1)
  
  if(!is.null(ekegg)) {
    ekegg_df <- as.data.frame(ekegg)
    if(nrow(ekegg_df) > 0) {
      write.xlsx(ekegg_df,
                 file = file.path(cluster_result_dir, "KEGG_enrichment.xlsx"),
                 rowNames = FALSE)
      cat(sprintf("KEGG analysis: Found %d enriched terms\n", nrow(ekegg_df)))
    }
  }
  
  # Reactome enrichment analysis
  ereactome <- enrichPathway(gene = unique(entrez_ids$ENTREZID),
                             organism = "human",
                             pvalueCutoff = 0.1,
                             qvalueCutoff = 0.1,
                             readable = TRUE)  
  
  if(!is.null(ereactome)) {
    ereactome_df <- as.data.frame(ereactome)
    if(nrow(ereactome_df) > 0) {
      write.xlsx(ereactome_df,
                 file = file.path(cluster_result_dir, "Reactome_enrichment.xlsx"),
                 rowNames = FALSE)
      cat(sprintf("Reactome analysis: Found %d enriched terms\n", nrow(ereactome_df)))
    }
  }
  
  cat(sprintf("Completed enrichment analysis for cluster %d\n", cluster_num))
}



library(ggplot2)
library(dplyr)
library(stringr)
library(openxlsx)
library(enrichplot) 
library(clusterProfiler)

# neural related pathways
extract_neural_pathways <- function(enrichment_df) {
  neural_keywords <- c("neur", "brain", "axon", "synap", "dendrit", 
                       "myelin", "glial", "nerve", "cerebr", "hippocamp",
                       "cortex", "retina")
  
  neural_paths <- enrichment_df %>%
    filter(grepl(paste(neural_keywords, collapse = "|"), 
                 Description, ignore.case = TRUE)) %>%
    filter(p.adjust < 0.1) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15)
  
  # Convert to enrichResult object format
  neural_paths <- new("enrichResult",
                      result = neural_paths,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.1)
  
  return(neural_paths)
}

# Process the results of each cluster
for(cluster_num in 1:6) {
  # GO
  go_file <- file.path(getwd(),
                       paste0("cluster_", cluster_num), 
                       "GO_enrichment.xlsx")
  
  if(file.exists(go_file)) {
    go_results <- read.xlsx(go_file)
    neural_go <- extract_neural_pathways(go_results)
    
    if(length(neural_go@result$Description) > 0) {
      p_go <- barplot(neural_go)
      
      ggsave(
        filename = file.path(getwd(), 
                            
                             paste0("cluster_", cluster_num), 
                             "neural_GO_pathways_barplot.pdf"),
        plot = p_go,
        width = 12,
        height = 8
      )
    }
  }
  
  # KEGG
  kegg_file <- file.path(getwd(),
                         
                         paste0("cluster_", cluster_num), 
                         "KEGG_enrichment.xlsx")
  
  if(file.exists(kegg_file)) {
    kegg_results <- read.xlsx(kegg_file)
    neural_kegg <- extract_neural_pathways(kegg_results)
    
    if(length(neural_kegg@result$Description) > 0) {
      p_kegg <- barplot(neural_kegg)
      
      ggsave(
        filename = file.path(getwd(), 
                             
                             paste0("cluster_", cluster_num), 
                             "neural_KEGG_pathways_barplot.pdf"),
        plot = p_kegg,
        width = 12,
        height = 8
      )
    }
  }
  
  # Reactome
  reactome_file <- file.path(getwd(), 
                             
                             paste0("cluster_", cluster_num), 
                             "Reactome_enrichment.xlsx")
  
  if(file.exists(reactome_file)) {
    reactome_results <- read.xlsx(reactome_file)
    neural_reactome <- extract_neural_pathways(reactome_results)
    
    if(length(neural_reactome@result$Description) > 0) {
      p_reactome <- barplot(neural_reactome)
      
      ggsave(
        filename = file.path(getwd(), 
                             
                             paste0("cluster_", cluster_num), 
                             "neural_Reactome_pathways_barplot.pdf"),
        plot = p_reactome,
        width = 12,
        height = 8
      )
    }
  }
}












