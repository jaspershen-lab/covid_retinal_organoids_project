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
                         "3_data_analysis/7-heatmap_analysis",  
                         paste0("cluster_", cluster_num, "_genes.csv"))
  
  if(!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  genes_df <- read.csv(file_path)
  return(genes_df)
}



######enrichement analysis
run_enrichment_analysis <- function() {
  # 创建基因ID映射
  symbol_to_ensembl <- as.data.frame(transcriptome_data@variable_info) %>%
    dplyr::select(SYMBOL, ENSEMBL) %>%
    na.omit() %>%
    distinct()
  
  # 结果目录
  results_dir <- getwd()
  
  for(cluster_num in 1:6) {
    cat(sprintf("\nProcessing cluster %d...\n", cluster_num))
    
    # 读取并处理基因数据
    genes_df <- read_cluster_genes(cluster_num)
    cat("Number of genes in cluster:", nrow(genes_df), "\n")
    
    # 将gene_id映射到ENSEMBL ID
    mapped_genes <- genes_df %>%
      left_join(symbol_to_ensembl, by = c("gene_id" = "SYMBOL"))
    
    # 检查映射结果
    valid_ensembl <- unique(na.omit(mapped_genes$ENSEMBL))
    cat(sprintf("Found %d unique genes with ENSEMBL IDs\n", length(valid_ensembl)))
    
    if(length(valid_ensembl) == 0) {
      cat("No valid ENSEMBL IDs found for enrichment analysis\n")
      next
    }
    
    # 创建cluster结果目录
    cluster_result_dir <- file.path(results_dir, paste0("cluster_", cluster_num))
    dir.create(cluster_result_dir, showWarnings = FALSE)
    
    # GO enrichment analysis
    ego <- enrichGO(gene = valid_ensembl,
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
        
        # 画GO富集分析图
        library(enrichplot)
        p1 <- barplot(ego, showCategory=20)
        ggsave(file.path(cluster_result_dir, "GO_enrichment_barplot.pdf"), p1, width = 10, height = 8)
      } else {
        cat("GO analysis: No significant enrichment found\n")
      }
    }
    
    # KEGG enrichment analysis
    entrez_ids <- bitr(valid_ensembl,
                       fromType = "ENSEMBL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    
    if(nrow(entrez_ids) > 0) {
      ekegg <- enrichKEGG(gene = unique(entrez_ids$ENTREZID),
                          organism = 'hsa',
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1)
      
      if(!is.null(ekegg)) {
        ekegg_df <- as.data.frame(ekegg)
        if(nrow(ekegg_df) > 0) {
          write.xlsx(ekegg_df,
                     file = file.path(cluster_result_dir, "KEGG_enrichment.xlsx"),
                     rowNames = FALSE)
          cat(sprintf("KEGG analysis: Found %d enriched terms\n", nrow(ekegg_df)))
          
          # 画KEGG富集分析图
          p2 <- barplot(ekegg, showCategory=20)
          ggsave(file.path(cluster_result_dir, "KEGG_enrichment_barplot.pdf"), p2, width = 10, height = 8)
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
          
          # 画Reactome富集分析图
          p3 <- barplot(ereactome, showCategory=20)
          ggsave(file.path(cluster_result_dir, "Reactome_enrichment_barplot.pdf"), p3, width = 10, height = 8)
        }
      }
    }
  }
}

# 运行分析
run_enrichment_analysis()




# neural related pathways
extract_neural_pathways <- function(enrichment_result) {
  neural_keywords <- c("neur", "brain", "axon", "synap", "dendrit", 
                       "myelin", "glial", "nerve", "cerebr", "hippocamp",
                       "cortex", "retina", "nervous", "neural crest", "neurite",
                       "thalamus", "hypothalamus", "amygdala", "ganglia", "neurotransmit")
  
  # Filter the result slot of enrichResult object
  neural_paths <- enrichment_result@result %>%
    filter(grepl(paste(neural_keywords, collapse = "|"), 
                 Description, ignore.case = TRUE)) %>%
    filter(p.adjust < 0.1) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15)
  
  if(nrow(neural_paths) > 0) {
    # Update the result slot while keeping the original enrichResult structure
    enrichment_result@result <- neural_paths
    return(enrichment_result)
  } else {
    return(NULL)
  }
}

# Process the results of each cluster
for(cluster_num in 1:6) {
  cat(sprintf("\nProcessing neural pathways for cluster %d...\n", cluster_num))
  cluster_dir <- file.path(getwd(), paste0("cluster_", cluster_num))
  
  # GO analysis
  go_file <- file.path(cluster_dir, "GO_enrichment.xlsx")
  if(file.exists(go_file)) {
    # Recreate enrichResult object for GO
    go_data <- read.xlsx(go_file)
    if(nrow(go_data) > 0) {
      ego <- new("enrichResult",
                 result = go_data,
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.1,
                 ontology = "ALL",
                 gene = unique(go_data$geneID),
                 keytype = "ENSEMBL")
      
      neural_go <- extract_neural_pathways(ego)
      if(!is.null(neural_go) && nrow(neural_go@result) > 0) {
        p_go <- barplot(neural_go, showCategory = 20) +
          ggtitle(paste("Cluster", cluster_num, "Neural GO Pathways")) +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(
          filename = file.path(cluster_dir, "neural_GO_pathways_barplot.pdf"),
          plot = p_go,
          width = 12,
          height = 8
        )
        
        # 保存结果到Excel
        write.xlsx(neural_go@result,
                   file = file.path(cluster_dir, "neural_GO_pathways.xlsx"),
                   rowNames = FALSE)
        
        cat("Generated neural GO pathway plot and Excel file\n")
      }
    }
  }
  
  # KEGG analysis
  kegg_file <- file.path(cluster_dir, "KEGG_enrichment.xlsx")
  if(file.exists(kegg_file)) {
    # Recreate enrichResult object for KEGG
    kegg_data <- read.xlsx(kegg_file)
    if(nrow(kegg_data) > 0) {
      ekegg <- new("enrichResult",
                   result = kegg_data,
                   pvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.1,
                   organism = "hsa",
                   gene = unique(kegg_data$geneID),
                   keytype = "ENTREZID")
      
      neural_kegg <- extract_neural_pathways(ekegg)
      if(!is.null(neural_kegg) && nrow(neural_kegg@result) > 0) {
        p_kegg <- barplot(neural_kegg, showCategory = 20) +
          ggtitle(paste("Cluster", cluster_num, "Neural KEGG Pathways")) +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(
          filename = file.path(cluster_dir, "neural_KEGG_pathways_barplot.pdf"),
          plot = p_kegg,
          width = 12,
          height = 8
        )
        
        # 保存结果到Excel
        write.xlsx(neural_kegg@result,
                   file = file.path(cluster_dir, "neural_KEGG_pathways.xlsx"),
                   rowNames = FALSE)
        
        cat("Generated neural KEGG pathway plot and Excel file\n")
      }
    }
  }
  
  # Reactome analysis
  reactome_file <- file.path(cluster_dir, "Reactome_enrichment.xlsx")
  if(file.exists(reactome_file)) {
    # Recreate enrichResult object for Reactome
    reactome_data <- read.xlsx(reactome_file)
    if(nrow(reactome_data) > 0) {
      ereactome <- new("enrichResult",
                       result = reactome_data,
                       pvalueCutoff = 0.1,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1,
                       organism = "human",
                       gene = unique(reactome_data$geneID),
                       keytype = "ENTREZID")
      
      neural_reactome <- extract_neural_pathways(ereactome)
      if(!is.null(neural_reactome) && nrow(neural_reactome@result) > 0) {
        p_reactome <- barplot(neural_reactome, showCategory = 20) +
          ggtitle(paste("Cluster", cluster_num, "Neural Reactome Pathways")) +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(
          filename = file.path(cluster_dir, "neural_Reactome_pathways_barplot.pdf"),
          plot = p_reactome,
          width = 12,
          height = 8
        )
        
        # 保存结果到Excel
        write.xlsx(neural_reactome@result,
                   file = file.path(cluster_dir, "neural_Reactome_pathways.xlsx"),
                   rowNames = FALSE)
        
        cat("Generated neural Reactome pathway plot and Excel file\n")
      }
    }
  }
}

# 创建一个汇总表格
create_neural_pathway_summary <- function() {
  summary_data <- data.frame()
  
  for(cluster_num in 1:6) {
    cluster_dir <- file.path(getwd(), paste0("cluster_", cluster_num))
    
    # 收集所有通路信息
    for(db_type in c("GO", "KEGG", "Reactome")) {
      file_path <- file.path(cluster_dir, paste0("neural_", db_type, "_pathways.xlsx"))
      if(file.exists(file_path)) {
        pathway_data <- read.xlsx(file_path)
        if(nrow(pathway_data) > 0) {
          pathway_data$Cluster <- cluster_num
          pathway_data$Database <- db_type
          summary_data <- rbind(summary_data, 
                                pathway_data %>% 
                                  select(Cluster, Database, Description, p.adjust, Count, GeneRatio))
        }
      }
    }
  }
  
  if(nrow(summary_data) > 0) {
    # 排序并保存
    summary_data <- summary_data %>%
      arrange(Cluster, Database, p.adjust)
    
    write.xlsx(summary_data, "neural_pathways_summary.xlsx", rowNames = FALSE)
    cat("\nCreated neural pathways summary file\n")
  }
}

# 运行汇总函数
create_neural_pathway_summary()





# eye related pathways
extract_eye_pathways <- function(enrichment_result) {
  eye_keywords <- c("eye", "ocul", "retina", "cornea", "lens", "vision",
                    "photoreceptor", "optic", "macula", "glaucoma", "cataract",
                    "uvea", "vitre", "iris", "choroid", "conjunctiv", "sclera",
                    "ophth", "blind", "visual", "refractive", "myopia", "presbyopia",
                    "astigmatism", "keratoconus", "strabismus", "amblyopia")
  
  # Filter the result slot of enrichResult object
  eye_paths <- enrichment_result@result %>%
    filter(grepl(paste(eye_keywords, collapse = "|"), 
                 Description, ignore.case = TRUE)) %>%
    filter(p.adjust < 0.1) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15)
  
  if(nrow(eye_paths) > 0) {
    # Update the result slot while keeping the original enrichResult structure
    enrichment_result@result <- eye_paths
    return(enrichment_result)
  } else {
    return(NULL)
  }
}

# Process the results of each cluster
for(cluster_num in 1:6) {
  cat(sprintf("\nProcessing eye-related pathways for cluster %d...\n", cluster_num))
  cluster_dir <- file.path(getwd(), paste0("cluster_", cluster_num))
  
  # GO analysis
  go_file <- file.path(cluster_dir, "GO_enrichment.xlsx")
  if(file.exists(go_file)) {
    # Recreate enrichResult object for GO
    go_data <- read.xlsx(go_file)
    if(nrow(go_data) > 0) {
      ego <- new("enrichResult",
                 result = go_data,
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.1,
                 ontology = "ALL",
                 gene = unique(go_data$geneID),
                 keytype = "ENSEMBL")
      
      eye_go <- extract_eye_pathways(ego)
      if(!is.null(eye_go) && nrow(eye_go@result) > 0) {
        p_go <- barplot(eye_go, showCategory = 20) +
          ggtitle(paste("Cluster", cluster_num, "Eye-related GO Pathways")) +
          theme(axis.text.y = element_text(size = 8))  # 调整y轴标签大小
        
        ggsave(
          filename = file.path(cluster_dir, "eye_GO_pathways_barplot.pdf"),
          plot = p_go,
          width = 12,
          height = 8
        )
        
        # 保存结果到Excel
        write.xlsx(eye_go@result,
                   file = file.path(cluster_dir, "eye_GO_pathways.xlsx"),
                   rowNames = FALSE)
        
        cat("Generated eye-related GO pathway plot and Excel file\n")
      }
    }
  }
  
  # KEGG analysis
  kegg_file <- file.path(cluster_dir, "KEGG_enrichment.xlsx")
  if(file.exists(kegg_file)) {
    # Recreate enrichResult object for KEGG
    kegg_data <- read.xlsx(kegg_file)
    if(nrow(kegg_data) > 0) {
      ekegg <- new("enrichResult",
                   result = kegg_data,
                   pvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.1,
                   organism = "hsa",
                   gene = unique(kegg_data$geneID),
                   keytype = "ENTREZID")
      
      eye_kegg <- extract_eye_pathways(ekegg)
      if(!is.null(eye_kegg) && nrow(eye_kegg@result) > 0) {
        p_kegg <- barplot(eye_kegg, showCategory = 20) +
          ggtitle(paste("Cluster", cluster_num, "Eye-related KEGG Pathways")) +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(
          filename = file.path(cluster_dir, "eye_KEGG_pathways_barplot.pdf"),
          plot = p_kegg,
          width = 12,
          height = 8
        )
        
        # 保存结果到Excel
        write.xlsx(eye_kegg@result,
                   file = file.path(cluster_dir, "eye_KEGG_pathways.xlsx"),
                   rowNames = FALSE)
        
        cat("Generated eye-related KEGG pathway plot and Excel file\n")
      }
    }
  }
  
  # Reactome analysis
  reactome_file <- file.path(cluster_dir, "Reactome_enrichment.xlsx")
  if(file.exists(reactome_file)) {
    # Recreate enrichResult object for Reactome
    reactome_data <- read.xlsx(reactome_file)
    if(nrow(reactome_data) > 0) {
      ereactome <- new("enrichResult",
                       result = reactome_data,
                       pvalueCutoff = 0.1,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1,
                       organism = "human",
                       gene = unique(reactome_data$geneID),
                       keytype = "ENTREZID")
      
      eye_reactome <- extract_eye_pathways(ereactome)
      if(!is.null(eye_reactome) && nrow(eye_reactome@result) > 0) {
        p_reactome <- barplot(eye_reactome, showCategory = 20) +
          ggtitle(paste("Cluster", cluster_num, "Eye-related Reactome Pathways")) +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(
          filename = file.path(cluster_dir, "eye_Reactome_pathways_barplot.pdf"),
          plot = p_reactome,
          width = 12,
          height = 8
        )
        
        # 保存结果到Excel
        write.xlsx(eye_reactome@result,
                   file = file.path(cluster_dir, "eye_Reactome_pathways.xlsx"),
                   rowNames = FALSE)
        
        cat("Generated eye-related Reactome pathway plot and Excel file\n")
      }
    }
  }
}

# 创建一个汇总表格
create_eye_pathway_summary <- function() {
  summary_data <- data.frame()
  
  for(cluster_num in 1:6) {
    cluster_dir <- file.path(getwd(), paste0("cluster_", cluster_num))
    
    # 收集所有通路信息
    for(db_type in c("GO", "KEGG", "Reactome")) {
      file_path <- file.path(cluster_dir, paste0("eye_", db_type, "_pathways.xlsx"))
      if(file.exists(file_path)) {
        pathway_data <- read.xlsx(file_path)
        if(nrow(pathway_data) > 0) {
          pathway_data$Cluster <- cluster_num
          pathway_data$Database <- db_type
          summary_data <- rbind(summary_data, 
                                pathway_data %>% 
                                  select(Cluster, Database, Description, p.adjust, Count, GeneRatio))
        }
      }
    }
  }
  
  if(nrow(summary_data) > 0) {
    # 排序并保存
    summary_data <- summary_data %>%
      arrange(Cluster, Database, p.adjust)
    
    write.xlsx(summary_data, "eye_pathways_summary.xlsx", rowNames = FALSE)
    cat("\nCreated eye pathways summary file\n")
  }
}





