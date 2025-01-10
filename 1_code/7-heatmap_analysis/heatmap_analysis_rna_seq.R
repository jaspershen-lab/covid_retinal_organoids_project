library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')


library(tidymass)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(Mfuzz)
library(tidyr)
library(ggplot2)


load("3_data_analysis/2-data_cleaning/1-transcriptome/transcriptome_data.RData")
load("2_data/pathway_human/pathway_GO.rda")
load("2_data/pathway_human/pathway_kegg.rda")
load("2_data/pathway_human/pathway_Reactome.rda")


dir.create(
  "3_data_analysis/7-heatmap_analysis/rna_seq",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/7-heatmap_analysis/rna_seq")



# 数据预处理函数
preprocess_data <- function(expression_data, min_value = 1) {
  # 对数转换
  if(max(expression_data, na.rm = TRUE) > 100) {
    expression_data <- log2(expression_data + min_value)
  }
  # 行标准化
  expression_data_scaled <- t(scale(t(expression_data)))
  return(expression_data_scaled)
}



# screen sample ID
all_sample_ids <- colnames(transcriptome_data@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "WT_1", "WT_2", "WT_3", "BA52_1", "BA52_2", "BA52_3")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)

# 获取共同基因
common_genes <- transcriptome_data@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids)) %>%
  na.omit() %>%
  rownames()

common_gene_info <- transcriptome_data@variable_info %>%
  dplyr::filter(variable_id %in% common_genes)

# 提取表达矩阵
expression_data <- transcriptome_data@expression_data[common_genes, filtered_sample_ids]
unique_ensembl <- make.unique(common_gene_info$ENSEMBL)
rownames(expression_data) <- unique_ensembl

# 计算pathway重叠
calculate_overlap <- function(pathway_data, gene_list) {
  pathway_overlap <- lapply(pathway_data$gene_list, function(pathway_genes) {
    common_genes <- intersect(pathway_genes$ID, gene_list)
    overlap_percentage <- length(common_genes) / length(pathway_genes$ID) * 100
    list(
      common_genes = common_genes,
      overlap_percentage = overlap_percentage
    )
  })
  return(pathway_overlap)
}

# 筛选pathway基因
common_ensembl_genes <- rownames(expression_data)
go_overlap <- calculate_overlap(go_database, common_ensembl_genes)
kegg_overlap <- calculate_overlap(kegg_database, common_ensembl_genes)
reactome_overlap <- calculate_overlap(reactome_database, common_ensembl_genes)

# 筛选overlap大于50%的pathway
filter_pathways <- function(overlap_results, threshold = 50) {
  Filter(function(x) x$overlap_percentage >= threshold, overlap_results)
}

go_filtered <- filter_pathways(go_overlap, 50)
kegg_filtered <- filter_pathways(kegg_overlap, 50)
reactome_filtered <- filter_pathways(reactome_overlap, 50)

# 合并筛选后的pathway基因
all_pathway_genes <- unique(c(
  unlist(lapply(go_filtered, function(x) x$common_genes)),
  unlist(lapply(kegg_filtered, function(x) x$common_genes)),
  unlist(lapply(reactome_filtered, function(x) x$common_genes))
))


# Count the number of genes screened in each database
go_gene_count <- length(unique(unlist(lapply(go_filtered, function(x) x$common_genes))))
kegg_gene_count <- length(unique(unlist(lapply(kegg_filtered, function(x) x$common_genes))))
reactome_gene_count <- length(unique(unlist(lapply(reactome_filtered, function(x) x$common_genes))))
total_gene_count <- length(all_pathway_genes)

cat("Number of genes mapped in GO pathways:", go_gene_count, "\n")
cat("Number of genes mapped in KEGG pathways:", kegg_gene_count, "\n")
cat("Number of genes mapped in Reactome pathways:", reactome_gene_count, "\n")
cat("Total unique genes mapped across all pathways:", total_gene_count, "\n")

# 提取热图矩阵并预处理
heatmap_matrix <- expression_data[all_pathway_genes, ]
heatmap_matrix <- preprocess_data(heatmap_matrix)

# 设置基因名称
rownames(heatmap_matrix) <- common_gene_info$SYMBOL[match(rownames(heatmap_matrix), common_gene_info$ENSEMBL)]
rownames(heatmap_matrix) <- make.unique(rownames(heatmap_matrix))

# 选择top变异基因
gene_vars <- apply(heatmap_matrix, 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(100, length(gene_vars))])


# Set heat map colors
color_scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Add annotation
annotation_col <- data.frame(
  Group = transcriptome_data@sample_info$group[match(filtered_sample_ids, transcriptome_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

annotation_row <- data.frame(
  Pathway = sample(c("Pathway1", "Pathway2", "Pathway3"), nrow(heatmap_matrix[top_var_genes, ]), replace = TRUE),
  stringsAsFactors = FALSE
)
rownames(annotation_row) <- rownames(heatmap_matrix[top_var_genes, ])

# plot heatmap
plot <- pheatmap(
  heatmap_matrix[top_var_genes, ], 
  color = color_scheme,                   
  cluster_rows = TRUE,                   
  cluster_cols = FALSE,               
  show_rownames = TRUE,                   
  show_colnames = TRUE,                   
  annotation_col = annotation_col,        
  # annotation_row = annotation_row,        
  main = "Heatmap of All Top Variable Genes",  
  fontsize_row = 6,                       
  fontsize_col = 10,
  cellwidth = 20,                         
  cellheight = 8,                        
  border_color = NA                       
)
plot

ggsave("heatmap_100.pdf", plot, width = 10, height = 18)



# 准备聚类数据
expression_matrix <- heatmap_matrix[top_var_genes, ]

# 创建时间点数据框
create_time_data <- function(sample_names) {
  data.frame(
    sample = sample_names,
    time = case_when(
      grepl("Ctrl", sample_names) ~ 1,
      grepl("WT", sample_names) ~ 2,
      grepl("BA52", sample_names) ~ 3
    )
  )
}

time_data <- create_time_data(colnames(expression_matrix))

# 将数据转换为Mfuzz所需的格式
expression_set <- ExpressionSet(assayData = as.matrix(expression_matrix))
data.s <- standardise(expression_set)

# 估计最佳的m参数
m1 <- mestimate(data.s)

# 确定最佳聚类数
Dmin(data.s, m = m1, crange = seq(2, 20, 2), repeats = 5, visu = TRUE)
ggsave("cluster_number_optimization.pdf", width = 8, height = 6)

# 执行聚类
cluster_number <- 6  # 可以根据上面的分析结果调整
c <- mfuzz(data.s, c = cluster_number, m = m1)

# 设置membership阈值
membership_cutoff <- 0.5

# 可视化每个cluster
for (idx in 1:cluster_number) {
  # 获取当前cluster的基因
  cluster_data <- data.frame(
    gene_id = names(c$cluster),
    membership = c$membership[,idx],
    cluster = c$cluster
  ) %>%
    dplyr::filter(cluster == idx, membership > membership_cutoff)
  
  # 提取基因表达数据
  gene_expr <- expression_matrix[cluster_data$gene_id,] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")
  
  # 转换为长格式
  gene_expr_long <- gene_expr %>%
    pivot_longer(-gene_id, 
                 names_to = "sample",
                 values_to = "expression") %>%
    left_join(time_data, by = "sample") %>%
    left_join(cluster_data, by = "gene_id")
  
  # 计算均值线数据
  mean_expr <- gene_expr_long %>%
    group_by(time) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE))
  
  # 创建cluster可视化
  plot <- ggplot() +
    # 显示个别基因表达曲线
    geom_line(data = gene_expr_long, 
              aes(x = time, y = expression, group = gene_id, alpha = membership),
              color = "grey50") +
    # 显示均值线
    geom_line(data = mean_expr,
              aes(x = time, y = mean_expression),
              color = "red", size = 1.5) +
    # 添加点以显示实际数据点
    geom_point(data = mean_expr,
               aes(x = time, y = mean_expression),
               color = "red", size = 3) +
    # 设置主题和标签
    theme_bw() +
    labs(title = paste("Cluster", idx, "(", nrow(cluster_data), " genes)"),
         x = "Time",
         y = "Expression") +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    scale_x_continuous(breaks = 1:3, labels = c("Ctrl", "WT", "BA52"))
  
  # 保存图片
  ggsave(
    file.path(paste0("cluster_", idx, "_improved.pdf")),
    plot,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(paste0("cluster_", idx, "_improved.png")),
    plot,
    width = 8,
    height = 6
  )
  
  # 保存cluster基因信息
  write.csv(
    cluster_data %>% 
      left_join(common_gene_info, by = c("gene_id" = "ENSEMBL")),
    file = paste0("cluster_", idx, "_genes.csv"),
    row.names = FALSE
  )
}

# 保存聚类结果
save(c, file = "clustering_results.RData")

# 分析cluster相关性
cluster_correlations <- cor(t(c$centers))
pdf("cluster_correlations.pdf", width = 8, height = 8)
corrplot::corrplot(
  cluster_correlations,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  diag = FALSE
)
dev.off()

# 保存最终的cluster信息
final_cluster_info <- data.frame(
  gene_id = names(c$cluster),
  cluster = c$cluster,
  max_membership = apply(c$membership, 1, max)
) %>%
  filter(max_membership >= membership_cutoff)

write.csv(final_cluster_info, "final_cluster_info.csv", row.names = FALSE)








#########neuro_related pathways
# neuro_related keywords
neuro_keywords <- c("neuro", "synapse", "axon", "nerve", "brain", "central nervous system", "CNS", "spinal cord")

# filter_neuro_pathways
filter_neuro_pathways <- function(pathway_data, keywords) {
  pathway_descriptions <- pathway_data$description
  selected_indices <- sapply(pathway_descriptions, function(desc) {
    any(grepl(paste(keywords, collapse = "|"), desc, ignore.case = TRUE))
  })
  pathway_data$gene_list[selected_indices]
}

# neuro_pathways
go_neuro_pathways <- filter_neuro_pathways(go_database, neuro_keywords)
kegg_neuro_pathways <- filter_neuro_pathways(kegg_database, neuro_keywords)
reactome_neuro_pathways <- filter_neuro_pathways(reactome_database, neuro_keywords)

# overlap gene
go_neuro_overlap <- calculate_overlap(list(gene_list = go_neuro_pathways), common_ensembl_genes)
kegg_neuro_overlap <- calculate_overlap(list(gene_list = kegg_neuro_pathways), common_ensembl_genes)
reactome_neuro_overlap <- calculate_overlap(list(gene_list = reactome_neuro_pathways), common_ensembl_genes)

# Screening threshold >= 50%
go_neuro_filtered <- filter_pathways(go_neuro_overlap, 50)
kegg_neuro_filtered <- filter_pathways(kegg_neuro_overlap, 50)
reactome_neuro_filtered <- filter_pathways(reactome_neuro_overlap, 50)

# neuro_pathway_genes
neuro_pathway_genes <- unique(c(
  unlist(lapply(go_neuro_filtered, function(x) x$common_genes)),
  unlist(lapply(kegg_neuro_filtered, function(x) x$common_genes)),
  unlist(lapply(reactome_neuro_filtered, function(x) x$common_genes))
))


# The total number of mapping genes of GO, KEGG and Reactome was calculated
go_mapped_genes <- unique(unlist(lapply(go_neuro_filtered, function(x) x$common_genes)))
kegg_mapped_genes <- unique(unlist(lapply(kegg_neuro_filtered, function(x) x$common_genes)))
reactome_mapped_genes <- unique(unlist(lapply(reactome_neuro_filtered, function(x) x$common_genes)))
all_mapped_genes <- unique(c(go_mapped_genes, kegg_mapped_genes, reactome_mapped_genes))

cat("GO Mapped Genes Count:", length(go_mapped_genes), "\n")
cat("KEGG Mapped Genes Count:", length(kegg_mapped_genes), "\n")
cat("Reactome Mapped Genes Count:", length(reactome_mapped_genes), "\n")
cat("Total Unique Mapped Genes Count:", length(all_mapped_genes), "\n")

# neuro_heatmap_matrix
neuro_heatmap_matrix <- expression_data[neuro_pathway_genes, ]
neuro_heatmap_matrix <- t(scale(t(neuro_heatmap_matrix)))  

# Z-score
neuro_heatmap_matrix <- t(scale(t(expression_data)))

# rownames SYMBOL
rownames(neuro_heatmap_matrix ) <- common_gene_info$SYMBOL[match(rownames(neuro_heatmap_matrix ), common_gene_info$ENSEMBL)]
rownames(neuro_heatmap_matrix ) <- make.unique(rownames(neuro_heatmap_matrix ))

#  Top 100 gene
gene_vars <- apply(neuro_heatmap_matrix, 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(100, length(gene_vars))])

# Set heat map colors
color_scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Add annotation
annotation_col <- data.frame(
  Group = transcriptome_data@sample_info$group[match(filtered_sample_ids, transcriptome_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

annotation_row <- data.frame(
  Pathway = sample(c("Pathway1", "Pathway2", "Pathway3"), nrow(neuro_heatmap_matrix[top_var_genes, ]), replace = TRUE),
  stringsAsFactors = FALSE
)
rownames(annotation_row) <- rownames(neuro_heatmap_matrix[top_var_genes, ])

# plot heatmap
plot <- pheatmap(
  neuro_heatmap_matrix[top_var_genes, ],
  color = color_scheme,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Heatmap of Neuro-Related Pathway Genes",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)

plot


ggsave("neuro_heatmap_100.pdf", plot, width = 10, height = 18)





######### eye disease-related pathways
# Define keywords for eye-related pathways
eye_keywords <- c("eye", "retina", "vision", "optic", "ocular", "cornea", "lens", "macula", "glaucoma", "cataract")

# Filter eye-related pathways
filter_eye_pathways <- function(pathway_data, keywords) {
  pathway_descriptions <- pathway_data$description
  selected_indices <- sapply(pathway_descriptions, function(desc) {
    any(grepl(paste(keywords, collapse = "|"), desc, ignore.case = TRUE))
  })
  pathway_data$gene_list[selected_indices]
}

# Screen eye-related pathways from GO, KEGG, Reactome
go_eye_pathways <- filter_eye_pathways(go_database, eye_keywords)
kegg_eye_pathways <- filter_eye_pathways(kegg_database, eye_keywords)
reactome_eye_pathways <- filter_eye_pathways(reactome_database, eye_keywords)

# Calculate overlapping genes
go_eye_overlap <- calculate_overlap(list(gene_list = go_eye_pathways), common_ensembl_genes)
kegg_eye_overlap <- calculate_overlap(list(gene_list = kegg_eye_pathways), common_ensembl_genes)
reactome_eye_overlap <- calculate_overlap(list(gene_list = reactome_eye_pathways), common_ensembl_genes)

# Filter pathways with a threshold >= 50%
go_eye_filtered <- filter_pathways(go_eye_overlap, 50)
kegg_eye_filtered <- filter_pathways(kegg_eye_overlap, 50)
reactome_eye_filtered <- filter_pathways(reactome_eye_overlap, 50)

# Combine genes from filtered pathways
eye_pathway_genes <- unique(c(
  unlist(lapply(go_eye_filtered, function(x) x$common_genes)),
  unlist(lapply(kegg_eye_filtered, function(x) x$common_genes)),
  unlist(lapply(reactome_eye_filtered, function(x) x$common_genes))
))


# The total number of mapping genes of GO, KEGG and Reactome was calculated
go_mapped_genes <- unique(unlist(lapply(go_eye_filtered, function(x) x$common_genes)))
kegg_mapped_genes <- unique(unlist(lapply(kegg_eye_filtered, function(x) x$common_genes)))
reactome_mapped_genes <- unique(unlist(lapply(reactome_eye_filtered, function(x) x$common_genes)))
all_mapped_genes <- unique(c(go_mapped_genes, kegg_mapped_genes, reactome_mapped_genes))

cat("GO Mapped Genes Count:", length(go_mapped_genes), "\n")
cat("KEGG Mapped Genes Count:", length(kegg_mapped_genes), "\n")
cat("Reactome Mapped Genes Count:", length(reactome_mapped_genes), "\n")
cat("Total Unique Mapped Genes Count:", length(all_mapped_genes), "\n")


# Extract the expression matrix for eye-related genes
eye_heatmap_matrix <- expression_data[eye_pathway_genes, ]
eye_heatmap_matrix <- t(scale(t(eye_heatmap_matrix)))  # Standardize by row (Z-score)

# Update rownames to SYMBOL
rownames(eye_heatmap_matrix) <- common_gene_info$SYMBOL[match(rownames(eye_heatmap_matrix), common_gene_info$ENSEMBL)]
rownames(eye_heatmap_matrix) <- make.unique(rownames(eye_heatmap_matrix))

# Select top 100 variable genes
gene_vars <- apply(eye_heatmap_matrix, 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(100, length(gene_vars))])

# Set heat map colors
color_scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Add column annotation
annotation_col <- data.frame(
  Group = transcriptome_data@sample_info$group[match(filtered_sample_ids, transcriptome_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

# Add row annotation
annotation_row <- data.frame(
  Pathway = sample(c("Pathway1", "Pathway2", "Pathway3"), nrow(eye_heatmap_matrix[top_var_genes, ]), replace = TRUE),
  stringsAsFactors = FALSE
)
rownames(annotation_row) <- rownames(eye_heatmap_matrix[top_var_genes, ])

# Plot heatmap
plot <- pheatmap(
  eye_heatmap_matrix[top_var_genes, ],
  color = color_scheme,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Heatmap of Eye Disease-Related Pathway Genes",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)

plot

# Save heatmap
ggsave("eye_heatmap_100.pdf", plot, width = 10, height = 18)
