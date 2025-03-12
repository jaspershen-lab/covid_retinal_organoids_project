library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')


library(tidymass)
library(dplyr)
library(Mfuzz)
library(ComplexHeatmap)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(corrplot)



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

# Set heat map colors
color_scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Add annotation
annotation_col <- data.frame(
  Group = transcriptome_data@sample_info$group[match(filtered_sample_ids, transcriptome_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

# Plot horizontal heatmap
plot <- pheatmap(
  t(heatmap_matrix[top_var_genes, ]), # 转置矩阵
  color = color_scheme,
  cluster_rows = FALSE,    
  cluster_cols = TRUE,     
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = annotation_col,  
  main = "Heatmap of All Top Variable Genes",
  fontsize_row = 10,     
  fontsize_col = 6,      
  cellwidth = 8,          
  cellheight = 20,       
  border_color = NA
)


ggsave("heatmap_100_horizontal.pdf", plot, width = 18, height = 10) 


#################### RNA-seq Clustering Analysis ####################
library(tidymass)
library(dplyr)
library(Mfuzz)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(pheatmap)

# Create output directory
dir.create("clustering_results", recursive = TRUE, showWarnings = FALSE)

######################## 1. Data Preprocessing ########################
# Screen sample IDs
# Screen sample IDs
all_sample_ids <- colnames(transcriptome_data@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", 
                   "WT_1", "WT_2", "WT_3", 
                   "BA52_1", "BA52_2", "BA52_3")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)

# Extract expression matrix
expression_matrix <- transcriptome_data@expression_data[, filtered_sample_ids]

# Remove rows with very low variance
row_vars <- apply(expression_matrix, 1, var)
expression_matrix <- expression_matrix[row_vars > 1e-10, ]

# Calculate group means
grouped_matrix <- matrix(0, nrow = nrow(expression_matrix), ncol = 3)
colnames(grouped_matrix) <- c("Ctrl", "WT", "BA52")
rownames(grouped_matrix) <- rownames(expression_matrix)

# Calculate means for each group
grouped_matrix[, "Ctrl"] <- rowMeans(expression_matrix[, grep("Ctrl", colnames(expression_matrix))])
grouped_matrix[, "WT"] <- rowMeans(expression_matrix[, grep("WT", colnames(expression_matrix))])
grouped_matrix[, "BA52"] <- rowMeans(expression_matrix[, grep("BA52", colnames(expression_matrix))])

# Scale the grouped data
grouped_matrix_scaled <- t(scale(t(grouped_matrix)))

# ######################## 2. Prepare Mfuzz Input ########################
# Create time points
temp_data <- rbind(
  time = 1:3,
  grouped_matrix_scaled
)

# Save temporary data
write.table(
  temp_data,
  file = "clustering_results/temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

# Create ExpressionSet
data <- table2eset(filename = "clustering_results/temp_data.txt")

######################## 3. Clustering Analysis ########################
# Estimate best m parameter
m1 <- mestimate(data)
print(paste("Estimated m parameter:", m1))

# Optimize cluster number
dmin_results <- Dmin(data, m = m1, crange = seq(2, 20, 2), repeats = 5, visu = TRUE)
ggsave("clustering_results/cluster_number_optimization.pdf", width = 8, height = 6)

# Perform clustering
cluster_number <- 4  # Based on previous analysis
# 设置随机种子以确保结果可重复
set.seed(123)
c <- mfuzz(data, c = cluster_number, m = m1)

######################## 4. Create Visualizations ########################
# 4.1 Correlation Heatmap
# Calculate cluster centers
n_clusters <- ncol(c$membership)
n_timepoints <- ncol(exprs(data))
centers <- matrix(0, nrow = n_clusters, ncol = n_timepoints)

for(i in 1:n_clusters) {
  memb <- c$membership[, i]
  valid_idx <- memb >= 0.5
  if(sum(valid_idx) > 0) {
    weights <- memb[valid_idx]
    expr_data <- exprs(data)[valid_idx, ]
    centers[i,] <- colSums(weights * expr_data) / sum(weights)
  }
}

rownames(centers) <- paste0("Cluster", 1:n_clusters)
colnames(centers) <- colnames(exprs(data))

# Show and save correlation plot
plot <- corrplot(
  cor(t(centers)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  col = colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))(100),
  number.cex = 0.7,
  addCoef.col = "black"
)

# Save plot
pdf("clustering_results/cluster_correlations.pdf", width = 10, height = 10)
corrplot(
  cor(t(centers)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  col = colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))(100),
  number.cex = 0.7,
  addCoef.col = "black"
)
dev.off()
# 4.2 Cluster Profiles
membership_cutoff <- 0.7

mfuzz.plot(
  eset = data,
  min.mem = membership_cutoff,
  cl = c,
  mfrow = c(2, 2),
  time.labels = 1:3,
  new.window = FALSE
)

# Plot cluster profiles
for(cluster_id in 1:cluster_number) {
  cluster_members <- which(c$cluster == cluster_id)
  cluster_data <- grouped_matrix_scaled[cluster_members, ]
  
  plot_data <- data.frame()
  for(i in 1:nrow(cluster_data)) {
    gene_data <- data.frame(
      Gene = paste0("Gene", i),
      Time = colnames(cluster_data),
      Value = as.numeric(cluster_data[i, ]),
      Membership = c$membership[cluster_members[i], cluster_id]
    )
    plot_data <- rbind(plot_data, gene_data)
  }
  
  plot_data$Time <- factor(plot_data$Time, levels = c("Ctrl", "WT", "BA52"))
  mean_profile <- aggregate(Value ~ Time, data = plot_data, FUN = mean)
  mean_profile$Time <- factor(mean_profile$Time, levels = c("Ctrl", "WT", "BA52"))
  
  p <- ggplot() +
    geom_line(data = plot_data, 
              aes(x = Time, y = Value, group = Gene,  
                  color = Membership),
              size = 0.8) +
    geom_line(data = mean_profile,
              aes(x = Time, y = Value, group = 1),
              color = "black",  
              size = 1.2) +
    scale_color_gradientn(
      colors = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                 "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"),   
      limits = c(0.2, 0.7),    # Adjusted from 0.9 to 0.7 since 3rd Qu. is 0.6336
      oob = scales::squish,    
      breaks = seq(0.2, 0.7, by = 0.1),  # Adjusted breaks
      name = "Membership"
    ) +
    theme_bw() +         
    labs(title = paste("Cluster", cluster_id, 
                       "(n=", length(cluster_members), ")"),
         x = "Group",
         y = "Expression") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      legend.position = "top",
      legend.key.width = unit(3, "cm")
    )
  
  ggsave(paste0("clustering_results/cluster_", cluster_id, "_profile.pdf"),
         p, width = 10, height = 7)
  
  
  # Save cluster gene information with symbols
  gene_info <- data.frame(
    gene_id = names(cluster_members),
    membership = c$membership[cluster_members, cluster_id],
    stringsAsFactors = FALSE
  ) %>%
    left_join(transcriptome_data@variable_info, by = c("gene_id" = "variable_id"))
  
  write.csv(
    gene_info,
    file = paste0("clustering_results/cluster_", cluster_id, "_genes.csv"),
    row.names = FALSE
  )
}

# Print final summary
print("Clustering Analysis Complete!")
print("Number of genes in each cluster:")
print(table(c$cluster))
print("\nMembership summary:")
print(summary(apply(c$membership, 1, max)))




# Create cluster info dataframe
cluster_info <- data.frame(
  gene_id = rownames(grouped_matrix_scaled),
  cluster = c$cluster
)

# Add membership info
cluster_info <- cbind(
  cluster_info,
  as.data.frame(c$membership, 
                stringsAsFactors = FALSE)
)
colnames(cluster_info)[3:ncol(cluster_info)] <- paste0("membership_cluster", 1:cluster_number)

cluster_info$max_membership <- apply(c$membership, 1, max)

# Now print distribution
print(table(cluster_info$cluster))


# Process clusters with membership cutoff
membership_cutoff <- 0.5
cluster_summary <- unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    cluster_info %>%
      dplyr::select(gene_id, paste0("membership_cluster", x), cluster) %>%
      dplyr::rename(membership = paste0("membership_cluster", x)) %>%
      dplyr::filter(membership >= membership_cutoff) %>%
      dplyr::mutate(
        cluster_raw = cluster,
        cluster = x
      )
  }) %>%
  dplyr::bind_rows()

# Count genes per cluster
print("\nTotal genes per cluster:")
print(cluster_summary %>% dplyr::count(cluster))

# Save results
final_cluster_info <- cluster_summary
save(final_cluster_info, file = "final_cluster_info.RData")

# Optional: Export to CSV for easier viewing
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
