library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(tidymass)
library(tidymass)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(org.Hs.eg.db)
library(clusterProfiler)

##read data
load("3_data_analysis/1-data_preparation/2-metabolome/gene_data.rda")
load("2_data/pathway_human/pathway_GO.rda")
load("2_data/pathway_human/pathway_kegg.rda")
load("2_data/pathway_human/pathway_Reactome.rda")

dir.create(
  "3_data_analysis/7-heatmap_analysis/gene_data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/7-heatmap_analysis/gene_data")



# Define sample IDs in desired order
desired_order <- c("Mock-A", "Mock-B", "Mock-C", 
                   "WT-A", "WT-B", "WT-C", 
                   "BA52-A", "BA52-B", "BA52-C")

# 找出不含NA的基因
complete_genes <- rownames(gene_data@expression_data)[complete.cases(gene_data@expression_data[, desired_order])]
print(paste("Number of genes without NA:", length(complete_genes)))

complete_gene_info <- gene_data@variable_info %>%
  dplyr::filter(variable_id %in% complete_genes)

# 创建variable_id到ENTREZID的映射
id_to_entrez <- setNames(complete_gene_info$ENTREZID, complete_gene_info$variable_id)

# 提取表达矩阵并确保它是matrix格式
expression_data <- as.matrix(gene_data@expression_data[complete_genes, desired_order])

# 将ENTREZID转换为ENSEMBL ID
id_mapping <- bitr(complete_gene_info$ENTREZID,
                   fromType = "ENTREZID",
                   toType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db)

# Calculate the overlap of pathways and genes
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

# Use ENSEMBL IDs for pathway analysis
ensembl_genes <- id_mapping$ENSEMBL

go_overlap <- calculate_overlap(go_database, ensembl_genes)
kegg_overlap <- calculate_overlap(kegg_database, ensembl_genes)
reactome_overlap <- calculate_overlap(reactome_database, ensembl_genes)

filter_pathways <- function(overlap_results, threshold = 10) {
  Filter(function(x) x$overlap_percentage >= threshold, overlap_results)
}

go_filtered <- filter_pathways(go_overlap)
kegg_filtered <- filter_pathways(kegg_overlap)
reactome_filtered <- filter_pathways(reactome_overlap)

# Get the corresponding ENTREZID for the pathway genes
get_entrez_for_ensembl <- function(ensembl_ids) {
  id_mapping$ENTREZID[match(ensembl_ids, id_mapping$ENSEMBL)]
}

# Combined screened pathway genes (convert back to ENTREZID)
pathway_entrez <- unique(c(
  get_entrez_for_ensembl(unlist(lapply(go_filtered, function(x) x$common_genes))),
  get_entrez_for_ensembl(unlist(lapply(kegg_filtered, function(x) x$common_genes))),
  get_entrez_for_ensembl(unlist(lapply(reactome_filtered, function(x) x$common_genes)))
))

# 将ENTREZID转换为variable_id
pathway_vars <- names(id_to_entrez)[id_to_entrez %in% pathway_entrez]

# Count the number of genes
go_gene_count <- length(unique(get_entrez_for_ensembl(unlist(lapply(go_filtered, function(x) x$common_genes)))))
kegg_gene_count <- length(unique(get_entrez_for_ensembl(unlist(lapply(kegg_filtered, function(x) x$common_genes)))))
reactome_gene_count <- length(unique(get_entrez_for_ensembl(unlist(lapply(reactome_filtered, function(x) x$common_genes)))))
total_gene_count <- length(pathway_vars)

cat("Number of genes mapped in GO pathways:", go_gene_count, "\n")
cat("Number of genes mapped in KEGG pathways:", kegg_gene_count, "\n")
cat("Number of genes mapped in Reactome pathways:", reactome_gene_count, "\n")
cat("Total unique genes mapped across all pathways:", total_gene_count, "\n")

# Extract the heat map matrix
heatmap_matrix <- expression_data[pathway_vars, ]

# 检查数据
print("Matrix dimensions:")
print(dim(heatmap_matrix))

# 对每一行进行Z-score标准化
heatmap_matrix <- t(scale(t(heatmap_matrix)))

# Set rownames to SYMBOL
gene_symbols <- complete_gene_info$SYMBOL[match(rownames(heatmap_matrix), complete_gene_info$variable_id)]
rownames(heatmap_matrix) <- make.unique(gene_symbols)

# Set heat map colors
color_scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Add annotation
annotation_col <- data.frame(
  Group = gene_data@sample_info$group[match(desired_order, gene_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- desired_order

# Plot heatmap
plot <- pheatmap(
  heatmap_matrix, 
  color = color_scheme,                   
  cluster_rows = TRUE,                   
  cluster_cols = FALSE,               
  show_rownames = TRUE,                   
  show_colnames = TRUE,                   
  annotation_col = annotation_col,        
  main = "Heatmap of All Pathway Genes",  
  fontsize_row = 6,                       
  fontsize_col = 10,
  cellwidth = 20,                         
  cellheight = 8,                        
  border_color = NA                       
)

plot

# Save the plot
pdf("heatmap_all_pathway_genes.pdf", width = 10, height = 18)
print(plot)
dev.off()






library(tidymass)
library(dplyr)
library(Mfuzz)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(pheatmap)


# Store original data in a different name to avoid conflicts
original_gene_data <- gene_data

# Create output directory
dir.create("clustering_results", recursive = TRUE, showWarnings = FALSE)

######################## 1. Data Preprocessing ########################
# Define sample IDs in desired order
desired_order <- c("Mock-A", "Mock-B", "Mock-C", 
                   "WT-A", "WT-B", "WT-C", 
                   "BA52-A", "BA52-B", "BA52-C")

# Extract expression matrix
expression_matrix <- gene_data@expression_data[, desired_order]

# Get complete genes
complete_genes <- rownames(expression_matrix)[complete.cases(expression_matrix)]
expression_matrix <- expression_matrix[complete_genes, ]

# Remove rows with very low variance
row_vars <- apply(expression_matrix, 1, var)
expression_matrix <- expression_matrix[row_vars > 1e-10, ]

# Calculate group means
grouped_matrix <- matrix(0, nrow = nrow(expression_matrix), ncol = 3)
colnames(grouped_matrix) <- c("Mock", "WT", "BA52")
rownames(grouped_matrix) <- rownames(expression_matrix)

# Calculate means for each group
grouped_matrix[, "Mock"] <- rowMeans(expression_matrix[, grep("Mock", colnames(expression_matrix))])
grouped_matrix[, "WT"] <- rowMeans(expression_matrix[, grep("WT", colnames(expression_matrix))])
grouped_matrix[, "BA52"] <- rowMeans(expression_matrix[, grep("BA52", colnames(expression_matrix))])

# Scale the grouped data
grouped_matrix_scaled <- t(scale(t(grouped_matrix)))

######################## 2. Prepare Mfuzz Input ########################
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
cluster_number <- 3
set.seed(123)
c <- mfuzz(data, c = cluster_number, m = m1)

######################## 4. Create Visualizations ########################
# 4.1 Correlation Heatmap
n_clusters <- ncol(c$membership)
n_timepoints <- ncol(exprs(data))
centers <- matrix(0, nrow = n_clusters, ncol = n_timepoints)

for(i in 1:n_clusters) {
  memb <- c$membership[, i]
  valid_idx <- memb >= 0.5
  if(sum(valid_idx) > 0) {
    weights <- memb[valid_idx]
    expr_data <- exprs(data)[valid_idx, , drop = FALSE]  # 添加 drop = FALSE
    centers[i,] <- colSums(weights * expr_data) / sum(weights)
  }
}

rownames(centers) <- paste0("Cluster", 1:n_clusters)
colnames(centers) <- colnames(exprs(data))


# Plot correlation heatmap
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
membership_cutoff <- 0.6

mfuzz.plot(
  eset = data,
  min.mem = membership_cutoff,
  cl = c,
  mfrow = c(2, 2),
  time.labels = 1:3,
  new.window = FALSE
)

# 4.2 Cluster Profiles
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
  
  plot_data$Time <- factor(plot_data$Time, levels = c("Mock", "WT", "BA52"))
  mean_profile <- aggregate(Value ~ Time, data = plot_data, FUN = mean)
  mean_profile$Time <- factor(mean_profile$Time, levels = c("Mock", "WT", "BA52"))
  
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
      limits = c(0.2, 0.7),
      oob = scales::squish,    
      breaks = seq(0.2, 0.7, by = 0.1),
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
  
  
  # Save cluster gene information using S4 object
  gene_info <- data.frame(
    gene_id = rownames(grouped_matrix_scaled)[cluster_members],
    membership = c$membership[cluster_members, cluster_id],
    SYMBOL = original_gene_data@variable_info$SYMBOL[match(
      rownames(grouped_matrix_scaled)[cluster_members], 
      original_gene_data@variable_info$variable_id
    )],
    ENTREZID = original_gene_data@variable_info$ENTREZID[match(
      rownames(grouped_matrix_scaled)[cluster_members], 
      original_gene_data@variable_info$variable_id
    )],
    stringsAsFactors = FALSE
  )
  
  write.csv(
    gene_info,
    file = paste0("clustering_results/cluster_", cluster_id, "_genes.csv"),
    row.names = FALSE
  )
}

# Create final cluster information
cluster_info <- data.frame(
  gene_id = rownames(grouped_matrix_scaled),
  cluster = c$cluster,
  stringsAsFactors = FALSE
)

cluster_info <- cbind(
  cluster_info,
  as.data.frame(c$membership, 
                stringsAsFactors = FALSE)
)
colnames(cluster_info)[3:ncol(cluster_info)] <- paste0("membership_cluster", 1:cluster_number)

cluster_info$max_membership <- apply(c$membership, 1, max)

# Add gene annotations using S4 object
cluster_info$SYMBOL <- original_gene_data@variable_info$SYMBOL[match(
  cluster_info$gene_id, 
  original_gene_data@variable_info$variable_id
)]
cluster_info$ENTREZID <- original_gene_data@variable_info$ENTREZID[match(
  cluster_info$gene_id, 
  original_gene_data@variable_info$variable_id
)]

# Calculate group means
mean_data <- data.frame(
  gene_id = rownames(expression_matrix),
  Mock_mean = rowMeans(expression_matrix[, grep("Mock", colnames(expression_matrix))]),
  WT_mean = rowMeans(expression_matrix[, grep("WT", colnames(expression_matrix))]),
  BA52_mean = rowMeans(expression_matrix[, grep("BA52", colnames(expression_matrix))])
)

# Combine all information
final_cluster_info <- cluster_info %>%
  left_join(mean_data, by = "gene_id")

# Save results
write.csv(
  final_cluster_info,
  file = "clustering_results/all_clusters_information.csv",
  row.names = FALSE
)

# Generate cluster summary
cluster_summary <- final_cluster_info %>%
  group_by(cluster) %>%
  summarise(
    n_genes = n(),
    mean_membership = mean(max_membership),
    sd_membership = sd(max_membership),
    high_confidence_genes = sum(max_membership > 0.7)
  )

write.csv(
  cluster_summary,
  file = "clustering_results/cluster_summary_statistics.csv",
  row.names = FALSE
)
