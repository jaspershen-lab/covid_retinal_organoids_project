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

# Load data
load("3_data_analysis/2-data_cleaning/2-metabolome/metabolome_data.rda")

# Create output directory
dir.create(
  "3_data_analysis/7-heatmap_analysis/metabolome_heatmap_analysis/metabolome_data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/7-heatmap_analysis/metabolome_heatmap_analysis/metabolome_data")

# Screen sample IDs
all_sample_ids <- colnames(metabolome_data@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4",
                   "WT_1", "WT_2", "WT_3", "WT_4",
                   "BA52_1", "BA52_2", "BA52_3", "BA52_4")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)

# Get common metabolites
common_metabolites <- metabolome_data@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids)) %>%
  na.omit() %>%
  rownames()

common_metabolite_info <- metabolome_data@variable_info %>%
  dplyr::filter(variable_id %in% common_metabolites)

# Extract and reorder expression matrix
expression_data <- metabolome_data@expression_data[common_metabolites, filtered_sample_ids]
unique_compounds <- make.unique(common_metabolite_info$Compound.name)
rownames(expression_data) <- unique_compounds

# Scale data
heatmap_matrix <- t(scale(t(expression_data)))

# Select top variable metabolites
metabolite_vars <- apply(heatmap_matrix, 1, var)
top_var_metabolites <- names(sort(metabolite_vars, decreasing = TRUE)[1:min(100, length(metabolite_vars))])
top_matrix <- heatmap_matrix[top_var_metabolites, ]

# Split POS and NEG metabolites
pos_metabolites <- grep("_POS", rownames(top_matrix), value = TRUE)
neg_metabolites <- grep("_NEG", rownames(top_matrix), value = TRUE)

pos_matrix <- top_matrix[pos_metabolites, ]
neg_matrix <- top_matrix[neg_metabolites, ]

# Create color scheme for heatmaps
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Prepare annotation
annotation_col <- data.frame(
  Group = metabolome_data@sample_info$group[match(filtered_sample_ids, metabolome_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

# Create color schemes
pos_colors <- colorRampPalette(c("white", "firebrick3"))(50)
neg_colors <- colorRampPalette(c("navy", "white"))(50)

# Plot overall heatmap
overall_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
overall_plot <- pheatmap(
  top_matrix,
  color = overall_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Overall Heatmap of Metabolite Changes",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)
ggsave("metabolome_heatmap_overall.pdf", overall_plot, width = 10, height = 18)

# Plot positive heatmap
pos_plot <- pheatmap(
  pos_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Heatmap of POS Mode Metabolites",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)
ggsave("metabolome_heatmap_pos_mode.pdf", pos_plot, width = 10, height = 18)

# Plot NEG metabolites heatmap
neg_plot <- pheatmap(
  neg_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Heatmap of NEG Mode Metabolites",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)
ggsave("metabolome_heatmap_neg_mode.pdf", neg_plot, width = 10, height = 18)


#################### Metabolome Clustering Analysis ####################
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
all_sample_ids <- colnames(metabolome_data@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4",
                   "WT_1", "WT_2", "WT_3", "WT_4",
                   "BA52_1", "BA52_2", "BA52_3", "BA52_4")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)

# Extract expression matrix
expression_matrix <- metabolome_data@expression_data[, filtered_sample_ids]

# Handle missing values
expression_matrix <- t(apply(expression_matrix, 1, function(x) {
  if(any(is.na(x))) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
  }
  return(x)
}))

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
m1 <- mestimate(data)
print(paste("Estimated m parameter:", m1))

# Optimize cluster number
dmin_results <- Dmin(data, m = m1, crange = seq(2, 20, 2), repeats = 5, visu = TRUE)
ggsave("clustering_results/cluster_number_optimization.pdf", width = 8, height = 6)

# Perform clustering
cluster_number <- 4 
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
    expr_data <- exprs(data)[valid_idx, , drop = FALSE]  # 添加 drop = FALSE
    centers[i,] <- colSums(weights * expr_data) / sum(weights)
  }
}

rownames(centers) <- paste0("Cluster", 1:n_clusters)
colnames(centers) <- colnames(exprs(data))

# Plot correlation heatmap
# 计算聚类中心时添加检查
centers <- matrix(0, nrow = n_clusters, ncol = n_timepoints)
rownames(centers) <- paste0("Cluster", 1:n_clusters)
colnames(centers) <- colnames(exprs(data))

for(i in 1:n_clusters) {
  memb <- c$membership[, i]
  valid_idx <- memb >= 0.5
  if(sum(valid_idx) > 0) {
    weights <- memb[valid_idx]
    expr_data <- exprs(data)[valid_idx, , drop = FALSE]
    centers[i,] <- colSums(weights * expr_data) / sum(weights)
  }
}

# 检查和处理NA值
if(any(is.na(centers))) {
  centers <- na.omit(centers)
}

# 绘制相关性热图
pdf("clustering_results/cluster_correlations.pdf", width = 10, height = 10)
if(nrow(centers) > 1) {
  cor_matrix <- cor(t(centers))
  corrplot(
    cor_matrix,
    type = "full",
    diag = TRUE,
    order = "original",  # 改用original排序
    col = colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))(100),
    number.cex = 0.7,
    addCoef.col = "black"
  )
}
dev.off()

# 可视化聚类结果
pdf("clustering_results/cluster_patterns.pdf", width = 10, height = 8)
mfuzz.plot(
  data,
  cl = c,
  mfrow = c(2,2),
  time.labels = c("Ctrl", "WT", "BA52")
)
dev.off()


# 4.2 Cluster Profiles
membership_cutoff <- 0.5

mfuzz.plot(
  eset = data,
  min.mem = membership_cutoff,
  cl = c,
  mfrow = c(2, 2),
  time.labels = 1:3,
  new.window = FALSE
)

# Cluster profile visualization with correct order
# Plot cluster profiles
for(cluster_id in 1:cluster_number) {
  cluster_members <- which(c$cluster == cluster_id)
  cluster_data <- grouped_matrix_scaled[cluster_members, ]
  
  plot_data <- data.frame()
  for(i in 1:nrow(cluster_data)) {
    metabolite_data <- data.frame(
      Metabolite = paste0("Metabolite", i),
      Time = colnames(cluster_data),
      Value = as.numeric(cluster_data[i, ]),
      Membership = c$membership[cluster_members[i], cluster_id]
    )
    plot_data <- rbind(plot_data, metabolite_data)
  }
  
  plot_data$Time <- factor(plot_data$Time, levels = c("Ctrl", "WT", "BA52"))
  mean_profile <- aggregate(Value ~ Time, data = plot_data, FUN = mean)
  mean_profile$Time <- factor(mean_profile$Time, levels = c("Ctrl", "WT", "BA52"))
  
  p <- ggplot() +
    geom_line(data = plot_data, 
              aes(x = Time, y = Value, group = Metabolite,  
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
}
  
######################## 保存聚类结果 ########################
# 1. 保存各聚类中心和相关性
write.csv(centers, "clustering_results/cluster_centers.csv")
write.csv(cor_matrix, "clustering_results/cluster_correlations.csv") 

# 2. 保存每个聚类的代谢物信息
for(cluster_id in 1:cluster_number) {
  cluster_members <- which(c$cluster == cluster_id)
  
  # 代谢物信息
  metabolite_info <- data.frame(
    metabolite_id = rownames(grouped_matrix_scaled)[cluster_members],
    membership = c$membership[cluster_members, cluster_id],
    cluster = cluster_id,
    expression_ctrl = grouped_matrix[cluster_members, "Ctrl"],
    expression_wt = grouped_matrix[cluster_members, "WT"], 
    expression_ba52 = grouped_matrix[cluster_members, "BA52"],
    stringsAsFactors = FALSE
  ) %>%
    left_join(metabolome_data@variable_info, 
              by = c("metabolite_id" = "variable_id")) %>%
    arrange(desc(membership))
  
  write.csv(metabolite_info,
            file = paste0("clustering_results/cluster_", 
                          cluster_id, "_metabolites.csv"),
            row.names = FALSE)
}

# 3. 保存聚类总结信息
cluster_summary <- data.frame(
  cluster = 1:cluster_number,
  size = sapply(1:cluster_number, function(i) sum(c$cluster == i)),
  mean_membership = sapply(1:cluster_number, function(i) 
    mean(c$membership[c$cluster == i, i]))
)
write.csv(cluster_summary, 
          "clustering_results/cluster_summary.csv",
          row.names = FALSE)

# 4. 保存分析设置
analysis_params <- list(
  cluster_number = cluster_number,
  m_parameter = m1,
  membership_cutoff = 0.5,
  data_scaling = "z-score"
)
saveRDS(analysis_params, "clustering_results/analysis_params.rds")

# 5. 保存完整聚类对象
save(c, file = "clustering_results/mfuzz_cluster_object.RData")

