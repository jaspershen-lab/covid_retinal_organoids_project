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
load("3_data_analysis/2-data_cleaning/3-lipidome/3-lipidome_all/lipidome_data_combined.RData")
load("3_data_analysis/2-data_cleaning/3-lipidome/3-lipidome_all/lipidome_data_pos_cleaned.RData")
load("3_data_analysis/2-data_cleaning/3-lipidome/3-lipidome_all/lipidome_data_neg_cleaned.RData")

# Create output directory
dir.create(
  "3_data_analysis/7-heatmap_analysis/lipidome_heatmap_analysis/lipidome_data_combined/clustering_results",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/7-heatmap_analysis/lipidome_heatmap_analysis/lipidome_data_combined/clustering_results")



# Screen sample IDs
all_sample_ids <- colnames(lipidome_data_combined@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4",
                   "WT_1", "WT_2", "WT_3", "WT_4",
                   "BA52_1", "BA52_2", "BA52_3", "BA52_4")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)



# Get common lipids
common_lipids <- lipidome_data_combined@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids)) %>%
  na.omit() %>%
  rownames()

common_lipid_info <- lipidome_data_combined@variable_info %>%
  dplyr::filter(variable_id %in% common_lipids)



# Extract and reorder expression matrix
expression_data <- lipidome_data_combined@expression_data[common_lipids, filtered_sample_ids]
unique_compounds <- make.unique(common_lipid_info$Compound.name)
rownames(expression_data) <- unique_compounds


# Scale data
heatmap_matrix <- t(scale(t(expression_data)))

# Select top variable metabolites
lipid_vars <- apply(heatmap_matrix, 1, var)
top_var_lipids <- names(sort(lipid_vars, decreasing = TRUE)[1:min(100, length(lipid_vars))])
top_matrix <- heatmap_matrix[top_var_lipids, ]

# Create color scheme for heatmaps
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)


# Prepare annotation
annotation_col <- data.frame(
  Group = lipidome_data_combined@sample_info$group[match(filtered_sample_ids, lipidome_data_combined@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids



overall_plot <- pheatmap(
  heatmap_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Overall Heatmap of Lipids Changes",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)

ggsave("lipidome_heatmap_overall_combined.pdf", overall_plot, width = 10, height = 18)


plot <- pheatmap(
  top_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Overall Heatmap of Lipids Changes",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)
ggsave("lipidome_heatmap_top100_combined.pdf", plot, width = 10, height = 18)


################################
#########miss pos na and combined
################################
# 1. First, create a mapping between compound names and variable IDs
name_to_id_map <- setNames(
  lipidome_data_combined@variable_info$variable_id,
  lipidome_data_combined@variable_info$Compound.name
)

# 2. Reset the original data filtering without na.omit
filtered_data <- lipidome_data_combined@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids))

# 3. Apply more lenient NA filtering (e.g., allow some NAs)
na_threshold <- 0.5  # Allow up to 50% NAs
keep_lipids <- rowSums(!is.na(filtered_data)) >= (ncol(filtered_data) * (1 - na_threshold))
filtered_data <- filtered_data[keep_lipids, ]

# 4. Scale the filtered data
scaled_data <- t(scale(t(filtered_data)))

# 5. Select top variable lipids by mode
n_per_mode <- 50  # Adjust as needed

# Get variances for non-NA values
row_variances <- apply(scaled_data, 1, function(x) var(x, na.rm = TRUE))

# Select top lipids by mode
top_lipids <- lipidome_data_combined@variable_info %>%
  filter(variable_id %in% rownames(filtered_data)) %>%
  mutate(variance = row_variances[Compound.name]) %>%
  group_by(ion_mode) %>%
  arrange(desc(variance), .by_group = TRUE) %>%
  slice_head(n = n_per_mode) %>%
  pull(variable_id)

# 6. Create final matrix
final_matrix <- scaled_data[top_lipids, ]

# 7. Create the heatmap using compound names for readability
rownames(final_matrix) <- lipidome_data_combined@variable_info$Compound.name[
  match(rownames(final_matrix), lipidome_data_combined@variable_info$variable_id)
]

# Create complete matrix with all lipids
all_matrix <- scaled_data  # 使用已经过NA过滤和标准化的数据

# Update rownames with compound names
rownames(all_matrix) <- make.unique(lipidome_data_combined@variable_info$Compound.name[
  match(rownames(all_matrix), lipidome_data_combined@variable_info$variable_id)
])

# Create complete heatmap
complete_plot <- pheatmap(
  all_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Complete Lipidome Heatmap (All pos/neg)",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA,
  na_col = "grey"
)

ggsave("lipidome_heatmap_overall_combined.pdf", complete_plot, width = 10, height = 30)

# Create the heatmap
overall_plot <- pheatmap(
  final_matrix,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Combined Lipidome Heatmap (Balanced pos/neg)",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA,
  na_col = "grey"  # Color for NA values
)

ggsave("lipidome_heatmap_top100_combined_balanced.pdf", overall_plot, width = 10, height = 18)


# 选择方差最大的top 100个脂质
row_variances <- apply(scaled_data, 1, function(x) var(x, na.rm = TRUE))
top_100_lipids <- names(sort(row_variances, decreasing = TRUE)[1:100])
final_matrix_top100 <- scaled_data[top_100_lipids, ]

# 获取对应的Compound.name
compound_names <- lipidome_data_combined@variable_info$Compound.name[
  match(rownames(final_matrix_top100), lipidome_data_combined@variable_info$variable_id)
]
rownames(final_matrix_top100) <- compound_names

# 创建热图
plot<-pheatmap(
  final_matrix_top100,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Combined Lipidome Heatmap (Top 100 most variable lipids)",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA,
  na_col = "grey"
)

ggsave("lipidome_heatmap_top100_combined.pdf", plot, width = 10, height = 18)

#################pos and neg#######
# Function to create heatmap for each mode
create_mode_heatmap <- function(data, mode) {
  # Get sample IDs excluding BQ11
  sample_ids <- colnames(data@expression_data)
  filtered_ids <- sample_ids[!grepl("BQ11", sample_ids)]
  
  # Define desired order
  desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4",
                     "WT_1", "WT_2", "WT_3", "WT_4",
                     "BA52_1", "BA52_2", "BA52_3", "BA52_4")
  filtered_ids <- intersect(desired_order, filtered_ids)
  
  # Get expression data
  expr_data <- data@expression_data[, filtered_ids]
  
  # Make unique compound names
  compounds <- make.unique(data@variable_info$Compound.name)
  rownames(expr_data) <- compounds
  
  # Scale data
  heatmap_matrix <- t(scale(t(expr_data)))
  
  # Get top variable compounds
  vars <- apply(heatmap_matrix, 1, var)
  top_vars <- names(sort(vars, decreasing = TRUE)[1:min(50, length(vars))])
  top_matrix <- heatmap_matrix[top_vars, ]
  
  # Create annotation
  annotation_col <- data.frame(
    Group = data@sample_info$group[match(filtered_ids, data@sample_info$sample_id)],
    stringsAsFactors = FALSE
  )
  rownames(annotation_col) <- filtered_ids
  
  # Color scheme
  colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
  
  # Create heatmaps
  # Overall heatmap
  overall <- pheatmap(
    heatmap_matrix,
    color = colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    main = paste("Overall Heatmap -", mode, "mode"),
    fontsize_row = 6,
    fontsize_col = 10,
    cellwidth = 20,
    cellheight = 8,
    border_color = NA
  )
  
  # Top 100 heatmap
  top50 <- pheatmap(
    top_matrix,
    color = colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    main = paste("Top 100 Variable Lipids -", mode, "mode"),
    fontsize_row = 6,
    fontsize_col = 10,
    cellwidth = 20,
    cellheight = 8,
    border_color = NA
  )
  
  # Save plots
  ggsave(paste0("lipidome_heatmap_overall_", mode, ".pdf"), 
         overall, width = 10, height = 18)
  ggsave(paste0("lipidome_heatmap_top100_", mode, ".pdf"), 
         top50, width = 10, height = 18)
}

# Create heatmaps for both modes
create_mode_heatmap(lipidome_data_combined_pos, "positive")
create_mode_heatmap(lipidome_data_combined_neg, "negative")


#######
# Get top 100 variable lipids matrix
filtered_sample_ids <- colnames(lipidome_data_combined_combined@expression_data)[!grepl("BQ11", colnames(lipidome_data_combined_combined@expression_data))]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4",
                   "WT_1", "WT_2", "WT_3", "WT_4",
                   "BA52_1", "BA52_2", "BA52_3", "BA52_4")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)

# Get common lipids without NA
common_lipids <- lipidome_data_combined_combined@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids)) %>%
  na.omit() %>%
  rownames()

# Get expression data
expression_data <- lipidome_data_combined_combined@expression_data[common_lipids, filtered_sample_ids]
unique_compounds <- make.unique(lipidome_data_combined_combined@variable_info$Compound.name[match(common_lipids, lipidome_data_combined_combined@variable_info$variable_id)])
rownames(expression_data) <- unique_compounds

# Scale data
heatmap_matrix <- t(scale(t(expression_data)))

# Get top 100 variable lipids
lipid_vars <- apply(heatmap_matrix, 1, var)
top_var_lipids <- names(sort(lipid_vars, decreasing = TRUE)[1:100])
top_matrix <- heatmap_matrix[top_var_lipids, ]

# Calculate mean expression for each group
ctrl_mean <- rowMeans(top_matrix[, grep("Ctrl", colnames(top_matrix))])
wt_mean <- rowMeans(top_matrix[, grep("WT", colnames(top_matrix))])
ba52_mean <- rowMeans(top_matrix[, grep("BA52", colnames(top_matrix))])

# Find downregulated lipids
down_lipids <- names(which(wt_mean < ctrl_mean & ba52_mean < ctrl_mean))
down_matrix <- top_matrix[down_lipids, ]

# Create annotation
annotation_col <- data.frame(
  Group = lipidome_data_combined_combined@sample_info$group[match(filtered_sample_ids, lipidome_data_combined_combined@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

# Create heatmap
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
plot <- pheatmap(
  down_matrix,
  color = colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Downregulated Lipids in WT and BA52",
  fontsize_row = 8,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 12,
  border_color = NA
)
plot

ggsave("lipidome_heatmap_combined_downregulated.pdf", 
       plot, 
       width = 12, 
       height = length(down_lipids) * 0.25 + 3)




# Get non-downregulated lipids
down_lipids <- names(which(wt_mean < ctrl_mean & ba52_mean < ctrl_mean))
non_down_lipids <- setdiff(top_var_lipids, down_lipids)

# Prepare data for circular heatmap
mat_for_circle <- top_matrix[non_down_lipids, ]

head(mat_for_circle)



library(ComplexHeatmap)
library(circlize)

# 将数据转换为矩阵并标准化
mat_scaled <- t(scale(t(mat_for_circle)))

# 设置颜色
# 更改为红色系配色
# 更改为粉红色系配色
col_fun = colorRamp2(c(-2, 0, 2), c("#393781", "white", "#f22942"))

pdf("heatmap_circle.pdf", width=15, height=15)

circos.par(gap.after = 15,  # Reduced gap
           start.degree = 90,
           track.height = 0.15,  # Increased height
           cell.padding = c(0, 0, 0, 0),
           points.overflow.warning = FALSE)

circos.heatmap(mat_scaled,
               col = col_fun,
               rownames.side = "outside",
               rownames.cex = 0.4,
               cluster = TRUE,
               dend.side = "inside",
               bg.border = NA,
               cell.border = NA,
               track.margin = c(0, 0))  # Removed margins

# 添加图例
lg = Legend(title = "Exp", 
            col_fun = col_fun,
            direction = "vertical",
            title_position = "topcenter")

draw(lg, x = unit(0.75, "npc"), y = unit(0.65, "npc"))

# 添加样本标签 
circos.track(track.index = get.current.track.index(), 
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 1) {
                 cn = colnames(mat_scaled)
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2] + 0.5, n),
                             (1:n) * 2.95 + 6,
                             cn,
                             cex = 0.9,
                             adj = c(0, 0.5),
                             facing = "inside")
               }
             }, bg.border = NA)

circos.clear()
# 绘图代码
dev.off()



#################### Lipidome Clustering Analysis ####################
library(tidymass)
library(dplyr)
library(Mfuzz)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(pheatmap)



######################## 1. Data Preprocessing ########################
# Screen sample IDs
all_sample_ids <- colnames(lipidome_data_combined_combined@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4",
                   "WT_1", "WT_2", "WT_3", "WT_4",
                   "BA52_1", "BA52_2", "BA52_3", "BA52_4")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids)

# Extract expression matrix
expression_matrix <- lipidome_data_combined@expression_data[, filtered_sample_ids]

# Handle missing values
expression_matrix <- t(apply(expression_matrix, 1, function(x) {
  if(any(is.na(x))) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
  }
  return(x)
}))

# Log transform if needed
if(max(expression_matrix, na.rm = TRUE) > 100) {
  expression_matrix <- log2(expression_matrix + 1)
}

# Calculate group means
grouped_matrix <- matrix(0, nrow = nrow(expression_matrix), ncol = 3)
colnames(grouped_matrix) <- c("Ctrl", "WT", "BA52")
rownames(grouped_matrix) <- rownames(expression_matrix)

grouped_matrix[, "Ctrl"] <- rowMeans(expression_matrix[, grep("Ctrl", colnames(expression_matrix))])
grouped_matrix[, "WT"] <- rowMeans(expression_matrix[, grep("WT", colnames(expression_matrix))])
grouped_matrix[, "BA52"] <- rowMeans(expression_matrix[, grep("BA52", colnames(expression_matrix))])

# Scale data
grouped_matrix_scaled <- t(scale(t(grouped_matrix)))


######################## 2. Prepare Mfuzz Input ########################
temp_data <- rbind(
  time = 1:3,
  grouped_matrix_scaled
)

write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

data <- table2eset(filename = "temp_data.txt")

######################## 3. Clustering Analysis ########################
m1 <- mestimate(data)
print(paste("Estimated m parameter:", m1))

dmin_results <- Dmin(data, m = m1, crange = seq(2, 20, 2), repeats = 5, visu = TRUE)
ggsave("cluster_number_optimization.pdf", width = 8, height = 6)

cluster_number <- 6
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

pdf("cluster_correlations.pdf", width = 10, height = 10)
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

for(cluster_id in 1:cluster_number) {
  cluster_members <- which(c$cluster == cluster_id)
  cluster_data <- grouped_matrix_scaled[cluster_members, ]
  
  plot_data <- data.frame()
  for(i in 1:nrow(cluster_data)) {
    gene_data <- data.frame(
      Lipid = paste0("Lipid", i),
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
              aes(x = Time, y = Value, group = Lipid,  
                  color = Membership),
              size = 0.8) +
    geom_line(data = mean_profile,
              aes(x = Time, y = Value, group = 1),
              color = "black",  
              size = 1.2) +
    scale_color_gradientn(
      colors = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                 "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"),   
      limits = c(0.2, 0.8),
      oob = scales::squish,    
      breaks = seq(0.2, 0.8, by = 0.2),
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
  
  ggsave(paste0("cluster_", cluster_id, "_profile.pdf"),
         p, width = 10, height = 7)
  

  # Save cluster lipid information
  cluster_df <- data.frame(
    lipid_id = rownames(expression_matrix)[cluster_members],
    lipid_name = lipidome_data_combined@variable_info$Compound.name[
      match(rownames(expression_matrix)[cluster_members], 
            lipidome_data_combined@variable_info$variable_id)
    ],
    membership = c$membership[cluster_members, cluster_id],
    stringsAsFactors = FALSE
  )
  
  write.csv(
    cluster_df,
    file = paste0("cluster_", cluster_id, "_lipids.csv"),
    row.names = FALSE
  )
}

# Print summary statistics
print("Clustering Analysis Complete!")
print("Number of lipids in each cluster:")
print(table(c$cluster))
print("\nMembership summary:")
print(summary(apply(c$membership, 1, max)))

