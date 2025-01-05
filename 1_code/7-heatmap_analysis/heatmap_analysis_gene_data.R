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




library(Mfuzz)

# Calculate mean expression for each group
calculate_group_means <- function(expression_data, sample_groups) {
  # Get sample indices for each group
  ctrl_idx <- grep("Mock", colnames(expression_data))
  wt_idx <- grep("WT", colnames(expression_data))
  ba52_idx <- grep("BA52", colnames(expression_data))
  
  # Calculate means
  group_means <- data.frame(
    "1" = rowMeans(expression_data[, ctrl_idx, drop = FALSE]),
    "2" = rowMeans(expression_data[, wt_idx, drop = FALSE]),
    "3" = rowMeans(expression_data[, ba52_idx, drop = FALSE])
  )
  
  return(group_means)
}

# Use all metabolites for clustering
temp_data <- calculate_group_means(heatmap_matrix, filtered_sample_ids)
expression_data <- temp_data

# Add time points
time <- colnames(temp_data)
temp_data <- rbind(time, temp_data)
row.names(temp_data)[1] <- "time"

# Write to file for Mfuzz input
write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

# Create expression set and standardize
data <- table2eset(filename = "temp_data.txt")
data.s <- standardise(data)
m1 <- mestimate(data.s)

# Get mfuzz center function
get_mfuzz_center <- function(data, c, membership_cutoff) {
  centers <- c$centers
  membership <- c$membership
  
  for (i in 1:nrow(centers)) {
    cluster_members <- which(c$cluster == i & 
                               apply(membership, 1, max) >= membership_cutoff)
    
    if (length(cluster_members) > 0) {
      centers[i, ] <- colMeans(exprs(data)[cluster_members, , drop = FALSE])
    }
  }
  
  return(centers)
}

# Perform clustering with fixed number
cluster_number <- 3
c <- mfuzz(data.s, c = cluster_number, m = m1)
save(c, file = "c")

# Analyze cluster correlations
membership_cutoff <- 0.5
center <- get_mfuzz_center(data = data.s,
                           c = c,
                           membership_cutoff = 0.5)

rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

# Plot correlation matrix
corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

# Plot clusters
mfuzz.plot(
  eset = data.s,
  min.mem = 0.5,
  cl = c,
  mfrow = c(2, 3),
  time.labels = time,
  new.window = FALSE
)

# Generate cluster information
cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)


# Plot individual clusters
for (idx in 1:cluster_number) {
  cat("Processing cluster", idx, "\n")
  
  cluster_data <-
    cluster_info %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  
  # Print cluster information for debugging
  print(paste("Number of metabolites in cluster", idx, ":", nrow(cluster_data)))
  
  path <- paste("cluster", idx, sep = "_")
  dir.create(path, showWarnings = FALSE)
  
  # Save cluster members
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  # Get center data
  temp_center <- data.frame(
    time = 1:3,  # explicitly set time points
    value = as.numeric(center[idx, ])
  )
  
  # Print center data for debugging
  print("Center values:")
  print(temp_center)
  
  # Get individual metabolite data
  if(nrow(cluster_data) > 0) {
    temp <- data.frame(
      time = rep(1:3, each = nrow(cluster_data)),
      value = as.vector(t(expression_data[cluster_data$variable_id, ])),
      metabolite = rep(cluster_data$variable_id, times = 3),
      membership = rep(cluster_data$membership, times = 3)
    )
    
    # Print first few rows of temp for debugging
    print("Sample of metabolite data:")
    print(head(temp))
    
    # Create plot with both center line and individual metabolites
    plot <- ggplot() +
      # Individual metabolite lines
      geom_line(data = temp,
                aes(x = time, y = value, group = metabolite),
                alpha = 0.7,
                color = "grey60") +
      # Center line
      geom_line(data = temp_center,
                aes(x = time, y = value),
                size = 1.5,
                color = "red") +
      # Add points to make the values more visible
      geom_point(data = temp,
                 aes(x = time, y = value),
                 alpha = 0.5,
                 size = 2) +
      # Customize the plot
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12)
      ) +
      labs(
        x = "Time point",
        y = "Z-score",
        title = paste("Cluster", idx, "(", nrow(cluster_data), "metabolites)")
      ) +
      scale_x_continuous(breaks = 1:3, 
                         labels = c("Ctrl", "WT", "BA52")) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")
    
    # Display plot
    print(plot)
    
    # Save plot
    ggsave(
      plot,
      filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
      width = 8,
      height = 6
    )
  } else {
    warning(paste("No metabolites in cluster", idx, "with membership >", membership_cutoff))
  }
}

# Save final cluster information
final_cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
    temp <-
      temp %>%
      dplyr::filter(membership >= membership_cutoff) %>%
      dplyr::mutate(cluster_raw = cluster,
                    cluster = x)
    temp
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

save(final_cluster_info, file = "final_cluster_info")
