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

##read data
load("3_data_analysis/2-data_cleaning/1-transcriptome/transcriptome_data.RData")
load("2_data/pathway_human/pathway_GO.rda")
load("2_data/pathway_human/pathway_kegg.rda")
load("2_data/pathway_human/pathway_Reactome.rda")


dir.create(
  "3_data_analysis/7-heatmap_analysis",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/7-heatmap_analysis")



# Screening sample ID
all_sample_ids <- colnames(transcriptome_data@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]  
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "WT_1", "WT_2", "WT_3", "BA52_1", "BA52_2", "BA52_3")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids) 



# common_genes
common_genes <- transcriptome_data@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids)) %>%
  na.omit() %>%
  rownames()

common_gene_info <- transcriptome_data@variable_info %>%
  dplyr::filter(variable_id %in% common_genes)

# Extract the representation matrix and reorder it
expression_data <- transcriptome_data@expression_data[common_genes, filtered_sample_ids]

unique_ensembl <- make.unique(common_gene_info$ENSEMBL)

rownames(expression_data) <- unique_ensembl

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

# Screening pathway gene
common_ensembl_genes <- rownames(expression_data)

go_overlap <- calculate_overlap(go_database, common_ensembl_genes)
kegg_overlap <- calculate_overlap(kegg_database, common_ensembl_genes)
reactome_overlap <- calculate_overlap(reactome_database, common_ensembl_genes)

filter_pathways <- function(overlap_results, threshold = 50) {
  Filter(function(x) x$overlap_percentage >= threshold, overlap_results)
}

go_filtered <- filter_pathways(go_overlap, 50)
kegg_filtered <- filter_pathways(kegg_overlap, 50)
reactome_filtered <- filter_pathways(reactome_overlap, 50)

# Combined screened pathway genes
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


# Extract the heat map matrix
heatmap_matrix <- expression_data[all_pathway_genes, ]
heatmap_matrix <- t(scale(t(heatmap_matrix)))  # 按行标准化

# rownames SYMBOL
rownames(heatmap_matrix) <- common_gene_info$SYMBOL[match(rownames(heatmap_matrix), common_gene_info$ENSEMBL)]
rownames(heatmap_matrix) <- make.unique(rownames(heatmap_matrix))

#  Top 100 gene
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


#######clustering
library(Mfuzz)

# Calculate mean expression for each group
calculate_group_means <- function(expression_data, sample_groups) {
  # Get sample indices for each group
  ctrl_idx <- grep("Ctrl", colnames(expression_data))
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

# Use the top genes from previous analysis
temp_data <- calculate_group_means(heatmap_matrix[top_var_genes, ], filtered_sample_ids)
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
data.s <- data
m1 <- mestimate(data.s)

# Determine optimal cluster number
plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 40, 2),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 40, 2)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  geom_segment(aes(
    x = k,
    y = 0,
    xend = k,
    yend = distance
  )) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "Cluster number",
       y = "Min. centroid distance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot,
       filename = "distance_k_number.pdf",
       width = 7,
       height = 7)

# Perform clustering
cluster_number <- 9
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
  mfrow = c(3, 4),
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
  cat(idx, " ")
  
  cluster_data <-
    cluster_info %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  temp_center <-
    center[idx, , drop = TRUE] %>%
    unlist() %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = as.numeric(time))
  
  temp <-
    expression_data[cluster_data$variable_id,] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time = as.numeric(time))
  
  if ("class" %in% colnames(variable_info)) {
    temp <- temp %>%
      dplyr::left_join(variable_info[, c("variable_id", "class")],
                       by = "variable_id")
  }
  
  plot <-
    temp %>%
    dplyr::arrange(membership, variable_id) %>%
    {if ("class" %in% colnames(.)) dplyr::arrange(., desc(class)) else .} %>%
    dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    {if ("class" %in% colnames(temp)) 
      geom_line(aes(color = class), alpha = 0.7)
      else 
        geom_line(alpha = 0.7)} +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.justification = c(0, 1),
      panel.grid = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste("Cluster ",
                    idx,
                    " (",
                    nrow(cluster_data),
                    " molecules)",
                    sep = "")
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    ) +
    geom_hline(yintercept = 0) +
    {if ("class" %in% colnames(temp) && exists("omics_color")) 
      scale_color_manual(values = omics_color)
      else 
        scale_color_viridis_d()}
  
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
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

# neuro_heatmap_matrix
neuro_heatmap_matrix <- expression_data[neuro_pathway_genes, ]
neuro_heatmap_matrix <- t(scale(t(neuro_heatmap_matrix)))  # 按行标准化

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

