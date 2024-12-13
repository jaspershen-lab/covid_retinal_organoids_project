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

