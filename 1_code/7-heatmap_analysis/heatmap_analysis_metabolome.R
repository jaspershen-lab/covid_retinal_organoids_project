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

##read data
load("3_data_analysis/2-data_cleaning/2-metabolome/metabolome_data.RData")

dir.create(
  "3_data_analysis/7-heatmap_analysis/metabolome_heatmap_analysis/metabolome_data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/7-heatmap_analysis/metabolome_heatmap_analysis/metabolome_data")

# Screening sample ID
all_sample_ids <- colnames(metabolome_data@expression_data)
filtered_sample_ids <- all_sample_ids[!grepl("BQ11", all_sample_ids)]  
desired_order <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "WT_1", "WT_2", "WT_3", "BA52_1", "BA52_2", "BA52_3")
filtered_sample_ids <- intersect(desired_order, filtered_sample_ids) 

# common metabolites
common_metabolites <- metabolome_data@expression_data %>%
  dplyr::select(all_of(filtered_sample_ids)) %>%
  na.omit() %>%
  rownames()

common_metabolite_info <- metabolome_data@variable_info %>%
  dplyr::filter(variable_id %in% common_metabolites)

# Extract the expression matrix and reorder it
expression_data <- metabolome_data@expression_data[common_metabolites, filtered_sample_ids]

# Make compound names unique
unique_compounds <- make.unique(common_metabolite_info$Compound.name)
rownames(expression_data) <- unique_compounds



# Initial heatmap of all metabolites
heatmap_matrix <- t(scale(t(expression_data)))  # Scale by row (Z-score)

# Select top 100 variable metabolites
metabolite_vars <- apply(heatmap_matrix, 1, var)
top_var_metabolites <- names(sort(metabolite_vars, decreasing = TRUE)[1:min(100, length(metabolite_vars))])

# Set heat map colors
color_scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Add annotation
annotation_col <- data.frame(
  Group = metabolome_data@sample_info$group[match(filtered_sample_ids, metabolome_data@sample_info$sample_id)],
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- filtered_sample_ids

# Plot heatmap
plot <- pheatmap(
  heatmap_matrix,
  color = color_scheme,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  main = "Heatmap of Metabolites",
  fontsize_row = 6,
  fontsize_col = 10,
  cellwidth = 20,
  cellheight = 8,
  border_color = NA
)
plot

ggsave("metabolome_heatmap.pdf", plot, width = 10, height = 18)



#######clustering analysis
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

# Use the top metabolites from previous analysis
temp_data <- calculate_group_means(heatmap_matrix[top_var_metabolites, ], filtered_sample_ids)
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

# Perform clustering
cluster_number <- 2  # can adjust this number based on data
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

