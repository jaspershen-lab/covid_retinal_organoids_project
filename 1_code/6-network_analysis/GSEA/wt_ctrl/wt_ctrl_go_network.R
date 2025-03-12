library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load data
load("3_data_analysis/4_differential_analysis/1-transcriptome/transcriptome_data")
data <- transcriptome_data@variable_info[, c("SYMBOL", "fc_wt_ctrl", "p_value_adjust_wt_ctrl")]
go_results <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/GO/GSEA_go_wt_ctrl.csv")

# Create output directory
output_dir <- "3_data_analysis/6-network_analysis/GSEA/wt_ctrl/wt_ctrl_go"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Select top pathways
top_pathways <- go_results %>%
  as.data.frame() %>%
  mutate(type = if_else(NES > 0, "positive", "negative")) %>%
  group_by(type) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(Description)

print(top_pathways)

# Process edge data
edge_data <- go_results %>%
  filter(Description %in% top_pathways) %>%
  mutate(genes = strsplit(core_enrichment, "/")) %>%
  unnest(genes) %>%
  group_by(Description) %>%
  summarise(
    genes_symbol = list(bitr(genes, 
                             fromType = "ENTREZID",
                             toType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db)$SYMBOL)
  ) %>%
  unnest(genes_symbol) %>%
  rename(from = genes_symbol) %>%
  mutate(
    to = Description,
    Correlation = go_results$NES[match(Description, go_results$Description)],
    fdr = -log10(go_results$pvalue[match(Description, go_results$Description)])
  )

# Process gene data
data_id <- bitr(data$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
data_all <- merge(data, data_id, by = "SYMBOL") %>%
  arrange(desc(fc_wt_ctrl))

# Get gene information
gene_info <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = unique(edge_data$from),
                                   columns = c("SYMBOL", "GENENAME"),
                                   keytype = "SYMBOL")

# Create node data
gene_nodes <- tibble(
  node = unique(edge_data$from)
) %>%
  mutate(Class = "Gene") %>%
  left_join(gene_info, by = c("node" = "SYMBOL")) %>%
  left_join(
    tibble(
      node = data_all$SYMBOL,
      size = -log10(data_all$p_value_adjust_wt_ctrl),
      color = log2(data_all$fc_wt_ctrl)
    ),
    by = "node"
  ) %>%
  dplyr::select(node, Class, size, color, GENENAME)

pathway_nodes <- data.frame(
  node = unique(edge_data$to),
  Class = "Pathway",
  size = -log10(go_results$pvalue[match(unique(edge_data$to), go_results$Description)]),
  color = go_results$NES[match(unique(edge_data$to), go_results$Description)],
  GENENAME = unique(edge_data$to),
  stringsAsFactors = FALSE
)

# Combine node data
node_data <- rbind(gene_nodes, pathway_nodes) %>%
  arrange(Class) %>%
  mutate(radius = if_else(Class == "Pathway", 0.7, 1))

# Create network graph
temp_data <- tidygraph::tbl_graph(
  nodes = node_data,
  edges = edge_data,
  directed = TRUE
) %>%
  mutate(Degree = centrality_degree(mode = 'all'))

# Create layout
g <- temp_data
V(g)$type <- bipartite_mapping(g)$type

# Use bipartite layout
coords <- layout_as_bipartite(g) %>%
  as.data.frame()
colnames(coords) <- c("x", "y")

coords$index <- 1:nrow(coords)
coords <- coords %>% 
  dplyr::arrange(x)

coords$y[coords$y == 0] <- 0.3

# Convert to a circular layout
coords <- coords %>% 
  dplyr::arrange(index) %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(theta = x / (max(x) + 1) * 2 * pi,
                r = y + 1,
                x = r * cos(theta),
                y = r * sin(theta))

# Create visualization
plot <- ggraph(temp_data,
               layout = 'manual',
               x = coords$x,
               y = coords$y) +
  geom_edge_diagonal(aes(color = Correlation),
                     alpha = 0.6,
                     strength = 1.0,
                     show.legend = TRUE) +
  geom_node_point(aes(fill = color,
                      size = size),
                  shape = 21,
                  alpha = 1,
                  show.legend = TRUE) +
  geom_node_text(
    aes(
      x = x * 1.03,
      y = y * 1.03,
      label = GENENAME,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      hjust = ifelse(Class == "Pathway", "outward", "outward"),
      colour = Class
    ),
    size = ifelse(V(temp_data)$Class == "Pathway", 0.5, 2.5),
    check_overlap = TRUE,
    show.legend = FALSE,
  ) +
  scale_fill_gradient2(
    low = "#00468BFF",
    mid = "white", 
    high = "#ED0000FF",
    midpoint = 0
  ) +
  scale_edge_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#E64B35FF",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(1, 7)) +
  guides(
    edge_color = guide_edge_colorbar(title = "NES"),
    fill = guide_colorbar(title = "Gene: log2(FC)\nPathway: NES"),
    size = guide_legend(title = "Gene & Pathway:\n-log10(P)")
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  ) +
  coord_fixed()
plot


# Save output
ggsave(
  plot,
  filename = "top10_NES_pathway_network_wt_ctrl_go.pdf",
  width = 8.5,
  height = 7,
  bg = "transparent",
  dpi = 300
)

ggsave(
  plot,
  filename = "top10_NES_pathway_network_wt_ctrl_go.png",
  width = 8.5,
  height = 7,
  bg = "transparent",
  dpi = 300
)

# Output pathway information
pathway_info <- go_results %>%
  as.data.frame() %>%
  filter(Description %in% top_pathways) %>%
  dplyr::select(Description, NES, pvalue, p.adjust) %>%  
  arrange(desc(abs(NES)))

print(pathway_info)
