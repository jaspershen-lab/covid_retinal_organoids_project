library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(wesanderson)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggridges)

##read data wt_ctrl
load("3_data_analysis/4-differential_analysis/1-transcriptome/transcriptome_data")
data <- transcriptome_data@variable_info[, c("SYMBOL", "fc_wt_ctrl")]
go_results<- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/GO/GSEA_GO_wt_ctrl.csv")

#####
if (!dir.exists("3_data_analysis/6-network_analysis/GSEA/wt_ctrl/wt_ctrl_go")) {
  dir.create(
    "3_data_analysis/6-network_analysis/GSEA/wt_ctrl/wt_ctrl_go",
    recursive = TRUE,
    showWarnings = FALSE
  )
}


setwd("3_data_analysis/6-network_analysis/GSEA/wt_ctrl/wt_ctrl_go")

#####GO
# Top 10 pathways
top_pathways <- go_results%>%
  as.data.frame() %>%
  dplyr::mutate(abs_NES = abs(NES)) %>%
  dplyr::arrange(desc(abs_NES)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(ID)

print(top_pathways)

# edge_data
  edge_data <- go_results %>%
    filter(ID %in% top_pathways) %>%
    mutate(genes = strsplit(core_enrichment, "/")) %>%
    unnest(genes) %>%
    group_by(ID) %>%
    summarise(
      genes_symbol = list(bitr(genes, 
                               fromType = "ENTREZID",
                               toType = "SYMBOL", 
                               OrgDb = org.Hs.eg.db)$SYMBOL)
    ) %>%
    unnest(genes_symbol) %>%
    rename(from = genes_symbol) %>%
    mutate(
      to = ID,
      Correlation = go_results$NES[match(ID, go_results$ID)],
      fdr = -log10(go_results$pvalue[match(ID, go_results$ID)])
    )

# node_data
# Use bitr to convert gene ids and merge
data_id <- bitr(data$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
data_all <- merge(data, data_id, by = "SYMBOL")
head(data_all)

data_all <-
  as.data.frame(data_all) %>%
  dplyr::arrange(desc(fc_wt_ctrl))

##gene_nodes
gene_nodes <- tibble(
  node = unique(edge_data$from)
) %>%
  mutate(Class = "Gene") %>%
  left_join(
    tibble(
      node = data_all$SYMBOL,
      size = abs(log2(data_all$fc_wt_ctrl)),
      color = log2(data_all$fc_wt_ctrl)
    ),
    by = "node"
  )


##pathway_node
pathway_nodes <- data.frame(
  node = unique(edge_data$to),
  Class = "Pathway",
  size = -log10(go_results$pvalue[match(unique(edge_data$to), go_results$ID)]),
  color = go_results$NES[match(unique(edge_data$to), go_results$ID)],
  stringsAsFactors = FALSE
)

node_data <- rbind(gene_nodes, pathway_nodes) %>%
  dplyr::arrange(Class)

# network data
temp_data <- tidygraph::tbl_graph(
  nodes = node_data,
  edges = edge_data,
  directed = TRUE
) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

# network plot
plot <- ggraph(temp_data, layout = 'linear', circular = TRUE) +
  # add edge
  geom_edge_arc(aes(color = Correlation),
                alpha = 0.6,  # transparency
                show.legend = TRUE) +
  # add node
  geom_node_point(aes(fill = color,
                      size = size),
                  shape = 21,
                  stroke = 0.5,  # borders
                  show.legend = TRUE) +
  # add text tag
  geom_node_text(
    aes(
      label = node,
      x = x * 1.15,  # add text distance
      y = y * 1.15,
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      colour = Class
    ),
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  # legend
  guides(
    edge_color = ggraph::guide_edge_colorbar(title = "NES"),
    fill = guide_colorbar(title = "log2(FC)"),
    size = guide_legend(title = "Gene: |log2(FC)|\nPathway: -log10(P)")
  ) +
  # colors and range
  ggraph::scale_edge_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#E64B35FF",
    midpoint = 0
  ) +
  scale_fill_gradient2(
    low = "#DC0000FF",
    mid = "white", 
    high = "#3C5488FF",
    midpoint = 0
  ) +
  ggraph::scale_edge_width(range = c(0.2, 1)) +
  scale_size_continuous(range = c(1, 7)) +  # node size range
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )

plot

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

###pathway_info
pathway_info <- go_results %>%
  as.data.frame() %>%
  dplyr::filter(ID %in% top_pathways) %>%
  dplyr::select(ID, NES, pvalue, p.adjust) %>%
  dplyr::arrange(desc(abs(NES)))

print(pathway_info)





