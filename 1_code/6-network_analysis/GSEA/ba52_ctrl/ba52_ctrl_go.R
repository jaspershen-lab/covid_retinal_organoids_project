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
library(ggnewscale) 

##read data ba52_ctrl
load("3_data_analysis/4-differential_analysis/1-transcriptome/transcriptome_data")
data <- transcriptome_data@variable_info[, c("SYMBOL", "fc_ba52_ctrl", "p_value_adjust_ba52_ctrl")]
go_results<- read_csv("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/GO/GSEA_GO_ba52_ctrl.csv")

#####
if (!dir.exists("3_data_analysis/6-network_analysis/GSEA/ba52_ctrl/ba52_ctrl_go")) {
  dir.create(
    "3_data_analysis/6-network_analysis/GSEA/ba52_ctrl/ba52_ctrl_go",
    recursive = TRUE,
    showWarnings = FALSE
  )
}


setwd("3_data_analysis/6-network_analysis/GSEA/ba52_ctrl/ba52_ctrl_go")

#####GO
# Top 10 pathways
top_pathways <- go_results %>%
  as.data.frame() %>%
  mutate(type = if_else(NES > 0, "positive", "negative")) %>%
  group_by(type) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(Description)

print(top_pathways)

# edge_data
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
  dplyr::arrange(desc(fc_ba52_ctrl))

##gene_nodes
gene_nodes <- tibble(
  node = unique(edge_data$from)
) %>%
  mutate(Class = "Gene") %>%
  left_join(
    tibble(
      node = data_all$SYMBOL,
      size = -log10(data_all$p_value_adjust_ba52_ctrl),  
      color = log2(data_all$fc_ba52_ctrl)               
    ),
    by = "node"
  )

##pathway_node
pathway_nodes <- data.frame(
  node = unique(edge_data$to),
  Class = "Pathway",
  size = -log10(go_results$pvalue[match(unique(edge_data$to), go_results$Description)]),
  color = go_results$NES[match(unique(edge_data$to), go_results$Description)],
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
  # 添加边
  geom_edge_arc(aes(color = Correlation),
                alpha = 0.6,
                show.legend = TRUE) +
  # 设置边的颜色
  scale_edge_color_gradient2(
    low = "#4DBBD5FF",    # 下调/负相关-蓝色
    mid = "white",        
    high = "#E64B35FF",   # 上调/正相关-红色
    midpoint = 0
  ) +
  
  # 首先添加Gene节点
  geom_node_point(data = function(x) x %>% filter(Class == "Gene"),
                  aes(fill = color,
                      size = size),
                  shape = 21,
                  stroke = 0.5,
                  show.legend = TRUE) +
  # Gene节点的颜色映射
  scale_fill_gradient2(
    name = "Gene: log2(FC)",
    low = "#00468BFF",    # 下调-蓝色
    mid = "white",
    high = "#ED0000FF",   # 上调-红色
    midpoint = 0
  ) +
  
  # 添加新的颜色映射
  new_scale_fill() +
  
  # 添加Pathway节点
  geom_node_point(data = function(x) x %>% filter(Class == "Pathway"),
                  aes(fill = color,
                      size = size),
                  shape = 21,
                  stroke = 0.5,
                  show.legend = TRUE) +
  # Pathway节点的颜色映射
  scale_fill_gradient2(
    name = "Pathway: NES",
    low = "#4DBBD5FF",    # 负向富集-深蓝
    mid = "white",
    high = "#E64B35FF",   # 正向富集-深红
    midpoint = 0
  ) +
  
  # 添加文本标签
  geom_node_text(
    aes(
      label = node,
      x = x * 1.15,
      y = y * 1.15,
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      colour = Class
    ),
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  # 设置文本颜色
  scale_colour_manual(
    values = c("Gene" = "black", "Pathway" = "#20854EFF")
  ) +
  
  # 设置大小
  scale_size_continuous(
    name = "-log10(P)",
    range = c(1, 7)
  ) +
  
  # 设置图例
  guides(
    edge_color = guide_edge_colorbar(title = "NES"),
    size = guide_legend(title = "-log10(P)")
  ) +
  
  # 主题设置
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

plot

ggsave(
  plot,
  filename = "top10_NES_pathway_network_ba52_ctrl_go.pdf",
  width = 8.5,
  height = 7,
  bg = "transparent",
  dpi = 300
)

ggsave(
  plot,
  filename = "top10_NES_pathway_network_ba52_ctrl_go.png",
  width = 8.5,
  height = 7,
  bg = "transparent",
  dpi = 300
)

###pathway_info
pathway_info <- go_results %>%
  as.data.frame() %>%
  dplyr::filter(Description %in% top_pathways) %>%
  dplyr::select(Description, NES, pvalue, p.adjust) %>%
  dplyr::arrange(desc(abs(NES)))

print(pathway_info)





