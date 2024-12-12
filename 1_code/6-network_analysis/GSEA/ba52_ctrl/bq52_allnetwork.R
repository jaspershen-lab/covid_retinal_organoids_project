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
load("3_data_analysis/4-differential_analysis/1-transcriptome/transcriptome_data")
data <- transcriptome_data@variable_info[, c("SYMBOL", "fc_wt_ctrl", "p_value_adjust_wt_ctrl")]
go_results <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/GO/GSEA_go_wt_ctrl.csv")
kegg_results <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/KEGG/GSEA_KEGG_wt_ctrl.csv")
reactome_results <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/REACTOME/GSEA_reactome_wt_ctrl.csv")

# Create output directory
output_dir <- "3_data_analysis/6-network_analysis/GSEA/wt_ctrl"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# top pathways
get_top_pathways <- function(results, n = 5) {
  results %>%
    as.data.frame() %>%
    mutate(type = if_else(NES > 0, "positive", "negative")) %>%
    group_by(type) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = n) %>%
    ungroup()
}

# 合并三个结果并添加来源标签
all_results <- bind_rows(
  get_top_pathways(go_results) %>% mutate(source = "GO"),
  get_top_pathways(kegg_results) %>% mutate(source = "KEGG"),
  get_top_pathways(reactome_results) %>% mutate(source = "Reactome")
) %>%
  # 分别选择NES最大和最小的5个通路
  group_by(NES > 0) %>%  # 按NES正负分组
  arrange(desc(abs(NES))) %>%  # 在每组内按绝对值排序
  slice_head(n = 5) %>%  # 每组取前5个
  ungroup()

# 获取最终的top pathways
top_pathways <- all_results %>% pull(Description)

# 2. 处理edge data
edge_data <- all_results %>%
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
    Correlation = all_results$NES[match(Description, all_results$Description)],
    fdr = -log10(all_results$pvalue[match(Description, all_results$Description)])
  ) %>%
  # 添加基因的fold change信息
  left_join(
    data %>% 
      select(SYMBOL, fc_wt_ctrl),
    by = c("from" = "SYMBOL")
  )

# 3. 处理gene data
data_id <- bitr(data$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
data_all <- merge(data, data_id, by = "SYMBOL") %>%
  arrange(desc(fc_wt_ctrl))

# 4. 获取基因信息
gene_info <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = unique(edge_data$from),
                                   columns = c("SYMBOL", "GENENAME"),
                                   keytype = "SYMBOL")

# 5. 创建node data
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

# 在pathway_nodes中添加source信息
pathway_nodes <- data.frame(
  node = unique(edge_data$to),
  Class = "Pathway",
  size = -log10(all_results$pvalue[match(unique(edge_data$to), all_results$Description)]),
  color = all_results$NES[match(unique(edge_data$to), all_results$Description)],
  source = all_results$source[match(unique(edge_data$to), all_results$Description)],
  GENENAME = unique(edge_data$to),
  stringsAsFactors = FALSE
)

# 在gene_nodes中添加source列（为NA）
gene_nodes <- gene_nodes %>%
  mutate(source = "Gene")


# 合并node data
node_data <- rbind(gene_nodes, pathway_nodes) %>%
  arrange(Class) %>%
  mutate(radius = if_else(Class == "Pathway", 0.7, 1))%>%
  mutate(color_group = case_when(
    Class == "Pathway" & color > 0 ~ "pathway_pos",
    Class == "Pathway" & color <= 0 ~ "pathway_neg",
    Class == "Gene" & color > 0 ~ "gene_pos",
    Class == "Gene" & color <= 0 ~ "gene_neg"
  ))

#筛选每个通路的top5基因
# filtered_edges <- edge_data %>%
#   group_by(to) %>%
#   # 直接按fold change的绝对值排序
#   arrange(desc(abs(fc_wt_ctrl))) %>%
#   slice_head(n = 5) %>%
#   ungroup()
# 
# genes_to_show <- unique(filtered_edges$from)

####根据NES筛选top5基因
filtered_edges <- edge_data %>%
  group_by(to) %>%
  mutate(pathway_nes = first(Correlation)) %>%
  group_modify(~{
    if(first(.x$pathway_nes) > 0) {
      arrange(.x, desc(fc_wt_ctrl))  # NES正，选择fold change最大的基因
    } else {
      arrange(.x, fc_wt_ctrl)        # NES负，选择fold change最小的基因
    }
  }) %>%
  slice_head(n = 3) %>%
  ungroup()

genes_to_show <- unique(filtered_edges$from)

# 8. 创建网络图
temp_data <- tidygraph::tbl_graph(
  nodes = node_data,  # node_data现在包含了color_group列
  edges = edge_data,
  directed = TRUE
) %>%
  mutate(Degree = centrality_degree(mode = 'all'))

# 9. 绘制网络图
library(ggnewscale)

plot <- ggraph(temp_data, 
               layout = "kk"
) +
  geom_edge_diagonal(aes(color = Correlation),
                     alpha = 0.3,
                     width = 0.5) +
  # Gene节点使用第一个填充色系统
  geom_node_point(data = function(x) x %>% filter(Class == "Gene"),
                  aes(fill = color,
                      size = size,
                      shape = source),
                  color = "black", 
                  stroke = 0.5) +
  scale_fill_gradient2(
    name = "Gene log2(FC)",
    low = "#00468BFF",
    high = "#ED0000FF",
    mid = "white", 
    midpoint = 0
  ) +
  # 添加新的填充色系统用于Pathway
  new_scale_fill() +
  # Pathway节点使用第二个填充色系统
  geom_node_point(data = function(x) x %>% filter(Class == "Pathway"),
                  aes(fill = color,
                      size = size,
                      shape = source),
                  color = "black", 
                  stroke = 0.5) +
  scale_fill_gradient2(
    name = "Pathway NES",
    low = "#00B2A9FF",
    high = "#C50084FF",
    mid = "white", 
    midpoint = 0
  ) +
  geom_node_text(aes(
    label = case_when(
      Class == "Pathway" ~ GENENAME,
      Class == "Gene" & node %in% genes_to_show ~ GENENAME,
      TRUE ~ ""
    ),
    color = color_group  
  ),
  repel = TRUE,
  max.overlaps = Inf,
  size = ifelse(V(temp_data)$Class == "Pathway", 4, 3),
  fontface = ifelse(V(temp_data)$Class == "Pathway", "bold", "plain"),
  check_overlap = TRUE
  ) +
  scale_edge_color_gradient2(
    low = "#00468BFF",
    mid = "white",
    high = "#ED0000FF",
    midpoint = 0
  ) +
  scale_color_manual(
    values = c(
      "pathway_pos" = "#C50084FF",
      "pathway_neg" = "#00B2A9FF",
      "gene_pos" = "#ED0000FF",
      "gene_neg" = "#00468BFF"
    ),
    guide = "none"
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_shape_manual(
    values = c(
      "GO" = 21,      
      "KEGG" = 22,    
      "Reactome" = 23 
    ),
    na.value = 21
  ) +
  guides(
    edge_color = guide_edge_colorbar(title = "NES",
                                     order = 1),
    size = guide_legend(title = "-log10(P)",
                        order = 2),
    shape = guide_legend(title = "Database",
                         override.aes = list(
                           size = 5,
                           fill = "black",  
                           color = "white",  
                           alpha = 1        
                         ),
                         order = 3)
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "white"),
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white", color = NA)
  )

plot

# Save output
ggsave(
  plot,
  filename = "top10_NES_pathway_network_wt_ctrl_all_top10.pdf",
  width = 10,
  height = 5,
  bg = "transparent",
  dpi = 300
)

ggsave(
  plot,
  filename = "top10_NES_pathway_network_wt_ctrl_all_top10.png",
  width = 10,
  height = 5,
  bg = "transparent",
  dpi = 300
)

######barplot
# 获取top 10 pathways (正负各5个)
top10_pathways <- bind_rows(
  get_top_pathways(go_results, n = 5) %>% mutate(source = "GO"),
  get_top_pathways(kegg_results, n = 5) %>% mutate(source = "KEGG"),
  get_top_pathways(reactome_results, n = 5) %>% mutate(source = "Reactome")
) %>%
  group_by(NES > 0) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  # 分组排序
  arrange(NES) %>%  # 先按NES排序
  # 创建一个简短的Description
  mutate(
    short_desc = str_wrap(Description, width = 35),  # 限制文本长度
    short_desc = factor(short_desc, levels = unique(short_desc))  # 保持顺序
  )

# 创建barplot
barplot <- ggplot(top10_pathways, 
                  aes(x = NES, 
                      y = short_desc,
                      fill = source)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, 
             linetype = "solid", 
             color = "black") +
  scale_fill_manual(values = c("GO" = "#F39B7FFF",
                               "KEGG" = "#8491B4FF",
                               "Reactome" = "#91D1C2FF")) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "",
       fill = "Database") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # 修改这一行，使用element_text()的自定义函数来设置颜色
    axis.text.y = element_text(
      size = 10,
      color = ifelse(top10_pathways$NES > 0, "#C50084FF", "#00B2A9FF")
    ),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  # 添加虚线网格
  geom_vline(xintercept = seq(-6, 6, by = 2),
             linetype = "dashed",
             color = "grey90") +
  # 设置x轴范围和刻度
  scale_x_continuous(
    limits = c(min(top10_pathways$NES) - 0.5, max(top10_pathways$NES) + 0.5),
    breaks = seq(-6, 6, by = 2)
  )

barplot

# Save output
ggsave(
  barplot,
  filename = "top10_NES_pathway_network_wt_ctrl_all_top10_barplot.pdf",
  width = 8.5,
  height = 7,
  bg = "transparent",
  dpi = 300
)

ggsave(
  barplot,
  filename = "top10_NES_pathway_network_wt_ctrl_all_top10_barplot.png",
  width = 8.5,
  height = 7,
  bg = "transparent",
  dpi = 300
)
