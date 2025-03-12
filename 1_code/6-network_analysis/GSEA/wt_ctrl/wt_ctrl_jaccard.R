library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggnewscale)

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

# top pathways 函数
get_top_pathways <- function(results, n = 15) {
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
  group_by(NES > 0) %>%  # 按NES正负分组
  arrange(desc(abs(NES))) %>%  # 在每组内按绝对值排序
  slice_head(n = 10) %>%  # 每组取前5个
  ungroup()

# 1. 获取每个pathway包含的基因集
pathway_genes <- all_results %>%
  mutate(genes = strsplit(core_enrichment, "/")) %>%
  select(Description, source, genes, NES, pvalue) %>%
  mutate(size = sapply(genes, length))  # 添加基因数量

# 2. 计算pathway之间的Jaccard index函数
calc_jaccard <- function(genes1, genes2) {
  intersection <- length(intersect(genes1, genes2))
  union <- length(union(genes1, genes2))
  return(intersection/union)
}

# 3. 创建pathway edge data
pathway_edges <- expand.grid(
  from = pathway_genes$Description,
  to = pathway_genes$Description,
  stringsAsFactors = FALSE
) %>%
  filter(from < to) %>%  # 避免重复连接
  rowwise() %>%
  mutate(
    jaccard = calc_jaccard(
      pathway_genes$genes[pathway_genes$Description == from][[1]],
      pathway_genes$genes[pathway_genes$Description == to][[1]]
    )
  ) %>%
  filter(jaccard > 0)  # 只保留有共享基因的pathway对

# 4. 创建node data
node_data <- pathway_genes %>%
  transmute(
    node = Description,
    short_name = str_wrap(Description, width = 20), 
    Class = "Pathway",
    size = size,  # 使用基因数量作为节点大小
    color = NES,
    source = source,
    GENENAME = Description
  ) %>%
  mutate(color_group = if_else(color > 0, "pathway_pos", "pathway_neg"))

# 5. 创建网络图对象
temp_data <- tidygraph::tbl_graph(
  nodes = node_data,
  edges = pathway_edges,
  directed = FALSE
)


# 6. 绘制网络图
plot <- ggraph(temp_data, 
               layout = 'linear',
               circular = TRUE , # 设置为环形布局
) +
  # 使用弧形边
  geom_edge_arc(aes(width = jaccard,
                    alpha = jaccard),
                edge_colour = "grey50",  
                start_cap = circle(2, "mm"),
                end_cap = circle(2, "mm"),
                fold = TRUE,  # 控制弧线的方向
                show.legend = TRUE) +
  # 节点
  geom_node_point(aes(fill = color,
                      size = size,
                      shape = source),
                  color = "black", 
                  stroke = 0.5) +
  # 文本标签
  geom_node_text(aes(
    label = short_name,
    color = color_group,
    x = x *1.03,  # 调整文本距离
    y = y *1.05,
    angle = -((-node_angle(x, y) + 90) %% 180) + 90,  # 调整文本角度
    hjust = 'outward'  # 文本向外对齐
  ),
  repel = FALSE,  # 环形布局不需要repel
  size = 12,
  fontface = "bold",
  show.legend = FALSE
  ) +
  scale_fill_gradient2(
    name = "Pathway NES",
    low = "#00B2A9FF",
    high = "#C50084FF",
    mid = "white", 
    midpoint = 0
  ) +
  scale_color_manual(
    values = c(
      "pathway_pos" = "#C50084FF",
      "pathway_neg" = "#00B2A9FF"
    ),
    guide = "none"
  ) +
  scale_edge_width(range = c(0.2, 0.5)) +
  scale_edge_alpha(range = c(0.3, 0.8)) +
  scale_size_continuous(
    name = "Gene Count",
    range = c(10, 20)
  ) +
  scale_shape_manual(
    name = "Database",
    values = c(
      "GO" = 21,      
      "KEGG" = 22,    
      "Reactome" = 24 
    )
  ) +
  guides(
    edge_width = guide_legend(title = "Jaccard Index",
                              order = 1),
    fill = guide_colorbar(title = "NES",
                          order = 2),
    size = guide_legend(title = "Gene Count",
                        order = 3),
    shape = guide_legend(title = "Database",
                         override.aes = list(
                           size = 5,
                           fill = "black",
                           color = "white"
                         ),
                         order = 4)
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

######渐变
plot <- ggraph(temp_data, 
               layout = 'linear',
               circular = TRUE) +
  
  # 使用弧形边，并调整 Jaccard Index 颜色和透明度
  geom_edge_arc(aes(width = jaccard,
                    alpha = jaccard,
                    color = jaccard),  # 让 Jaccard Index 颜色有梯度变化
                start_cap = circle(2, "mm"),
                end_cap = circle(2, "mm"),
                fold = TRUE,
                show.legend = TRUE) +
  
  # **Jaccard Index 颜色渐变**
  scale_edge_color_gradientn(
    name = "Jaccard Index",
    colors = c("#E0E0E0", "#A0A0A0", "#505050", "#000000"),  # 由浅到深
    limits = c(min(temp_data$jaccard, na.rm = TRUE), max(temp_data$jaccard, na.rm = TRUE)),
    guide = guide_colorbar(order = 1)
  ) +
  
  # **Jaccard Index 线条透明度渐变**
  scale_edge_alpha(
    name = "Jaccard Index",
    range = c(0.4, 1),  # 让最低值更透明，最高值不透明
    guide = guide_legend(
      override.aes = list(edge_width = c(0.1, 0.2, 0.3, 0.4))  # **放大图例线段**
    )
  ) +
  
  # **Jaccard Index 线条宽度调整**
  scale_edge_width(
    name = "Jaccard Index",
    range = c(0.5, 3),  # 让 Jaccard Index 低值更细，高值更粗
    guide = guide_legend(order = 1)
  ) +
  
  # 绿色节点 (负值)
  geom_node_point(data = . %>% filter(color < 0),
                  aes(fill = color,
                      size = size,
                      shape = source),
                  color = "black",
                  stroke = 0.5) +
  
  scale_fill_gradientn(
    name = "Negative NES",
    colors = c("#00B2A9", "#4CCCC5", "#99E5E1", "#CCEEED"),
    limits = c(-1.9, -1.7),
    breaks = seq(-1.9, -1.7, 0.1),
    guide = guide_colorbar(order = 2)
  ) +
  
  # 添加新的填充色标
  new_scale_fill() +
  
  # 粉色节点 (正值)
  geom_node_point(data = . %>% filter(color >= 0),
                  aes(fill = color,
                      size = size,
                      shape = source),
                  color = "black",
                  stroke = 0.5) +
  
  scale_fill_gradientn(
    name = "Positive NES",
    colors = c("#FFE4F0", "#FFD0E6", "#FFB6D8", "#E57CAE", "#C50084FF"),
    limits = c(2.0, 2.9),
    breaks = seq(2.0, 2.9, 0.2),
    guide = guide_colorbar(order = 3)
  ) +
  
  geom_node_text(aes(
    label = short_name,
    color = color_group,
    x = x *1.03,
    y = y *1.05,
    angle = -((-node_angle(x, y) + 90) %% 180) + 90,
    hjust = 'outward'
  ),
  repel = FALSE,  
  size = 12,
  fontface = "bold",
  show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "pathway_pos" = "#C50084FF",
      "pathway_neg" = "#00B2A9FF"
    ),
    guide = "none"
  ) +
  
  # **放大 Gene Count**
  scale_size_continuous(
    name = "Gene Count",
    range = c(12, 25)  # 让节点更明显
  ) +
  
  scale_shape_manual(
    name = "Database",
    values = c(
      "GO" = 21,      
      "KEGG" = 22,    
      "Reactome" = 24 
    )
  ) +
  
  guides(
    edge_width = guide_legend(title = "Jaccard Index",
                              order = 1),
    size = guide_legend(title = "Gene Count",
                        order = 3),
    shape = guide_legend(title = "Database",
                         override.aes = list(
                           size = 5,
                           fill = "black",
                           color = "white"
                         ),
                         order = 4)
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


# 保存图片
ggsave(
  plot,
  filename = "pathway_network_jaccard_wt_ctrl.pdf",
  width = 14.5,
  height = 12,
  bg = "white",
  dpi = 300
)

ggsave(
  plot,
  filename = "pathway_network_jaccard_wt_ctrl.png",
  width = 14.5,
  height = 12,
  bg = "white",
  dpi = 300
)


########NES max min
# 先看看每个数据源的top 10 pathways的NES分布
# GO
top_go <- get_top_pathways(go_results, n = 10)
cat("GO NES range:\n")
cat("Positive NES:", range(top_go$NES[top_go$NES > 0]), "\n")
cat("Negative NES:", range(top_go$NES[top_go$NES < 0]), "\n\n")

# KEGG
top_kegg <- get_top_pathways(kegg_results, n = 10)
cat("KEGG NES range:\n")
cat("Positive NES:", range(top_kegg$NES[top_kegg$NES > 0]), "\n")
cat("Negative NES:", range(top_kegg$NES[top_kegg$NES < 0]), "\n\n")

# Reactome
top_reactome <- get_top_pathways(reactome_results, n = 10)
cat("Reactome NES range:\n")
cat("Positive NES:", range(top_reactome$NES[top_reactome$NES > 0]), "\n")
cat("Negative NES:", range(top_reactome$NES[top_reactome$NES < 0]), "\n\n")

# 合并后的整体分布
cat("Overall distribution after combining and filtering:\n")
print(summary(all_results$NES))
cat("\nRange by type:\n")
all_results %>%
  group_by(type) %>%
  summarise(
    min = min(NES),
    max = max(NES),
    mean = mean(NES),
    n = n()
  ) %>%
  print()
