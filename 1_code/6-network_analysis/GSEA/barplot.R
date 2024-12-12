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
go_results_wt <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/GO/GSEA_go_wt_ctrl.csv")
kegg_results_wt <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/KEGG/GSEA_KEGG_wt_ctrl.csv")
reactome_results_wt <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/wt_ctrl/REACTOME/GSEA_reactome_wt_ctrl.csv")


go_results_ba52 <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/GO/GSEA_go_ba52_ctrl.csv")
kegg_results_ba52 <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/KEGG/GSEA_KEGG_ba52_ctrl.csv")
reactome_results_ba52 <- read_csv("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/REACTOME/GSEA_reactome_ba52_ctrl.csv")


# Create output directory
output_dir <- "3_data_analysis/6-network_analysis/GSEA/barplot"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)


# Function to count positive and negative pathways
count_pathways <- function(results) {
  results %>%
    as.data.frame() %>%
    summarise(
      positive = sum(NES > 0),
      negative = sum(NES < 0)
    )
}

# Count pathways for each database and comparison
wt_counts <- bind_rows(
  count_pathways(go_results_wt) %>% mutate(database = "GO", comparison = "wt_vs_ctrl"),
  count_pathways(kegg_results_wt) %>% mutate(database = "KEGG", comparison = "wt_vs_ctrl"),
  count_pathways(reactome_results_wt) %>% mutate(database = "Reactome", comparison = "wt_vs_ctrl")
)

ba52_counts <- bind_rows(
  count_pathways(go_results_ba52) %>% mutate(database = "GO", comparison = "ba52_vs_ctrl"),
  count_pathways(kegg_results_ba52) %>% mutate(database = "KEGG", comparison = "ba52_vs_ctrl"),
  count_pathways(reactome_results_ba52) %>% mutate(database = "Reactome", comparison = "ba52_vs_ctrl")
)

# Combine all counts
all_counts <- bind_rows(wt_counts, ba52_counts) %>%
  pivot_longer(
    cols = c(positive, negative),
    names_to = "direction",
    values_to = "count"
  ) %>%
  mutate(
    plot_count = if_else(direction == "negative", -count, count),
    database = factor(database, levels = c("GO", "KEGG", "Reactome")),
    comparison = factor(comparison, 
                        levels = c("wt_vs_ctrl", "ba52_vs_ctrl"),
                        labels = c("WT", "BA52"))
  )

# Find max count for y-axis scaling
max_count <- max(abs(all_counts$plot_count))
# Round up to nearest 50
max_axis <- ceiling(max_count/50) * 50

# Create the plot
pathway_counts_plot <- ggplot(all_counts, 
                              aes(x = database, 
                                  y = plot_count, 
                                  fill = database,
                                  group = interaction(database, comparison))) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9),  # 减小柱子间距
           width = 0.8) +  # 减小柱子宽度
  geom_hline(yintercept = 0, 
             color = "black", 
             linewidth = 0.5) +
  # # Add value labels on the bars
  # geom_text(aes(label = abs(plot_count)),
  #           position = position_dodge(width = 0.6),  # 匹配柱子间距
  #           vjust = ifelse(all_counts$plot_count >= 0, -0.5, 1.5),
  #           size = 3) +
  # # Add WT/BA52 labels
  # geom_text(aes(label = comparison),
  #           position = position_dodge(width = 0.6),  # 匹配柱子间距
  #           vjust = ifelse(all_counts$plot_count >= 0, -2, 3),
  #           size = 3) +
  scale_fill_manual(
    name = "Database",
    values = c(
      "GO" = "#F39B7FFF",
      "KEGG" = "#8491B4FF",
      "Reactome" = "#91D1C2FF"
    )
  ) +
  labs(x = "Database",
       y = "Number of Enriched Pathways",
       title = "Pathway Enrichment Analysis Results") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.line.y = element_line(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # Set y-axis breaks and labels to show actual numbers
  scale_y_continuous(
    # breaks = seq(-max_axis, max_axis, by = 50),
    # labels = function(x) abs(x),
    limits = c(-80, 20)
  )
pathway_counts_plot

# Save the plot
ggsave(
  pathway_counts_plot,
  filename = "pathway_counts_comparison.pdf",
  width = 3,
  height = 8,
  bg = "white",
  dpi = 300
)

ggsave(
  pathway_counts_plot,
  filename = "pathway_counts_comparison.png",
  width = 10,
  height = 8,
  bg = "white",
  dpi = 300
)



####pie chart
# Create data frames for both conditions
ba52_data <- data.frame(
  category = c("Neuro-related", "Other"),
  count = c(40, 60),
  percentage = c(40, 60)
)
ba52_data$label <- sprintf("%d (%g%%)", ba52_data$count, ba52_data$percentage)

wt_data <- data.frame(
  category = c("Neuro-related", "Other"),
  count = c(28, 35),
  percentage = c(round(28/63*100, 1), round(35/63*100, 1))
)
wt_data$label <- sprintf("%d (%g%%)", wt_data$count, wt_data$percentage)

# 创建BA52的饼图
p1 <- ggplot(ba52_data, aes(x = 2, y = percentage, fill = category)) + 
  geom_bar(width = 1, stat = "identity") +
  xlim(1, 2.5) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_polar(theta = "y", start = 0) +
  labs(x = NULL, y = NULL, title = "BA52\n(100)") +
  scale_fill_manual(
    values = c("Neuro-related" = "#F0AB00FF", "Other" = "#2A6EBBFF"),
    name = "Pathway Type"
  )
p1

# 创建WT的饼图
p2 <- ggplot(wt_data, aes(x = 2, y = percentage, fill = category)) + 
  geom_bar(width = 1, stat = "identity") +
  xlim(1, 2.5) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_polar(theta = "y", start = 0) +
  labs(x = NULL, y = NULL, title = "WT\n(63)") +
  scale_fill_manual(
    values = c("Neuro-related" = "#F0AB00FF", "Other" = "#2A6EBBFF"),
    name = "Pathway Type"
  )

p2


# Save the plots
ggsave(
  p1,
  filename = "neuro_pathway_proportions_ba52.pdf",
  width = 8,
  height = 4,
  bg = "white",
  dpi = 300
)

ggsave(
  p2,
  filename = "neuro_pathway_proportions_wt.pdf",
  width = 8,
  height = 4,
  bg = "white",
  dpi = 300
)
