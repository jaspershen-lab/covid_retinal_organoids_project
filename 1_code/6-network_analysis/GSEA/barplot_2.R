library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(readxl)
library(ggplot2)
library(dplyr)


# Load data
data <- read_excel("2_data/ba52_wt_vs_ctrl.xlsx")

head(data)

# Create output directory
output_dir <- "3_data_analysis/6-network_analysis/GSEA/barplot"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

###sort data
data <- data %>%
  arrange(desc(abs(NES_BA52)))

data_long <- data %>%
  pivot_longer(cols = c(NES_BA52, NES_WT), names_to = "Group", values_to = "NES")

# plot
plot<-ggplot(data_long, aes(x = reorder(Description_x, abs(NES)), y = NES, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  coord_flip() +  # 翻转坐标轴，水平条形图
  scale_fill_manual(values = c("NES_BA52" = "skyblue", "NES_WT" = "salmon")) +
  labs(
    title = "Comparison of NES for Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

plot

ggsave(
  plot,
  filename = "ba5_wt_vs_ctrl_shared_pathway_NES.pdf",
  width = 12,
  height = 9,
  bg = "white",
  dpi = 300
)
