library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(tidymass)

##read data
load("3_data_analysis/2-data_cleaning/2-metabolome/metabolome_data.RData")

dir.create(
  "3_data_analysis/4-differential_analysis/2-metabolome_zj",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/4-differential_analysis/2-metabolome_zj")

metabolome_data@sample_info
metabolome_data@variable_info
metabolome_data@expression_data

metabolome_data

ba52_id <-
  metabolome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BA52") %>%
  pull(sample_id)

bq11_id <-
  metabolome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BQ11") %>%
  pull(sample_id)

wt_id <-
  metabolome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "WT") %>%
  pull(sample_id)


ctrl_id <-
  metabolome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Ctrl") %>%
  pull(sample_id)

#####ba52 VS wt
metabolome_data <-
  metabolome_data %>%
  mutate_fc(
    control_sample_id = wt_id,
    case_sample_id = ba52_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  mutate_p_value(
    control_sample_id = wt_id,
    case_sample_id = ba52_id,
    method = "t.test",
    p_adjust_methods = "BH",
    return_mass_dataset = TRUE
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::rename(
    fc_ba52_wt = fc,
    p_value_ba52_wt = p_value,
    p_value_adjust_ba52_wt = p_value_adjust
  )

head(metabolome_data@variable_info)
total_genes <- nrow(metabolome_data@variable_info)
total_genes


sig_genes <- metabolome_data@variable_info %>%
  filter(p_value_ba52_wt < 0.05)

upregulated_genes <- sig_genes %>%
  filter(fc_ba52_wt > 2) %>%
  nrow()

downregulated_genes <- sig_genes %>%
  filter(fc_ba52_wt < 0.5) %>%
  nrow()

cat("Number of upregulated genes: ", upregulated_genes, "\n")
cat("Number of downregulated genes: ", downregulated_genes, "\n")


#####bq11 VS wt
metabolome_data <-
  metabolome_data %>%
  mutate_fc(
    control_sample_id = wt_id,
    case_sample_id = bq11_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  mutate_p_value(
    control_sample_id = wt_id,
    case_sample_id = bq11_id,
    method = "t.test",
    p_adjust_methods = "BH",
    return_mass_dataset = TRUE
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::rename(
    fc_bq11_wt = fc,
    p_value_bq11_wt = p_value,
    p_value_adjust_bq11_wt = p_value_adjust
  )

head(metabolome_data@variable_info)

sig_genes <- metabolome_data@variable_info %>%
  filter(p_value_adjust_bq11_wt < 0.05)

upregulated_genes <- sig_genes %>%
  filter(fc_bq11_wt > 2) %>%
  nrow()

downregulated_genes <- sig_genes %>%
  filter(fc_bq11_wt < 0.5) %>%
  nrow()

cat("Number of upregulated genes: ", upregulated_genes, "\n")
cat("Number of downregulated genes: ", downregulated_genes, "\n")



#####wt VS ctrl
metabolome_data <-
  metabolome_data %>%
  mutate_fc(
    control_sample_id = ctrl_id,
    case_sample_id = wt_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  mutate_p_value(
    control_sample_id = ctrl_id,
    case_sample_id = wt_id,
    method = "t.test",
    p_adjust_methods = "BH",
    return_mass_dataset = TRUE
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::rename(
    fc_wt_ctrl = fc,
    p_value_wt_ctrl = p_value,
    p_value_adjust_wt_ctrl = p_value_adjust
  )

head(metabolome_data@variable_info)
total_genes <- nrow(metabolome_data@variable_info)
total_genes


sig_genes <- metabolome_data@variable_info %>%
  filter(p_value_adjust_wt_ctrl < 0.05)

upregulated_genes <- sig_genes %>%
  filter(fc_wt_ctrl > 2) %>%
  nrow()

downregulated_genes <- sig_genes %>%
  filter(fc_wt_ctrl < 0.5) %>%
  nrow()

cat("Number of upregulated genes: ", upregulated_genes, "\n")
cat("Number of downregulated genes: ", downregulated_genes, "\n")



#####ba52 VS ctrl
metabolome_data <-
  metabolome_data %>%
  mutate_fc(
    control_sample_id = ctrl_id,
    case_sample_id = ba52_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  mutate_p_value(
    control_sample_id = ctrl_id,
    case_sample_id = ba52_id,
    method = "t.test",
    p_adjust_methods = "BH",
    return_mass_dataset = TRUE
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::rename(
    fc_ba52_ctrl = fc,
    p_value_ba52_ctrl = p_value,
    p_value_adjust_ba52_ctrl = p_value_adjust
  )

head(metabolome_data@variable_info)
total_genes <- nrow(metabolome_data@variable_info)
total_genes


sig_genes <- metabolome_data@variable_info %>%
  filter(p_value_adjust_ba52_ctrl < 0.05)

upregulated_genes <- sig_genes %>%
  filter(fc_ba52_ctrl > 2) %>%
  nrow()

downregulated_genes <- sig_genes %>%
  filter(fc_ba52_ctrl < 0.5) %>%
  nrow()

cat("Number of upregulated genes: ", upregulated_genes, "\n")
cat("Number of downregulated genes: ", downregulated_genes, "\n")



#####bq11 VS ctrl
metabolome_data <-
  metabolome_data %>%
  mutate_fc(
    control_sample_id = ctrl_id,
    case_sample_id = bq11_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  mutate_p_value(
    control_sample_id = ctrl_id,
    case_sample_id = bq11_id,
    method = "t.test",
    p_adjust_methods = "BH",
    return_mass_dataset = TRUE
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::rename(
    fc_bq11_ctrl = fc,
    p_value_bq11_ctrl = p_value,
    p_value_adjust_bq11_ctrl = p_value_adjust
  )

head(metabolome_data@variable_info)
total_genes <- nrow(metabolome_data@variable_info)
total_genes


sig_genes <- metabolome_data@variable_info %>%
  filter(p_value_adjust_bq11_ctrl < 0.05)

upregulated_genes <- sig_genes %>%
  filter(fc_bq11_ctrl > 2) %>%
  nrow()

downregulated_genes <- sig_genes %>%
  filter(fc_bq11_ctrl < 0.5) %>%
  nrow()

cat("Number of upregulated genes: ", upregulated_genes, "\n")
cat("Number of downregulated genes: ", downregulated_genes, "\n")


###Volcano plot
###ba52 vs WT
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ba52_wt",
    p_value_column_name = "p_value_ba52_wt",
    labs_x = "log2(Fold change, BA52/WT)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_ba52_wt.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_ba52_wt.pdf",
       width = 7,
       height = 6)

plot_adjust <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ba52_wt",
    p_value_column_name = "p_value_adjust_ba52_wt",
    labs_x = "log2(Fold change, BA52/WT)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot_adjust

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_ba52_wt.png",
       width = 7,
       height = 6)

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_ba52_wt.pdf",
       width = 7,
       height = 6)



###bq11 vs WT
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_bq11_wt",
    p_value_column_name = "p_value_bq11_wt",
    labs_x = "log2(Fold change, BQ11/WT)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_bq11_wt.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_bq11_wt.pdf",
       width = 7,
       height = 6)


plot_adjust <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_bq11_wt",
    p_value_column_name = "p_value_adjust_bq11_wt",
    labs_x = "log2(Fold change, BQ11/WT)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot_adjust

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_bq11_wt.png",
       width = 7,
       height = 6)

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_bq11_wt.pdf",
       width = 7,
       height = 6)



###WT vs ctrl
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_wt_ctrl",
    p_value_column_name = "p_value_wt_ctrl",
    labs_x = "log2(Fold change, WT/Ctrl)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_wt_ctrl.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_wt_ctrl.pdf",
       width = 7,
       height = 6)

plot_adjust <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_wt_ctrl",
    p_value_column_name = "p_value_adjust_wt_ctrl",
    labs_x = "log2(Fold change, WT/Ctrl)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot_adjust

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_wt_ctrl.png",
       width = 7,
       height = 6)

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_wt_ctrl.pdf",
       width = 7,
       height = 6)



###ba52 vs ctrl
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ba52_ctrl",
    p_value_column_name = "p_value_ba52_ctrl",
    labs_x = "log2(Fold change, BA52/Ctrl)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_ba52_ctrl.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_ba52_ctrl.pdf",
       width = 7,
       height = 6)

plot_adjust <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ba52_ctrl",
    p_value_column_name = "p_value_adjust_ba52_ctrl",
    labs_x = "log2(Fold change, BA52/Ctrl)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot_adjust

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_ba52_ctrl.png",
       width = 7,
       height = 6)

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_ba52_ctrl.pdf",
       width = 7,
       height = 6)




###bq11 vs ctrl
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_bq11_ctrl",
    p_value_column_name = "p_value_bq11_ctrl",
    labs_x = "log2(Fold change, BQ11/Ctrl)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_bq11_ctrl.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_bq11_ctrl.pdf",
       width = 7,
       height = 6)

plot_adjust <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_bq11_ctrl",
    p_value_column_name = "p_value_adjust_bq11_ctrl",
    labs_x = "log2(Fold change, BQ11/Ctrl)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot_adjust

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_bq11_ctrl.png",
       width = 7,
       height = 6)

ggsave(plot_adjust,
       filename = "volcano_plot_adjust_bq11_ctrl.pdf",
       width = 7,
       height = 6)

save(metabolome_data, file = "metabolome_data")




#####venn plot#######
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

library(VennDiagram)
library(dplyr)

# extract genelist
#wt_ctrl
sig_genes_wt_ctrl <- metabolome_data@variable_info %>%
  filter(p_value_adjust_wt_ctrl < 0.05)

upregulated_genes_wt_ctrl <- sig_genes_wt_ctrl %>%
  filter(fc_wt_ctrl > 2) %>%
  nrow()

downregulated_genes_wt_ctrl <- sig_genes %>%
  filter(fc_wt_ctrl < 0.5) %>%
  nrow()

#ba52_ctrl
sig_genes_ba52_ctrl <- metabolome_data@variable_info %>%
  filter(p_value_adjust_ba52_ctrl < 0.05)

upregulated_genes_ba52_ctrl <- sig_genes_ba52_ctrl %>%
  filter(fc_ba52_ctrl > 2) %>%
  nrow()

downregulated_genes_ba52_ctrl <- sig_genes %>%
  filter(fc_ba52_ctrl < 0.5) %>%
  nrow()


#bq11_ctrl
sig_genes_bq11_ctrl <- metabolome_data@variable_info %>%
  filter(p_value_adjust_bq11_ctrl < 0.05)

upregulated_genes_bq11_ctrl <- sig_genes_bq11_ctrl %>%
  filter(fc_bq11_ctrl > 2) %>%
  nrow()

downregulated_genes_bq11_ctrl <- sig_genes %>%
  filter(fc_bq11_ctrl < 0.5) %>%
  nrow()

# 提取基因 ID 列表
# wt_ctrl 上调和下调基因
upregulated_genes_wt_ctrl <- sig_genes_wt_ctrl %>%
  filter(fc_wt_ctrl > 2) %>%
  pull(SYMBOL)

downregulated_genes_wt_ctrl <- sig_genes_wt_ctrl %>%
  filter(fc_wt_ctrl < 0.5) %>%
  pull(SYMBOL)

# ba52_ctrl 上调和下调基因
upregulated_genes_ba52_ctrl <- sig_genes_ba52_ctrl %>%
  filter(fc_ba52_ctrl > 2) %>%
  pull(SYMBOL)

downregulated_genes_ba52_ctrl <- sig_genes_ba52_ctrl %>%
  filter(fc_ba52_ctrl < 0.5) %>%
  pull(SYMBOL)

# bq11_ctrl 上调和下调基因
upregulated_genes_bq11_ctrl <- sig_genes_bq11_ctrl %>%
  filter(fc_bq11_ctrl > 2) %>%
  pull(SYMBOL)

downregulated_genes_bq11_ctrl <- sig_genes_bq11_ctrl %>%
  filter(fc_bq11_ctrl < 0.5) %>%
  pull(SYMBOL)


###veen
library(VennDiagram)
library(grid)

get_common_genes <- function(set1, set2) {
  intersect(set1, set2)
}

plot_venn <- function(set1, set2, set1_name, set2_name, filename) {
  # 绘制 Venn 图对象
  venn.plot <- venn.diagram(
    x = list(set1 = set1, set2 = set2),
    category.names = c(set1_name, set2_name),
    filename = NULL,  # 不直接保存，稍后用 grid 绘制
    col = "black",
    fill = c("#4682B4", "#DB7093"),  # 自定义颜色
    alpha = 0.5,  # 设置透明度
    cat.col = c("#4682B4", "#DB7093"),  # 类别名称颜色
    cat.cex = 1.5,  # 类别名称字体大小
    cex = 1.5,  # Venn 图数字字体大小
    cat.fontface = "bold",  # 类别名称字体加粗
    fontface = "bold",  # 数字字体加粗
    margin = 0.1
  )
  
  # 保存为 PDF 格式
  pdf(file = filename, width = 8, height = 8)  # 指定宽高
  grid.draw(venn.plot)  # 使用 grid 绘制图形
  dev.off()  # 关闭设备
  
  # 返回公共基因集
  common_genes <- get_common_genes(set1, set2)
  return(common_genes)
}

# 上调基因的 Venn 图（WT vs Ctrl 和 BA52 vs Ctrl）
common_upregulated_genes_wt_ba52 <- plot_venn(
  set1 = upregulated_genes_wt_ctrl,
  set2 = upregulated_genes_ba52_ctrl,
  set1_name = "WT vs Ctrl Upregulated",
  set2_name = "BA52 vs Ctrl Upregulated",
  filename = "wt_ba52_upregulated_venn.pdf"
)

# 将公共上调基因保存到文件
write.csv(common_upregulated_genes_wt_ba52, "common_upregulated_genes_wt_ba52.csv", row.names = FALSE)

# 下调基因的 Venn 图（WT vs Ctrl 和 BA52 vs Ctrl）
common_downregulated_genes_wt_ba52 <- plot_venn(
  set1 = downregulated_genes_wt_ctrl,
  set2 = downregulated_genes_ba52_ctrl,
  set1_name = "WT vs Ctrl Downregulated",
  set2_name = "BA52 vs Ctrl Downregulated",
  filename = "wt_ba52_downregulated_venn.pdf"
)

# 将公共下调基因保存到文件
write.csv(common_downregulated_genes_wt_ba52, "common_downregulated_genes_wt_ba52.csv", row.names = FALSE)


# 上调基因的Venn图（WT vs Ctrl 和 BQ11 vs Ctrl）
common_upregulated_genes_wt_bq11 <- plot_venn(
  set1 = upregulated_genes_wt_ctrl,
  set2 = upregulated_genes_bq11_ctrl,
  set1_name = "WT vs Ctrl Upregulated",
  set2_name = "BQ11 vs Ctrl Upregulated",
  filename = "wt_bq11_upregulated_venn.pdf"
)

# 将公共上调基因保存到文件
write.csv(common_upregulated_genes_wt_bq11, "common_upregulated_genes_wt_bq11.csv", row.names = FALSE)

# 下调基因的Venn图（WT vs Ctrl 和 BQ11 vs Ctrl）
common_downregulated_genes_wt_bq11 <- plot_venn(
  set1 = downregulated_genes_wt_ctrl,
  set2 = downregulated_genes_bq11_ctrl,
  set1_name = "WT vs Ctrl Downregulated",
  set2_name = "BQ11 vs Ctrl Downregulated",
  filename = "wt_bq11_downregulated_venn.pdf"
)

# 将公共下调基因保存到文件
write.csv(common_downregulated_genes_wt_bq11, "common_downregulated_genes_wt_bq11.csv", row.names = FALSE)

