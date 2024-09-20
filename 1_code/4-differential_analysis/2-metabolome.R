library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

##read data
load("3-data_analysis/2-data-cleaning/2-metabolome/metabolome_data.RData")

dir.create(
  "3-data_analysis/4-differential_analysis/2-metabolome",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3-data_analysis/4-differential_analysis/2-metabolome")

metabolome_data@sample_info

ctrl_id <-
  metabolome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Ctrl") %>%
  pull(sample_id)

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

##Ctrl vs WT
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
    fc_ctrl_wt = fc,
    p_value_ctrl_wt = p_value,
    p_value_adjust_ctrl_wt = p_value_adjust
  )


##Ctrl vs BA52
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
    fc_ctrl_ba52 = fc,
    p_value_ctrl_ba52 = p_value,
    p_value_adjust_ctrl_ba52 = p_value_adjust
  )


##Ctrl vs BQ11
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
    fc_ctrl_bq11 = fc,
    p_value_ctrl_bq11 = p_value,
    p_value_adjust_ctrl_bq11 = p_value_adjust
  )


###Volcano plot
###Ctrl vs WT
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ctrl_wt",
    p_value_column_name = "p_value_adjust_ctrl_wt",
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
       filename = "volcano_plot_ctrl_wt.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_ctrl_wt.pdf",
       width = 7,
       height = 6)

###Ctrl vs BA52
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ctrl_ba52",
    p_value_column_name = "p_value_adjust_ctrl_ba52",
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
       filename = "volcano_plot_ctrl_ba52.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_ctrl_ba52.pdf",
       width = 7,
       height = 6)

###Ctrl vs BQ11
plot <-
  volcano_plot(
    metabolome_data,
    fc_column_name = "fc_ctrl_bq11",
    p_value_column_name = "p_value_adjust_ctrl_bq11",
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
       filename = "volcano_plot_ctrl_bq11.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_ctrl_bq11.pdf",
       width = 7,
       height = 6)


###Heatmap to show the DEGs in different comparisons
temp <-
  metabolome_data %>%
  extract_variable_info() %>%
  dplyr::filter(
    p_value_adjust_ctrl_wt < 0.05 |
      p_value_adjust_ctrl_ba52 < 0.05 |
      p_value_adjust_ctrl_bq11 < 0.05
  )

p_value_adjust <-
  temp %>%
  dplyr::select(contains("p_value_adjust"))

fc <-
  temp %>%
  dplyr::select(contains("fc"))

temp <-
  log(fc, 2)

temp[which(temp > 0, arr.ind = TRUE)] <- 1
temp[which(temp < 0, arr.ind = TRUE)] <- -1

temp <-
  apply(temp, 2, function(x) {
    as.character(x)
  })

temp[which(p_value_adjust > 0.05, arr.ind = TRUE)] <- NA

###Heatmap
library(ComplexHeatmap)

colors = structure(c("blue", "red"), names = c("-1", "1"))

plot <-
  Heatmap(
    matrix = t(temp),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = colors,
    show_row_names = TRUE,
    show_column_names = FALSE,
    name = "Marker",
    row_title = "",
    column_title = "Gene",
    border = TRUE
  )


plot <-
  ggplotify::as.ggplot(plot)

plot

ggsave(plot,
       filename = "marker_overlap.pdf",
       width = 10,
       height = 3)

ggsave(plot,
       filename = "marker_overlap.png",
       width = 10,
       height = 3)

final_result <-
  temp %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    x[is.na(x)] <- "0"
    data.frame(class = c("No", "Down", "UP"),
               number = c(sum(x == "0"), sum(x == "-1"), sum(x == "1")))
  })

####genes are not consistant
final_result %>%
  purrr::map(function(x) {
    x$number[2] > 0 & x$number[3] > 0
  }) %>%
  unlist() %>%
  sum()

final_result %>%
  purrr::map(function(x) {
    x$number[2] >= 2 & x$number[3] == 0
  }) %>%
  unlist() %>%
  sum()

final_result %>%
  purrr::map(function(x) {
    x$number[2] == 0 & x$number[3] >= 2
  }) %>%
  unlist() %>%
  sum()

final_result %>%
  purrr::map(function(x) {
    x$number[1] == 2
  }) %>%
  unlist() %>%
  sum()
