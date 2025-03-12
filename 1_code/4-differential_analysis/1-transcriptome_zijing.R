library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(tidymass)

##read data
load("3_data_analysis/2-data-cleaning/1-transcriptome/transcriptome_data.RData")

dir.create(
  "3_data_analysis/4_differential_analysis/1-transcriptome",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/4_differential_analysis/1-transcriptome")

transcriptome_data@sample_info


transcriptome_data

ba52_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BA52") %>%
  pull(sample_id)

bq11_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BQ11") %>%
  pull(sample_id)

wt_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "WT") %>%
  pull(sample_id)

#####ba52 VS wt
transcriptome_data <-
  transcriptome_data %>%
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

head(transcriptome_data@variable_info)



#####bq11 VS wt
transcriptome_data <-
  transcriptome_data %>%
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

head(transcriptome_data@variable_info)


###Volcano plot
###ba52 vs WT
plot <-
  volcano_plot(
    transcriptome_data,
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

###bq11 vs WT
plot <-
  volcano_plot(
    transcriptome_data,
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

save(transcriptome_data, file = "transcriptome_data")

