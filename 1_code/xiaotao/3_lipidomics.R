library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(tidyverse)
library(tidymass)

##read data
load("3_data_analysis/2_data_cleaning/3_lipidome/lipidome_data.RData")

dir.create(
  "3_data_analysis/xiaotao/2_lipidomics",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/xiaotao/2_lipidomics")

####remove some lipids
lipidome_data <-
  lipidome_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(lipid_class)) %>%
  dplyr::filter(!lipid_class %in% c("maybe", "RIKEN"))

temp_data <-
  unique(lipidome_data@variable_info$lipid_class) %>%
  purrr::map(function(x) {
    idx <-
      which(lipidome_data@variable_info$lipid_class == x)
    lipidome_data[idx, ] %>%
      scale_data(center = FALSE) %>%
      extract_expression_data() %>%
      colMeans()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(temp_data) <-
  unique(lipidome_data@variable_info$lipid_class)

temp_data <-
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

range(temp_data)


library(ComplexHeatmap)

# colors = structure(c("blue", "red"), names = c("-3", "3"))

library(circlize)

# colors <-
#   circlize::colorRamp2(c(-3, 3), c("blue", "red"))

plot <-
  Heatmap(
    matrix = temp_data,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    # col = colors,
    show_row_names = TRUE,
    show_column_names = TRUE,
    name = "Marker",
    row_title = "",
    column_title = "",
    border = TRUE
  )

plot <-
  ggplotify::as.ggplot(plot)

plot
