library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

##read data
load("3-data_analysis/2-data-cleaning/2-metabolome/metabolome_data.RData")

dir.create(
  "3-data_analysis/3-data-overview/2-metabolome",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3-data_analysis/3-data-overview/2-metabolome")

##PCA analysis
library(tidymass)

dim(metabolome_data)

metabolome_data <-
  metabolome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class != "Blank")

pca_object <-
  metabolome_data %>%
  scale_data() %>%
  massstat::run_pca()

plot <-
  pca_score_plot(
    object = metabolome_data,
    pca_object = pca_object,
    color_by = "group",
    frame = FALSE
  ) +
  scale_color_manual(values = group_color) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3)

plot

ggsave(plot,
       file = "pca_score_plot.pdf",
       width = 7,
       height = 6)
ggsave(plot,
       file = "pca_score_plot.png",
       width = 7,
       height = 6)
