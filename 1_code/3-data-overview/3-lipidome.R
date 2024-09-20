library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

##read data
load("3-data_analysis/2-data-cleaning/3-lipidome/lipidome_data.RData")

dir.create(
  "3-data_analysis/3-data-overview/3-lipidome",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3-data_analysis/3-data-overview/3-lipidome")

##PCA analysis
library(tidymass)

dim(lipidome_data)

lipidome_data <-
  lipidome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class != "Blank")
# filter(class != "QC")

pca_object <-
  lipidome_data %>%
  scale_data() %>%
  massstat::run_pca()

plot <-
  pca_score_plot(
    object = lipidome_data,
    pca_object = pca_object,
    color_by = "group",
    frame = FALSE
  ) +
  scale_color_manual(values = group_color) +
  scale_fill_manual(values = group_color) +
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
