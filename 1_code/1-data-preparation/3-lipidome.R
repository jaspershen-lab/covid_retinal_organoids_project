library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

###read data
data <- readxl::read_xlsx("2-data/2. ye organoid lipid all useful features list.xlsx") %>%
  as.data.frame()

dir.create("3-data_analysis/1-data-preparation/3-lipidome/",
           showWarnings = "FALSE")

setwd("3-data_analysis/1-data-preparation/3-lipidome/")

###remove some genes
expression_data <-
  data %>%
  dplyr::select(`Mock-1`:`QC-7`) %>%
  as.data.frame()

variable_id <-
  paste("lipid", 1:nrow(data), sep = "_")

rownames(expression_data) <-
  variable_id

variable_info <-
  data  %>%
  dplyr::select(`Alignment ID`:`Post curation result`) %>%
  dplyr::mutate(variable_id = variable_id) %>%
  dplyr::rename(
    rt = `Average Rt(min)`,
    mz = `Average Mz`,
    Compound.name = `Metabolite name`,
    Adduct = `Adduct type`,
    Post_curation_result = `Post curation result`
  ) %>%
  dplyr::select(variable_id, everything()) %>%
  dplyr::mutate(rt = rt * 60)


colnames(expression_data) <-
  colnames(expression_data) %>%
  stringr::str_replace("-", "_") %>%
  stringr::str_replace("Mock", "Ctrl")


sample_info <-
  data.frame(sample_id = colnames(expression_data)) %>%
  dplyr::mutate(class = "Subject") %>%
  dplyr::mutate(group = stringr::str_replace(sample_id, "_[0-9]{1,2}", ""))

sample_info$class[sample_info$group == "QC"] <- "QC"

library(massdataset)
library(tidymass)

lipidome_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )

BA52_id <-
  lipidome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BA52") %>%
  dplyr::pull(sample_id)

BQ11_id <-
  lipidome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BQ11") %>%
  dplyr::pull(sample_id)

WT_id <-
  lipidome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "WT") %>%
  dplyr::pull(sample_id)

Ctrl_id <-
  lipidome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Ctrl") %>%
  dplyr::pull(sample_id)

lipidome_data <-
  lipidome_data %>%
  mutate_variable_zero_freq() %>%
  mutate_variable_zero_freq(according_to_samples = BA52_id) %>%
  mutate_variable_zero_freq(according_to_samples = BQ11_id) %>%
  mutate_variable_zero_freq(according_to_samples = WT_id) %>%
  mutate_variable_zero_freq(according_to_samples = Ctrl_id) %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(zero_freq.1 < 0.5 |
           zero_freq.2 < 0.5 |
           zero_freq.3 < 0.5 | zero_freq.4 < 0.5)

variable_info <-
  extract_variable_info(lipidome_data)

variable_info$note <-
  variable_info$Compound.name %>%
  stringr::str_extract("\\[MS2 confirm\\]|w/o MS2")

variable_info$Compound.name <-
  variable_info$Compound.name %>%
  stringr::str_replace("\\[MS2 confirm\\]|w/o MS2:", "") %>%
  stringr::str_trim(side = "both")

lipidome_data@variable_info <-
  variable_info

save(lipidome_data, file = "lipidome_data.RData")
