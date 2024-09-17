library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

##read data

load("3-data_analysis/1-data-preparation/3-lipidome/lipidome_data.RData")

dir.create("3-data_analysis/2-data-cleaning/3-lipidome", showWarnings = FALSE)

setwd("3-data_analysis/2-data-cleaning/3-lipidome/")

variable_info <-
  extract_variable_info(lipidome_data)

variable_info <-
  variable_info %>%
  dplyr::mutate(All_identifications = Compound.name)

variable_info$note[grep("RIKEN", variable_info$Compound.name)]
grep("RIKEN", variable_info$Compound.name, value = TRUE)

variable_info$All_identifications <-
  variable_info$All_identifications %>%
  stringr::str_replace_all("TG", "TAG")

Compound.name <-
  stringr::str_split(variable_info$All_identifications, ";|[|]") %>%
  lapply(function(x) {
    x[1]
  }) %>%
  unlist()

variable_info$Compound.name <- Compound.name

lipid_class <-
  variable_info$Compound.name %>%
  stringr::str_split(" ") %>%
  lapply(function(x) {
    x[1] %>%
      stringr::str_extract("^[A-Za-z]+")
  }) %>%
  unlist()

variable_info$lipid_class <- lipid_class

lipidome_data@variable_info <- variable_info

###data normalization
lipidome_data <-
  lipidome_data %>% 
  normalize_data(method = "median")

save(lipidome_data, file = "lipidome_data.RData")
