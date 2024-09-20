library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

##read data

load("3-data_analysis/1-data-preparation/2-metabolome/metabolome_data.RData")

dir.create("3-data_analysis/2-data-cleaning/2-metabolome",
           showWarnings = FALSE)

setwd("3-data_analysis/2-data-cleaning/2-metabolome/")

variable_info <-
  extract_variable_info(metabolome_data)

library(tidymass)

KEGG_ID <-
  c(
    "guanosine" = "C00387",
    "Uracil" = "C00106",
    "D-glucose-6-phosphate" = "C00092",
    "mannose" = "C00159",
    "glucose" = "C00031",
    "D-Fructose-6-phosphate" = "C05345",
    "Uridine" = "C00299",
    "2-deoxyguanosine" = "C00330",
    "Thymidine" = "C00214",
    "Guanine" = "C00242",
    "Aspartic acid" = "C00049",
    "Cystine" = "C00491",
    "O-phospho-L-threonine" = "C12147",
    "2-deoxycytidine" = "C00881",
    "Cholesterol sulphate" = "C18043",
    "Phospho-choline" = "C00588"
  )

names(KEGG_ID) == variable_info$Compound.name

variable_info$KEGG_ID <- KEGG_ID

###HMDB ID
library(tidymass)

HMDB_ID <-
  variable_info$KEGG_ID %>%
  lapply(function(x) {
    trans_ID(x, "KEGG", "Human Metabolome Database", server = "cts.fiehnlab")
  })

HMDB_ID <-
  HMDB_ID %>%
  do.call(rbind, .) %>%
  as.data.frame()

variable_info <-
  variable_info %>%
  dplyr::left_join(HMDB_ID, by = c("KEGG_ID" = "KEGG")) %>%
  dplyr::rename("HMDB_ID" = "Human Metabolome Database")


metabolome_data@variable_info <- variable_info

###data normalization
metabolome_data <-
metabolome_data %>% 
  normalize_data(method = "median")

save(metabolome_data, file = "metabolome_data.RData")
