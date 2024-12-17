library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

###read data
data <- readxl::read_xlsx("2_data/3. Target metabolomics for eye organoids_report_v1.xlsx") %>%
  as.data.frame()

dir.create("3_data_analysis/1_data_preparation/2-metabolome/",
           showWarnings = "FALSE")

setwd("3_data_analysis/1_data_preparation/2-metabolome/")

data <-
  data %>%
  dplyr::filter(!is.na(`Sample Name`)) %>%
  dplyr::filter(stringr::str_detect(`Sample Name`, "Eye"))

data$`Sample Name` <-
  data$`Sample Name` %>%
  stringr::str_replace("Eye Organoid-", "") %>%
  stringr::str_replace("Eye organoid-", "") %>%
  stringr::str_replace("-", "_") %>%
  stringr::str_replace("B5", "BA52") %>%
  stringr::str_replace("BQ", "BQ11") %>%
  stringr::str_replace("M", "Ctrl") %>%
  stringr::str_replace("BLK smaple", "BLK")

###remove some genes
expression_data <-
  data %>%
  dplyr::select(-`Sample Name`) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  apply(2, as.numeric) %>%
  as.data.frame()

rownames(expression_data) <- colnames(data)[-1]

colnames(expression_data) <- data$`Sample Name`

variable_id <-
  paste("metabolite", 1:nrow(expression_data), sep = "_")

variable_info <-
  data.frame(variable_id = variable_id,
             Compound.name = rownames(expression_data))

rownames(expression_data) <-
  variable_id

sample_info <-
  data.frame(sample_id = colnames(expression_data)) %>%
  dplyr::mutate(class = "Subject") %>%
  dplyr::mutate(group = stringr::str_replace(sample_id, "_[0-9]{1,2}", ""))

sample_info$group[sample_info$group == "BLK"] <- "Blank"

sample_info$class[sample_info$group == "Blank"] <- "Blank"
sample_info$class[sample_info$group == "Blank"] <- "Blank"

expression_data[which(is.na(expression_data), arr.ind = TRUE)] <- 0

library(massdataset)
library(tidymass)

variable_info

library(masstools)

variable_info$HMDB <-
  c("HMDB0000133", "HMDB0000300", "HMDB0001401",
    "HMDB0000169", "HMDB0000122", "HMDB0000124",
    "HMDB0000296", "HMDB0000085", "HMDB0000273",
    "HMDB0000132", "HMDB0006483", "HMDB0000192",
    "HMDB0011185", "HMDB0000014", "HMDB0000653",
    "HMDB0001565")

variable_info$KEGG <-
  c("C00387", "C00106", "C00092",
    "C00936", "C00221", "C00085",
    "C00299", "C00330", "C00214",
    "C00242", "C00402", "C00491",
    "C12147", "C00881", "C18043",
    "C00588")

metabolome_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )

# variable_info <-
#   extract_variable_info(metabolome_data)
# 
# variable_info$note <-
#   variable_info$Compound.name %>%
#   stringr::str_extract("\\[MS2 confirm\\]|w/o MS2")
# 
# variable_info$Compound.name <-
#   variable_info$Compound.name %>%
#   stringr::str_replace("\\[MS2 confirm\\]|w/o MS2:", "") %>%
#   stringr::str_trim(side = "both")
# 
# metabolome_data@variable_info <-
#   variable_info

save(metabolome_data, file = "metabolome_data.rda")
