library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

process_lipid_data <- function(sheet_num, ion_mode) {
  original_dir <- getwd()  # 保存原始目录
  on.exit(setwd(original_dir))
  
  data <- readxl::read_xlsx("2_data/2. ye organoid lipid all useful features list.xlsx",
                            sheet = sheet_num) %>%
    as.data.frame()
  
  dir_name <- paste0("3_data_analysis/1-data_preparation/3-lipidome/3-lipidome_", ion_mode, "/")
  dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
  setwd(dir_name)
  
  expression_data <-
    data %>%
    dplyr::select(`Mock-1`:`QC-7`) %>%
    as.data.frame()
  
  variable_id <-
    paste("lipid", ion_mode, 1:nrow(data), sep = "_")
  rownames(expression_data) <- variable_id
  
  variable_info <-
    data %>%
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
  
  # Add ion mode suffix to Compound.name
  variable_info$Compound.name <- paste0(variable_info$Compound.name, "_", ion_mode)
  
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
  
  # Filter based on sample groups
  BA52_id <- lipidome_data %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::filter(group == "BA52") %>%
    dplyr::pull(sample_id)
  
  BQ11_id <- lipidome_data %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::filter(group == "BQ11") %>%
    dplyr::pull(sample_id)
  
  WT_id <- lipidome_data %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::filter(group == "WT") %>%
    dplyr::pull(sample_id)
  
  Ctrl_id <- lipidome_data %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::filter(group == "Ctrl") %>%
    dplyr::pull(sample_id)
  
  lipidome_data <- lipidome_data %>%
    mutate_variable_zero_freq() %>%
    mutate_variable_zero_freq(according_to_samples = BA52_id) %>%
    mutate_variable_zero_freq(according_to_samples = BQ11_id) %>%
    mutate_variable_zero_freq(according_to_samples = WT_id) %>%
    mutate_variable_zero_freq(according_to_samples = Ctrl_id) %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(zero_freq.1 < 0.5 |
             zero_freq.2 < 0.5 |
             zero_freq.3 < 0.5 |
             zero_freq.4 < 0.5)
  
  variable_info <- extract_variable_info(lipidome_data)
  variable_info$note <- variable_info$Compound.name %>%
    stringr::str_extract("\\[MS2 confirm\\]|w/o MS2")
  
  # Remove ion mode suffix for processing but keep original name
  variable_info$original_name <- variable_info$Compound.name
  variable_info$Compound.name <- variable_info$Compound.name %>%
    stringr::str_remove("_positive|_negative") %>%
    stringr::str_replace("\\[MS2 confirm\\]|w/o MS2:", "") %>%
    stringr::str_trim(side = "both")
  
  lipidome_data@variable_info <- variable_info
  
  save_name <- paste0("lipidome_data_", ion_mode, ".RData")
  save(lipidome_data, file = save_name)
  
  setwd("../../../")
  
  return(lipidome_data)
}

# Process both modes
lipidome_data_pos <- process_lipid_data(sheet_num = 1, ion_mode = "positive")
lipidome_data_neg <- process_lipid_data(sheet_num = 2, ion_mode = "negative")

# Save combined results
dir.create("3_data_analysis/1-data_preparation/3-lipidome/3-lipidome_combined/",
           showWarnings = FALSE, recursive = TRUE)
setwd("3_data_analysis/1-data_preparation/3-lipidome/3-lipidome_combined/")

# Get common compounds
pos_info <- extract_variable_info(lipidome_data_pos)
neg_info <- extract_variable_info(lipidome_data_neg)
common_compounds <- intersect(
  stringr::str_remove(pos_info$Compound.name, "_positive"),
  stringr::str_remove(neg_info$Compound.name, "_negative")
)

# Mark duplicates in variable info
pos_info$duplicate <- pos_info$Compound.name %in% common_compounds
neg_info$duplicate <- neg_info$Compound.name %in% common_compounds

# Update datasets
lipidome_data_pos@variable_info <- pos_info
lipidome_data_neg@variable_info <- neg_info

save(lipidome_data_pos, lipidome_data_neg,
     file = "lipidome_data_combined.RData")
