library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

clean_lipid_data <- function(ion_mode) {
  load_path <- paste0("3_data_analysis/1-data_preparation/3-lipidome/3-lipidome_", ion_mode, "/lipidome_data_", ion_mode, ".RData")
  load(load_path)
  
  variable_info <- extract_variable_info(lipidome_data)
  variable_info <- variable_info %>%
    dplyr::mutate(All_identifications = original_name)
  
  variable_info$note[grep("RIKEN", variable_info$original_name)]
  grep("RIKEN", variable_info$original_name, value = TRUE)
  
  variable_info$All_identifications <-
    variable_info$All_identifications %>%
    stringr::str_replace_all("TG", "TAG")
  
  Compound.name <-
    stringr::str_split(variable_info$All_identifications, ";|[|]") %>%
    lapply(function(x) {
      x[1] %>%
        stringr::str_remove_all("_positive|_negative") %>%  # 移除已有后缀
        paste0("_", ifelse(ion_mode == "positive", "pos", "neg"))  # 添加新后缀
    }) %>%
    unlist()
  
  variable_info$Compound.name <- Compound.name
  
  lipid_class <-
    variable_info$Compound.name %>%
    stringr::str_remove(paste0("_", ion_mode)) %>%  # Remove suffix for class extraction
    stringr::str_split(" ") %>%
    lapply(function(x) {
      x[1] %>%
        stringr::str_extract("^[A-Za-z]+")
    }) %>%
    unlist()
  
  variable_info$lipid_class <- lipid_class
  variable_info$ion_mode <- ion_mode  # Add ion mode column
  
  lipidome_data@variable_info <- variable_info
  lipidome_data <- lipidome_data %>% 
    normalize_data(method = "median")
 
  
  return(lipidome_data)
}

# Clean both modes
lipidome_data_pos <- clean_lipid_data("positive")
lipidome_data_neg <- clean_lipid_data("negative")

dir.create(paste0("3_data_analysis/2-data_cleaning/3-lipidome/3-lipidome_all"),
           showWarnings = FALSE,
           recursive = TRUE)
setwd(paste0("3_data_analysis/2-data_cleaning/3-lipidome/3-lipidome_all"))



save(lipidome_data_pos, file =  "lipidome_data_pos_cleaned.RData")
save(lipidome_data_neg, file =  "lipidome_data_neg_cleaned.RData")




# Update duplicate information
pos_info <- extract_variable_info(lipidome_data_pos)
neg_info <- extract_variable_info(lipidome_data_neg)

common_compounds <- intersect(
  stringr::str_remove(pos_info$Compound.name, "_pos"),
  stringr::str_remove(neg_info$Compound.name, "_neg")
)

pos_info$duplicate <- pos_info$Compound.name %>%
  stringr::str_remove("_positive") %>%
  (`%in%`)(common_compounds)

neg_info$duplicate <- neg_info$Compound.name %>%
  stringr::str_remove("_negative") %>%
  (`%in%`)(common_compounds)

lipidome_data_pos@variable_info <- pos_info
lipidome_data_neg@variable_info <- neg_info

# Get expression data from both modes
pos_expr <- lipidome_data_pos@expression_data
neg_expr <- lipidome_data_neg@expression_data

# Add missing columns to positive mode data
missing_cols <- c("WT_4", "BQ11_3")
pos_expr[, missing_cols] <- NA

# Ensure column order matches between datasets
pos_expr <- pos_expr[, colnames(neg_expr)]

# Combine variable info
combined_info <- rbind(
  lipidome_data_pos@variable_info,
  lipidome_data_neg@variable_info
)

# Create combined massFeaturesObject
lipidome_data_combined <- create_mass_dataset(
  expression_data = rbind(pos_expr, neg_expr),
  sample_info = lipidome_data_neg@sample_info,
  variable_info = combined_info
)

# Save combined data
save(lipidome_data_combined, file = "lipidome_data_combined.RData")

