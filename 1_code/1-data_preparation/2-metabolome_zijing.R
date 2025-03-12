library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

###read data from both sheets
data_neg <- readxl::read_xlsx("2_data/3. Target metabolomics for eye organoids_report_v1.xlsx", 
                              sheet = 1) %>%
  as.data.frame()

data_pos <- readxl::read_xlsx("2_data/3. Target metabolomics for eye organoids_report_v1.xlsx", 
                              sheet = 2) %>%
  as.data.frame()

dir.create("3_data_analysis/1-data_preparation/2-metabolome/",
           recursive = TRUE,
           showWarnings = "FALSE")
setwd("3_data_analysis/1-data_preparation/2-metabolome/")

# Function to process sample names
process_sample_names <- function(names) {
  names %>%
    stringr::str_replace("Eye Organoid-", "") %>%
    stringr::str_replace("Eye organoid-", "") %>%
    stringr::str_replace("-", "_") %>%
    stringr::str_replace("B5", "BA52") %>%
    stringr::str_replace("BQ", "BQ11") %>%
    stringr::str_replace("M", "Ctrl") %>%
    stringr::str_replace("BLK smaple", "BLK")
}

# Process NEG data
data_neg <- data_neg %>%
  dplyr::filter(!is.na(`Sample Name`)) %>%
  dplyr::filter(stringr::str_detect(`Sample Name`, "Eye"))
data_neg$`Sample Name` <- process_sample_names(data_neg$`Sample Name`)

# Process POS data
data_pos <- data_pos %>%
  dplyr::filter(!is.na(`Sample Name`)) %>%
  dplyr::filter(stringr::str_detect(`Sample Name`, "Eye"))
data_pos$`Sample Name` <- process_sample_names(data_pos$`Sample Name`)

# Function to create expression data
create_expression_data <- function(data) {
  data %>%
    dplyr::select(-`Sample Name`) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    apply(2, as.numeric) %>%
    as.data.frame()
}

# Create expression data for both modes
expression_data_neg <- create_expression_data(data_neg)
expression_data_pos <- create_expression_data(data_pos)

rownames(expression_data_neg) <- colnames(data_neg)[-1]
rownames(expression_data_pos) <- colnames(data_pos)[-1]
colnames(expression_data_neg) <- data_neg$`Sample Name`
colnames(expression_data_pos) <- data_pos$`Sample Name`

# Create combined variable info with mode suffixes
variable_id_neg <- paste("metabolite", 1:nrow(expression_data_neg), "_NEG", sep = "_")
variable_id_pos <- paste("metabolite", 1:nrow(expression_data_pos), "_POS", sep = "_")

# Combine expression data
expression_data_combined <- rbind(expression_data_neg, expression_data_pos)

# Create combined variable info
variable_info_neg <- data.frame(
  variable_id = variable_id_neg,
  Compound.name = paste0(rownames(expression_data_neg), "_NEG")
)

variable_info_pos <- data.frame(
  variable_id = variable_id_pos,
  Compound.name = paste0(rownames(expression_data_pos), "_POS")
)

variable_info_combined <- rbind(variable_info_neg, variable_info_pos)
rownames(expression_data_combined) <- variable_info_combined$variable_id

# Create sample info
sample_info <- data.frame(sample_id = colnames(expression_data_combined)) %>%
  dplyr::mutate(class = "Subject") %>%
  dplyr::mutate(group = stringr::str_replace(sample_id, "_[0-9]{1,2}", ""))

sample_info$group[sample_info$group == "BLK"] <- "Blank"
sample_info$class[sample_info$group == "Blank"] <- "Blank"

# Replace NA with 0
expression_data_combined[which(is.na(expression_data_combined), arr.ind = TRUE)] <- 0

# Add HMDB and KEGG IDs for both modes
variable_info_combined$HMDB <- c(
  # NEG mode IDs
  c("HMDB0000133", "HMDB0000300", "HMDB0001401",
    "HMDB0000169", "HMDB0000122", "HMDB0000124",
    "HMDB0000296", "HMDB0000085", "HMDB0000273",
    "HMDB0000132", "HMDB0006483", "HMDB0000192",
    "HMDB0011185", "HMDB0000014", "HMDB0000653",
    "HMDB0001565"),
  # Add POS mode IDs here
  # You'll need to add the correct HMDB IDs for POS mode metabolites
  rep(NA, nrow(expression_data_pos))
)

variable_info_combined$KEGG <- c(
  # NEG mode IDs
  c("C00387", "C00106", "C00092",
    "C00936", "C00221", "C00085",
    "C00299", "C00330", "C00214",
    "C00242", "C00402", "C00491",
    "C12147", "C00881", "C18043",
    "C00588"),
  # Add POS mode IDs here
  # You'll need to add the correct KEGG IDs for POS mode metabolites
  rep(NA, nrow(expression_data_pos))
)

# Create mass dataset
metabolome_data <- create_mass_dataset(
  expression_data = expression_data_combined,
  variable_info = variable_info_combined,
  sample_info = sample_info
)

save(metabolome_data, file = "metabolome_data.rda")
