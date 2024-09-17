library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

##read data

load("3-data_analysis/1-data-preparation/1-transcriptome/transcriptome_data.RData")

dir.create("3-data_analysis/2-data-cleaning/1-transcriptome",
           showWarnings = FALSE)

setwd("3-data_analysis/2-data-cleaning/1-transcriptome/")

###remove low abundance genes
remain_idx <-
  apply(transcriptome_data@expression_data, 1, function(x) {
    sum(x > 10)
  }) %>%
  `>=`(1) %>%
  which()

transcriptome_data <-
  transcriptome_data[remain_idx, ]

###remove the genes who don't have all the IDs (ENTREZID, UNIPROT)
transcriptome_data <-
  transcriptome_data %>%
  activate_mass_dataset(what  = "variable_info") %>%
  dplyr::filter(!is.na(ENTREZID) | !is.na(UNIPROT))


library(biomaRt)

# Use Ensembl to get gene lengths
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene lengths (in base pairs)
gene_info <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "gene_biotype",
    "start_position",
    "end_position"
  ),
  filters = "ensembl_gene_id",
  values = transcriptome_data@variable_info$ENSEMBL,
  mart = mart
)

# Calculate gene lengths
gene_info <- gene_info %>%
  mutate(length = end_position - start_position + 1)

variable_info <-
  extract_variable_info(transcriptome_data)

variable_info <-
  variable_info %>%
  dplyr::left_join(gene_info[, c("ensembl_gene_id", "length")], by = c("ENSEMBL" = "ensembl_gene_id"))

transcriptome_data@variable_info <-
  variable_info

transcriptome_data <-
  transcriptome_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(length))


# Convert counts to FPKM
expression_data <-
  extract_expression_data(transcriptome_data)
variable_info <-
  extract_variable_info(transcriptome_data)

fpkm <- sweep(expression_data, 1, variable_info$length / 1000, FUN = "/") # Normalize by gene length (kb)
fpkm <- sweep(fpkm, 2, colSums(expression_data) / 1e6, FUN = "/") # Normalize by total counts (in millions)

# Convert FPKM to TPM
tpm <- sweep(fpkm, 2, colSums(fpkm) / 1e6, FUN = "/")

plot(as.numeric(expression_data[1, ]), as.numeric(tpm[1, ]))

plot(as.numeric(expression_data[2, ]), as.numeric(tpm[2, ]))

expression_data <-
  tpm

transcriptome_data@expression_data <-
  expression_data


expression_data <-
  extract_expression_data(transcriptome_data)

sample_info <-
  extract_sample_info(transcriptome_data)

colnames(expression_data) <-
  colnames(expression_data) %>%
  stringr::str_replace_all("\\.", "")

sample_info <-
  sample_info %>%
  dplyr::mutate(sample_id = sample_id %>% stringr::str_replace_all("\\.", "")) %>%
  dplyr::mutate(group = group %>% stringr::str_replace_all("\\.", ""))

transcriptome_data@expression_data <-
  expression_data

transcriptome_data@sample_info <-
  sample_info

save(transcriptome_data, file = "transcriptome_data.RData")
