library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

###read data
data <- readr::read_csv(
  "2-data/1. RO infected with SARS-CoV-2 strains  RNA-seq data all_counts_with_symbol_new.csv"
) %>%
  as.data.frame()

data <-
  data %>%
  dplyr::filter(!is.na(symbol))

setwd("3-data_analysis/1-data-preparation/1-transcriptome/")

###remove some genes
expression_data <-
  data %>%
  dplyr::select(-symbol) %>%
  as.data.frame()

variable_id <-
  paste("gene", 1:nrow(data), sep = "_")

rownames(expression_data) <-
  variable_id

variable_info <-
  data  %>%
  dplyr::select(symbol) %>%
  dplyr::mutate(variable_id = variable_id) %>%
  dplyr::mutate(ENSEMBL = case_when(
    stringr::str_detect(symbol, "ENSG") ~ symbol,
    TRUE ~ NA_character_
  )) %>%
  dplyr::mutate(SYMBOL = case_when(
    stringr::str_detect(symbol, "ENSG") ~ NA_character_,
    TRUE ~ symbol
  )) %>%
  dplyr::select(variable_id, SYMBOL, ENSEMBL)

sample_info <-
  data.frame(sample_id = colnames(expression_data)) %>%
  dplyr::mutate(class = "Subject") %>%
  dplyr::mutate(group = stringr::str_replace(sample_id, "_[0-9]{1,2}", ""))

library(massdataset)
library(tidymass)

transcriptome_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )

BA.5.2_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BA.5.2") %>%
  dplyr::pull(sample_id)

BQ.1.1_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "BQ.1.1") %>%
  dplyr::pull(sample_id)

WT_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "WT") %>%
  dplyr::pull(sample_id)

Ctrl_id <-
  transcriptome_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Ctrl") %>%
  dplyr::pull(sample_id)

transcriptome_data <-
  transcriptome_data %>%
  mutate_variable_zero_freq() %>%
  mutate_variable_zero_freq(according_to_samples = BA.5.2_id) %>%
  mutate_variable_zero_freq(according_to_samples = BQ.1.1_id) %>%
  mutate_variable_zero_freq(according_to_samples = WT_id) %>%
  mutate_variable_zero_freq(according_to_samples = Ctrl_id) %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(zero_freq.1 < 0.5 |
           zero_freq.2 < 0.5 |
           zero_freq.3 < 0.5 | zero_freq.4 < 0.5)


variable_info <-
  extract_variable_info(transcriptome_data)

idx1 <- which(is.na(variable_info$SYMBOL))
idx2 <- which(is.na(variable_info$ENSEMBL))

intersect(idx1, idx2)

###give more information to the genes
library(clusterProfiler)
library(org.Hs.eg.db)

ENSEMBL <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(ENSEMBL, by = "SYMBOL") %>%
  dplyr::mutate(ENSEMBL = case_when(is.na(ENSEMBL.x) ~ ENSEMBL.y, TRUE ~ ENSEMBL.x)) %>%
  dplyr::select(-c(ENSEMBL.x, ENSEMBL.y))

SYMBOL <-
  clusterProfiler::bitr(
    geneID = variable_info$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(SYMBOL, by = "ENSEMBL") %>%
  dplyr::mutate(SYMBOL = case_when(is.na(SYMBOL.x) ~ SYMBOL.y, TRUE ~ SYMBOL.x)) %>%
  dplyr::select(-c(SYMBOL.x, SYMBOL.y))


###ENTREZID
ENTREZID1 <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

ENTREZID2 <-
  clusterProfiler::bitr(
    geneID = variable_info$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(ENTREZID1, by = "SYMBOL") %>%
  dplyr::left_join(ENTREZID2, by = "ENSEMBL") %>%
  dplyr::mutate(ENTREZID = case_when(is.na(ENTREZID.x) ~ ENTREZID.y, TRUE ~ ENTREZID.x)) %>%
  dplyr::select(-c(ENTREZID.x, ENTREZID.y))


###GENENAME
GENENAME1 <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "GENENAME",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

GENENAME2 <-
  clusterProfiler::bitr(
    geneID = variable_info$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "GENENAME",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(GENENAME1, by = "SYMBOL") %>%
  dplyr::left_join(GENENAME2, by = "ENSEMBL") %>%
  dplyr::mutate(GENENAME = case_when(is.na(GENENAME.x) ~ GENENAME.y, TRUE ~ GENENAME.x)) %>%
  dplyr::select(-c(GENENAME.x, GENENAME.y))





###GENETYPE
GENETYPE1 <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "GENETYPE",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

GENETYPE2 <-
  clusterProfiler::bitr(
    geneID = variable_info$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "GENETYPE",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(GENETYPE1, by = "SYMBOL") %>%
  dplyr::left_join(GENETYPE2, by = "ENSEMBL") %>%
  dplyr::mutate(GENETYPE = case_when(is.na(GENETYPE.x) ~ GENETYPE.y, TRUE ~ GENETYPE.x)) %>%
  dplyr::select(-c(GENETYPE.x, GENETYPE.y))







###UNIPROT
UNIPROT1 <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "UNIPROT",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

UNIPROT2 <-
  clusterProfiler::bitr(
    geneID = variable_info$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "UNIPROT",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(UNIPROT1, by = "SYMBOL") %>%
  dplyr::left_join(UNIPROT2, by = "ENSEMBL") %>%
  dplyr::mutate(UNIPROT = case_when(is.na(UNIPROT.x) ~ UNIPROT.y, TRUE ~ UNIPROT.x)) %>%
  dplyr::select(-c(UNIPROT.x, UNIPROT.y))


transcriptome_data@variable_info <-
  variable_info

save(transcriptome_data, file = "transcriptome_data.RData")
