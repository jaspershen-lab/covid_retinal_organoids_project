library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

###metabolomics data
load("3_data_analysis/1_data_preparation/2-metabolome/metabolome_data.rda")

# ##qpcr data
# gene_data <-
#   readxl::read_xlsx("2_data/qpcr_data/2024-12-17 RO-CDNA-qPCR -WT and BA.5.2 (n=3) Yuanchen-2.xlsx")
#
# gene_data <-
#   gene_data[, c(1, 2, 9)]
#
# colnames(gene_data) <-
#   c("sample_name", "gene_name", "intensity")
#
# gene_data <-
#   gene_data %>%
#   dplyr::filter(!is.na(sample_name)) %>%
#   dplyr::mutate(intensity = as.numeric(intensity))
#
# gene_data <-
#   gene_data %>%
#   dplyr::mutate(
#     sample_id = case_when(
#       stringr::str_detect(sample_name, "Mock-A") ~ "Mock-A",
#       stringr::str_detect(sample_name, "Mock-B") ~ "Mock-B",
#       stringr::str_detect(sample_name, "Mock-C") ~ "Mock-C",
#       stringr::str_detect(sample_name, "Mock-1") ~ "Mock-A",
#       stringr::str_detect(sample_name, "Mock-2") ~ "Mock-B",
#       stringr::str_detect(sample_name, "Mock-3") ~ "Mock-C",
#       stringr::str_detect(sample_name, "WT-A") ~ "WT-A",
#       stringr::str_detect(sample_name, "WT-B") ~ "WT-B",
#       stringr::str_detect(sample_name, "WT-C") ~ "WT-C",
#       stringr::str_detect(sample_name, "BA.5.2-A") ~ "BA52-A",
#       stringr::str_detect(sample_name, "BA.5.2-B") ~ "BA52-B",
#       stringr::str_detect(sample_name, "BA.5.2-C") ~ "BA52-C",
#       stringr::str_detect(sample_name, "BQ.1.1-A") ~ "BQ11-A",
#       stringr::str_detect(sample_name, "BQ.1.1-B") ~ "BQ11-B",
#       stringr::str_detect(sample_name, "BQ.1.1-C") ~ "BQ11-C",
#       TRUE ~ "NA"
#     )
#   )  %>%
#   dplyr::arrange(sample_id, gene_name)
#
# gene_data <-
#   gene_data %>%
#   dplyr::mutate(temp_id = paste(sample_id, gene_name, sep = "_"))
#
# library(plyr)
#
# gene_data <-
#   gene_data %>%
#   plyr::dlply(.(temp_id), function(x) {
#     if (nrow(x) == 1) {
#       return(x)
#     }
#     x <-
#       x %>%
#       dplyr::filter(!is.na(intensity)) %>%
#       dplyr::filter(intensity == min(intensity))
#     return(x)
#   })  %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# gene_data <-
#   gene_data %>%
#   dplyr::select(-c(temp_id, sample_name)) %>%
#   tidyr::pivot_wider(names_from = sample_id, values_from = intensity)
#
# gene_data <-
#   gene_data %>%
#   tibble::column_to_rownames(var = "gene_name")
#
# sample_info <-
#   data.frame(sample_id = colnames(gene_data)) %>%
#   dplyr::mutate(class = "Subject") %>%
#   dplyr::mutate(
#     group = case_when(
#       stringr::str_detect(sample_id, "Mock") ~ "control",
#       stringr::str_detect(sample_id, "WT") ~ "WT",
#       stringr::str_detect(sample_id, "BA52") ~ "BA52",
#       TRUE ~ "NA"
#     )
#   )
#
# variable_info <-
#   data.frame(variable_id = rownames(gene_data))
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
#
# variable_info <-
#   variable_info %>%
#   dplyr::mutate(SYMBOL = stringr::str_replace(variable_id, "-", ""))
#
# variable_info$SYMBOL[variable_info$SYMBOL == "CathepsinB"] <- "CTSB"
# variable_info$SYMBOL[variable_info$SYMBOL == "CathepsinL"] <- "CTSL"
# variable_info$SYMBOL[variable_info$SYMBOL == "IFNa"] <- "IFNA1"
# variable_info$SYMBOL[variable_info$SYMBOL == "IFNB"] <- "IFNB1"
# variable_info$SYMBOL[variable_info$SYMBOL == "IFNGamma"] <- "IFNG"
# variable_info$SYMBOL[variable_info$SYMBOL == "IFNLambda1"] <- "IFNL1"
# variable_info$SYMBOL[variable_info$SYMBOL == "IL12"] <- "IL12A"
# variable_info$SYMBOL[variable_info$SYMBOL == "IL8"] <- "CXCL8"
# variable_info$SYMBOL[variable_info$SYMBOL == "IP10"] <- "CXCL10"
# variable_info$SYMBOL[variable_info$SYMBOL == "MCP1"] <- "CCL2"
# variable_info$SYMBOL[variable_info$SYMBOL == "MKi67"] <- "MKI67"
# variable_info$SYMBOL[variable_info$SYMBOL == "Pou4f2"] <- "POU4F2"
# variable_info$SYMBOL[variable_info$SYMBOL == "SAP97"] <- "DLG1"
# variable_info$SYMBOL[variable_info$SYMBOL == "TNFa"] <- "TNF"
#
# variable_info$ENTREZID <-
#   bitr(
#     variable_info$SYMBOL,
#     fromType = "SYMBOL",
#     toType = "ENTREZID",
#     OrgDb = org.Hs.eg.db,
#     drop = FALSE
#   ) %>%
#   pull(ENTREZID)
#
# expression_data <-
#   gene_data
#
# library(massdataset)
#
# gene_data <-
#   create_mass_dataset(
#     expression_data = expression_data,
#     variable_info = variable_info,
#     sample_info = sample_info
#   )

###kegg pathway
load("2_data/qpcr_data/pathway.rda")
pathway

dir.create("3_data_analysis/xiaotao/1_qpcr_data", recursive = TRUE)
setwd("3_data_analysis/xiaotao/1_qpcr_data")

# save(gene_data, file = "gene_data.rda")
load("gene_data.rda")


##map gene to pathways
mapped_result_gene <-
  pathway@gene_list %>%
  purrr::map(function(x) {
    if (nrow(x) == 0) {
      return(c(0, 0))
    }
    
    mapped_number <-
      x %>%
      dplyr::filter(KEGG.ID %in% gene_data@variable_info$ENTREZID) %>%
      nrow()
    
    c(mapped_number, mapped_number / nrow(x))
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(mapped_result_gene) <-
  c("number", "rate")

mapped_result_gene$pathway_id <-
  pathway@pathway_id

mapped_result_gene$pathway_name <-
  pathway@pathway_name

mapped_result_gene$pathway_class <-
  pathway@pathway_class

mapped_result_gene %>%
  dplyr::filter(number > 3) %>%
  dplyr::arrange(desc(number), desc(rate))

# ###map metabolites to pathways
# intersect(
#   metabolome_data@variable_info$KEGG,
#   pathway@compound_list %>%
#     purrr::map(~ .$KEGG.ID) %>%
#     unlist()
# )
#
#
# mapped_result_metabolite <-
#   pathway@compound_list %>%
#   purrr::map(function(x) {
#     if (nrow(x) == 0) {
#       return(c(0, 0))
#     }
#
#     mapped_number <-
#       x %>%
#       dplyr::filter(KEGG.ID %in% metabolome_data@variable_info$KEGG) %>%
#       nrow()
#
#     c(mapped_number, mapped_number / nrow(x))
#
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# colnames(mapped_result_metabolite) <-
#   c("number", "rate")
#
# mapped_result_metabolite$pathway_id <-
#   pathway@pathway_id
#
# mapped_result_metabolite$pathway_name <-
#   pathway@pathway_name
#
# mapped_result_metabolite$pathway_class <-
#   pathway@pathway_class
#
# mapped_result_metabolite %>%
#   dplyr::filter(number > 0) %>%
#   dplyr::arrange(desc(rate), desc(number))
#
#
# mapped_result_gene %>%
#   dplyr::left_join(mapped_result_metabolite,
#                    by = c("pathway_id", "pathway_name", "pathway_class")) %>%
#   dplyr::select(pathway_id,
#                 pathway_name,
#                 pathway_class,
#                 number.x,
#                 rate.x,
#                 number.y,
#                 rate.y) %>%
#   dplyr::rename(
#     gene_number = number.x,
#     gene_rate = rate.x,
#     metabolite_number = number.y,
#     metabolite_rate = rate.y
#   ) %>%
#   dplyr::filter(gene_number > 0 & metabolite_number > 0)
#