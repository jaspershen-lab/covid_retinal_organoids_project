library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

###metabolomics data
load("3_data_analysis/1_data_preparation/2-metabolome/metabolome_data.rda")

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