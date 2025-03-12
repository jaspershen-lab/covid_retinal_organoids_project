library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(tidymass)
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggridges)



#devtools::install_github("davidsjoberg/ggsankey")
#devtools::install_github("junjunlab/GseaVis")
library(GseaVis)


##read data
load("3_data_analysis/4_differential_analysis/1-transcriptome/transcriptome_data")


if (!dir.exists("3_data_analysis/5-pathway_analysis/GSEA/ba52")) {
  dir.create(
    "3_data_analysis/5-pathway_analysis/GSEA/ba52",
    recursive = TRUE,
    showWarnings = FALSE
  )
}



setwd("3_data_analysis/5-pathway_analysis/GSEA/bq11")


###data preparation####
# Extract SYMBOL and log2FoldChange
data <- transcriptome_data@variable_info[, c("SYMBOL", "fc_bq11_wt")]
head(data)


transcriptome_data %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  dplyr::filter(fc_ba52_wt > 300) %>% 
  extract_expression_data() %>% 
  view()

transcriptome_data %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  dplyr::filter(fc_ba52_wt == 0) %>% 
  extract_expression_data() %>% 
  view()

# Use bitr to convert gene ids and merge
data_id <- bitr(data$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
data_all <- merge(data, data_id, by = "SYMBOL")
head(data_all)

data_all <-
  data_all %>% 
  dplyr::arrange(desc(fc_bq11_wt))

# # Create the named vector geneList, delete NA and sort
# filtered_data <- data_all[abs(data_all$fc_ba52_wt) > 2, ]
# head(filtered_data )


# # Display duplicate values
# duplicate_fc_filtered <- filtered_data %>%
#   group_by(fc_ba52_wt) %>%
#   filter(n() > 1) %>%
#   summarize(count = n())
# 
# if(nrow(duplicate_fc_filtered) > 0) {
#   print("There are duplicate fc_ba52_wt values in filtered_data:")
#   print(duplicate_fc_filtered, n = Inf)  # Use n = Inf to display all rows
# } else {
#   print("There are no duplicate fc_ba52_wt values in filtered_data.")
# }

geneList <- setNames(data_all$fc_bq11_wt, data_all$ENTREZID)
# geneList <- sort(geneList, decreasing = TRUE)
summary(geneList)
head(geneList)

geneList[geneList == 0] <- 
  min(geneList[geneList != 0])
summary(geneList)

#####GO#####
#plot(log(geneList, 2))
set.seed(1234) 
go_results <- gseGO(geneList = log(geneList, 2),
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE,
                    keyType = "ENTREZID")

head(go_results)
go_results@result %>% dim()
go_results@result$ID


if (!dir.exists("GO")) {
  dir.create("GO")
}

write.csv(go_results, file = c('GO/GSEA_GO_BQ.1.1.csv'))

go_ids <- go_results@result$ID
for (go_id in go_ids) {
  plot <- gseaNb(object = go_results, geneSetID = go_id)
  png_filename <- paste0("GO/GSEA_", go_id, "_BQ.1.1.png")
  pdf_filename <- paste0("GO/GSEA_", go_id, "_BQ.1.1.pdf")
  ggsave(plot, filename = png_filename, width = 5, height = 4)
  ggsave(plot, filename = pdf_filename, width = 5, height = 4)
}


#gseaplot(x = go_results, geneSetID = "GO:0045211")

#p1<- gseaplot2(go_results, 
#               geneSetID= "GO:0045211", 
#                color= 'red',
#                rel_heights= c(1.5, 0.5, 1), 
#                subplots= 1:3,
#                pvalue_table= F, 
#                title= go_results$Description[18],
#                ES_geom= 'line') #'dot'
#p1




#####KEGG#####
kegg_results <- gseKEGG(geneList     = geneList,
                        organism     = 'hsa',
                        minGSSize    = 3,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)

head(kegg_results)

kegg_results@result %>% dim()
kegg_results@result$ID

if (!dir.exists("KEGG")) {
  dir.create("KEGG")
}

write.csv(kegg_results, file = c('KEGG/GSEA_KEGG_BQ.1.1.csv'))


kegg_ids <- kegg_results@result$ID
for (kegg_id in kegg_ids) {
  plot <- gseaNb(object = kegg_results, geneSetID = kegg_id)
  png_filename <- paste0("KEGG/GSEA_", kegg_id, "_BQ.1.1.png")
  pdf_filename <- paste0("KEGG/GSEA_", kegg_id, "_BQ.1.1.pdf")
  ggsave(plot, filename = png_filename, width = 5, height = 4)
  ggsave(plot, filename = pdf_filename, width = 5, height = 4)
}




#######Reactome#####
library(ReactomePA)
library(pathview)

rec_results <- gsePathway(geneList, 
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", 
                          verbose = FALSE)
head(rec_results)

if (!dir.exists("REACTOME")) {
  dir.create("REACTOME")
}

rec_results@result %>% dim()

write.csv(go_results, file = c('REACTOME/GSEA_Reactome_BQ.1.1.csv'))

rec_ids <- rec_results@result$ID
for (rec_id in rec_ids) {
  plot <- gseaNb(object = rec_results, geneSetID = rec_id)
  png_filename <- paste0("REACTOME/GSEA_", rec_id, "_BQ.1.1.png")
  pdf_filename <- paste0("REACTOME/GSEA_", rec_id, "_BQ.1.1.pdf")
  ggsave(plot, filename = png_filename, width = 5, height = 4)
  ggsave(plot, filename = pdf_filename, width = 5, height = 4)
}
