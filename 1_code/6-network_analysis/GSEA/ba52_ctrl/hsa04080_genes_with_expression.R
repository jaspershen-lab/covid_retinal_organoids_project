# 读取数据
library(tidyverse)
library(AnnotationDbi) 
setwd(get_project_wd())
# 读取GSEA结果
gsea_results <- read.csv("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/KEGG/GSEA_KEGG_ba52_ctrl.csv")

# 获取hsa04080通路的基因
hsa04080_genes <- gsea_results %>%
  filter(ID == "hsa04080") %>%
  pull(core_enrichment) %>%
  str_split("/") %>%
  unlist()

# 从transcriptome_data中提取表达量信息
expression_data <- transcriptome_data@variable_info %>%
  dplyr::select(SYMBOL, fc_ba52_ctrl) 

# 转换基因ID
gene_id_mapping <- bitr(hsa04080_genes,
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db)

# 合并基因信息和表达量
result_df <- gene_id_mapping %>%
  left_join(expression_data, by = "SYMBOL") %>%
  arrange(desc(fc_ba52_ctrl))

# 添加基因注释信息
gene_info <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = result_df$ENTREZID,
                                   columns = c("SYMBOL", "GENENAME"),
                                   keytype = "ENTREZID")

final_result <- result_df %>%
  left_join(gene_info, by = c("ENTREZID", "SYMBOL"))

# 保存结果
if (!dir.exists("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/KEGG")) {
  dir.create(
    "3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/KEGG",
    recursive = TRUE,
    showWarnings = FALSE
  )
}


setwd("3_data_analysis/5-pathway_analysis/GSEA/ba52_ctrl/KEGG")
write.csv(final_result, "hsa04080_genes_with_expression.csv", row.names = FALSE)

# 创建表达量可视化
ggplot(final_result, aes(x = reorder(SYMBOL, fc_ba52_ctrl), y = fc_ba52_ctrl)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Changes in hsa04080 Pathway",
       x = "Gene Symbol",
       y = "Fold Change (BA.5.2 vs Control)")

# 输出基本统计信息
cat("基因数量:", nrow(final_result), "\n")
cat("上调基因数量:", sum(final_result$fc_ba52_ctrl > 0), "\n")
cat("下调基因数量:", sum(final_result$fc_ba52_ctrl < 0), "\n")

