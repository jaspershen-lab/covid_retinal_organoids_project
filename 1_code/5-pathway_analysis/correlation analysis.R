library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(tidymass)
library(dplyr)

##read data
load("3_data_analysis/2-data-cleaning/1-transcriptome/transcriptome_data.RData")

amd <- read.table("2_data/RNA-Seq/RNA-Seq.Disease.Age-related_macular_degeneration.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(amd)

dr<-read.table("2_data/RNA-Seq/RNA-Seq.Disease.Diabetic_retinopathy.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(dr)

kt<-read.table("2_data/RNA-Seq/RNA-Seq.Disease.Keratoconus.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(kt)

poag<-read.table("2_data/RNA-Seq/RNA-Seq.Disease.Primary_open-angle_glaucoma.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(poag)

rp<-read.table("2_data/RNA-Seq/RNA-Seq.Disease.Retinitis_pigmentosa.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rp)

rb<-read.table("2_data/RNA-Seq/RNA-Seq.Disease.Retinoblastoma.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rb)

ctrl<- read.table("2_data/RNA-Seq/RNA-Seq.Normal.Retina.TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(ctrl)


#####Mean expression value
colnames(amd)
amd_avg <- amd %>%
  rowwise() %>%
  mutate(amd_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, amd_avg_expr)

colnames(dr)
dr_avg <- dr %>%
  rowwise() %>%
  mutate(dr_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, dr_avg_expr)

colnames(kt)
kt_avg <- kt %>%
  rowwise() %>%
  mutate(kt_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, kt_avg_expr)

colnames(poag)
poag_avg <- poag %>%
  rowwise() %>%
  mutate(poag_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, poag_avg_expr)

colnames(rp)
rp_avg <- rp %>%
  rowwise() %>%
  mutate(rp_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, rp_avg_expr)

colnames(rb)
rb_avg <- rb %>%
  rowwise() %>%
  mutate(rb_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, rb_avg_expr)

colnames(ctrl)
ctrl_avg <- ctrl %>%
  rowwise() %>%
  mutate(ctrl_avg_expr = mean(c_across(starts_with("GSM")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Symbol, ctrl_avg_expr)




####differential genes
dir.create(
  "3_data_analysis/5-pathway_analysis/correlation_analysis",
  recursive = TRUE,
  showWarnings = FALSE
)


setwd("3_data_analysis/5-pathway_analysis/correlation_analysis")

#### AMD
merged_data_amd <- merge(amd_avg, ctrl_avg, by = "Symbol")

#fold change & log2 fold change
#merged_data <- merged_data %>%
 # mutate(fold_change = amd_avg_expr / ctrl_avg_expr,
 #        log2_fold_change = log2(fold_change))

merged_data_amd <- merged_data_amd %>%
  mutate(fold_change = (amd_avg_expr + 1e-6) / (ctrl_avg_expr + 1e-6),
         log2_fold_change = log2(fold_change))

differential_genes_amd <- merged_data_amd %>%
  filter(abs(log2_fold_change) > 1)  


write.csv(differential_genes, file = c('amd_differential_genes.csv'))


####DR
merged_data_dr <- merge(dr_avg, ctrl_avg, by = "Symbol")

merged_data_dr <- merged_data_dr %>%
  mutate(fold_change = (dr_avg_expr + 1e-6) / (ctrl_avg_expr + 1e-6),
         log2_fold_change = log2(fold_change))

differential_genes_dr <- merged_data_dr %>%
  filter(abs(log2_fold_change) > 1)  

write.csv(differential_genes, file = c('dr_differential_genes.csv'))


####KT
merged_data_kt <- merge(kt_avg, ctrl_avg, by = "Symbol")

merged_data_kt <- merged_data_kt %>%
  mutate(fold_change = (kt_avg_expr + 1e-6) / (ctrl_avg_expr + 1e-6),
         log2_fold_change = log2(fold_change))

differential_genes_kt <- merged_data_kt %>%
  filter(abs(log2_fold_change) > 1)  

write.csv(differential_genes, file = c('kt_differential_genes.csv'))


####POAG
merged_data_poag <- merge(poag_avg, ctrl_avg, by = "Symbol")

merged_data_poag <- merged_data_poag %>%
  mutate(fold_change = (poag_avg_expr + 1e-6) / (ctrl_avg_expr + 1e-6),
         log2_fold_change = log2(fold_change))

differential_genes_poag <- merged_data_poag %>%
  filter(abs(log2_fold_change) > 1)  

write.csv(differential_genes, file = c('poag_differential_genes.csv'))


####rp
merged_data_rp <- merge(rp_avg, ctrl_avg, by = "Symbol")

merged_data_rp <- merged_data_rp %>%
  mutate(fold_change = (rp_avg_expr + 1e-6) / (ctrl_avg_expr + 1e-6),
         log2_fold_change = log2(fold_change))

differential_genes_rp <- merged_data_rp %>%
  filter(abs(log2_fold_change) > 1)  

write.csv(differential_genes, file = c('rp_differential_genes.csv'))


####rb
merged_data_rb <- merge(rb_avg, ctrl_avg, by = "Symbol")

merged_data_rb <- merged_data_rb %>%
  mutate(fold_change = (rb_avg_expr + 1e-6) / (ctrl_avg_expr + 1e-6),
         log2_fold_change = log2(fold_change))

differential_genes_rb <- merged_data_rb %>%
  filter(abs(log2_fold_change) > 1)  

write.csv(differential_genes, file = c('rb_differential_genes.csv'))





# 提取 variable_info 中需要的列
variable_info_subset <- transcriptome_data@variable_info[, c("SYMBOL", "ENSEMBL", "ENTREZID")]

expression_data <- transcriptome_data@expression_data

data <- cbind(variable_info_subset, expression_data)

head(data)

#####correlation analysis#####3
# 导入必要的库
library(dplyr)      # 用于数据操作
library(tidyr)      # 用于数据整理
library(ggplot2)    # 用于绘图
library(reshape2)   # 用于数据重塑
library(corrplot)   # 用于绘制相关性热图

# 获取上调的top30
top30_upregulated <- differential_genes_amd %>%
  arrange(desc(log2_fold_change)) %>%
  slice_head(n = 30)

# 获取下调的top30
top30_downregulated <- differential_genes_amd %>%
  arrange(log2_fold_change) %>%
  slice_head(n = 30)


# 提取共同基因
common_genes <- intersect(differential_genes$Symbol, data$SYMBOL) # 替换为实际基因列名

# 计算 AMD 数据集中共同基因的平均表达值
amd_avg <- amd %>%
  filter(Symbol %in% common_genes) %>%
  summarize(
    amd_avg_expr = mean(c_across(GSM3985312:GSM3191480), na.rm = TRUE),
    .by = "Symbol"
  )


# 计算 data 数据集中共同基因的平均表达值
data_avg <- data %>%
  filter(gene_symbol %in% common_genes) %>%
  group_by(gene_symbol) %>%
  summarise(
    WT_avg_expr = rowMeans(select(., starts_with("WT_")), na.rm = TRUE),
    BA52_avg_expr = rowMeans(select(., starts_with("BA52_")), na.rm = TRUE),
    BQ11_avg_expr = rowMeans(select(., starts_with("BQ11_")), na.rm = TRUE)
  )

# 合并 AMD 和 data 的平均表达值
combined_data <- left_join(amd_avg, data_avg, by = "gene_symbol")

# 计算相关性
correlation_matrix <- cor(combined_data[, -1], use = "pairwise.complete.obs") # 去掉基因名列

# 绘制相关性热图
corrplot(correlation_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)
