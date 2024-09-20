library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100-tools.R')

library(DESeq2)
library(dplyr)
library(ggplot2)
if(!require(remotes)){
  install.packages("remotes")
}
remotes::install_gitlab("tidymass/tidymass")
library(tidymass)



###data clean
data <- read.csv("2_data/1. RO infected with SARS-CoV-2 strains  RNA-seq data all_counts_with_symbol_new.csv")
# 去除symbol列中的NA值
data <- data[!is.na(data$symbol), ]
library(dplyr)


# 按symbol分组，计算平均值
data_clean <- data %>%
  group_by(symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = 'drop')

head(data_clean)

# 保存为clean_csv文件
setwd("/Volumes/sirius/Study/covid_retinal_organoids_project/3_data_analysis/5-pathway_analysis")
write.csv(data_clean, "clean_csv.csv", row.names = FALSE)


####DEG analysis###
DEG_data <- read.csv("clean_csv.csv", header=TRUE, row.names=1)
head(DEG_data)

DEG_data<- round(DEG_data)  # 四舍五入

# 构建样本信息
sample_info <- data.frame(
  row.names = colnames(DEG_data),
  condition = factor(c(rep("WT", 3), rep("BA.5.2", 3), rep("BQ.1.1", 3), rep("Ctrl", 3)))
)
head(sample_info)

dds <- DESeqDataSetFromMatrix(countData = DEG_data,
                              colData = sample_info,
                              design = ~ condition)
dds

dds <- DESeq(dds)

# 提取结果（以 BA.5.2 相对于 WT 为例）
results_ba52 <- results(dds, contrast = c("condition", "BA.5.2", "WT"))
head(results_ba52)
summary(results_ba52)
write.csv(as.data.frame(results_ba52), file = "DEG_results_BA.5.2_vs_WT.csv")

# 提取结果（以 BQ.1.1 相对于 WT 为例）
results_bq11 <- results(dds, contrast = c("condition", "BQ.1.1", "WT"))
head(results_bq11)
summary(results_bq11)
write.csv(as.data.frame(results_bq11), file = "DEG_results_BQ.1.1_vs_WT.csv")


###Volcano plot
###WT vs BA.5.2
results_ba52_df <- as.data.frame(results_ba52)
head(results_ba52_df)
plot <-
  volcano_plot(
    results_ba52_df,
    fc_column_name = "log2FoldChange",
    p_value_column_name = "padj",
    labs_x = "Log2 Fold Change",
    labs_y = "-Log10 P-value",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.5,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_wt_ba52.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_wt_ba52.pdf",
       width = 7,
       height = 6)


###WT vs BQ.1.1
results_bq11_df <- as.data.frame(results_bq11)
head(results_bq11_df)
plot <-
  volcano_plot(
    results_bq11_df,
    fc_column_name = "log2FoldChange",
    p_value_column_name = "padj",
    labs_x = "Log2 Fold Change",
    labs_y = "-Log10 P-value",
    fc_up_cutoff = 2,
    fc_down_cutoff = 0.5,
    p_value_cutoff = 0.05,
    line_color = "red",
    up_color = "#EE0000FF",
    down_color = "#3B4992FF",
    no_color = "#808180FF",
    point_size = 2,
    point_alpha = 0.5,
    point_size_scale = "log10_p",
    line_type = 1,
    add_text = FALSE
  ) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_x_continuous(limits = c(-10, 7)) +
  scale_y_continuous(limits = c(0, 3))

plot

ggsave(plot,
       filename = "volcano_plot_wt_bq11.png",
       width = 7,
       height = 6)

ggsave(plot,
       filename = "volcano_plot_wt_bq11.pdf",
       width = 7,
       height = 6)



##pathway_analysis
library(devtools)
remotes::install_github("YuLab-SMU/GOSemSim")
install_github("YuLab-SMU/clusterProfiler")

library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)


###BA52
####差异基因提取
results_ba52$group <- case_when(
  results_ba52$log2FoldChange > 2 & results_ba52$pvalue < 0.05 ~ "up",
  results_ba52$log2FoldChange < -2 & results_ba52$pvalue < 0.05 ~ "down",
  abs(results_ba52$log2FoldChange) <= 2 ~ "none",
  results_ba52$pvalue >= 0.05 ~ "none"
)
head(results_ba52)
rownames(results_ba52)
table(results_ba52$group)

#根据研究目的筛选差异基因(仅上调、下调或者全部)：
up <- rownames(results_ba52)[results_ba52$group=="up"]#差异上调
down <- rownames(results_ba52)[results_ba52$group=="down"]#差异下调
diff <- c(up,down)#所有差异基因
head(up)
head(down)

#ID转换：
#查看可转换的ID类型：
columns(org.Hs.eg.db)

##使用clusterProfiler包自带ID转换函数bitr(基于org.Hs.eg.db)：
#up：
up_entrez <- bitr(up,
                  fromType= "SYMBOL",#现有的ID类型
                  toType= "ENTREZID",#需转换的ID类型
                  OrgDb= "org.Hs.eg.db")
#down：
down_entrez<- bitr(down,
                    fromType= "SYMBOL",
                    toType= "ENTREZID",
                    OrgDb= "org.Hs.eg.db")
#diff：
diff_entrez<- bitr(diff,
                    fromType= "SYMBOL",
                    toType= "ENTREZID",
                    OrgDb= "org.Hs.eg.db")
head(diff_entrez)

#GO(Gene Ontology)富集分析：
##BP(我们以总差异基因的GO富集为例):
GO_BP_diff<- enrichGO(gene = diff_entrez$ENTREZID, #用来富集的差异基因
                       OrgDb= org.Hs.eg.db, #指定包含该物种注释信息的org包
                       ont= "BP", #可以三选一分别富集,或者"ALL"合并
                       pAdjustMethod= "BH", #多重假设检验矫正方法
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.05,
                       readable= TRUE) #是否将gene ID映射到gene name
#提取结果表格：
GO_BP_result<- GO_BP_diff@result
View(GO_BP_result)

#CC:
GO_CC_diff<- enrichGO(gene = diff_entrez$ENTREZID,
                       OrgDb= org.Hs.eg.db,
                       ont= "CC",
                       pAdjustMethod= "BH",
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.05,
                       readable= TRUE)
#提取结果表格：
GO_CC_result<- GO_CC_diff@result

#MF:
GO_MF_diff<- enrichGO(gene = diff_entrez$ENTREZID,
                       OrgDb= org.Hs.eg.db,
                       ont= "BP",
                       pAdjustMethod= "BH",
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.05,
                       readable= TRUE)
#提取结果表格：
GO_MF_result<- GO_BP_diff@result

#MF、CC、BP三合一：
GO_all_diff<- enrichGO(gene = diff_entrez$ENTREZID,
                        OrgDb= org.Hs.eg.db,
                        ont= "ALL", #三合一选择“ALL”
                        pAdjustMethod= "BH",
                        pvalueCutoff= 0.05,
                        qvalueCutoff= 0.05,
                        readable= TRUE)
#提取结果表格：
GO_all_result<- GO_all_diff@result

##保存GO富集结果：
save(GO_MF_diff,GO_CC_diff,GO_BP_diff,GO_all_diff,file = c("GO_diff.Rdata"))

#####快速可视化探索(以CC为例)：
BiocManager::install("topGO")
install.packages("SparseM")
library(topGO)
library(enrichplot)

#GO有向无环图绘制(两种形式)：
goplot(GO_CC_diff) #来自enrichplot
plotGOgraph(GO_CC_diff) #来自topGO

#GO富集气泡图：
dotplot(
  GO_BP_diff,
  x= "GeneRatio",
  color= "p.adjust",
  title= "Top 20 of GO CC terms Enrichment",
  showCategory= 20,
  label_format= 30
)


#富集网络图:
edo <- pairwise_termsim(GO_CC_diff)
emapplot(edo,
         layout = "kk", #布局形式
         showCategory = 30) #展示GO terms的数量


####使用ggplot2进行可视化:
#取前top20，并简化命名：
MF<- GO_MF_result[1:20,]
CC<- GO_CC_result[1:20,]
BP<- GO_BP_result[1:20,]

#在MF的Description中存在过长字符串，我们将长度超过50的部分用...代替：
MF2<- MF
MF2$Description <- str_trunc(MF$Description,width = 50,side = "right")
MF2$Description

CC2<- CC
CC2$Description <- str_trunc(CC$Description,width = 50,side = "right")
CC2$Description

BP2<- BP
BP2$Description <- str_trunc(BP$Description,width = 50,side = "right")
BP2$Description

#自定义主题
mytheme<- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust= 0.5,
                                           face= "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

# 创建唯一的描述
MF2$unique_description <- paste(MF2$Description, seq_along(MF2$Description), sep = "_")
BP2$unique_description <- paste(BP2$Description, seq_along(BP2$Description), sep = "_")
# 设置因子，确保水平是唯一的
MF2$term <- factor(MF2$unique_description, levels = rev(MF2$unique_description))

CC2$term <- factor(CC2$Description,levels = rev(CC2$Description))
BP2$term <- factor(BP2$unique_description, levels = rev(BP2$unique_description))



#GO富集柱形图：
GO_bar<- function(x){
  y<- get(x)
  ggplot(data = y,
         aes(x = Count,
             y= term,
             fill= -log10(pvalue))) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 50) ) + #label换行，部分term描述太长
    geom_bar(stat = "identity",width = 0.8) +
    labs(x = "Gene Number",
         y= "Description",
         title= paste0(x," of GO enrichment barplot")) +
    theme_bw() +
    mytheme
}

#MF:
p1<- GO_bar("MF2")+scale_fill_distiller(palette = "Blues",direction = 1)
p1
#BP：
p2<- GO_bar("CC2")+scale_fill_distiller(palette = "Reds",direction = 1)
p2
#BP:
p3 <- GO_bar("BP2")+scale_fill_distiller(palette = "Oranges",direction = 1)
p3


#将三个ontology拉通取top30(按照p值排序)绘图：
all_result<- arrange(GO_all_result,pvalue) #默认升序

#取top30：
all<- all_result[1:30,]

#指定绘图顺序（转换为因子）：
all$term <- factor(all$Description,levels = rev(all$Description))
#自定义y轴标签颜色(区分不同ontology)：
col_function<- function(x){
  col<- rep("black", length(x))
  BP<- which(x %in% c("BP"))
  CC<- which(x %in% c("CC"))
  MF<- which(x %in% c("MF"))
  col[BP] <- "#fc4d26"
  col[CC] <- "#1792c1"
  col[MF] <- "#3fad5d"
  col
}

y_text_color<- col_function(all$ONTOLOGY)

#绘制富集柱形图(代码同上，这里不再赘述)：
ppp1 <- GO_bar("all")+scale_fill_distiller(palette = "Blues",direction = 1)
ppp1
