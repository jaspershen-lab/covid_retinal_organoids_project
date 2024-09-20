#KEGG富集分析(超几何分布检验)：
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa", #物种Homo sapiens (智人)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize= 500)

#将ENTREZ重转为symbol：
KEGG_diff <- setReadable(KEGG_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")
View(KEGG_diff@result)

colnames(KEGG_diff@result) #解释不同列名含义

parse_ratio <- function(ratio) {
  if (is.character(ratio)) {
    parts <- strsplit(ratio, "/")[[1]]
    return(as.numeric(parts[1]) / as.numeric(parts[2]))
  } else {
    return(as.numeric(ratio))
  }
}


#计算Rich Factor（富集因子）：
KEGG_diff2<- mutate(KEGG_diff,
                     RichFactor= Count / as.numeric(sub("/\\d+", "", BgRatio)))

#计算Fold Enrichment（富集倍数）：
KEGG_diff2<- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]


#提取KEGG富集分析结果表：
KEGG_result <- KEGG_diff2@result

#保存富集结果到本地：
setwd("/Volumes/sirius/Study/covid_retinal_organoids_project/3_data_analysis/5-pathway_analysis")
save(KEGG_diff2, KEGG_result, file = c("KEGG_diff.Rdata"))
write.csv(KEGG_result, file = c('KEGG_diff_BA.5.2.csv'))


#富集条形图绘制：
library(enrichplot)
?enrichplot::barplot.enrichResult #查看barplot函数说明
barplot(
  KEGG_diff2, 
  x= "Count", #or "GeneRatio"
  color= "pvalue", #or "p.adjust", "qvalue"
  showCategory= 20, #显示pathway的数量
  font.size = 12, #字号
  title= "KEGG enrichment barplot", #标题
  label_format= 30 #pathway标签长度超过30个字符串换行
)

#以富集结果表Top20为例：
KEGG_top20<- KEGG_result[1:20,]

#指定绘图顺序（转换为因子）：
KEGG_top20$pathway <- factor(KEGG_top20$Description, levels = rev(KEGG_top20$Description))

#Top20富集数目条形图：
mytheme<- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 7), 
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
                 legend.title = element_text(size = 13), 
                 legend.text = element_text(size = 11)) #自定义主题

#将pathway顺序按照富集的gene number数排列：
p3<- ggplot(data = KEGG_top20, 
             aes(x = Count, y = reorder(pathway, Count))) + 
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_color_distiller(palette = "Spectral",direction = 1) +
  labs(x = "Gene Number", 
       y= "",
       title= "Dotplot of Enriched KEGG Pathways",
       size= "gene number") +
  theme_bw() + 
  mytheme
p3

ggsave(p3,
       filename = "KEGG_wt_ba52.png",
       width = 7,
       height = 6)

ggsave(p3,
       filename = "KEGG_wt_ba52.pdf",
       width = 7,
       height = 6)

#富集倍数版显著性气泡图
p4<- ggplot(data = KEGG_top20, 
             aes(x = FoldEnrichment, y = pathway)) + 
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_color_distiller(palette = "Spectral", direction = -1) +
  labs(x = "Fold Enrichment", 
       y= "",
       title= "Dotplot of Enriched KEGG Pathways",
       size= "Gene Number") +
  theme_bw() +
  mytheme
p4
