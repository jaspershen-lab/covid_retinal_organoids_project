BiocManager::install('ReactomePA')
library(ReactomePA)

#富集分析
eReac<-enrichPathway(gene=diff_entrez$ENTREZID,
                       organism='human',
                       pvalueCutoff=0.05)

#通路网络图

eReacx<-setReadable(eReac,'org.Hs.eg.db','ENTREZID')


cnetplot<-cnetplot(eReacx,
                      #circular=TRUE,
                      colorEdge=TRUE,
                      categorySize="pvalue",
                      showCategory=5,
                      layout='kk')
cnetplot
