library(readxl)
library(Seurat)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(wordcloud2)
library(dplyr)
library(jiebaR)
library(tm)
library(htmlwidgets)
library(webshot)
library(patchwork)
set.seed(111)

data1 <- read_excel("C57_CLP_DF.xls", sheet = 1)
data2 <- read_excel("C57_CLPLPS_DF.xls", sheet = 1)
data3 <- read_excel("C57_C57LPS_DF.xls", sheet = 1)
data4 <- read_excel("C57_CLP2DG_DF.xls", sheet = 1)
data5 <- read_excel("C57_CLP2DGLPS_DF.xls", sheet = 1)

data1$sample <- "C57_CLP"
data2$sample <- "C57_CLPLPS"
data3$sample <- "C57_C57LPS"
data4$sample <- "C57_CLP2DG"
data5$sample <- "C57_CLP2DGLPS"

data <- rbind(data1[,c(1:3,6,9)], data2[,c(1:3,6,9)], data3[,c(1:3,6,9)], data4[,c(1:3,6,9)], data5[,c(1:3,6,9)])
data_wide <- dcast(data, `Diseases or Functions Annotation`~data$sample,value.var = 'Activation z-score')
data_wide[is.na(data_wide)] <- 0
rownames(data_wide) <- data_wide[,1]

##################################################################################################
scRNA <- CreateSeuratObject(counts = t(data_wide[,-1]))
anno <- unique(data[match(colnames(scRNA), data$`Diseases or Functions Annotation`),1:3])
rownames(anno) <- anno$`Diseases or Functions Annotation`
scRNA <- AddMetaData(scRNA,metadata = anno)
data_matrix <- GetAssayData(scRNA, assay = "RNA", slot = "data")
cells_data <- t(as.matrix(data_matrix))
duplicate_cells <- duplicated(cells_data)
scRNA_filtered <- scRNA[, !duplicate_cells]
scRNA_filtered <- RunTSNE(scRNA_filtered, features = rownames(scRNA_filtered), dims = NULL)
DimPlot(scRNA_filtered, reduction = "tsne")

genes <- rownames(scRNA)
plots <- lapply(genes, function(gene) {
  FeaturePlot(scRNA_filtered, features = gene,pt.size = 3, reduction = "tsne") +
    scale_color_gradient2(low = '#193E8F', high = '#E53528', mid = "#ebebeb")
})
pdf("FeaturePlot_sample.pdf",width=16,height=10)
wrap_plots(plots, ncol = 3)
dev.off()

###########################################################################################

p <- FeaturePlot(scRNA_filtered, features = rownames(scRNA)[1],pt.size = .5) +
  scale_color_gradient2(low = '#193E8F', high = '#E53528', mid = "#ebebeb")
select.cells <- CellSelector(p)
cluster1 <- select.cells

p <- FeaturePlot(scRNA_filtered, features = rownames(scRNA)[1],pt.size = .5) +
  scale_color_gradient2(low = '#193E8F', high = '#E53528', mid = "#ebebeb")
select.cells <- CellSelector(p)
cluster1 <- c(cluster1,select.cells)
write.csv(anno[cluster1,],"cluster1_DF.csv")

engine <- worker()
segment <- segment(anno[cluster1,]$Functions, engine)
wordfreqs <- freq(segment)
# wordfreqs <- as.data.frame(table(anno[select.cells,2]))
wordf <- arrange(wordfreqs, -freq)
write.csv(wordf,"cluster1_wordf.csv")
wordf <- read_xlsx("cluster1_wordf.xlsx",sheet = 2)
wordcloud <- wordcloud2(data=wordf,shuffle = F,rotateRatio = 0)
saveWidget(wordcloud,"cluster1_wordcloud.html")

##########################################################################################

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.1 <- CellSelector(p)
cluster2 <- select.cells.1

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.2 <- CellSelector(p)
cluster2 <- c(cluster2,select.cells.2)
write.csv(anno[cluster2,],"cluster2_DF.csv")

engine <- worker()
segment <- segment(anno[cluster2,]$Functions, engine)
wordfreqs <- freq(segment)
wordf <- arrange(wordfreqs, -freq)
write.csv(wordf,"cluster2_wordf.csv")
wordf <- read_xlsx("cluster2_wordf.xlsx",sheet = 2)
wordcloud <- wordcloud2(data=wordf,shuffle = F,rotateRatio = 0)
saveWidget(wordcloud,"cluster2_wordcloud.html")

###########################################################################################

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.1 <- CellSelector(p)
cluster3 <- select.cells.1

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.2 <- CellSelector(p)
cluster3 <- c(cluster3,select.cells.2)
write.csv(anno[cluster3,],"cluster3_DF.csv")

###########################################################################################

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.1 <- CellSelector(p)
cluster4 <- select.cells.1

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.2 <- CellSelector(p)
cluster4 <- unique(c(cluster4,select.cells.2))
write.csv(anno[cluster4,],"cluster4_DF.csv")

c4_df <- read.csv("cluster4_DF.csv")
engine <- worker()
segment <- segment(c4_df$Functions, engine)
wordfreqs <- freq(segment)
wordf <- arrange(wordfreqs, -freq) 
write.csv(wordf,"cluster4_wordf.csv")
wordf <- read_xlsx("cluster4_wordf.xlsx",sheet = 1)
wordcloud <- wordcloud2(data=wordf,shuffle = F,rotateRatio = 0)
saveWidget(wordcloud,"cluster4_wordcloud.html")

###########################################################################################

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.1 <- CellSelector(p)
cluster5 <- select.cells.1

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.2 <- CellSelector(p)
cluster5 <- unique(c(cluster5,select.cells.2))

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.3 <- CellSelector(p)
cluster5 <- unique(c(cluster5,select.cells.3))

write.csv(anno[cluster5,],"cluster5_DF.csv")

###########################################################################################

p <- DimPlot(scRNA_filtered, reduction = "tsne")
select.cells.1 <- CellSelector(p)
cluster6 <- select.cells.1

write.csv(anno[cluster6,],"cluster6_DF.csv")

###########################################################################################

df1 <- read.csv("cluster1_DF.csv")
df2 <- read.csv("cluster2_DF.csv")
df3 <- read.csv("cluster3_DF.csv")
df4 <- read.csv("cluster4_DF.csv")
df5 <- read.csv("cluster5_DF.csv")
df6 <- read.csv("cluster6_DF.csv")

df1$function_cluster <- "1"
df2$function_cluster <- "2"
df3$function_cluster <- "3"
df4$function_cluster <- "4"
df5$function_cluster <- "5"
df6$function_cluster <- "6"

df <- rbind(df1,df2,df3,df4,df5,df6)
rownames(df) <- df$Diseases.or.Functions.Annotation
cluster_anno <- df[scRNA_filtered$Diseases.or.Functions.Annotation,]
scRNA_filtered <- AddMetaData(scRNA_filtered,cluster_anno$function_cluster,col.name = "function_cluster")

pdf("DF_dimplot.pdf",width = 5.5,height = 5)
DimPlot(scRNA_filtered, reduction = "tsne",group.by = "function_cluster",pt.size = 3,label = T,label.size = 6,
        cols = c("#EA8379","#299D8F","#FFB77F","#CB95BB","#EBCC96","#AED0DF"))+
  labs(x = "tSNE1", y = "tSNE2",title = NULL) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size=4)))

dev.off()



                                   