library(ggplot2)
library(reshape2)
library(pheatmap)

datfpkm <- read.table("./FPKM_coding.txt", header=T, sep="\t")
datfpkm=datfpkm[,c(-1,-2)]
datfpkm[is.na(datfpkm)]= 0
custom_order <- c("C57.1", "C57.2", "C57.3","CLP.1","CLP.2","CLP.3","CLP.2DG.1","CLP.2DG.2","CLP.2DG.3","C57LPS1h.1","C57LPS1h.2","C57LPS1h.3","newCLP.LPS1h.1","newCLP.LPS1h.2","newCLP.LPS1h.3","CLP.2DG.LPS1h.1","CLP.2DG.LPS1h.2","CLP.2DG.LPS1h.3") 
datfpkm=datfpkm[,custom_order]

calculate_means <- function(df) {
  n <- ncol(df)
  result <- data.frame(matrix(nrow = nrow(df), ncol = n / 3))
  for (i in seq(1, n, by = 3)) {
    result[, (i + 2) / 3] <- rowMeans(df[, i:(i + 2)])
  }
  return(result)
}

datfpkm_mean <- calculate_means(datfpkm)
datfpkm_mean_2 <- datfpkm_mean[,c("C57","C57LPS","CLP","CLPLPS","CLP2DG","CLP2DGLPS")]
colnames(datfpkm_mean_2) <- c("Sham","Sham+LPS","CLP","CLP+LPS","CLP+2-DG","CLP+2-DG+LPS")
p=pheatmap(datfpkm_mean_2, scale="row", show_rownames=F, clustering_method="ward.D2",cluster_cols = F,fontsize_col = 15)
save_pheatmap_pdf(p, "./Fig2B.pdf",7,7)