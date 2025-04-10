library(factoextra)
library(ggforce)

fpkm <- read.table("FPKM_coding_filt.txt",sep ="\t",header = T)
pca_data <- fpkm[,-c(1,2,6:8,12:14)]
colnames(pca_data) <- c("CLPLPS-1","CLPLPS-2","CLPLPS-3","CLP-1","CLP-2","CLP-3","CLP2DGLPS-1","CLP2DGLPS-2","CLP2DGLPS-3","CLP2DG-1","CLP2DG-2","CLP2DG-3")
pca <- t(pca_data)
pca_result <- prcomp(pca[ , which(apply(pca, 2, var) != 0)], scale. = T)
fviz_screeplot(pca_result, addlabels = TRUE, ylim = c(0, 40))
group = c(rep("CLPLPS",3),rep("CLP",3),rep("CLP2DGLPS",3),rep("CLP2DG",3))
names(group) <- colnames(pca_data)
group <- factor(group, levels = c("CLP","CLP2DG","CLPLPS","CLP2DGLPS"))

fviz_pca_ind(pca_result, 
             mean.point=F,
             geom.ind = c('point','text'),
             habillage = group,
             pointsize = 4,
             #addEllipses = T,
             #ellipse.level = c(0.5),
             palette = c("#55B7E6","#B395BD","#EA8379","#F09739"),
             repel = T,
             pointshape = 16) +
  theme_bw((base_size=14))+
  theme(text = element_text(size = 12),
        legend.margin = margin(-10),
        axis.text = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))+
  ggtitle('')+
ggforce::geom_mark_ellipse(aes(fill = group,color = group))+coord_cartesian(clip = "off")
ggsave("Fig3A-pca.pdf",width = 8,height=7)
