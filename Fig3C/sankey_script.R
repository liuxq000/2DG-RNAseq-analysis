library(ggplot2)
library(reshape2)
library(ggalluvial)

dg <- read.delim("./mouse.txt", header=T)

data1 <- read.csv("CLP_vs_CLPLPS_res1.csv", header=T)
data2 <- read.csv("CLPLPS_vs_CLP2DGLPS_res1.csv", header=T)

data1m <- merge(data1, dg, by.x="X", by.y="ENSEMBLE", all.x=T)
data1m <- data1m[order(data1m$padj),]
data1m$log10p <- -log(data1m$padj,10)

data2m <- merge(data2, dg, by.x="X", by.y="ENSEMBLE", all.x=T)
data2m <- data2m[order(data2m$padj),]
data2m$log10p <- -log(data2m$padj,10)

df1 <- subset(data1m, padj<0.05 & log2FoldChange < -1)
df2 <- subset(data1m, padj<0.05 & log2FoldChange > 1)
df3 <- subset(data2m, padj<0.05 & log2FoldChange < -1)
df4 <- subset(data2m, padj<0.05 & log2FoldChange > 1)

df1$flag <- "CLP_vs_CLPLPS_dn"
df2$flag <- "CLP_vs_CLPLPS_up"
df3$flag <- "CLPLPS_vs_CLP2DGLPS_dn"
df4$flag <- "CLPLPS_vs_CLP2DGLPS_up"

newdf <- rbind(df1, df2, df3, df4)
newdf <- newdf[,c(1,15)]
newdfd <- dcast(newdf, X~flag)
newdfd[is.na(newdfd)] <- "flat" # 4129 5
# colnames(newdfd) <- c("X","I2DGdn","I2DGup","LPSdn","LPSup")

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="CLP_vs_CLPLPS_dn" & CLP_vs_CLPLPS_up=="flat" & CLPLPS_vs_CLP2DGLPS_dn=="flat" & CLPLPS_vs_CLP2DGLPS_up=="flat")
dim(dx)
ddd <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="dn", CLPLPS_vs_CLP2DGLPS="flat")

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="flat" & CLP_vs_CLPLPS_up=="CLP_vs_CLPLPS_up" & CLPLPS_vs_CLP2DGLPS_dn=="flat" & CLPLPS_vs_CLP2DGLPS_up=="flat")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="up", CLPLPS_vs_CLP2DGLPS="flat")
ddd <- rbind(ddd, dy)

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="flat" & CLP_vs_CLPLPS_up=="flat" & CLPLPS_vs_CLP2DGLPS_dn=="CLPLPS_vs_CLP2DGLPS_dn" & CLPLPS_vs_CLP2DGLPS_up=="flat")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="flat", CLPLPS_vs_CLP2DGLPS="dn")
ddd <- rbind(ddd, dy)

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="flat" & CLP_vs_CLPLPS_up=="flat" & CLPLPS_vs_CLP2DGLPS_dn=="flat" & CLPLPS_vs_CLP2DGLPS_up=="CLPLPS_vs_CLP2DGLPS_up")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="flat", CLPLPS_vs_CLP2DGLPS="up")
ddd <- rbind(ddd, dy)

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="CLP_vs_CLPLPS_dn" & CLP_vs_CLPLPS_up=="flat" & CLPLPS_vs_CLP2DGLPS_dn=="CLPLPS_vs_CLP2DGLPS_dn" & CLPLPS_vs_CLP2DGLPS_up=="flat")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="dn", CLPLPS_vs_CLP2DGLPS="dn")
ddd <- rbind(ddd, dy)

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="CLP_vs_CLPLPS_dn" & CLP_vs_CLPLPS_up=="flat" & CLPLPS_vs_CLP2DGLPS_dn=="flat" & CLPLPS_vs_CLP2DGLPS_up=="CLPLPS_vs_CLP2DGLPS_up")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="dn", CLPLPS_vs_CLP2DGLPS="up")
ddd <- rbind(ddd, dy)

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="flat" & CLP_vs_CLPLPS_up=="CLP_vs_CLPLPS_up" & CLPLPS_vs_CLP2DGLPS_dn=="flat" & CLPLPS_vs_CLP2DGLPS_up=="CLPLPS_vs_CLP2DGLPS_up")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="up", CLPLPS_vs_CLP2DGLPS="up")
ddd <- rbind(ddd, dy)

dx <- subset(newdfd, CLP_vs_CLPLPS_dn=="flat" & CLP_vs_CLPLPS_up=="CLP_vs_CLPLPS_up" & CLPLPS_vs_CLP2DGLPS_dn=="CLPLPS_vs_CLP2DGLPS_dn" & CLPLPS_vs_CLP2DGLPS_up=="flat")
dim(dx)
dy <- data.frame(ensemble=dx$X, CLP_vs_CLPLPS="up", CLPLPS_vs_CLP2DGLPS="dn")
ddd <- rbind(ddd, dy)

dm <- data.frame(table(paste0(ddd$CLP_vs_CLPLPS, "-", ddd$CLPLPS_vs_CLP2DGLPS)))
dmm <- data.frame(CLP_vs_CLPLPS=unlist(lapply(strsplit(as.vector(dm$Var1), "-"), function(a) a[1])),
                  CLPLPS_vs_CLP2DGLPS=unlist(lapply(strsplit(as.vector(dm$Var1), "-"), function(a) a[2])),
                  Freq=dm$Freq)
dmm$color <- paste0(dmm$CLP_vs_CLPLPS, "-", dmm$CLPLPS_vs_CLP2DGLPS)

pdf("CLP_vs_CLPLPS--CLPLPS_vs_CLP2DGLPS--sankey.pdf")
ggplot(dmm, aes(axis1=CLP_vs_CLPLPS, axis2=CLPLPS_vs_CLP2DGLPS, y=Freq)) +
  scale_x_discrete(limits = c("CLP_vs_CLPLPS", "CLPLPS_vs_CLP2DGLPS"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = color)) +
  #scale_fill_manual(values=c("#f7acbc", "#94d6da", "#dec674", "#f2eada")) +
  geom_stratum() + geom_text(stat="stratum",aes(label=after_stat(stratum))) +
  theme_bw()
dev.off()



d1 <- subset(ddd, CLP_vs_CLPLPS=="up" & CLPLPS_vs_CLP2DGLPS=="dn")
d2 <- subset(ddd, CLP_vs_CLPLPS=="dn" & CLPLPS_vs_CLP2DGLPS=="up")
d <- rbind(d1,d2)
positive_effect_geneset <- data2m[data2m$X %in% d$ensemble,]
write.csv(positive_effect_geneset,"positive_effect_geneset.csv")

d1 <- subset(ddd, CLP_vs_CLPLPS=="flat" & CLPLPS_vs_CLP2DGLPS=="dn")
d2 <- subset(ddd, CLP_vs_CLPLPS=="flat" & CLPLPS_vs_CLP2DGLPS=="up")
d <- rbind(d1,d2)
side_effect_geneset <- data2m[data2m$X %in% d$ensemble,]
write.csv(side_effect_geneset,"side_effect_geneset.csv")

d1 <- subset(ddd, CLP_vs_CLPLPS=="dn" & CLPLPS_vs_CLP2DGLPS=="flat")
d2 <- subset(ddd, CLP_vs_CLPLPS=="up" & CLPLPS_vs_CLP2DGLPS=="flat")
d <- rbind(d1,d2)
no_effect_geneset <- data1m[data1m$X %in% d$ensemble,]
write.csv(no_effect_geneset,"no_effect_geneset.csv")

