#########
#
# TU1 = read.table("TU1_scimpute.txt",header=T,row.names=1,sep=" ")
# dim(TU1)
# TU4 = read.table("TU4_scimpute.txt",header=T,row.names=1,sep=" ")
# dim(TU4)
#########
library(Seurat); library(plyr); library(dplyr); library(GSVA); library(gplots); library(RColorBrewer); library(car)

source("~/Research/scripts/r_scripts/useful_functions.R")
source("~/Research/scripts/r_scripts/plotfns.R")

TU1 = read.delim("TU1.tsv",header=T,row.names=1)
TU4 = read.delim("TU4.tsv",header=T,row.names=1)

### Batch correction using mutual nearest neighbors method
# library(scater); library(scran)
# TU1_sce = SingleCellExperiment(assays = list(counts = as.matrix(TU1)), colData = data.frame("Cell"=colnames(TU1)))
# TU1_sce = normalizeSCE(TU1_sce)
# TU4_sce = SingleCellExperiment(assays = list(counts = as.matrix(TU4)), colData = data.frame("Cell"=colnames(TU4)))
# TU4_sce = normalizeSCE(TU4_sce)
# TU_new = mnnCorrect(exprs(TU1_sce),exprs(TU4_sce))
#
# TU_comb = do.call("cbind",TU_new$corrected)
# colnames(TU_comb) = c(colnames(TU1),colnames(TU4))
#
# ## Seurat
# library(Seurat)
# pdata = data.frame("Batch"=c(rep("TU1",ncol(TU1)),rep("TU4",ncol(TU4))))
# rownames(pdata) = colnames(TU_comb)
# seu = CreateSeuratObject(raw.data = TU_comb,project = "TU1_TU4",orig.ident=NULL,meta.data=pdata)
# VlnPlot(object=seu,features.plot=c("nGene", "nUMI"),nCol=2)
# GenePlot(object = seu, gene1 = "nUMI", gene2 = "nGene")
###
pdf("TU1_TU4_plots.pdf",width=11,height=7)
seu1 = CreateSeuratObject(raw.data=TU1)
seu2 = CreateSeuratObject(raw.data=TU4)
seu_new = MergeSeurat(seu1,seu2)
seu_new@meta.data$Batch = c(rep("TU1",ncol(TU1)),rep("TU4",ncol(TU4)))
VlnPlot(object=seu_new,features.plot=c("nGene", "nUMI"),nCol=2)
GenePlot(object = seu_new, gene1 = "nUMI", gene2 = "nGene")
quantile(seu_new@meta.data$nGene)
quantile(seu_new@meta.data$nUMI)
seu_new = FilterCells(object = seu_new, subset.names = c("nGene"), low.thresholds = c(100), high.thresholds = c(Inf))
seu_new = NormalizeData(object = seu_new, normalization.method = "LogNormalize", scale.factor = 10000)
seu_new = FindVariableGenes(seu_new, do.plot = T, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.01, x.high.cutoff = 3, y.cutoff = 0.5)
seu_new = ScaleData(object = seu_new,vars.to.regress = c("nUMI", "Batch"))
seu_new = RunPCA(seu_new, pcs.print = 0, pc.genes = seu_new@var.genes)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2)
set.seed(123)
# seu_new = RunUMAP(seu_new, dims.use = 1:10)
# DimPlot(seu_new,reduction.use = 'umap',group.by="Batch",cols.use=c("red","blue"))
seu_new = RunTSNE(seu_new, dims.use = 1:10,do.fast=T)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 2.5,group.by="Batch")

my_TU1 = grep("TU1", row.names(seu_new@meta.data), value = TRUE)
my_TU4 = grep("TU4", row.names(seu_new@meta.data), value = TRUE)
my_TU1_tsne = TSNEPlot(seu_new,cells.use = my_TU1,do.return=TRUE,colors.use = "red")
my_TU4_tsne = TSNEPlot(object = seu_new, cells.use = my_TU4,do.return = TRUE,colors.use = "blue")
plot_grid(my_TU1_tsne,my_TU4_tsne)

seu_new = FindClusters(object = seu_new, reduction.type = "pca", resolution = c(0.4, 0.8, 1.2), dims.use = 1:10, save.SNN = TRUE)
table(seu_new@meta.data$res.0.4,seu_new@meta.data$Batch)
table(seu_new@meta.data$res.0.8,seu_new@meta.data$Batch)
table(seu_new@meta.data$res.1.2,seu_new@meta.data$Batch)
seu_new = SetAllIdent(seu_new, id='res.0.4')
TSNEPlot(seu_new, do.label = TRUE, pt.size = 2.5)
PCHeatmap(object = seu_new, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE,cexRow=0.5)
FeaturePlot(seu_new, features.plot = c("Krt5","Krt14","Krt6a"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("Krt8","Krt18","Krt19"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

markers <- FindAllMarkers(object = seu_new,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seu_new, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 5,group.cex = 5,col.low="blue",col.mid="white",col.high="red")

############
library(GSVA)
# cp = read.gmt.file("~/Research/pathways/mouse_pathways/Mouse_Human_MSigdb_February_01_2019_symbol.gmt")
# hm = read.gmt.file("~/Research/pathways/h.all.v5.2.symbols.gmt")
ach = read.gmt.file("~/Research/pathways/mouse_pathways/acharya_genesets_mouse.gmt")

ssgsea_out = gsva(as.matrix(seu_new@data),gset.idx.list=ach$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
FM = ssgsea_out
rownames(FM)= gsub(" ","_",rownames(FM))

seu_ssgsea = CreateSeuratObject(raw.data=FM)
seu_ssgsea@meta.data$Batch = seu_ssgsea@meta.data$orig.ident
seu_ssgsea = ScaleData(object = seu_ssgsea,vars.to.regress = c("Batch"))
seu_ssgsea@ident = seu_new@ident
DoHeatmap(object = seu_ssgsea, genes.use = rownames(FM), slim.col.label = TRUE, remove.key = TRUE,cex.row = 5,group.cex = 8,col.low="blue",col.mid="white",col.high="red")
col_breaks = c( seq(-1,0,length=100),               # for red
                seq(0.01,0.8,length=100),           # for yellow
                seq(0.81,1,length=100))             # for green
my_palette <- rev(colorRampPalette(c("red", "white", "blue"))(n = 299))

for(i in 0:(nlevels(seu_ssgsea@ident)-1)){
  print(i)
  seu_ssgsea_sub = SubsetData(seu_ssgsea,ident.use=i)
  ColSideColors = c("green","orange")[seu_ssgsea_sub@meta.data$Batch]
  heatmap.2(as.matrix(seu_ssgsea_sub@scale.data),main=paste("Gene-set enrichment within the cluster ", i, "(TU1/TU4)",sep=""),col=my_palette,notecol="black",breaks=col_breaks,scale="none",key=F,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.8,cexCol=0.8,srtCol=90,dendrogram=c("none"),Rowv=NULL,Colv=NULL,trace="none",margins=c(9,15),ColSideColors=ColSideColors)
  legend("bottomleft",legend = c("TU1","TU4"),fill = c("green","orange"),bty="o",cex=0.6)
}

# seu_ssgsea = RunPCA(seu_ssgsea, pcs.print = 0, pc.genes = rownames(FM))
# PCAPlot(object = seu_ssgsea, dim.1 = 1, dim.2 = 2)
# seu_ssgsea = RunUMAP(seu_ssgsea, dims.use = 1:10,n_neighbors = 30L,min_dist = 0.05,seed.use=1234)
# DimPlot(seu_ssgsea,reduction.use = 'umap')
# seu_ssgsea = FindClusters(object = seu_ssgsea, reduction.type = "pca", resolution = c(0.4, 0.8, 1.2), dims.use = 1:10, save.SNN = TRUE)
# table(seu_ssgsea@meta.data$res.0.4,seu_ssgsea@meta.data$Batch)
# table(seu_ssgsea@meta.data$res.0.8,seu_ssgsea@meta.data$Batch)
# table(seu_ssgsea@meta.data$res.1.2,seu_ssgsea@meta.data$Batch)
# seu_ssgsea = SetAllIdent(seu_ssgsea, id='res.0.8')
# DimPlot(seu_ssgsea,reduction.use = 'umap',do.label=T)

# path_markers <- FindAllMarkers(object = seu_ssgsea,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
# top10_path = path_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# DoHeatmap(object = seu_ssgsea, genes.use = top10_path$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 5,group.cex = 5)
# PCHeatmap(object = seu_ssgsea, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE,cexRow=0.5)

# ColSideColors = genespring.colors(nlevels(FM_new1$Cluster))[FM_new1$Cluster]
# RowSideColors = c(rep("gray", 14),rep("green", 15),rep("black", 9),rep("orange",5))
#
# FM_new1 = t(FM_new1[,1:nrow(FM)])
# rownames(FM_new1)[24] = rownames(FM)[24]
# FM_new1 = FM_new1[match(rownames(FM),rownames(FM_new1)),]
# col_breaks = c( seq(-1,0,length=100),               # for red
#                 seq(0.01,0.8,length=100),           # for yellow
#                 seq(0.81,1,length=100))             # for green
# my_palette <- rev(colorRampPalette(c("red", "white", "blue"))(n = 299))
# # heatmap.2(as.matrix(FM_new1),main="",Rowv=NULL,Colv=NULL,col=rev(brewer.pal(11,"RdBu")),scale="row",key=TRUE,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.8,labCol="",dendrogram=c("none"),trace="none",margins=c(8,15),ColSideColors=ColSideColors)
# heatmap.2(as.matrix(FM_new1),main="Gene-set enrichment within the clusters",Rowv=NULL,Colv=NULL,col=my_palette,notecol="black",breaks=col_breaks,scale="none",key=TRUE,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.8,labCol="",dendrogram=c("none"),trace="none",margins=c(9,15),ColSideColors=ColSideColors,RowSideColors=RowSideColors)
# legend("bottomleft",legend = c("Immune signaling", "Oncogenic signaling", "Tumor microenvironment","DNA repair"),fill = c("gray", "green", "black","orange"),bty="o",cex=0.6)
# dev.off()

FM_new = data.frame(t(seu_ssgsea@scale.data),"Cluster"=as.factor(as.numeric(seu_new@ident)))
FM_new1 = FM_new[order(FM_new$Cluster),]
FM_new$Cluster_new = FM_new$Cluster
p_out = list()
for(i in 1:nlevels(FM_new$Cluster)){
  print(i)
  FM_new$Cluster_new = ifelse(FM_new$Cluster==i,1,2)
  p = p.adjust(apply(FM_new[,1:nrow(FM)],2,function(x) Anova(lm(x ~ FM_new$Cluster_new))$"Pr(>F)"[1]))
  p = data.frame("FDR"=p[p<=0.01],"Cluster"=i)
  p = data.frame("Pathways"=rownames(p),p[,1:2],row.names=NULL)
  p = p[order(p$FDR),]
  p_out[[i]] = p
}
p_out1 = do.call("rbind",llply(1:length(p_out),function(i) head(p_out[[i]],10)))
write.table(p_out1,"TU1_TU4_genesets.txt",sep="\t",col.names=NA,quote=F)

#############################################
#

save.image("final_analysis.rds")
