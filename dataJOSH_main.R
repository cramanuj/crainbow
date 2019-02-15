library(Seurat); library(plyr); library(dplyr); library(GSVA); library(gplots); library(RColorBrewer); library(car)

source("~/Research/scripts/r_scripts/useful_functions.R")
source("~/Research/scripts/r_scripts/plotfns.R")
ach = read.gmt.file("~/Research/pathways/mouse_pathways/acharya_genesets_mouse.gmt")

run_analysis = function(dat,str){
    pdf(paste(str,"_plots.pdf",sep=""),height=10,width=12)
    seu = CreateSeuratObject(raw.data=dat);
    VlnPlot(object=seu,features.plot=c("nGene", "nUMI"),nCol=2)
    GenePlot(object = seu, gene1 = "nUMI", gene2 = "nGene")
    seu = FilterCells(object = seu, subset.names = c("nGene", "nUMI"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 9000))
    print(dim(seu@data))
    seu = NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
    seu = FindVariableGenes(seu, do.plot = T, mean.function = ExpMean, dispersion.function = LogVMR)
    print(length(seu@var.genes))
    seu = ScaleData(object = seu,vars.to.regress = c("nUMI"))
    seu = RunPCA(seu, pcs.print = 0, pc.genes = seu@var.genes)
    PCAPlot(object = seu, dim.1 = 1, dim.2 = 2)
    set.seed(123)
    seu = RunTSNE(seu, dims.use = 1:10,do.fast=T)
    seu = FindClusters(object = seu, reduction.type = "pca", resolution = 0.8, dims.use = 1:10, save.SNN = T,print.output=F)
    TSNEPlot(seu_new, do.label = TRUE, pt.size = 2.5)
    PCHeatmap(object = seu, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE,cexRow=0.5)
    FeaturePlot(seu, features.plot = c("Krt5","Krt14","Krt6a"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
    FeaturePlot(seu, features.plot = c("Krt8","Krt18","Krt19"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
    markers <- FindAllMarkers(object = seu,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
    top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
    print(DoHeatmap(object = seu, genes.use = top10$gene, slim.col.label = TRUE, remove.key = T,cex.row = 5,group.cex = 8,col.low="blue",col.mid="white",col.high="red"))

    FM = gsva(as.matrix(seu@data),gset.idx.list=ach$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
    seu_ssgsea = CreateSeuratObject(raw.data=FM)
    seu_ssgsea@meta.data$Batch = seu_ssgsea@meta.data$orig.ident
    seu_ssgsea = ScaleData(object = seu_ssgsea)
    seu_ssgsea@ident = seu@ident
    print(DoHeatmap(object = seu_ssgsea, genes.use = rownames(FM), slim.col.label = TRUE, remove.key = T,cex.row = 5,group.cex = 8,col.low="blue",col.mid="white",col.high="red"))
    # ColSideColors = genespring.colors(nlevels(FM_new1$Cluster))[FM_new1$Cluster]
    # RowSideColors = c(rep("gray", 14),rep("green", 15),rep("black", 9),rep("orange",5))
    #
    # FM_new1 = t(FM_new1[,1:nrow(FM)])
    # rownames(FM_new1)[24] = rownames(FM)[24]
    # FM_new1 = FM_new1[match(rownames(FM),rownames(FM_new1)),]
    # col_breaks = c( seq(-1,0,length=50),               # for red
    #                 seq(0.01,0.8,length=50),           # for yellow
    #                 seq(0.81,1,length=50))             # for blue
    # my_palette <- rev(colorRampPalette(c("red", "yellow", "blue"))(n = 149))
    # heatmap.2(as.matrix(FM_new1),main="Gene-set enrichment within the clusters",Colv=NULL,col=my_palette,breaks=col_breaks,scale="none",key=TRUE,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.8,labCol="",dendrogram=c("none"),trace="none",margins=c(9,15),ColSideColors=ColSideColors,RowSideColors=RowSideColors)
    # legend("bottomleft",legend = c("Immune signaling", "Oncogenic signaling", "Tumor microenvironment","DNA repair"),fill = c("gray", "green", "black","orange"),bty="o",cex=0.6)
    for(i in 0:(nlevels(seu_ssgsea@ident)-1)){
      print(i)
      seu_ssgsea_sub = SubsetData(seu_ssgsea,ident.use=i)
      heatmap.2(as.matrix(seu_ssgsea_sub@scale.data),main=paste("Gene-set enrichment within the cluster ", i," (",str,")",sep=""),col=my_palette,notecol="black",breaks=col_breaks,scale="none",key=F,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.8,cexCol=0.8,srtCol=90,dendrogram=c("none"),Rowv=NULL,Colv=NULL,trace="none",margins=c(9,15))
    }
    dev.off()

    FM_new = data.frame(t(seu_ssgsea@scale.data),"Cluster"=as.factor(as.numeric(seu@ident)))
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
    return(list("ssgsea"=FM,"pval"=p_out1,"data"=seu@data,"genes"=seu@var.genes))
}

TU2 = read.delim("TU2.tsv",header=T,row.names=1)
dim(TU2)
TU2_out = run_analysis(TU2,"TU2")

TU3 = read.delim("TU3.tsv",header=T,row.names=1)
dim(TU3)
TU3_out = run_analysis(TU3,"TU3")

MT1 = read.delim("MT1.tsv",header=T,row.names=1)
dim(MT1)
MT1_out = run_analysis(MT1,"MT1")

MT2 = read.delim("MT2.tsv",header=T,row.names=1)
dim(MT2)
MT2_out = run_analysis(MT2,"MT2")

MT3 = read.delim("MT3.tsv",header=T,row.names=1)
dim(MT3)
MT3_out = run_analysis(MT3,"MT3")

MT4 = read.delim("MT4.tsv",header=T,row.names=1)
dim(MT4)
MT4_out = run_analysis(MT4,"MT4")

d16 = read.delim("d16.tsv",header=T,row.names=1)
dim(d16)
d16_out = run_analysis(d16,"d16")

####
seu_ssgsea_TU2 = CreateSeuratObject(raw.data=TU2_out$ssgsea)
seu_ssgsea_TU3 = CreateSeuratObject(raw.data=TU3_out$ssgsea)
seu_ssgsea_MT1 = CreateSeuratObject(raw.data=MT1_out$ssgsea)
seu_ssgsea_MT2 = CreateSeuratObject(raw.data=MT2_out$ssgsea)
seu_ssgsea_MT3 = CreateSeuratObject(raw.data=MT3_out$ssgsea)
seu_ssgsea_MT4 = CreateSeuratObject(raw.data=MT4_out$ssgsea)
seu_ssgsea_d16 = CreateSeuratObject(raw.data=d16_out$ssgsea)

seu_ssgsea@meta.data$Batch = seu_ssgsea@meta.data$orig.ident
seu_ssgsea = ScaleData(object = seu_ssgsea)
