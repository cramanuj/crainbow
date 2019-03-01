
rm(list=ls())

####################################################################
## Code to run multivariate CCA on various mouse tumors (scRNA-seq)
##
## Author: Chaitanya Acharya
## Updated: Feb 27, 2019
####################################################################
## Set your working directory
setwd("/users/jcs48/Repository/Her2/") ## setwd("your directory")

## Install BiocManager if not yet installed
if(!requireNamespace("BiocManager")) install.packages("BiocManager")

## Load the required R libraries
## If the libraries are not installed, use BiocManager
lib.list=c("Seurat","plyr","ggpubr","dplyr","ggthemes","data.table","GSVA","biomaRt")

for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			BiocManager::install(lib.list[i])
			library(lib.list[i],character.only=T)
    };
}

## R scripts with utility functions
source("/users/jcs48/Repository/Her2/crainbow/useful_functions.R")
source("/users/jcs48/Repository/Her2/crainbow/plotfns.R")

# Assuming that these R scripts are in the current working dir.
# source("useful_functions.R")
# source("plotfns.R")

### Read the gene-set file
ach = read.gmt.file("/users/jcs48/Repository/Her2/crainbow/acharya_genesets_mouse.gmt")
pam50 = scan("/users/jcs48/Repository/Her2/crainbow/pam50_genes.txt",what="")
pam50_mm = convertHumanGeneList(pam50)
ach$geneset.descriptions[51]="PAM50_mouse_homologs"
ach$geneset.names[51]="PAM50"
ach$genesets[[51]]=pam50_mm
names(ach$genesets)[51]="PAM50"

## Cell cycle genes
cc_genes = scan("/users/jcs48/regev_lab_cell_cycle_genes.txt",what="")

# We can segregate this list into markers of G2/M phase and markers of S phase and convert them from human to mouse markers
s_genes = convertHumanGeneList(cc_genes[1:43])
g2m_genes = convertHumanGeneList(cc_genes[44:97])

### Read all the possible .tsv files and identity select files
### Create a list object containing gene expression matrices of the select files
files = list.files(pattern=".tsv")
files = files[files %in% c("MT1.tsv","TU1.tsv","TU2.tsv","TU3.tsv","TU4.tsv")]
data_list = llply(1:length(files),.progress="time",function(i) read.delim(files[i],header=T,row.names=1))

### Create a SEURAT object from each gene expression matrix within the list
### Each SEURAT object is preprocessed separately and highly variable genes are computed
seu_list = gene_list = vector(mode="list",length=length(data_list))
for(i in 1:length(seu_list)){
  print(i)
  seu = CreateSeuratObject(raw.data=data_list[[i]],min.cells=5)
  seu@meta.data$group = gsub(".tsv","",files[i])
  seu = FilterCells(seu, subset.names = "nGene", low.thresholds = round(quantile(seu@meta.data$nGene,0.1)), high.thresholds = Inf)
  seu = NormalizeData(seu,display.progress = F)
  seu = CellCycleScoring(seu, s.genes = s_genes, g2m.genes = g2m_genes, set.ident = F)
  seu@meta.data$CC.Difference = seu@meta.data$S.Score - seu@meta.data$G2M.Score
  seu = ScaleData(seu, display.progress = F,vars.to.regress = c("nUMI"))
  seu = FindVariableGenes(seu,display.progress = F, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR)
  genes_use = head(rownames(seu@hvg.info), 2000)
  seu_list[[i]] = seu
  gene_list[[i]] = genes_use
}

#######################################################################################
## Gene selection for CCA
## We find genes that are highly variable in at least two datasets
#######################################################################################

genes_use <- names(which(table(unlist(gene_list)) > 1))
for (i in 1:length(seu_list)) {
  genes_use <- genes_use[genes_use %in% rownames(seu_list[[i]]@scale.data)]
}

#######################################################################################
### Run multivariate Canonical Correlation Analysis (CCA) on the common genes across all the genesets
### Calculate the ratio of total variance explained by PPCA vs total variance explained by CCA,
###       and filter cells based on these values
### CCA Align the matrices and compute the alignment metric score
### Map a t-SNE plot and find clusters
#######################################################################################

pdf("CCA_plots.pdf",width=10,height=7)
cca_out = RunMultiCCA(seu_list,genes.use=genes_use,num.ccs = 20)
DimPlot(object = cca_out, reduction.use = "cca", group.by = "group", pt.size = 1, do.return = F)
VlnPlot(object = cca_out, features.plot = "CC1", group.by = "group", do.return = F)
MetageneBicorPlot(cca_out, grouping.var = "group", dims.eval = 1:20, display.progress = FALSE)
DimHeatmap(object = cca_out, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
cca_out = CalcVarExpRatio(cca_out,reduction.type = "pca", grouping.var = "group", dims.use = 1:20)
cca_out = SubsetData(cca_out, subset.name = "var.ratio.pca",accept.low = 0.40)
cca_out = AlignSubspace(cca_out, reduction.type = "cca", grouping.var = "group", dims.align = 1:10,num.genes=50)
metric = rep(0,10)
for(i in 1:10){
  metric[i]=CalcAlignmentMetric(cca_out,reduction.use = "cca.aligned",dims.use = 1:i, grouping.var =  "group")
  cat("Number of PCs: ",i," Alignment metric: ",metric[i],"\n")
}
dim_max = metric[3]
dev.off()

pdf("TSNE_plots.pdf",width=10,height=7)
cca_out = RunTSNE(cca_out, reduction.use = "cca.aligned", dims.use = 1:3, do.fast = T, check_duplicates=F)
TSNEPlot(object = cca_out, do.label = F,group.by="group")
TSNEPlot(object = cca_out, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5,group.by="Phase")
cca_out = FindClusters(cca_out, reduction.type = "cca.aligned", resolution = c(0.6), dims.use = 1:3,print.output = F)
TSNEPlot(object = cca_out, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5)
FeaturePlot(cca_out, features.plot = c("Krt5","Krt14","Krt6a"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Krt8","Krt18","Krt19"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Il6","Il12b","Serpinb2","Cxcl3"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Arg1","Mrc1","Egr2","Cd83"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out, features.plot = c("Ptprc","Foxp3","Cd3d","Cd4","Cd14"), reduction.use="tsne",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)

### Compute single sample GSEA (ssGSEA) scores
library(GSVA); library(ggpubr); library(dplyr); library(ggthemes); library(data.table);
ssgsea_out = gsva(as.matrix(cca_out@data),gset.idx.list=ach$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
FM = ssgsea_out
FM = as.matrix(Matrix::t(scale(Matrix::t(FM))))
FM = melt(data.frame(t(FM),"Cluster"=cca_out@ident,"Group"=cca_out@meta.data$group))
colnames(FM)[3:4]=c("GeneSets","EnrichmentScore")

### Heatmap of cluster-wise gene-set enrichment scores for every gene-set
ES_mean = FM %>% dplyr::group_by(Cluster,GeneSets,Group) %>% summarize_all(list(mean))
ggplot(ES_mean, aes(Cluster, GeneSets)) + geom_tile(aes(fill = EnrichmentScore), colour = "white") + scale_fill_gradient2(low ="dark blue", high ="dark red", mid ="white",space = "Lab") +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + facet_grid(~Group)

### Heatmap of cluster-specific gene-set enrichment scores
ES1_mean = FM[,-1] %>% dplyr::group_by(GeneSets,Group) %>% summarize_all(list(mean))
ggplot(ES1_mean, aes(Group, GeneSets)) + geom_tile(aes(fill = EnrichmentScore), colour = "white") + scale_fill_gradient2(low ="dark blue", high ="dark red", mid ="white",space = "Lab") +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90))

### Gene-set feature plots
cca_out@meta.data = cbind(cca_out@meta.data,t(as.matrix(Matrix::t(scale(Matrix::t(ssgsea_out))))))
FeaturePlot(cca_out,features.plot=c("BASAL","LUMINAL","MASC","STROMAL","LUMINAL PROGENITORS","EMT"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("MACROPHAGE_ACTIVITY","Treg","INFLAMMATORY RESPONSE","NK_CELLS","TCELL_ACTIVATION","G2M CHECKPOINT"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("HYPOXIA","INVASIVENESS GENE SIGNATURE","ANGIOGENESIS","ECM","FIBROBLAST GROWTH FACTOR RECEPTOR","EPIGENETIC STEM CELL"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("WNT BETA CATENIN SIGNALING","TGF BETA SIGNALING","E2F TARGETS","P53 PATHWAY","HEDGEHOG SIGNALING","PI3K-AKT-MTOR SIGNALING"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("BASE EXCISION REPAIR","MISMATCH EXCISION REPAIR","NUCLEOTIDE EXCISION REPAIR","HOMOLOGOUS RECOMBINATION","NONHOMOLOGOUS END JOINING"),reduction.use="tsne",nCol=3,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
FeaturePlot(cca_out,features.plot=c("PAM50"),reduction.use="tsne",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "red"), pt.size = 1)
###

dev.off()
