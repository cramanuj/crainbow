
## 1) Do PCA
doPCA = function(dat,scale=T,center=T,numPC=10){
	data_pca = prcomp(t(dat),center=center,scale=scale)
	pca = as.data.frame((t(dat) %*% data_pca$rotation))
	pca_var = data.frame("Variance"=round(data_pca$sdev^2/sum(data_pca$sdev^2) * 100 ,2)[1:numPC],"PC"=1:numPC);
	pca_var$PC = as.character(pca_var$PC); pca_var$PC=factor(pca_var$PC,levels=unique(pca_var$PC))
	return(list("PC"=pca,"Var"=pca_var))
}

## 2) Read GMT files

read.gmt.file = function(x){
	require(GSA)
	xgs = GSA.read.gmt(x)
	names(xgs$genesets) = xgs$geneset.names
	for(i in 1:length(xgs$genesets)){
		a = xgs$genesets[[i]]
		a = a[a!=""]
		xgs$genesets[[i]] = a
	}
	return(xgs)
}

## 3) Flatten correlation matrix

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

# usage: flatternCorrMatrix(res$r,res$P)

## 4) Pull out p-values from a linear regression

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

## 6) Scale data
normalize = function(D){
	norm_factors = sqrt(apply(D^2,1,sum))
	return(D / norm_factors);
}

## 7) CALC LD
calc_LD <- function( geno, inds=1:nrow(geno), get.D=T, get.Dprime=F, get.rsq=T, get.chisq=T, get.chisq_prime=F ) {
### Eva KF Chan
### 23 Feb 2009
### Last Modified: 29 Nov 2013
###
### Calculates D, D', r, chisq, chisq'
### Given locus A with allele frequencies pA & pa and locus B with allele frequencies pB and pb
### Let pAB be the allele frequencies of allele A/B.  As the AB/ab is indistinguishable from Ab/aB, we assume equal probability for either assortment; i.e. For individuals with Aa at locus A and Bb at locus b, we assume p(AB/ab)=p(Ab/aB)=0.5.  NOTE that this is assumption is part of the reason why this function is relatively fast, compare to, for example, the LD() function in R/genetics which estimates p(AB) using a maximum likelihood approach.
### D = pAB - pApB
### D' = { D/Dmax for D>=0 , where Dmax = min( pApb, papB )
###    = { D/Dmin for D<0  , where Dmin = max(-pApB,-papb )
### r = D / sqrt( pApapBpb )
### chi2 = (2ND^2) / (pApapBpb)
### chi2'= chisq / ( 2N(l-1) ) where l=min(k,m)
###                                  N = # individuals
###                                  k & m = # alelles in locus A & B
###
### Arguments:
### 	geno:            m x n matrix of genotypes {0,1,2,NA} where m=number of markers, n=number of individuals
###	inds:            integer vector of marker indices (rows of geno) for subseting markers for LD calculation
###	get.D:           {T,F} Boolean value to indicate whether the D measure is to be calculated
###	get.Dprime:      {T,F} Boolean value to indicate whether the D' measure is to be calculated
###	get.rsq:         {T,F} Boolean value to indicate whether the r^2 measure is to be calculated
###	get.chisq:       {T,F} Boolean value to indicate whether the chi2 measure is to be calculated
###	get.chisq_prime: {T,F} Boolean value to indicate whether the chi2' measure is to be calculated


	if( all(!get.D, !get.Dprime, !get.rsq, !get.chisq, !get.chisq_prime) ) { stop('Must request at least one LD statistic.\n') }
	D_prime <- rsq <- chisq <- chisq_prime <- df <- NULL
	D <- matrix(NA, nrow=nrow(geno), ncol=length(inds))
	if( get.Dprime ) { D_prime <- matrix(NA, nrow=nrow(geno), ncol=length(inds)) }
	if( get.rsq ) { rsq <- matrix(NA, nrow=nrow(geno), ncol=length(inds)) }
	if( get.chisq | get.chisq_prime ) {
		chisq <- matrix(NA, nrow=nrow(geno), ncol=length(inds))
		df <- matrix(NA, nrow=nrow(geno), ncol=length(inds))
		if( get.chisq_prime ) { chisq_prime <- matrix(NA, nrow=nrow(geno), ncol=length(inds)) }
	}

	if( all(as.logical(!is.na(geno))) ) {	#no missing data
		tmp.geno <- geno	## genotypes at locus A
		N <- ncol(tmp.geno)	#number of individuals (diploid is assumed)
		pA <- ((2*apply(tmp.geno==0,1,sum,na.rm=T))+apply(tmp.geno==1,1,sum,na.rm=T)) / (2*N)
		pa <- 1-pA
		for(i in 1:length(inds)) {
			tmp.Bgeno <- matrix(tmp.geno[inds[i],],nrow=nrow(tmp.geno),ncol=ncol(tmp.geno),byrow=T)	## genotypes at locus B
			pB <- ((2*apply(tmp.Bgeno==0,1,sum,na.rm=T))+apply(tmp.Bgeno==1,1,sum,na.rm=T)) / (2*N)
			pb <- 1-pB
			pAB <- ((apply(tmp.geno==0 & tmp.Bgeno==0, 1, sum,na.rm=T)*2) + (apply(tmp.geno==1 & tmp.Bgeno==0, 1, sum,na.rm=T)) + (apply(tmp.geno==0 & tmp.Bgeno==1, 1, sum,na.rm=T)) + (apply(tmp.geno==1 & tmp.Bgeno==1, 1, sum,na.rm=T)*0.5)) / (2*N)
			D[,i] <- pAB-(pA*pB)
			if( get.Dprime ) {
				Dmax <- pmin(pA*pb, pa*pB)
				Dmin <- pmax(-pA*pB, -pa*pb)
				pos <- (D[,i]>=0)
				D_prime[which(pos),i] <- D[which(pos),i] / Dmax[which(pos)]
				D_prime[which(!pos),i] <- D[which(!pos),i] / Dmin[which(!pos)]
			}
			if( get.rsq ) {
				rsq[,i] <- (D[,i]*D[,i]) / (pA*pa*pB*pb)
			}
			if( get.chisq | get.chisq_prime ) {
				chisq[,i] <- (2*N*D[,i]*D[,i]) / (pA*pa*pB*pb)
				if( get.chisq_prime ) {
					k=2-as.integer(pA==0|pa==0)
					m=2-as.integer(pB==0|pb==0)
					#df[,i] <- (k-1)*(m-1)
					chisq_prime[,i] <- chisq[,i] / (2*N*pmin(k,m))
				}
			}
		}
	} else {	#at least one missing data point in geno
		for(i in 1:length(inds)) {
			tmp.geno <- geno[,!is.na(geno[inds[i],])]	## genotypes at locus A; i.e. all loci, but excluding samples with missing data at lcous B (i)
			tmp.Bgeno <- matrix(tmp.geno[inds[i],],nrow=nrow(tmp.geno),ncol=ncol(tmp.geno),byrow=T)	## genotypes at locus B (i.e. i-th locus); pulling from tmp.geno, so samples with missing data at i-th locus (B) will also be excluded
			tmp.Bgeno[is.na(tmp.geno)] <- NA 	#anytime where locus A (i.e. all non i-th locus) is missing, set as missing
			N <- rowSums(!is.na(tmp.geno))
			pA <- ((2*apply(tmp.geno==0,1,sum,na.rm=T))+apply(tmp.geno==1,1,sum,na.rm=T)) / (2*N)
			pB <- ((2*apply(tmp.Bgeno==0,1,sum,na.rm=T))+apply(tmp.Bgeno==1,1,sum,na.rm=T)) / (2*N)
			pa <- 1-pA
			pb <- 1-pB
			pAB <- ((apply(tmp.geno==0 & tmp.Bgeno==0, 1, sum,na.rm=T)*2) + (apply(tmp.geno==1 & tmp.Bgeno==0, 1, sum,na.rm=T)) + (apply(tmp.geno==0 & tmp.Bgeno==1, 1, sum,na.rm=T)) + (apply(tmp.geno==1 & tmp.Bgeno==1, 1, sum,na.rm=T)*0.5)) / (2*N)
			D[,i] <- pAB-(pA*pB)
			if( get.Dprime ) {
				Dmax <- pmin(pA*pb, pa*pB)
				Dmin <- pmax(-pA*pB, -pa*pb)
				pos <- (D[,i]>=0)
				D_prime[which(pos),i] <- D[which(pos),i] / Dmax[which(pos)]
				D_prime[which(!pos),i] <- D[which(!pos),i] / Dmin[which(!pos)]
			}
			if( get.rsq ) {
				rsq[,i] <- (D[,i]*D[,i]) / (pA*pa*pB*pb)
			}
			if( get.chisq | get.chisq_prime ) {
				chisq[,i] <- (2*N*D[,i]*D[,i]) / (pA*pa*pB*pb)
				k=2-as.integer(pA==0|pa==0)
				m=2-as.integer(pB==0|pb==0)
				df[,i] <- (k-1)*(m-1)
				if( get.chisq_prime ) {
					chisq_prime[,i] <- chisq[,i] / (2*N*pmin(k,m))
				}
			}
		}
	}
	if( !get.D ) { D <- NULL }
	if( !get.chisq ) { chisq <- NULL }
	return(list(D=D, Dprime=D_prime, rsq=rsq, chisq=chisq, chisq_prime=chisq_prime, chisq_df=df))
}

### 8) preprocess genotype data

preprocess_genotype = function(data){

# Step 1: Subset samples with <5% NAs across all SNPs
	cat("Filter by missingness in SNPs (<5% NAs across the rows) \n")
	pNA.row=rowSums(is.na(data))/ncol(data)
	if(any(pNA.row>0.05)) data = data[which(pNA.row<0.05),]
	cat("SNPs: ",dim(data)[1]," Samples: ",dim(data)[2],"\n")

# Step 2: Subset samples with <10% NAs across all samples
	cat("Filter by missingness in individuals (<10% NAs across all samples) \n")
	pNA.col=colSums(is.na(data))/dim(data)[1]
	if(any(pNA.col>0.10)) data = data[,which(pNA.col<0.10)]
	cat("SNPs: ",dim(data)[1]," Samples: ",dim(data)[2],"\n")

# Step 3: Filter by minor allele frequency
	cat("Filtering SNPs with MAF <0.05 \n")
	n0 = apply(data==0,1,sum,na.rm=T)
	n1 = apply(data==1,1,sum,na.rm=T)
	n2 = apply(data==2,1,sum,na.rm=T)
	n = n0 + n1 + n2
	p = ((2*n0)+n1)/(2*n)
	q = 1 - p
	maf = pmin(p, q)
	if(any(maf<0.05)) data = data[which(maf>0.05),]
	cat("SNPs: ",dim(data)[1]," Samples: ",dim(data)[2],"\n")

## Step 4: HWE chi-sq test
	cat("Filtering SNPs by HWE chi-sq-pval<0.001 \n")
	n0 = apply(data==0,1,sum,na.rm=T)
	n1 = apply(data==1,1,sum,na.rm=T)
	n2 = apply(data==2,1,sum,na.rm=T)
	n = n0 + n1 + n2
	p = ((2*n0)+n1)/(2*n)
	q = 1 - p
	obs = cbind(n0=n0,n1=n1,n2=n2)
	exp = cbind(p*p, 2*p*q, q*q)
	exp = exp*n
	chisq = (obs-exp)
	chisq = (chisq*chisq) /exp
	hwe.chisq = apply(chisq,1,sum)
	hwe.chisq.p = 1-pchisq(hwe.chisq,df=1)
	if(any(hwe.chisq.p<0.001)) data = data[which(hwe.chisq.p<0.001),]
	cat("SNPs: ",dim(data)[1]," Samples: ",dim(data)[2],"\n")

## Finished preprocessing data
	return(data)
}

## 9) EIGENSTRAT
eigenstrat<-function(geno,maxMis=0,minMaf=0.01){
## geno: snp x ind matrix of genotypes \in 0,1,2
##maxMis maximum allowed missing genotypes for a site

  nMis <- rowSums(is.na(geno))
  freq <- rowMeans(geno,na.rm=T)/2               # get allele frequency
  keep <- freq>minMaf&freq<(1-minMaf) & nMis<=maxMis         # remove sites with non-polymorphic data
  freq<-freq[keep]
  geno<-geno[keep,]
  snp<-nrow(geno)                           #number of snps used in analysis
  ind<-ncol(geno)                           #number of individuals used in analuysis
  M <- (geno-freq*2)/sqrt(freq*(1-freq))       #normalize the genotype matrix
  M[is.na(M)]<-0
  X<-t(M)%*%M                               #get the (almost) covariance matrix
  X<-X/(sum(diag(X))/(snp-1))
  E<-eigen(X)

  mu<-(sqrt(snp-1)+sqrt(ind))^2/snp         #for testing significance (assuming no LD!)
  sigma<-(sqrt(snp-1)+sqrt(ind))/snp*(1/sqrt(snp-1)+1/sqrt(ind))^(1/3)
  E$TW<-(E$values[1]*ind/sum(E$values)-mu)/sigma
  E$mu<-mu
  E$sigma<-sigma
  E$nSNP <- nrow(geno)
  E$nInd <- ncol(geno)
  class(E)<-"eigenstrat"
  E
}
plot.eigenstrat<-function(x,col=1,...)
  plot(x$vectors[,1:2],col=col,...)

print.eigenstrat<-function(x)
  cat("statistic",x$TW,"\n")

## 9) Raw counts to TPM
countToTpm <- function(counts, effLen)
	{
	    rate <- log(counts) - log(effLen)
	    denom <- log(sum(exp(rate)))
	    exp(rate - denom + log(1e6))
	}

## 10) Raw counts to FPKM
countToFpkm <- function(counts, effLen)
	{
	    N <- sum(counts)
	    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
	}

## 11) FPKM to TPM
fpkmToTpm <- function(fpkm)
	{
	    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}

## 12) Raw counts to effective counts
countToEffCounts <- function(counts, len, effLen)
	{
	    counts * (len / effLen)
	}

## 13) Function to tSNE dim reduction and density clustering
cluster = function(FM, FMclin, scale=T, max_components = 2, num_dim = 3, num_clusters = 4, perp = 40, plot=F, seed=123){

  if(scale) FM = as.matrix(Matrix::t(scale(Matrix::t(FM))))
  set.seed(seed)
  irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
  irlba_pca_res <- irlba_res$x
  topDim_pca <- irlba_pca_res
  set.seed(seed)
  tsne_res <- Rtsne(as.matrix(topDim_pca), dims = max_components,pca = F,perplexity=perp)
  tsne_data <- tsne_res$Y[, 1:max_components]
  colnames(tsne_data) = c("Component_1", "Component_2")
  rownames(tsne_data) = colnames(FM)
  ## Clustering
  dataDist <- dist(tsne_data)
  set.seed(seed)
  dataClust <- densityClust(dataDist, gaussian = T)
  delta_rho_df <- data.frame(delta = dataClust$delta, rho = dataClust$rho)
  rho_threshold <- 0
  delta_threshold <- sort(delta_rho_df$delta, decreasing = T)[num_clusters] - .Machine$double.eps
  dataClust <- densityClust::findClusters(dataClust, rho = rho_threshold, delta = delta_threshold)
  tsne_data = data.frame(tsne_data)
  tsne_data$Cluster = as.factor(dataClust$cluster)
  tsne_data$Cluster = as.factor(dataClust$cluster)
  # tcenters = tsne_data[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))
  FMclin$Cluster = as.factor(dataClust$cluster)
  if(plot){
    ggplot(data.frame(tsne_data),aes(x=Component_1,y=Component_2,colour=Cluster,label=Cluster)) + geom_point(size=1) + geom_point(data = tcenters, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) +
          geom_text(data=tcenters,mapping = aes(label = Cluster), colour="black",size = 6) + theme_pander(boxes=T)
    # ggsave("tSNE_densClust_plot.pdf",width=11.1,height=6.64,units="in")
  }
  return(FMclin)
}

## 14) Function to CV filter
CVfilter = function(dat,plot=F,threshold=0.75){
  mean_values=apply(dat,MARGIN=1,FUN="mean")
  sd_values=apply(dat,MARGIN=1,FUN="sd")
  CV_values=sd_values/mean_values
  quantile_cut_off=quantile(CV_values,probs=threshold)
  cv_dat = data.frame("Means"=mean_values,"Sigma"=sd_values,"Variance"=apply(dat,MARGIN=1,FUN="var"),"CV"=CV_values)
  cv_dat$Dispersion = cv_dat$Variance/cv_dat$Means
  if(plot) ggplot(aes(x=Means, y=Variance, col=CV, size=CV), data=cv_dat) + geom_point() + theme_pander()
  return(dat[CV_values>quantile_cut_off,])
}

## 15) Basic function to convert mouse to human gene names
#musGenes <- c("Hmmr", "Tlx3", "Cpeb4")
##
convertMouseGeneList <- function(x){
	require("biomaRt")
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

	genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
	humanx <- unique(genesV2[, 2])

	# Print the first 6 genes found to the screen
	return(humanx)
}

## 16) Basic function to convert human to mouse gene names

convertHumanGeneList <- function(x){
	require("biomaRt")
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	humanx <- unique(genesV2[, 2])

	# Print the first 6 genes found to the screen
	return(humanx)
}

## 17) Permutation test for PCA

pca_eigenperm<- function(data, nperm = 1000){
        pca_out<- prcomp(data, scale. = T)
        eigenperm<- data.frame(matrix(NA, nperm, ncol(data)))
        n<- ncol(data)
        data_i<- data.frame(matrix(NA, nrow(data), ncol(data)))
        for (j in 1: nperm){
        for (i in 1:n){
                data_i[,i]<- sample(data[,i], replace = F)
        }
        pca.perm<- prcomp(data_i, scale. = T)
        eigenperm[j,]<- pca.perm$sdev^2
        }
        colnames(eigenperm)<- colnames(pca_out$rotation)
        eigenperm

}
