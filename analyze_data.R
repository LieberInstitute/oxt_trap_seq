###########################
## analysis of oxt data  ##
## andrew jaffe ###########

## libraries
library(jaffelab)
library(edgeR)
library(SummarizedExperiment)
library(recount)

## load in data
load("rseObjs_oxtMerge_n18_4features.Rdata")

## filter lowly exprssed genes by RPKM overall
geneIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.1
rse_gene = rse_gene[geneIndex,]

## pca
pca = prcomp(t(log2(getRPKM(rse_gene, "Length")+1)))

####################################
## start w/ gene level, ovation
rse_gene_ovation = rse_gene[,rse_gene$Kit == "Ovation"]
rse_gene_ovation$Fraction[rse_gene_ovation$Fraction == "Input"] = "IP"

## pca

## modeling
mod_ovation = model.matrix(~Fraction + totalAssignedGene, 
	data=colData(rse_gene_ovation))

##### DGE ######
dge_ovation = DGEList(counts = assays(rse_gene_ovation)$counts, 
	genes = rowData(rse_gene_ovation))
dge_ovation = calcNormFactors(dge_ovation)
vGene_ovation = voom(dge_ovation,mod_ovation,plot=TRUE)

fitGene_ovation = lmFit(vGene_ovation)
eBGene_ovation = eBayes(fitGene_ovation)
outGene_ovation = topTable(eBGene_ovation,coef=2,
	p.value = 0.05,number=nrow(rse_gene_ovation))
ebGene_ovation = ebayes(fitGene_ovation)
