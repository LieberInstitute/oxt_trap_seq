###
library(SummarizedExperiment)
library(recount)
library(jaffelab)

## load data
load("preprocessed_data/Bdnf_SoLo_Oxt/rse_gene_Bdnf_SoLo_Oxt_MiSeq_n6.Rdata")
rse_gene_bdnf = rse_gene
load("rseObjs_oxtMerge_n18_4features.Rdata")
rse_gene_ip = rse_gene

## merge
geneRpkm = cbind(getRPKM(rse_gene_bdnf, "Length"), getRPKM(rse_gene_ip, "Length"))
x = c(c("Oxt","Bdnf", "Oxt","Oxt","Bdnf","Bdnf"), rse_gene_ip$CellType)
experiment = c(rep("Bdnf", 6), rse_gene_ip$Kit)

pca = prcomp(t(log2(geneRpkm+1)))
pcaVars = getPcaVars(pca)

palette(brewer.pal(5,"Dark2"))
plot(pca$x, pch = 20 + as.numeric(factor(experiment)),
	bg= factor(x),cex=2,
	xlab=paste0("PC1: ", pcaVars[1],"% Var Expl"),
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"))
legend("bottomleft", levels(factor(x)), col = 1:4, pch = 15,cex=2)
legend("bottomright", levels(factor(experiment)), pch = 21:24, 
	pt.bg = "black",pt.cex=2,cex=1.5)

table(experiment)