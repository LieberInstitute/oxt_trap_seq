###########################
## analysis of oxt data  ##
## andrew jaffe ###########

## libraries
library(jaffelab) # for helper functions
library(edgeR) # for limma+DGEList
library(SummarizedExperiment) # for RSEs
library(recount) # for getRPKM()
library(clusterProfiler) # gene set tests
library(org.Mm.eg.db) # gene set test annotation
library(RColorBrewer) # for plotting

## load in data
load("preprocessed_data/Bdnf_SoLo_Oxt/rse_gene_Bdnf_SoLo_Oxt_MiSeq_n6.Rdata")
load("preprocessed_data/Bdnf_SoLo_Oxt/rse_exon_Bdnf_SoLo_Oxt_MiSeq_n6.Rdata")
load("preprocessed_data/Bdnf_SoLo_Oxt/rse_jx_Bdnf_SoLo_Oxt_MiSeq_n6.Rdata")

## filter lowly exprssed genes by RPKM in ovation
geneIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.1
rse_gene = rse_gene[geneIndex,]

## phenotype
colData(rse_gene) = colData(rse_gene)
colData(rse_gene)$Genotype = factor(c("WT","BDNF", "WT","WT","BDNF","BDNF"))
colData(rse_gene)$Genotype = relevel(colData(rse_gene)$Genotype, "WT")
colData(rse_gene)$Run = c(1,1,2,2,2,2)

## pca
pca = prcomp(t(log2(getRPKM(rse_gene, "Length")+1)))
plot(pca$x, pch = 21, cex=1.5, bg = factor(rse_gene$Genotype))

####################################
## start w/ gene level

## modeling
oxt = log2(getRPKM(rse_gene, "Length")["ENSMUSG00000027301.7",]+1)
mod = model.matrix(~Genotype + totalAssignedGene + oxt, 
	data=colData(rse_gene))

##### DGE ######
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)

## do analysis
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)

## top table
eBGene = eBayes(fitGene)
outGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene$sigColor = as.numeric(outGene$adj.P.Val < 0.01)+1
sigGene = outGene[outGene$adj.P.Val < 0.1,]
write.csv(sigGene[,-1], "tables/bdnf_genotype_effect_fdr10.csv",
	row.names=FALSE)
write.csv(outGene[,-1], "tables/bdnf_genotype_effect_allGenes.csv",
	row.names=FALSE)

## no oxt adjustment
mod0 = model.matrix(~Genotype + totalAssignedGene, 
	data=colData(rse_gene))
vGene0 = voom(dge,mod0,plot=TRUE)
fitGene0 = lmFit(vGene0)
ebGene0 = ebayes(fitGene0)

## top table
eBGene0 = eBayes(fitGene0)
outGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
sigGene0= outGene0[outGene0$adj.P.Val < 0.1,]
write.csv(sigGene0[,-1], "tables/bdnf_genotype_effect_fdr10_noOxtAdj.csv",
	row.names=FALSE)

################
## plots #######

## ma-plot
par(mfrow = c(1,2))
palette(brewer.pal(5, "Dark2"))
plot(logFC ~ AveExpr, pch = 21, bg=sigColor, 
	data = outGene,ylab="IP vs Input log2FC")
	
## volanco plot
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene, xlab = "IP vs Input log2FC")
	
################
## gene set ####
sigGene_marginal = outGene[outGene$P.Value < 0.005,]
sigGeneList = split(sigGene_marginal$EntrezID,
	sign(sigGene_marginal$logFC))
sigGeneList = lapply(sigGeneList, 
			function(x) x[!is.na(x)])
geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

goBP <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)
goMF <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)		
goCC <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)

## write out				
goBP_DF = as.data.frame(goBP)
goBP_DF$Ontology = "BP"
goMF_DF = as.data.frame(goMF)
goMF_DF$Ontology = "MF"
goCC_DF = as.data.frame(goCC)
goCC_DF$Ontology = "CC"

goOut = rbind(goBP_DF, goMF_DF, goCC_DF)
colnames(goOut)[1] = "Direction"
goOut = data.frame(BioProc = goBP$Description[1:15],
	MolFunc = goMF$Description[1:15],
	CellComp = goCC$Description[1:15],
	stringsAsFactors=FALSE)
goOut

##############################
## compare to trap inputs ####
outGene_ovation = outGene_ovation[match(rownames(outGene), rownames(outGene_ovation)),]
plot(outGene$t, outGene_ovation$t)
plot(outGene$logFC, outGene_ovation$logFC,
	xlab = "BDNF Effect", ylab = "IP vs Input Effect",
	main = "log2FC",pch=21,bg="grey")
	
plot(outGene0$t, outGene_ovation$t)
plot(outGene0$logFC, outGene_ovation$logFC,
	xlab = "BDNF Effect", ylab = "IP vs Input Effect",
	main = "log2FC",pch=21,bg="grey")
	
