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

## load in data
load("rseObjs_oxtMerge_n18_4features.Rdata")

## filter lowly exprssed genes by RPKM overall
geneIndex = rowMeans(getRPKM(rse_gene, "Length")) > 0.1
rse_gene = rse_gene[geneIndex,]

## pca
pca = prcomp(t(log2(getRPKM(rse_gene, "Length")+1)))
plot(pca$x, pch = 21, cex=1.5, bg = factor(rse_gene$CellType))

####################################
## start w/ gene level, ovation
rse_gene_ovation = rse_gene[,rse_gene$Kit == "Ovation"]

## modeling
mod_ovation = model.matrix(~Fraction + totalAssignedGene, 
	data=colData(rse_gene_ovation))

##### DGE ######
dge_ovation = DGEList(counts = assays(rse_gene_ovation)$counts, 
	genes = rowData(rse_gene_ovation))
dge_ovation = calcNormFactors(dge_ovation)

## mean-variance
vGene_ovation = voom(dge_ovation,mod_ovation,plot=TRUE)

## do analysis
fitGene_ovation = lmFit(vGene_ovation)
ebGene_ovation = ebayes(fitGene_ovation)

## top table
eBGene_ovation = eBayes(fitGene_ovation)
outGene_ovation = topTable(eBGene_ovation,coef=2,
	p.value = 1,number=nrow(rse_gene_ovation))
outGene_ovation$sigColor = as.numeric(outGene_ovation$adj.P.Val < 0.01)+1
sigGene_ovation = outGene_ovation[outGene_ovation$adj.P.Val < 0.01,]

## FDR < 1% leads to better replication w/ solo, 95% versus 80%

################
## plots #######

## ma-plot
plot(logFC ~ AveExpr, pch = 21, bg=sigColor, 
	data = outGene_ovation)
	
## volanco plot
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene_ovation)
	
################
## gene set ####
sigGene_ = split(sigGene_ovation$EntrezID,
	sign(sigGene_ovation$logFC))
sigGeneList_ovation = lapply(sigGeneList_ovation, 
			function(x) x[!is.na(x)])
geneUniverse = as.character(outGene_ovation$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

goBP_ovation <- enrichGO(gene = as.character(sigGeneList_ovation[["1"]]),
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)
goMF_ovation <- enrichGO(gene = as.character(sigGeneList_ovation[["1"]]),
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)
goCC_ovation <- enrichGO(gene = as.character(sigGeneList_ovation[["1"]]),
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)
			
###########################
## replication with solo ##
rse_gene_solo = rse_gene[,rse_gene$Kit == "SoLo"]

## modeling
mod_solo = model.matrix(~Fraction + totalAssignedGene, 
	data=colData(rse_gene_solo))

##### DGE ######
dge_solo = DGEList(counts = assays(rse_gene_solo)$counts, 
	genes = rowData(rse_gene_solo))
dge_solo = calcNormFactors(dge_solo)

## mean-variance
vGene_solo = voom(dge_solo,mod_solo,plot=TRUE)

## do analysis
fitGene_solo = lmFit(vGene_solo)
ebGene_solo = ebayes(fitGene_solo)

################
## compare #####

## top table
eBGene_solo = eBayes(fitGene_solo)
outGene_solo = topTable(eBGene_solo,coef=2,
	p.value = 1,number=nrow(rse_gene_solo))

## match up stats
replicationStats = outGene_solo[rownames(sigGene_ovation),]
replicationStats$sameDir = sign(replicationStats$logFC) == 
	sign(sigGene_ovation$logFC)
prop.table(table(replicationStats$sameDir &
	replicationStats$P.Value < 0.05))
	
## plots
plot(sigGene_ovation$logFC, replicationStats$logFC,
	xlab = "Ovation", ylab = "SoLo", main = "log2 Fold Changes")
plot(sigGene_ovation$t, replicationStats$t)

####################################3
### CST comparison to ribo null? ####

rse_gene_clone = rse_gene[,rse_gene$Kit == "Clonetech"]

## modeling
mod_clone = model.matrix(~Fraction + totalAssignedGene, 
	data=colData(rse_gene_clone))

##### DGE ######
dge_clone = DGEList(counts = assays(rse_gene_clone)$counts, 
	genes = rowData(rse_gene_clone))
dge_clone = calcNormFactors(dge_clone)

## mean-variance
vGene_clone = voom(dge_clone,mod_clone,plot=TRUE)

## do analysis
fitGene_clone = lmFit(vGene_clone)
ebGene_clone = ebayes(fitGene_clone)

## top table
eBGene_clone = eBayes(fitGene_clone)
outGene_clone = topTable(eBGene_clone,coef=2,
	p.value = 1,number=nrow(rse_gene_clone))

## match up stats
compareStats = outGene_clone[rownames(sigGene_ovation),]
compareStats$sameDir = sign(compareStats$logFC) == 
	sign(sigGene_ovation$logFC)
prop.table(table(compareStats$sameDir &
	compareStats$P.Value < 0.01))
	
## plots
plot(sigGene_ovation$logFC, compareStats$logFC,
	xlab = "Ovation", ylab = "SoLo", main = "log2 Fold Changes")
plot(sigGene_ovation$t, compareStats$t)

## check GO again
cleanGene = sigGene_ovation$EntrezID[sigGene_ovation$logFC > 0 &
	!(compareStats$sameDir & compareStats$P < 0.05)]
cleanGene = as.character(cleanGene[!is.na(cleanGene)])

goBP_clean <- enrichGO(gene = cleanGene,
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)
goMF_clean <- enrichGO(gene = cleanGene,
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)
goCC_clean <- enrichGO(gene = cleanGene,
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
				readable= TRUE)