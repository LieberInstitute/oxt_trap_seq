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
library(TeachingDemos) # for shadow text

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
outGene$sigColor = as.numeric(outGene$adj.P.Val < 0.1)+1
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
	
## volanco plot, thresholded 
outGeneThresh = outGene
outGeneThresh$P.Value[outGeneThresh$P.Value < 2.2e-16] = 2.2e-16
outGeneThresh$sigColor[outGeneThresh$sigColor==2] = 3

# highlight genes
g = c("Gabra2", "Cckar")
m = match(g, outGeneThresh$Symbol)

pdf("plots/figure5c_volcano_plot_bdnfEffect.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGeneThresh, xlab = "BDNF Mut vs WT log2FC")
shadowtext(outGeneThresh$logFC[m], -log10(outGeneThresh$P.Value[m]),
	letters[16:17],font=2,cex=2,col="grey")
dev.off()

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
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)
goMF <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)		
goCC <- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)

## write out				
goBP_DF = as.data.frame(goBP)
goMF_DF = as.data.frame(goMF)
goCC_DF = as.data.frame(goCC)

goOut = rbind(goBP_DF, goMF_DF, goCC_DF)
goOut$Ontology = rep(c("BP", "MF", "CC"), 
	times = c(nrow(goBP), nrow(goMF), nrow(goCC)))
goOut = goOut[goOut$p.adjust < 0.05,]
colnames(goOut)[1] = "Direction"
goOut = goOut[order(goOut$p.adjust),]
write.csv(goOut, file = "tables/GO_bdnf_enrichment_DE_at_p005.csv",row.names=FALSE)

## example GO plots
goOut = read.csv("tables/GO_bdnf_enrichment_DE_at_p005.csv",as.is=TRUE)

setsToPlot = c("GO:0098793", "GO:0060077", "GO:0005184",
		"GO:0043679", "GO:0021761", "GO:0003954",
		"GO:0042773", "GO:0019199", "GO:0007270", "GO:0050808")
goUp = goOut[goOut$Direction == 1,]
goUp = goUp[match(setsToPlot[6:10], goUp$ID),]	
goDown = goOut[goOut$Direction == -1,]
goDown = goDown[match(setsToPlot[1:5], goDown$ID),]		
goExample = rbind(goUp, goDown)

pdf("plots/go_figure5_barplot.pdf",h=4,w=8)
par(mar=c(5,23,2,2),cex.axis=1.2,cex.lab=1.5)
barplot(-log10(goExample$qvalue),width=0.5,
	names = goExample$Description,horiz=TRUE,
	xlab="-log10(FDR)",las=1)
abline(h=3.05,lwd=1.5,lty=2)
abline(v=-log10(0.05), col="blue")
dev.off()
		
##############################
## compare to trap inputs ####
outGene_ovation = read.csv("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/tables/TRAP_discovery_differential_expression.csv.gz",
	row.names=1)
outGene_ovation = outGene_ovation[match(rownames(outGene), rownames(outGene_ovation)),]

pdf("plots/suppFigure_Bdnf_versus_IP_effects.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(outGene$logFC, outGene_ovation$logFC,
	xlab = "BDNF Effect", ylab = "Oxt IP vs Input Effect",
	main = "log2FC",pch=21,bg="grey")
plot(outGene$t, outGene_ovation$t,
	xlab = "BDNF Effect", ylab = "Oxt IP vs Input Effect",
	main = "T-statistics",pch=21,bg="grey")
dev.off()


