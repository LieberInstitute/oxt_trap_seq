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
load("rseObjs_oxtMerge_n18_4features.Rdata")

## filter lowly exprssed genes by RPKM in ovation
geneIndex = rowMeans(getRPKM(rse_gene[,rse_gene$Kit=="Ovation"], "Length")) > 0.1
rse_gene = rse_gene[geneIndex,]

## pca
pca = prcomp(t(log2(getRPKM(rse_gene, "Length")+1)))
plot(pca$x, pch = 21, cex=1.5, bg = factor(rse_gene$CellType))

####################################
## start w/ gene level, ovation
rse_gene_ovation = rse_gene[,rse_gene$Kit == "Ovation"]

t.test(rse_gene_ovation$totalAssignedGene ~ rse_gene_ovation$Fraction)
tapply(rse_gene_ovation$numMapped, rse_gene_ovation$Fraction,mean)

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
dim(sigGene_ovation)

save(outGene_ovation, file= "TRAP_discovery_differential_expression.rda")

## FDR < 1% leads to better replication w/ solo, 95% versus 80%

################
## plots #######

## ma-plot
palette(brewer.pal(5, "Dark2"))
plot(logFC ~ AveExpr, pch = 21, bg=sigColor, 
	data = outGene_ovation,ylab="IP vs Input log2FC")
	
## volanco plot
pdf("plots/figure4c_volcano_plot_trapSeq.pdf")
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene_ovation, xlab = "IP vs Input log2FC")
abline(v=c(-1,1), lty=2,lwd=2)
dev.off()

## genes to talk about in text
outGene_ovation[match(c("Oxt", "Agrp", "Cartpt", "Pmch"), outGene_ovation$Symbol),]
2^outGene_ovation$logFC[match(c("Oxt", "Agrp", "Cartpt", "Pmch"), outGene_ovation$Symbol)]
1/2^outGene_ovation$logFC[match(c("Oxt", "Agrp", "Cartpt", "Pmch"), outGene_ovation$Symbol)]

################
## gene set ####
sigGeneList_ovation = split(sigGene_ovation$EntrezID,
	sign(sigGene_ovation$logFC))
sigGeneList_ovation = lapply(sigGeneList_ovation, 
			function(x) x[!is.na(x)])
geneUniverse = as.character(outGene_ovation$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

goBP <- compareCluster(sigGeneList_ovation, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)
goMF <- compareCluster(sigGeneList_ovation, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)		
goCC <- compareCluster(sigGeneList_ovation, fun = "enrichGO",
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
colnames(goOut)[1] = "Direction"
goOut = goOut[order(goOut$p.adjust),]
goOut[which(goOut$Direction == "1")[1:10],]
write.csv(goOut, file = "tables/trap_GO_analysis_DE_FDR01.csv",row.names=FALSE)

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
	
## check proportions
prop.table(table(replicationStats$sameDir))
prop.table(table(replicationStats$sameDir &
	replicationStats$P.Value < 0.05))
	
## plots
plot(sigGene_ovation$logFC, replicationStats$logFC,bg=2, pch=21,
	xlab = "Ovation", ylab = "SoLo", main = "log2 Fold Changes")
abline(h=0,v=0,lty=2)
plot(sigGene_ovation$t, replicationStats$t,bg=2, pch=21,
	xlab = "Ovation", ylab = "SoLo", main = "T-statistics")
abline(h=0,v=0,lty=2)


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
compareStatsClone = outGene_clone[rownames(sigGene_ovation),]
compareStatsClone$sameDir = sign(compareStatsClone$logFC) == 
	sign(sigGene_ovation$logFC)

	
###################################
#### other TRAP data ##############
load("preprocessed_data/Nectow_TRAP/rse_gene_Nectow_TRAP_Single_n4.Rdata")
pdTrap = read.delim("preprocessed_data/Nectow_TRAP/Nectow_SraRunTable.txt",as.is=TRUE)
rse_gene$GEO_ID = pdTrap$Sample_Name_s[match(colnames(rse_gene), pdTrap$Run_s)]
rse_gene_trap = rse_gene[rownames(rse_gene_ovation),]
mod_trap = model.matrix(~c("IP","IP","Input","Input") + rse_gene_trap$totalAssignedGene)
colnames(mod_trap)[2] = "Fraction"

dge_trap = DGEList(counts = assays(rse_gene_trap)$counts, 
	genes = rowData(rse_gene_ovation))
dge_trap = calcNormFactors(dge_trap)

vGene_trap = voom(dge_trap,mod_trap,plot=TRUE)
fitGene_trap = lmFit(vGene_trap)
ebGene_trap = ebayes(fitGene_trap)
eBGene_trap = eBayes(fitGene_trap)
outGene_trap = topTable(eBGene_trap,coef=2,
	p.value = 1,number=nrow(rse_gene_trap))
compareStatsTrap = outGene_trap[rownames(sigGene_ovation),]

###################
# compare plots ###
###################

prop.table(table(compareStatsClone$sameDir))
prop.table(table(compareStatsClone$sameDir &
	compareStatsClone$P.Value < 0.05))
	
## Cst
plot(sigGene_ovation$logFC, compareStatsClone$logFC,bg=2, pch=21,
	xlab = "Oxt", ylab = "Cst", main = "log2 Fold Changes")
abline(h=0,v=0,lty=2)

plot(sigGene_ovation$t, compareStatsClone$t,bg=2, pch=21,
	xlab = "Oxt", ylab = "Cst", main = "T-statistics")
abline(h=0,v=0,lty=2)
abline(h=10,lty=3)

## retina
plot(sigGene_ovation$logFC, compareStatsTrap$logFC,bg=2, pch=21,
	xlab = "Oxt vs Input", ylab = "Ntsr1 vs Input", main = "log2 Fold Changes")
abline(h=0,v=0,lty=2)

plot(sigGene_ovation$t, compareStatsTrap$t,bg=2, pch=21,
	xlab = "Oxt vs Input", ylab = "Ntsr1 vs Input", main = "T-statistics")
abline(h=0,v=0,lty=2)
abline(h=10,lty=3)

##### write output ######
sigGene_ovation$t_SoLo = replicationStats$t
sigGene_ovation$logFC_SoLo = replicationStats$logFC
sigGene_ovation$P.Value_SoLo = replicationStats$P.Value

sigGene_ovation$t_CST = compareStatsClone$t
sigGene_ovation$logFC_CST = compareStatsClone$logFC
sigGene_ovation$P.Value_CST = compareStatsClone$P.Value

sigGene_ovation$t_Ntsr1 = compareStatsTrap$t
sigGene_ovation$logFC_Ntsr1 = compareStatsTrap$logFC
sigGene_ovation$P.Value_Ntsr1 = compareStatsTrap$P.Value

outStats = sigGene_ovation
outStats$sigColor = NULL
outStats$chr_mm10 = as.character(seqnames(rse_gene_ovation[rownames(outStats)]))
outStats$start_mm10 = start(rse_gene_ovation[rownames(outStats)])
outStats$end_mm10 = end(rse_gene_ovation[rownames(outStats)])
outStats$strand_mm10 = as.character(strand(rse_gene_ovation[rownames(outStats)]))
outStats$bed_mm10 = paste0(outStats$chr_mm10, ":", 
	outStats$start_mm10, "-", outStats$end_mm10)
outStats$coord_mm10 = paste0(outStats$chr_mm10, ":", 
	outStats$start_mm10, "-", outStats$end_mm10, "(", outStats$strand_mm10,")")

## filter
outStats = outStats[,c(2, 5, 8:13, 28, 1, 3:4, 6, 14:22, 27)]
write.csv(outStats, file= "tables/TRAP_stats_allTests.csv", row.names=FALSE)

 	write.csv(sigGene_ovation, file = 
compareStatsClone$Symbol[compareStatsClone$t > 10]


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