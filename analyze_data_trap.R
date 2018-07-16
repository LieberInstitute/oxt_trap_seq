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
load("rseObjs_oxtMerge_n18_4features.Rdata")

## filter lowly exprssed genes by RPKM in ovation
geneIndex = rowMeans(getRPKM(rse_gene[,rse_gene$Kit=="Ovation"], "Length")) > 0.1
rse_gene = rse_gene[geneIndex,]
dim(rse_gene)

## pca
pca = prcomp(t(log2(getRPKM(rse_gene, "Length")+1)))
plot(pca$x, pch = 21, cex=1.5, bg = factor(rse_gene$CellType))

####################################
## start w/ gene level, ovation
rse_gene_ovation = rse_gene[,rse_gene$Kit == "Ovation"]
rse_gene_all = rse_gene # save for below

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
statsOut = outGene_ovation[,c(2,5, 8,10:13,3,1,4,6,9)]
write.csv(statsOut, file = gzfile("tables/TRAP_discovery_differential_expression.csv.gz"),
	row.names=FALSE)
## FDR < 1% leads to better replication w/ solo, 95% versus 80%

################
## plots #######

## ma-plot
palette(brewer.pal(5, "Dark2"))
plot(logFC ~ AveExpr, pch = 21, bg=sigColor, 
	data = outGene_ovation,ylab="IP vs Input log2FC")
	
## volanco plot
g = c("Oxt", "Agrp", "Cartpt", "Pmch")
m = match(g, outGene_ovation$Symbol)

pdf("plots/figure4c_volcano_plot_trapSeq.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene_ovation, xlab = "IP vs Input log2FC")
shadowtext(outGene_ovation$logFC[m], -log10(outGene_ovation$P.Value[m]),
	letters[21:24],font=2,cex=2,col="grey")
abline(v=c(-1,1), lty=2,lwd=2)
dev.off()

## genes to talk about in text
outGene_ovation[m,]
2^outGene_ovation$logFC[m]
1/2^outGene_ovation$logFC[m]

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
goOut = goOut[goOut$p.adjust < 0.01,]
colnames(goOut)[1] = "Direction"
goOut = goOut[order(goOut$p.adjust),]
goOut[which(goOut$Direction == "1")[1:10],]
write.csv(goOut, file = "tables/trap_GO_analysis_DE_FDR01.csv",row.names=FALSE)

## plots
goOut = read.csv("tables/trap_GO_analysis_DE_FDR01.csv",as.is=TRUE)
setsToPlot = c("GO:0043209", "GO:0044455", "GO:0042063",
	"GO:0010001", "	GO:0048709",
	"GO:0003779", "GO:0005516","GO:0098794",
	"GO:0030424", "GO:0000226")
	
goUp = goOut[goOut$Direction == 1,]
goUp = goUp[match(setsToPlot[6:10], goUp$ID),]	
goDown = goOut[goOut$Direction == -1,]
goDown = goDown[match(setsToPlot[1:5], goDown$ID),]		
goExample = rbind(goUp, goDown)

goExample$qvalue[goExample$qvalue < 2.2e-16] = 2.2e-16

pdf("plots/go_figure4_barplot.pdf",h=3.5,w=7)
par(mar=c(5,20,2,2),cex.axis=1.2,cex.lab=1.5)
barplot(-log10(goExample$qvalue),width=0.5,
	names = goExample$Description,horiz=TRUE,
	xlab="-log10(FDR)",las=1)
abline(h=3.05,lwd=1.5,lty=2)
abline(v=-log10(0.05), col="blue")
dev.off()
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
signif(prop.table(table(replicationStats$sameDir)),3)
signif(prop.table(table(replicationStats$sameDir &
	replicationStats$P.Value < 0.05)),3)
	
## plots
pdf("plots/suppFigure_ovation_versus_solo_replication.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(sigGene_ovation$logFC, replicationStats$logFC,bg=2, pch=21,
	xlab = "Ovation (Discovery)", ylab = "SoLo (Replication)", 
	main = "log2 Fold Changes") # i can change labels
abline(h=0,v=0,lty=2)
plot(sigGene_ovation$t, replicationStats$t,bg=2, pch=21,
	xlab = "Ovation (Discovery)", ylab = "SoLo (Replication)", 
	main = "T-statistics")
abline(h=0,v=0,lty=2)
dev.off()


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

## plots
pdf("plots/suppFigure_oxt_versus_cst.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(sigGene_ovation$logFC, compareStatsClone$logFC,bg=2, pch=21,
	xlab = "Oxt (Discovery)", ylab = "CST (Replication)", 
	main = "log2 Fold Changes") # i can change labels
abline(h=0,v=0,lty=2)
plot(sigGene_ovation$t, compareStatsClone$t,bg=2, pch=21,
	xlab = "Oxt (Discovery)", ylab = "CST", 
	main = "T-statistics")
abline(h=0,v=0,lty=2)
dev.off()


	
###################################
#### other TRAP data ##############
load("preprocessed_data/Nectow_TRAP/rse_gene_Nectow_TRAP_Single_n4.Rdata")
pdTrap = read.delim("preprocessed_data/Nectow_TRAP/Nectow_SraRunTable.txt",as.is=TRUE)
rse_gene$GEO_ID = pdTrap$Sample_Name_s[match(colnames(rse_gene), pdTrap$Run_s)]
rse_gene_trap = rse_gene[rownames(rse_gene_ovation),]
rse_gene_trap$CellType = rep(c("Ntsr1", "Bulk"), each=2)
mod_trap = model.matrix(~CellType + totalAssignedGene,
	data = colData(rse_gene_trap))

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


## plots
pdf("plots/suppFigure_oxt_versus_ntsr1.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(sigGene_ovation$logFC, compareStatsTrap$logFC,bg=2, pch=21,
	xlab = "Oxt (Discovery)", ylab = "Ntsr1", 
	main = "log2 Fold Changes") # i can change labels
abline(h=0,v=0,lty=2)
plot(sigGene_ovation$t, compareStatsTrap$t,bg=2, pch=21,
	xlab = "Oxt (Discovery)", ylab = "Ntsr1", 
	main = "T-statistics")
abline(h=0,v=0,lty=2)
dev.off()

pdf("plots/suppFigure_cst_versus_ntsr1_onlyOxtGenes.pdf",
	useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(compareStatsClone$logFC, compareStatsTrap$logFC,bg=2, pch=21,
	xlab = "Cst", ylab = "Ntsr1", 
	main = "log2 Fold Changes") # i can change labels
abline(h=0,v=0,lty=2)
plot(compareStatsClone$t, compareStatsTrap$t,bg=2, pch=21,
	xlab = "Cst", ylab = "Ntsr1", 
	main = "T-statistics")
abline(h=0,v=0,lty=2)
dev.off()

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

######
## add stats for Oxt specificity
rse_gene_ip = rse_gene_all[,rse_gene_all$Fraction == "IP"]
rse_gene_ntsr1 = rse_gene_trap[,rse_gene_trap$CellType == "Ntsr1"]

n = intersect(colnames(colData(rse_gene_ip)),
	colnames(colData(rse_gene_trap)))
colData(rse_gene_ip) = colData(rse_gene_ip)[,n]
colData(rse_gene_ntsr1) = colData(rse_gene_ntsr1)[,n]
counts = cbind(assays(rse_gene_ip)$counts, 
		assays(rse_gene_ntsr1)$counts)
colnames(counts) = c(colnames(rse_gene_ip), colnames(rse_gene_ntsr1))
rse_gene_ip = SummarizedExperiment(assays = list(counts = counts), 
	colData = rbind(colData(rse_gene_ip), colData(rse_gene_ntsr1)),
		rowRanges = rowRanges(rse_gene_ip))

rse_gene_ip$Oxt = ifelse(rse_gene_ip$CellType == "Oxt", 1, 0)

##### DGE ######
mod_ip = model.matrix(~Oxt + totalAssignedGene, data=colData(rse_gene_ip))
dge_ip = DGEList(counts = assays(rse_gene_ip)$counts, 
	genes = rowData(rse_gene_ip))
dge_ip = calcNormFactors(dge_ip)
vGene_ip = voom(dge_ip,mod_ip,plot=TRUE)
fitGene_ip = lmFit(vGene_ip)
ebGene_ip = ebayes(fitGene_ip)
eBGene_ip = eBayes(fitGene_ip)
outGene_ip = topTable(eBGene_ip,coef=2,
	p.value = 1,number=nrow(rse_gene_ip))
	
## add to stats
sigGene_ip = outGene_ip[outStats$gencodeID,]
outStats$OxtEnr_logFC = sigGene_ip$logFC
outStats$OxtEnr_t = sigGene_ip$t
outStats$OxtEnr_pval = sigGene_ip$P.Value

## write out
write.csv(outStats, file= "tables/TRAP_stats_allTests.csv", row.names=FALSE)

## oxt specificity
outGene_ip = outGene_ip[,c(colnames(outStats)[1:8], "EntrezID")]
outGene_ip = outGene_ip[order(outGene_ip$P.Value),]

write.csv(outGene_ip, file = "tables/Oxt_Specificity_geneLevel.csv",
	row.names=FALSE)

## do GO
sigGene_ip = outGene_ip[outGene_ip$adj.P.Val < 0.01,]	
sigGeneList_ip = split(sigGene_ip$EntrezID,
	sign(sigGene_ip$logFC))
sigGeneList_ip = lapply(sigGeneList_ip, 
			function(x) x[!is.na(x)])

goBP_ip <- compareCluster(sigGeneList_ip, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)
goMF_ip <- compareCluster(sigGeneList_ip, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)		
goCC_ip <- compareCluster(sigGeneList_ip, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)

## write out				
goBP_DF_ip = as.data.frame(goBP_ip)
goMF_DF_ip = as.data.frame(goMF_ip)
goCC_DF_ip = as.data.frame(goCC_ip)

goOut_ip = rbind(goBP_DF_ip, goMF_DF_ip, goCC_DF_ip)
goOut_ip$Ontology = rep(c("BP", "MF", "CC"), 
	times = c(nrow(goBP_ip), nrow(goMF_ip), nrow(goCC_ip)))
goOut_ip = goOut_ip[goOut_ip$p.adjust < 0.01,]
colnames(goOut_ip)[1] = "Direction"
goOut_ip = goOut_ip[order(goOut_ip$p.adjust),]
write.csv(goOut_ip, file = "tables/trap_GO_analysis_OxtEnr_FDR01.csv",row.names=FALSE)


### how many?
otherCells = outStats[outStats$P.Value_CST < 0.01 & outStats$P.Value_Ntsr1 < 0.01,]
split(otherCells$Symbol, sign(otherCells$logFC))

otherCells = outStats[which(outStats$P.Value_CST < 0.01 & 
	outStats$P.Value_Ntsr1 < 0.01 & 
	sign(outStats$logFC) == sign(outStats$logFC_CST) & 
	sign(outStats$logFC) == sign(outStats$logFC_Ntsr1)),]
	
split(otherCells$Symbol, sign(otherCells$logFC))