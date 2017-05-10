##
library(jaffelab)
library(SummarizedExperiment)
library(readxl)
library(stringr)

# read in phenotype
pdList = lapply(1:3, read_excel, 
	path = "kristen_analysis_manifests.xlsx")
names(pdList) = c("SoLo_Oxt","Ovation_Oxt", "Clonetech_CST")
pdList = lapply(pdList, function(x) x[,1:4])

## merge
pd = as.data.frame(do.call("rbind", pdList))
pd$Kit = str_trim(pd$Kit)

#####################
## load count data ##

################
# clonetech ## 
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Clonetech_CST/rpkmCounts_Clonetech_CST_HiSeq_n6.rda")
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Clonetech_CST/rawCounts_Clonetech_CST_HiSeq_n6.rda")

txTpm_clone = txTpm
exonCounts_clone= exonCounts
geneCounts_clone = geneCounts
jCounts_clone = as.matrix(as.data.frame(jCounts))
jMap_clone = jMap
metrics_clone = metrics

###########3
# solo
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/SoLo_Oxt/rpkmCounts_SoLo_Oxt_MiSeq_n6.rda")
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/SoLo_Oxt/rawCounts_SoLo_Oxt_MiSeq_n6.rda")

txTpm_solo = txTpm
exonCounts_solo= exonCounts
geneCounts_solo = geneCounts
jCounts_solo = as.matrix(as.data.frame(jCounts))
jMap_solo = jMap
metrics_solo = metrics

###############
# ovation ######
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Ovation_Oxt/rawCounts_Ovation_Oxt_HiSeq_n6.rda")
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Ovation_Oxt/rpkmCounts_Ovation_Oxt_HiSeq_n6.rda")

txTpm_ovation = txTpm
exonCounts_ovation= exonCounts
geneCounts_ovation = geneCounts
jCounts_ovation = as.matrix(as.data.frame(jCounts))
jMap_ovation = jMap
metrics_ovation = metrics

########################3
### combine ########
metrics = rbind(metrics_solo, 
	metrics_ovation[,names(metrics_solo)],
	metrics_clone[,names(metrics_solo)])
pd = cbind(pd, metrics[pd$SampleID,-1])
pd$ReadType = rep(c("Single","Paired"), c(6,12))

## gene
geneCounts = cbind(geneCounts_solo, geneCounts_ovation,
	geneCounts_clone)
geneCounts = geneCounts[,pd$SampleID]

## exon
exonCounts = cbind(exonCounts_solo, exonCounts_ovation,
	exonCounts_clone)
exonCounts = exonCounts[,pd$SampleID]

## transcripts
txTpm = cbind(txTpm_solo, txTpm_ovation, txTpm_clone)

## junction
jMap_clone$meanExprs = jMap_ovation$meanExprs = jMap_solo$meanExprs = NULL 
jMap = sort(unique(c(jMap_solo, jMap_clone, jMap_ovation)))

jCounts = cbind(jCounts_solo[match(names(jMap),names(jMap_solo)),],
	jCounts_ovation[match(names(jMap),names(jMap_ovation)),],
	jCounts_clone[match(names(jMap),names(jMap_clone)),])
rownames(jCounts) = names(jMap)
jCounts[is.na(jCounts)] = 0
colnames(jCounts) = gsub(".", "-", colnames(jCounts), fixed=TRUE)

### match order ####
geneCounts = geneCounts[,rownames(metrics)]
exonCounts = exonCounts[,rownames(metrics)]
jCounts = jCounts[,rownames(metrics)]
txTpm = txTpm[,rownames(metrics)]

#########################################
####### export merged rse objects #######
## Create gene,exon RangedSummarizedExperiment objects

## gene
gr_genes <- GRanges(seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_gene <- SummarizedExperiment(assays = list('counts' = geneCounts),
    rowRanges = gr_genes, colData = metrics)

## exon
gr_exons <- GRanges(seqnames = exonMap$Chr,
    IRanges(exonMap$Start, exonMap$End), strand = exonMap$Strand)
names(gr_exons) <- rownames(exonMap)
mcols(gr_exons) <- DataFrame(exonMap[, - which(colnames(exonMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_exon <- SummarizedExperiment(assays = list('counts' = exonCounts),
    rowRanges = gr_exons, colData = metrics)

## jxn
jIndex = which(rowSums(jCounts > 0) > 2)
rse_jx <- SummarizedExperiment(
	assays = list('counts' = jCounts[jIndex,]),
    rowRanges = jMap[jIndex,], colData = metrics)

## tx
txMap_coord = geneMap[txMap$gencodeID,]
gr_txs <- GRanges(seqnames = txMap_coord$Chr,
    IRanges(txMap_coord$Start, txMap_coord$End), 
	strand = txMap_coord$Strand)
names(gr_txs) <- rownames(txMap)
mcols(gr_txs) <- DataFrame(txMap)

rse_tx <- SummarizedExperiment(assays = list('tpm' = txTpm),
    rowRanges = gr_txs, colData = metrics)

save(rse_gene, rse_exon, rse_jx, rse_tx,
	file = "rseObjs_oxtMerge_n18_4features.Rdata")
