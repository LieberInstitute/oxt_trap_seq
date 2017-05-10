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
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Clonetech_CST/rawCounts_Clonetech_CST_HiSeq_n6.rda")

erccTPM_clone = erccTPM
exonCounts_clone= exonCounts
geneCounts_clone = geneCounts
jCounts_clone = as.matrix(as.data.frame(jCounts))
jMap_clone = jMap
metrics_clone = metrics

###########3
# solo
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/SoLo_Oxt/rawCounts_SoLo_Oxt_MiSeq_n6.rda")

erccTPM_solo = erccTPM
exonCounts_solo= exonCounts
geneCounts_solo = geneCounts
jCounts_solo = as.matrix(as.data.frame(jCounts))
jMap_solo = jMap
metrics_solo = metrics

###############
# ovation ######
load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Ovation_Oxt/rawCounts_Ovation_Oxt_HiSeq_n6.rda")

erccTPM_ovation = erccTPM
exonCounts_ovation= exonCounts
geneCounts_ovation = geneCounts
jCounts_ovation = as.matrix(as.data.frame(jCounts))
jMap_ovation = jMap
metrics_ovation = metrics

########################3
### combine ########
metrics = rbind(metrics_solo, metrics_ovation[,-c(2,4,9)],
	metrics_clone[,-c(2,4,9)])
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

## junction
u_jxn = unique(c(jC
jCounts = cbind(jCounts_solo, jCounts_ovation,jCounts_clone)
jCounts = jCounts[,pd$SampleID]