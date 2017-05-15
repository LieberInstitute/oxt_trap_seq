###########################
## analysis of oxt data  ##
## andrew jaffe ###########

## libraries
library(jaffelab)
library(edgeR)
library(SummarizedExperiment)
library(derfinder)

## load in data
load("rseObjs_oxtMerge_n18_4features.Rdata")

## start w/ gene level, ovation
rse_gene_ovation = rse_gene[,rse_gene$
