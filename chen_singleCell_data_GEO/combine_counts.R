
library(jaffelab)
library(parallel)

# read in data
fn = list.files("geo_data", full.names=TRUE)
names(fn) = ss(list.files("geo_data"), "_")
datList = mclapply(fn, read.delim, row.names=1,mc.cores=6)

## fix column names
for(i in seq(along=datList)) colnames(datList[[i]]) = paste0(
	colnames(datList[[i]]), "_", names(datList)[i])

## get unique genes
g = sort(unique(unlist(sapply(datList,rownames))))
sampleNames = do.call("c", sapply(datList, colnames))
geneCounts = matrix(0, nr = length(g), nc = length(sampleNames),
	dimnames = list(g, sampleNames))

for(i in seq(along=datList)) {
	cat(".")
	geneCounts[rownames(datList[[i]]),colnames(datList[[i]])] = as.matrix(datList[[i]])
}

## load Oxt stats
load("../TRAP_discovery_differential_expression.rda")
g = unique(outGene_ovation$Symbol[outGene_ovation$Symbol %in% rownames(geneCounts)])
outGene_ovation = outGene_ovation[match(g, outGene_ovation$Symbol),]
geneCounts2 = geneCounts[g,]

## check Oxt
yGene = log2(geneCounts2+1)

## check out different expression levels
library(genefilter)
minExprsCutoffs = c(2,4,6,8,10,15,20,25,50)
ttList = lapply(minExprsCutoffs, function(x) {
	cat(".")
	group = rep(1, ncol(yGene))
	ind = which(geneCounts2["Oxt",] >= x)
	group[ind] = 0 # sign flip
	tt = rowttests(yGene, factor(group))
	tt = as.data.frame(tt)
	tt$Oxt = rowMeans(yGene[,group==0])
	return(tt)
})

## correlation of enrichment
ttStats = sapply(ttList, "[[", "statistic")
dmStats = sapply(ttList, "[[", "dm")
rownames(ttStats) =rownames(dmStats) = rownames(yGene)
## expressed?
expIndex = which(rowSums(geneCounts2) > 1000)
cor(outGene_ovation$t[expIndex], ttStats[expIndex,])
cor(outGene_ovation$logFC[expIndex], dmStats[expIndex,])

plot(outGene_ovation$t[expIndex], ttStats[expIndex,ncol(ttStats)] ,
	xlim = c(-10,10), ylim = c(-10,10))