
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
geneCounts = do.call("cbind", datList)