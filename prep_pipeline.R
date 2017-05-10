##

library(jaffelab)
library(readxl)

# read in phenotype
pdList = lapply(1:3, read_excel, 
	path = "kristen_analysis_manifests.xlsx")
names(pdList) = c("SoLo_Oxt","Ovation_Oxt", "Clonetech_CST")

## add md5
for(i in seq(along=pdList)) pdList[[i]]$md5=0

#### SoLo Oxt ####
dir.create("preprocessed_data/SoLo_Oxt")
man_SoLo = pdList$SoLo_Oxt[,c("singleRead","md5","SampleID")] 
write.table(man_SoLo, "preprocessed_data/SoLo_Oxt/samples.manifest",
	sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#### Ovation Oxt ####
dir.create("preprocessed_data/Ovation_Oxt")
man_Ovation = pdList$Ovation_Oxt[,c("leftRead", 
	"md5", "rightRead","md5", "SampleID")] 
write.table(man_Ovation, "preprocessed_data/Ovation_Oxt/samples.manifest",
	sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#### Ovation Oxt ####

dir.create("preprocessed_data/Clonetech_CST")
fqFiles = list.files("/dcl01/lieber/ajaffe/Keri/trapPilot/FASTQ",
	recur = TRUE, pattern = ".fastq.gz$", full.names=TRUE)
names(fqFiles) = ss(fqFiles, "/", 8)
fqList = split(fqFiles, names(fqFiles))
fqList = lapply(fqList, function(x) {
	xList = split(x, ss(x, "_", 9))
	do.call("cbind", xList)
})
fqMat = do.call("rbind", fqList)
man_cst = data.frame(leftRead = fqMat[,1], leftMd5 = 0,
	rightRead = fqMat[,2], rightMd5 = 0, 
	SampleID = rownames(fqMat))
write.table(man_cst, "preprocessed_data/Clonetech_CST/samples.manifest",
	sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)