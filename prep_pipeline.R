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
	
####### 
## BDNF experiment
dir.create("preprocessed_data/Bdnf_SoLo_Oxt")

## latest batch
fqFiles1 = list.files("/dcl01/ajaffe/data/Nina/Keri/SoLo_MiSeq_052617",
	recur = TRUE, pattern = ".fastq.gz$", full.names=TRUE)
fqFiles1 = fqFiles1[grep("_R1_", fqFiles1)]
names(fqFiles1) = ss(ss(fqFiles1, "/", 8), "_")
fqFiles1 = fqFiles1[c("J-IP", "L-IP", "N-IP", "O-IP")]

## previous data
fqFiles2 = list.files("/dcl01/ajaffe/data/Nina/Keri/SoLo_MiSeq_041117",
	recur = TRUE, pattern = ".fastq.gz$", full.names=TRUE)
fqFiles2 = fqFiles2[grep("_R1_", fqFiles2)]
names(fqFiles2) = ss(ss(fqFiles2, "/", 8), "_")
fqFiles2 = fqFiles2[c("IP-B", "IP-I")]
names(fqFiles2) = paste0(ss(names(fqFiles2), "-", 2), "-IP")

## previous data
fqFiles3 = list.files("/dcl01/ajaffe/data/Nina/Keri/SoLo_MiSeq_032217",
	recur = TRUE, pattern = ".fastq.gz$", full.names=TRUE)
fqFiles3 = fqFiles3[grep("_R1_", fqFiles3)]
names(fqFiles3) = ss(ss(fqFiles3, "/", 8), "_")
fqFiles3 = fqFiles3[c("IP-B", "IP-I")]
names(fqFiles3) = paste0(ss(names(fqFiles3), "-", 2), "-IP")

## combine
fqFiles = c(fqFiles1,fqFiles2,fqFiles3)

man_bdnf = data.frame(read = fqFiles, md5 = 0,
	SampleID = names(fqFiles))
write.table(man_bdnf, "preprocessed_data/Bdnf_SoLo_Oxt/samples.manifest",
	sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	