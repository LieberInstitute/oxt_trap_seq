## Check Oxt dataset DEGs with human asd and animal model genes from SFARI

## via badoi
lookfor= function(this,inThat) { #gene name matching
  this = toupper(this); inThat = toupper(inThat);
  tmp = sapply(this,function(x) grep(paste0('^',x,'$'),inThat))
  return(sapply(tmp,function(x) ifelse(length(x)==0,NA,x[1])))}

library(jaffelab)
library(biomaRt)
library(WriteXLS)

options(stringsAsFactors = FALSE)


##################################
# get Ensembl mouse to human genes
ensembl = useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
               mart = ensembl)


##############################
# load the Oxt results ####
load('TRAP_discovery_differential_expression.rda')
#############################

outGene_ovation$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(outGene_ovation$ensemblID,MMtoHG$ensembl_gene_id)]
outGene_ovation$sigColor = NULL 

# keep only expressed genes (DESeq2 will not adjust p-value if not expressed enough)
#outGene_ovation = outGene_ovation[!is.na(outGene_ovation$padj),] #take only genes w/ human homolog
outGene_ovation = outGene_ovation[grepl('ENSG',outGene_ovation$hsapien_homolog),] #take only genes w/ human homolog
dim(outGene_ovation) #14803 background genes
sigGene = outGene_ovation[which(outGene_ovation$adj.P.Val < 0.01),]
dim(sigGene) #1545 DEGs with human homolog

########################
# load SFARI human genes
humanSFARI = read.csv('/users/bphan/oxt/PTHS_mouse/oxt_mouse/tables/gene-score.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(-grepl('S',Score),ss(Score,'S')),])
humanSFARI = cbind(humanSFARI,outGene_ovation[lookfor(humanSFARI$Gene.Symbol,outGene_ovation$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 558 expressed in mouse oxt dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('/users/bphan/oxt/PTHS_mouse/oxt_mouse/tables/sfari_genetic-animal-model-summary.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[Model.Species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene_ovation[lookfor(mouseSFARI$Model.Gene.symbol,outGene_ovation$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 176 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
# 32 genes differentially expressed and in SFARI 
outGene_ovation$inMouseSFARI = outGene_ovation$Symbol %in% mouseSFARI$Symbol
(t1 = with(outGene_ovation,table(inMouseSFARI,inDEG = adj.P.Val < 0.01)))
fisher.test(t1) # OR = 1.9, p-value = 0.0019
ind1 = with(sigGene, which(mouseSFARI$Symbol %in% Symbol))

#######################################
# list of DEG in scored human SFARI list
# 131 human SFARI ASD-linked genes in our list of DEGs
outGene_ovation$inHumanSFARI =  outGene_ovation$Symbol %in% humanSFARI$Symbol
(t2 = with(outGene_ovation,table(inHumanSFARI,inDEG = adj.P.Val < 0.01)))
fisher.test(t2) # OR = 2.24, p-value = 5.76e-12
ind2 = with(sigGene, which(humanSFARI$Symbol %in% Symbol))

## venn diagrams
library(limma)
vennCount = vennCounts(data.frame(DEG = outGene_ovation$adj.P.Val < 0.01,
	`Hs SFARI` = outGene_ovation$inHumanSFARI,
	`Mm SFARI`= outGene_ovation$inMouseSFARI))

pdf("plots/vennDiagram_asd_oxt_enrichment.pdf")
vennDiagram(vennCount,cex=1.4)
dev.off()


####################
# save the two lists
sfariList = list(Scored_Human_SFARI_Genes = humanSFARI[ind2,],
                 Mouse_Model_Genes = mouseSFARI[ind1,])
WriteXLS(sfariList,ExcelFileName = 'tables/Oxt_SFARI_asd_gene_list.xls')


############################
##### Bdnf effects #########
############################
bdnfStats = read.csv("tables/bdnf_genotype_effect_allGenes.csv",as.is=TRUE)
bdnfStats$sigColor = NULL
bdnfStats$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(bdnfStats$ensemblID,MMtoHG$ensembl_gene_id)]
bdnfStats = bdnfStats[grepl('ENSG',bdnfStats$hsapien_homolog),] #take only genes w/ human homolog
dim(bdnfStats) #14567 background genes
sigGene = bdnfStats[which(bdnfStats$adj.P.Val < 0.1),]
dim(sigGene) # 71 DEGs with human homolog


########################
# load SFARI human genes
humanSFARI = read.csv('/users/bphan/oxt/PTHS_mouse/oxt_mouse/tables/gene-score.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(-grepl('S',Score),ss(Score,'S')),])
humanSFARI = cbind(humanSFARI,bdnfStats[lookfor(humanSFARI$Gene.Symbol,bdnfStats$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 556 expressed in mouse oxt dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('/users/bphan/oxt/PTHS_mouse/oxt_mouse/tables/sfari_genetic-animal-model-summary.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[Model.Species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,bdnfStats[lookfor(mouseSFARI$Model.Gene.symbol,bdnfStats$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 179 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
# 2 genes differentially expressed and in SFARI 
(t1 = with(bdnfStats,table(inMouseSFARI = Symbol %in% mouseSFARI$Symbol,inDEG = adj.P.Val < 0.1)))
fisher.test(t1) # OR = 6.23, p-value = 0.0018
ind1 = with(sigGene, which(mouseSFARI$Symbol %in% Symbol))

#######################################
# list of DEG in scored human SFARI list
# 131 human SFARI ASD-linked genes in our list of DEGs
(t2 = with(bdnfStats,table(inHumanSFARI = Symbol %in% humanSFARI$Symbol,inDEG = adj.P.Val < 0.1)))
fisher.test(t2) # OR = 2.79, p-value = 0.018
ind2 = with(sigGene, which(humanSFARI$Symbol %in% Symbol))

####################
# save the two lists
sfariList = list(Scored_Human_SFARI_Genes = humanSFARI[ind2,],
                 Mouse_Model_Genes = mouseSFARI[ind1,])
WriteXLS(sfariList,ExcelFileName = 'tables/BdnfEffect_Oxt_SFARI_asd_gene_list.xls')

