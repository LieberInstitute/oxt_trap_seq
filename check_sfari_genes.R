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
sigGene = outGene_ovation[outGene_ovation$adj.P.value < 0.01,]

########################
# load SFARI human genes
humanSFARI = read.csv('SFARI-Gene_genes_export03-11-2017.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(-grepl('S',gene.score),ss(gene.score,'S')),])
humanSFARI = cbind(humanSFARI,outGene_ovation[lookfor(humanSFARI$gene.symbol,outGene_ovation$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 842 expressed in mouse oxt dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('SFARI-Gene_animal-genes_export03-11-2017.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene_ovation[lookfor(mouseSFARI$gene.symbol,outGene_ovation$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 224 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
# 24 genes in Oxt expressed and in SFARI 
outGene_ovation$inMouseSFARI = outGene_ovation$Symbol %in% mouseSFARI$Symbol
(t1 = with(outGene_ovation,table(inMouseSFARI,inDEG = adj.P.Val < 0.01 & logFC >0)))
fisher.test(t1) # OR = 4.244335 p-value = 6.617e-10
ind1 = with(sigGene, which(mouseSFARI$Symbol %in% Symbol))

#######################################
# list of DEG in scored human SFARI list
# 131 human SFARI ASD-linked genes in our list of DEGs
outGene_ovation$inHumanSFARI =  outGene_ovation$Symbol %in% humanSFARI$Symbol
outGene_ovation_hs = outGene_ovation[grep("^ENSG", outGene_ovation$hsapien_homolog),]
(t2 = with(outGene_ovation_hs,table(inHumanSFARI,inDEG = adj.P.Val < 0.01& logFC >0)))
fisher.test(t2) # OR = 4.00, p-value <2.2e-16
ind2 = with(sigGene, which(humanSFARI$Symbol %in% Symbol))

## venn diagrams
library(limma)
vennCount = vennCounts(data.frame(OxtEnr = outGene_ovation$adj.P.Val < 0.01 & outGene_ovation$logFC>0,
	`Hs SFARI` = outGene_ovation$inHumanSFARI,
	`Mm SFARI`= outGene_ovation$inMouseSFARI))
vennCounts_hs = with(outGene_ovation[grep("^ENSG", outGene_ovation$hsapien_homolog),],
	vennCounts(data.frame(OxtEnr = adj.P.Val < 0.01 & logFC>0,
	`Hs SFARI` = inHumanSFARI,
	`Mm SFARI`= inMouseSFARI)))

pdf("plots/vennDiagram_asd_oxt_enrichment.pdf")
vennDiagram(vennCount,cex=1.4,main="Mouse Expressed")
vennDiagram(vennCounts_hs,cex=1.4, main = "Mouse Exprs (Hom)")
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
sigGene = bdnfStats[which(bdnfStats$adj.P.Val < 0.1),]

########################
# load SFARI human genes
humanSFARI = read.csv('SFARI-Gene_genes_export03-11-2017.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(-grepl('S',gene.score),ss(gene.score,'S')),])
humanSFARI = cbind(humanSFARI,bdnfStats[lookfor(humanSFARI$gene.symbol,bdnfStats$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 843 expressed in mouse bdnf dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('SFARI-Gene_animal-genes_export03-11-2017.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,bdnfStats[lookfor(mouseSFARI$gene.symbol,bdnfStats$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 227 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
# 5 genes differentially expressed and in SFARI 
(t1 = with(bdnfStats,table(inMouseSFARI = Symbol %in% mouseSFARI$Symbol,inDEG = adj.P.Val < 0.1)))
fisher.test(t1) # OR = 5.2, p-value = 0.0036
ind1 = with(sigGene, which(mouseSFARI$Symbol %in% Symbol))

#######################################
# list of DEG in scored human SFARI list
# 131 human SFARI ASD-linked genes in our list of DEGs
bdnfStats_hs = bdnfStats[grepl('ENSG',bdnfStats$hsapien_homolog),] #take only genes w/ human homolog

(t2 = with(bdnfStats_hs,table(inHumanSFARI = Symbol %in% humanSFARI$Symbol,inDEG = adj.P.Val < 0.1)))
fisher.test(t2) # OR = 2.05, p-value = 0.018
ind2 = with(sigGene, which(humanSFARI$Symbol %in% Symbol))

####################
# save the two lists
sfariList = list(Scored_Human_SFARI_Genes = humanSFARI[ind2,],
                 Mouse_Model_Genes = mouseSFARI[ind1,])
WriteXLS(sfariList,ExcelFileName = 'tables/BdnfEffect_Oxt_SFARI_asd_gene_list.xls')

