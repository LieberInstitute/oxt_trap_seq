##
library(WriteXLS)
library(readxl)

### tabled
s1 = read.csv("tables/TRAP_discovery_differential_expression.csv.gz", as.is=TRUE)

# DE stats for Oxt IP vs Input (FDR controlled) with other datasets
s2 = read.csv("tables/TRAP_stats_allTests.csv", as.is=TRUE)
s2 = s2[,c(1:8, 14:16, 9:13, 17:ncol(s2))]
s2$bed_mm10 = NULL

## GO for Oxt: IP vs input 
s3 = read.csv("tables/trap_GO_analysis_DE_FDR01.csv", as.is=TRUE)
s3$Direction = ifelse(s3$Direction == 1, "OXT_enrich", "OXT_deplete")

## DE Stats for Oxt IP vs CST IP & Ntsr1 IP
s4 = read.csv("tables/Oxt_Specificity_geneLevel.csv", as.is=TRUE)

## GO for Oxt IP vs CST IP & Ntsr1 IP
s5 = read.csv("tables/trap_GO_analysis_OxtEnr_FDR01.csv", as.is=TRUE)

s6 = as.data.frame(read_excel("tables/Oxt_SFARI_asd_gene_list.xls", sheet=2))
s7 = as.data.frame(read_excel("tables/Oxt_SFARI_asd_gene_list.xls", sheet=1))

s8 = read.csv("tables/bdnf_genotype_effect_allGenes.csv", as.is=TRUE)
s8$sigColor = NULL
s8$FDR_sig = s8$adj.P.Val < 0.1

s9 = read.csv("tables/GO_bdnf_enrichment_DE_at_p005.csv", as.is=TRUE)
s9$Direction = ifelse(s9$Direction == 1, "Bdnf_e1>Control", "Bdnf_e1<Control")

s10 = as.data.frame(read_excel("tables/BdnfEffect_Oxt_SFARI_asd_gene_list.xls",sheet=2))

WriteXLS(paste0("s", 1:10), ExcelFileName = "maynard_suppTables_OxtPaper.xlsx",
	SheetNames = paste0("Table S", 1:10), BoldHeaderRow=TRUE)