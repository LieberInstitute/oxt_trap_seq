##
library(WriteXLS)
library(readxl)

### tabled
s1 = read.csv("tables/TRAP_discovery_differential_expression.csv.gz", as.is=TRUE)

s2 = read.csv("tables/TRAP_stats_allTests.csv", as.is=TRUE)
s2 = s2[,c(1:8, 14:16, 9:13, 17:ncol(s2))]
s2$bed_mm10 = NULL

s3 = read.csv("tables/trap_GO_analysis_DE_FDR01.csv", as.is=TRUE)
s3$Direction = ifelse(s3$Direction == 1, "OXT_enrich", "OXT_deplete")

s4 = as.data.frame(read_excel("tables/Oxt_SFARI_asd_gene_list.xls", sheet=1))
s5 = as.data.frame(read_excel("tables/Oxt_SFARI_asd_gene_list.xls", sheet=2))

s6 = read.csv("tables/bdnf_genotype_effect_allGenes.csv", as.is=TRUE)
s6$sigColor = NULL
s6$FDR_sig = s6$adj.P.Val < 0.1

s7 = read.csv("tables/GO_bdnf_enrichment_DE_at_p005.csv", as.is=TRUE)
s7$Direction = ifelse(s7$Direction == 1, "Bdnf_e1>Control", "Bdnf_e1<Control")

s8 = as.data.frame(read_excel("tables/BdnfEffect_Oxt_SFARI_asd_gene_list.xls",sheet=2))

WriteXLS(paste0("s", 1:8), ExcelFileName = "maynard_suppTables_OxtPaper.xlsx",
	SheetNames = paste0("Table S", 1:8), BoldHeaderRow=TRUE)