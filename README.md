## BDNF-TrkB signaling in oxytocin neurons contributes to maternal behavior

### Authors

Kristen R. Maynard, PhD, John W. Hobbs, BA, Badoi N. Phan, BS, Amolika Gupta, Sumita Rajpurohit, Courtney Williams, BA, Nina Rajpurohit, BA, Joo Heon Shin, PhD, Andrew E. Jaffe, PhD, Keri Martinowich, PhD

### Abstract

Brain-derived neurotrophic factor (Bdnf) transcription is controlled by several promoters, which drive expression of multiple transcripts encoding an identical protein. We previously reported that BDNF derived from promoters I and II is highly expressed in hypothalamus and is critical for regulating aggression in male mice. Here we report that BDNF loss from these promoters causes reduced sexual receptivity and impaired maternal care in female mice, which is concomitant with decreased oxytocin (Oxt) expression during development. We identify a novellink between BDNF signaling, oxytocin, and maternal behavior by demonstrating that ablation of TrkB selectively in OXT neurons partially recapitulates maternal care impairments observed in BDNF-deficient females. Using translating ribosome affinity purification and RNA-sequencingwe define a molecular profile for OXT neurons and delineate how BDNF signaling impacts gene pathways critical for structural and functional plasticity. Our findings highlight BDNF as a modulator of sexually-dimorphic hypothalamic circuits that govern female-typical behaviors.

### Data and Code Availability

The [`Oxt_trap_seq` GitHub repository](https://github.com/LieberInstitute/oxt_trap_seq/edit/master/README.md) contains code and summarized data for RNA sequencing of Oxt neurons in the mouse. 

Raw gene counts, as a [RangedSummarizedExperiment object](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) can be downloaded from [here](rseObjs_oxtMerge_n18_4features.Rdata)

Code for various analyses described in the paper can be found below:
- [Comparing IP to Input samples](analyze_data_trap.R)
- [Comparing Bdnf-e1 mutant IP to control IP](analyze_data_bdnf.R)

The raw sequencing reads will be made available through [SRA](https://www.ncbi.nlm.nih.gov/sra) shortly

