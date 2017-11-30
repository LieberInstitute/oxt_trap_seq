## BDNF-TrkB signaling in oxytocin neurons contributes to maternal behavior

### Authors

Kristen R. Maynard, PhD, John W. Hobbs, BA, Badoi N. Phan, BS, Amolika Gupta, Sumita Rajpurohit, Courtney Williams, BA, Nina Rajpurohit, BA, Joo Heon Shin, PhD, Andrew E. Jaffe, PhD, Keri Martinowich, PhD

### Abstract

Brain-derived neurotrophic factor (BDNF) and its receptor tropomyosin receptor kinase B (TrkB) play a critical role in male-typical social behavior. However, it is unknown whether BDNF-TrkB signaling contributes to female-typical social behavior.  Furthermore, the sexually-dimorphic neural circuits mediating effects of BDNF on sex-specific social behaviors remain unexplored.  Transcription of Bdnf is controlled by several promoters, which drive expression of multiple transcripts encoding an identical protein. We previously reported that BDNF derived from promoters I and II is highly expressed in the hypothalamus and is critical for the regulation of aggression in male mice. Here we report that loss of promoter I or II-derived BDNF in Bdnf-e1 -/- and Bdnf-e2 -/- female mice causes reduced sexual receptivity and impaired maternal care.  We show that disruption of BDNF generated from promoters I and II, but not IV and VI, leads to decreased oxytocin (Oxt) gene expression during development.   We further demonstrate that ablation of TrkB selectively in OXT neurons recapitulates impairments in maternal care observed in BDNF-deficient females.  To determine how BDNF signaling impacts OXT neuron function, we use translating ribosome affinity purification (TRAP) in combination with RNA-sequencing to define a molecular profile for OXT-expressing neurons and identify how perturbations in BDNF signaling impact gene pathways critical for structural and functional plasticity.  These findings highlight BDNF as a key molecular player in modulating sexually-dimorphic hypothalamic circuits that govern complex female-typical behaviors, including sexual receptivity and maternal care.     

### Data and Code Availability

The [`Oxt_trap_seq` GitHub repository](https://github.com/LieberInstitute/oxt_trap_seq/edit/master/README.md) contains code and summarized data for RNA sequencing of Oxt neurons in the mouse. 

Raw gene counts, as a [RangedSummarizedExperiment object](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) can be downloaded from [here](rseObjs_oxtMerge_n18_4features.Rdata)

Code for various analyses described in the paper can be found below:
- [Comparing IP to Input samples](analyze_data_trap.R)
- [Comparing Bdnf hets to wild type samples](analyze_data_bdnf.R)

The raw sequencing reads will be made available through [SRA](https://www.ncbi.nlm.nih.gov/sra) shortly

