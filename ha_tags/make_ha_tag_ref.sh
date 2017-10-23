###
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=55G,h_fsize=100G
#$ -N salmon_build
#$ -pe local 1
#$ -o build_indexes_log
#$ -e build_indexes_log
#$ -m a
echo "**** Job starts ****"
date

MAINDIR=/dcl01/lieber/ajaffe/lab/oxt_trap_seq/ha_tags
fasta=/dcl01/lieber/ajaffe/lab/oxt_trap_seq/ha_tags/ha.fa

### hg38

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.8.2_linux_x86_64/bin/salmon index \
-t ${fasta} -i ${MAINDIR}/salmon_0.8.2_index_ha_tags -p 1 --type quasi -k 15

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon index \
-t ${fasta} -i ${MAINDIR}/salmon_0.7.2_index_ha_tags -p 1 --type quasi -k 15

echo "**** Job ends ****"
date
