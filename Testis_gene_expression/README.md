1. The inputs for Feature_count_on_bams.Rmd are sorted bam files after alignment.

Below are results from FeatureCount and save out as the rds format:
Cattle: counts_Cattle_UOA_Wagyu_1_Y.rda
Human: counts_Human_T2T-CHM13v2.0.rda
Chimpanzee: counts_Chimpanzee_T2Tv2.0.rda
Gorilla:counts_PRJEB48233_Gorilla.T2Tv2.0.rda
Bornean orangutan:counts_B_orangutan_T2Tv2.0.gtf.rda

Below are annotations downloaded from NCBI and saved as rda format in R.
Cattle: annotation_Cattle_UOA_Wagyu_1_Y.rda
Human: annotation_Human_T2T-CHM13v2.0.rda
Chimpanzee: annotation_Chimpanzee_T2Tv2.0.rda
Gorilla:annotation_PRJEB48233_Gorilla.T2Tv2.0.rda
Bornean orangutan:annotation_B_orangutan_T2Tv2.0.gtf.rda

2. The Conservatively_expressed_genes_testis.Rmd was run after script Feature_count_on_bams.Rmd.
The inputs are given in the folder of inputs.
However, as the title shown they also can be found as the sub tables of the paper.