# xylemGWAS
This is a repository of the data and code for Klein et al. (2024) "Integrating GWAS with a gene co-expression network better prioritizes candidate genes associated with root metaxylem phenes in maize." In brief, we used GWAS and a gene co-expression network to identify candidate genes associated with root metaxylem phenes using the Wisconsin Diversity Panel. 

File summary:
GAPITcode.R - R script for running GWAS using a mixed linear model in GAPIT. The SNP panel used to conduct studies is from Mazaheri et al. (2019).

root_coexpression_network.R - R script for constructing a gene co-expression network using WGCNA. Transcriptome data is the maize gene expression atlas (Stelpflug et al. 2016) available at www.maizegdb.com

WiDiv_root_phenotype_data_2015 - Root metaxylem phenotypes observed in the field under well-watered and water stress conditions in the Wisconsin Diversity Panel in 2015. This set served as the "corroboration set."

WiDiv_root_phenotype_data_2016 - Root metaxylem phenotypes observed in the field under well-watered and water stress conditions in the Wisconsin Diversity Panel in 2016. This set served as the "pilot set."
