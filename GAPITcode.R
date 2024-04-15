

######################################################################
########  GWAS  - WIDIV  REGULAR MIXED LINEAR MODEL Q + K ############
######################################################################

## Bioconductor
library(BiocManager)
BiocManager::install("multtest")
install.packages("gplots")
install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d

library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt") #EMMAx
source("http://www.zzlab.net/GAPIT/gapit_functions.txt") # GAPIT

### Installing required packages
install.packages("LDheatmap")
install.packages("genetics")

### Loading the required packages

library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library(data.table)
library(bigmemory)
library(biganalytics)

### Reading the phenotype data

myY  <- read.csv("widiv16_MXVphenes_gwas.csv", head = TRUE)


###### Reading SNP, Kinship and PC for GWAS analysis #####

myG <- fread("widiv_942g_899784SNPs_imputed_filteredGenos_noRTA_AGPv4.hmp.txt", header=F, sep="\t", stringsAsFactors=F, na.strings=c("NA", "N"))

#myKI <- read.csv("GAPIT.Kin.VanRaden.csv", head=FALSE)

#myCV <- read.csv("GAPIT.PCA.csv", head=TRUE)

#myCV <- read.table("C:/Users/German/Documents/WIDIV_GWAS/RNA-seq/GWAS_SNPs_Q_K_Matrices/GAPIT.Qmatrix.10.PC.txt", head = TRUE)\

#Running GAPIT
## Root traits code for all traits using Loiselle kinship
myGAPIT <- GAPIT(Y=myY, G=myG, group.from=413, group.to=413, PCA.total=3, SNP.impute="Major", Model.selection=TRUE)



