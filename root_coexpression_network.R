install.packages("BiocManager")
BiocManager::install("WGCNA")

###################################### Clean Up Data ##################################################

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/steph/OneDrive/Documents/Lab/Experiments/Coexpression Networks/Coexpression_Networks";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the root data set
rootData = read.csv("Stelpflug_RNAseq_root2.csv") #testing with subset of genes where >12 samples have expression data
######### For exercise choose ~12k genes ###############
#library(dplyr)
#rootData2 <- sample_n(rootData, 10000, replace = TRUE)


# Take a quick look at what is in the data set:
dim(rootData);
names(rootData);
# remove auxiliary data and transpose the expression data
datExpr0 = as.data.frame(t(rootData[, -c(1:5)]));
names(datExpr0) = rootData$gene;
rownames(datExpr0) = names(rootData)[-c(1:5)];

# check for genes with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# This section only necessary if not all genes have passed the cut
# To remove the offending genes and samples
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  #if (sum(!gsg$goodSamples)>0) 
    #printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}



# cluster the samples (not genes) to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#outliers can be removed by hand or automatic approach using a height cut
## cortical parenchyma_3d looks like an outlier
# Plot a line to show the cut
abline(h = 60000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 65000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

############################ Building the Co-Expression Network Step by step ####################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, corOptions = list(use = 'p', method = 'pearson'))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# calculate adjacencies
softPower = 6;
adjacency = adjacency(datExpr, power = softPower)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function using TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


# branch cutting
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# plot module assignments under gene dendrogram
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
module_assign <- cbind(allgenes, dynamicMods, dynamicColors) # first draft of assignments prior to merge


# merge modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# choose height of 0.25 whic corresponds to correlation of 0.75 to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
# plot dendrogram again with original and merged modules
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
module_assign <- cbind(allgenes, dynamicMods, dynamicColors, mergedColors) # new to include new module colors

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs


#################################### Relate modules to phenotypes ############################################################

# load trait data and match the samples for which they were measured to the expression samples
traitData = read.csv("root_sample_info.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
#allTraits = traitData[, -c(31, 16)];
#allTraits = allTraits[, c(2, 11:36) ];
allTraits <- traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
rootSamples = rownames(datExpr)
traitRows = match(rootSamples, allTraits$Tissue)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage();


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#visualize correlations
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



############################ Interfacing with GO and Functional Annotations #################################################

# Read in the probe annotation
annot = read.csv(file = "Zm_B73_5b_annotation_matrix_all_models_withFuncAnnot.csv");
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$Maize.Locus.Identifier)
# Get the corresponding Locus Link IDs
allgenes = annot$Maize.Locus.Identifier[probes2annot];
# $ Choose interesting modules
intModules = c("steelblue", "lightgreen", "turquoise", "salmon")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modgenes = allgenes[modGenes];
  # Write them into a file
  fileName = paste("Maize.Locus.Identifier-", module, ".txt", sep="");
  write.table(as.data.frame(allgenes), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("Maize.Locus.Identifier-all.txt", sep="");
write.table(as.data.frame(allgenes), file = fileName,
            row.names = FALSE, col.names = FALSE)

# enrichment analysis
GO <- read.csv("LawrenceDill_Zea_mays_GO_AGPV3.0.csv")
# filter all genes set for only intModules
int_modules0 <- as.data.frame(module_assign) %>% 
  filter(mergedColors %in% intModules) %>% 
  select(allgenes, mergedColors)
colnames(int_modules0) <- c("maize_gene", "moduleColor")
# merge GO terms with modules of interest
int_modules <- left_join(int_modules0, GO, by = c("maize_gene" = "db_object_id"))

# look for top represented GO terms of each module
GO_overrep <- int_modules %>% 
  select(moduleColor, term_accession) %>% 
  group_by(moduleColor, term_accession) %>%
  summarise(count = n())
# top 10 GO term expressions in each module
GO_overrep_top10 <- GO_overrep %>%
  arrange(desc(count)) %>% 
  group_by(moduleColor) %>%
  top_n(10, count) %>% 
  arrange(moduleColor)
#write table - can look up GO term annotations online until spreadsheet is made
write.table(GO_overrep_top10, "intModules_over-repGO.txt", row.names = F)



################################################# Visualizing the Network ##############################################################


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


# reduce number of genes to 400 to increase plotting speed
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

# visualize network of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
stele = as.data.frame(datTraits$Stele.relevant);
names(stele) = "stele"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, stele))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


# identify groups of correlated eigengenes
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)