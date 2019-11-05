#Compiled by Vy Tran, email: vtran21@jhmi.edu
#Date: 11/05/2019

#Most codes were modifed from ARCHS4 https://amp.pharm.mssm.edu/archs4/help.html (Ma'ayan et al.), and from the Weight Gene 
#Correlation Network Analysis (WGCNA) tutorial https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/ 
#(Peter Langelder and Steve Horvath)
#################################################################################################

## 0. Load libraries and data:

#Check working directory:
getwd()

#Load required packages:
if(!require("BiocManager", "plotrix")){
  install.packages("BiocManager")
}

BiocManager::install("WGCNA")
install.packages("tidyverse")
library(tidyverse)
library(WGCNA)
library(plotrix)
library(tidyverse)

# Important setting:
options(stringsAsFactors=FALSE) 

#Load RNA-seq data:
load(file = "Data Inputs/PRAD.Rfile") #Loaded as PRAD_DATA
load(file = "Data Inputs/COAD.Rfile") #Loaded as COAD_DATA
load(file = "Data Inputs/GBMLGG.Rfile") #Loaded as GBMLGG_DATA

#Load PubMed publications, unfavorable genes, and disease index data:
load(file = "Data Inputs/PubMed_ID_for_all_TCGA_genes.Rdata") #Loaded as PubMed
load(file = "Data Inputs/Prognostic_unfavorable_genes.Rdata") #Loaded as unfavorable_genes
load(file = "Data Inputs/Gleasonscore_PRAD.Rdata") #Loaded as Gleason
load(file = "Data Inputs/Aggressiveness score_COAD.Rdata")
load(file = "Data Inputs/Cindex_GBMLGG.Rdata")
#################################################################################################

## 1. Create Figure 1 for PubMed distribution of prognostic unfavorable genes:

# Create a data frame including both Pubmed publications and prognostic unfavorable genes:
colnames(unfavorable_genes) = c("genename")
Pubmed_unfavorable = merge(unfavorable_genes, PubMed, by = "genename")
colnames(Pubmed_unfavorable)

#Create Figure 1a:

#Open a jpeg file
jpeg("Plots/Figure_1.jpeg", width = 3200, height = 3200, res = 300) 

layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE),
       widths=c(1.5,1.5), heights=c(1.3,1,1))


#Plot Figure 1a:
par(cex.axis = 1.5)
density_unfavorable = density(Pubmed_unfavorable$PubMed)
gap.plot(x = density_unfavorable$x, y = density_unfavorable$y, gap = c(250, 7650), gap.axis = "x", type = "l", xlim = c(0, 7750),
         ylab = "Density", xlab = "Number of PubMed publications", lwd = 1.9, cex.lab = 1.5,
         xtics = c(0, 50, 100, 150, 300, 6500, 7700))
axis.break(1, 250, breakcol="white", style="gap") #to cover the lines created by the gap.plot fucntion
axis.break(axis = 1, 250, style = "slash") #to add a pair of slash markers on the x axis
axis.break(axis = 3, 250, style = "slash")
abline(v = 50, col = "red")
text(x = 50, y = 0.015,"cut-off",pos=4, col = "red", cex = 1.7)
title(main ="a", adj=0, line=2, font=2, cex.main = 2)


# Plot figures 1b-e:
par(cex.axis = 1)
slices_enigmatic_autosomal = c(2.5, 7.4, 90.1)
lbls = c("Autosomal dominant", "Autosomal recessive", "Other")
lbls = paste(lbls, slices_enigmatic_autosomal) # add percents to labels 
lbls = paste(lbls,"%",sep="") # add % to labels 
pie(slices_enigmatic_autosomal, radius = 0.84, labels = lbls, cex = 1.3, main = "Functionally enigmatic genes",
    cex.main = 1.5, col = c("darkgoldenrod1", "cornflowerblue", "darkgray"))
title(main ="b", adj=0, line=2, font=2, cex.main = 2)

slices_wellstudied_autosomal = c(14.7,15.0, 70.3)
lbls = c("Autosomal dominant", "Autosomal recessive", "Other")
lbls = paste(lbls, slices_wellstudied_autosomal) # add percents to labels 
lbls = paste(lbls,"%",sep="") # add % to labels 
pie(slices_wellstudied_autosomal, radius = 0.84,labels = lbls, cex = 1.3, cex.main = 1.5, main="Well-studied genes", col = c("darkgoldenrod1", "cornflowerblue", "darkgray"))
title(main ="c", adj=0, line=2, font=2, cex.main = 2)

slices_enigmatic_species = c(753,102, 3507)
lbls = c("Eukaryote", "Primate", "Other")
pct = round(slices_enigmatic_species/sum(slices_enigmatic_species)*100, digits = 1)
lbls = paste(lbls, pct) # add percents to labels 
lbls = paste(lbls,"%",sep="") # add % to labels 
pie(slices_enigmatic_species, radius = 0.84, labels = lbls, cex = 1.3, cex.main = 1.5, main = "Functionally enigmatic genes", col = c( "aquamarine4", "hotpink", "darkgray"))
title(main ="d", adj=0, line=2, font=2, cex.main = 2)

slices_wellstudied_species = c(555,19,1607)
lbls = c("Eukaryote", "Primate", "Other")
pct = round(slices_wellstudied_species/sum(slices_wellstudied_species)*100, digits = 1)
lbls = paste(lbls, pct) # add percents to labels 
lbls = paste(lbls,"%",sep="") # add % to labels 
pie(slices_wellstudied_species, radius = 0.84, cex = 1.3, cex.main = 1.5, main="Well-studied genes", labels = lbls, col = c("aquamarine4","hotpink", "darkgray"))
title(main ="e", adj=0, line=2, font=2, cex.main = 2)

dev.off()
#################################################################################################


## 2. WGCNA Analysis for prostate adenocarcinoma (PRAD) data set

#Remove uninformative data:
PRAD = PRAD_DATA
row.names(PRAD) = PRAD_DATA[,1]; colnames(PRAD) = PRAD_DATA[1,]
PRAD = PRAD[-c(1,2),-1]

#Convert data frame to numeric matrix:
PRAD1 = as.matrix(sapply(PRAD, as.numeric))  
colnames(PRAD1) = colnames(PRAD)
row.names(PRAD1) = row.names(PRAD)
class(PRAD1)
is.numeric(PRAD1)

#Transpose matrix and filter for 10,000 most variant genes:
PRADdata = t(PRAD1[order(apply(PRAD1,1,mad), decreasing = T)[1:10000],])

#Check if PRADdata have many missing values:
gsg = goodSamplesGenes(PRADdata, verbose = 3)
gsg$allOK

#The command returned "TRUE", so all genes have passed the cuts.

# Re-cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers:
sampleTree = hclust(dist(PRADdata), method = "average")
# Plot the sample tree:
sizeGrWindow(12,9)
pdf(file = "Plots/Sample clustering to detect outliers for PRAD.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers for PRAD", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()
#It appears there is no obvious outlier. The variable PRADdata now contains the expression data ready for WGCNA analysis.
#################################################################################################

#Choosing a soft-thresholding power. Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression
#similarity is raised to calculate adjacency. The choice of the soft thresholding power is based on the criterion of approximate scale-free topology:

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(PRADdata, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf(file = "Plots/Soft-thresholding power for PRAD.pdf", width = 12, height = 9)

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

dev.off()
#Based on the scale-free topology graph, the soft-thresholding power of 6 was chosen.
#################################################################################################

# Constructing the gene network and identifying modules for PRAD:
net_PRAD = blockwiseModules(PRADdata, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PRADTOM",
                       verbose = 3)

# To see how many modules were identified and what the module sizes are, one can use table(net$colors).
table(net_PRAD$colors)

#Now we can visualize the modules.

# Convert labels to colors for plotting
mergedColors_PRAD = labels2colors(net_PRAD$colors)

# Plot the dendrogram and the module colors underneath
par(mfrow = c(1,1))
pdf(file = "Plots/Cluster dendrogram for PRAD.pdf", width = 12, height = 9)

plotDendroAndColors(net_PRAD$dendrograms[[1]], mergedColors_PRAD[net_PRAD$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Cluster dendrogram for PRAD")
 dev.off()
 
# We now save the module assignment and module eigengene information necessary for subsequent analysis:
moduleLabels_PRAD = net_PRAD$colors
moduleColors_PRAD = labels2colors(net_PRAD$colors)
MEs_PRAD = net_PRAD$MEs;
geneTree_PRAD = net_PRAD$dendrograms[[1]];
save(MEs_PRAD, moduleLabels_PRAD, moduleColors_PRAD, geneTree_PRAD,
     file = "Data Outputs/PRADnetwork_modulecolor_and_label.RData")
#################################################################################################

#We now create a data frame holding the following information for all genes: gene names, 
#module color,and module membership and p-values in all modules:

# Extract eigengenes:
# Define numbers of genes and samples
nGenes = ncol(PRADdata)
nSamples = nrow(PRADdata)

# Recalculate MEs with color labels
MEs0_PRAD = moduleEigengenes(PRADdata, moduleColors_PRAD)$eigengenes
MEs_PRAD = orderMEs(MEs0_PRAD)

# Calculate module membership to identify important genes. Module membership is defined as the correlation of the module eigengene 
# with gene expression profile of the gene.

# names (colors) of the modules
modNames_PRAD = substring(names(MEs_PRAD), 3)

geneModuleMembership_PRAD = as.data.frame(cor(PRADdata, MEs_PRAD, use = "p"));
MMPvalue_PRAD = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_PRAD), nSamples));

names(geneModuleMembership_PRAD) = paste("MM", modNames_PRAD, sep="");
names(MMPvalue_PRAD) = paste("p.MM", modNames_PRAD, sep="");

# Create the starting data frame
geneInfoPRAD = data.frame(Genename = colnames(PRADdata),
                          moduleColor = moduleColors_PRAD,
                          geneModuleMembership_PRAD,
                          MMPvalue_PRAD)

#Order the genes in the geneInfo variable by module color:
geneOrder_PRAD = order(geneInfoPRAD$moduleColor)
geneInfoPRAD_1 = geneInfoPRAD[geneOrder_PRAD, ]

#Save the data frame into a text-format spreadsheet:
write.csv(geneInfoPRAD_1, file = "Data Outputs/PRAD_geneMM.csv")
#################################################################################################

#Now we calculate scaled connectivity of genes in the PRAD network:

#Create a TOM matrix:
tom_PRAD = TOMsimilarityFromExpr(PRADdata)
fun = fundamentalNetworkConcepts(tom_PRAD, GS = NULL) #May take > 30 minutes to compute

#We want the scaled connectivity k = connectivity/max(connectivity), which is an indication of hub gene significance.
connectivity_PRAD = as.data.frame(fun$ScaledConnectivity)
row.names(connectivity_PRAD) = colnames(PRADdata)

#Next step is to merge the gene significance data frame with TCGA PubMed ID data.

#First, we separate the row name into two columns (official gene symbol and Entrez ID) before merging with PubMed info:
connectivity_PRAD = tibble::rownames_to_column(connectivity_PRAD, "gene") #To make the rowname column into a new column, using tidiverse.

#Separate the "gene" column:
connectivity_PRAD = separate(connectivity_PRAD, 1, into = c("genename", "entrez"), sep = "([|])", remove = FALSE)

#Reformat the genesig_PRAD data frame:
row.names(connectivity_PRAD) = connectivity_PRAD$gene
connectivity_PRAD = connectivity_PRAD[,-1]

#Now we can merge genesig_PRAD and PubMed:
connectivity_PubMed_PRAD = merge(PubMed, connectivity_PRAD, by = "entrez")
connectivity_PubMed_PRAD = connectivity_PubMed_PRAD[, -4]
colnames(connectivity_PubMed_PRAD) = c("entrez", "genename", "PubMed", "scaledconnectivity")
head(connectivity_PubMed_PRAD, 3)

#Save the data:
save(connectivity_PubMed_PRAD, file = "Data Outputs/Scaled connectivity and PubMed ID for PRAD.Rdata")
#################################################################################################

#Now we calculate correlation between scaled connectivity and PubMed publications for PRAD genes:

#Calcualte Kendall correlation between scaled connectivity and PubMed number:
corAndPvalue(connectivity_PubMed_PRAD$PubMed, connectivity_PubMed_PRAD$scaledconnectivity, method = "kendall")

#Plot a scatterplot to show the correlation between scaled connectivity and publications:
pdf(file = "Plots/Figure_2a.pdf")
plot(data = connectivity_PubMed_PRAD, PubMed ~ scaledconnectivity, xlim = c(0, 1.19),
     xlab = "Scaled connectivity of genes", ylab = "Number of PubMed IDs", pch=19, main = "PRAD",
     col=ifelse(PubMed > 7000 | scaledconnectivity == max(scaledconnectivity), "red", "black"))
with(data = connectivity_PubMed_PRAD,
     text(PubMed ~ scaledconnectivity, pos = 4, cex = 0.80,
          labels=ifelse(PubMed > 7000 | scaledconnectivity == max(scaledconnectivity), 
                        as.character(connectivity_PubMed_PRAD$genename), "")))
dev.off()
#################################################################################################

#We also correlate scaled connectivity with Gleason score, a metric for prostate cancer severity:
colnames(Gleason) = c("genename", "entrez", "Q", "GleasonScore")

#Add Gleason score to connectivity_PubMed_PRAD dataframe:
Gleason_connectivity_PubMed_PRAD = merge(connectivity_PubMed_PRAD, Gleason, by = "entrez")
Gleason_connectivity_PubMed_PRAD = Gleason_connectivity_PubMed_PRAD[, -c(2)]
dim(Gleason_connectivity_PubMed_PRAD)
class(Gleason_connectivity_PubMed_PRAD)
colnames(Gleason_connectivity_PubMed_PRAD)[colnames(Gleason_connectivity_PubMed_PRAD) == "genename.y"] = "gemename"

#Save the data:
save(Gleason_connectivity_PubMed_PRAD, file = "Data Outputs/Scaled connectivity, PubMed number, and Gleason score for PRAD.RData")
#################################################################################################

# Calculate Kendall correlation between Gleason score and number of PubMed publications:
corAndPvalue(Gleason_connectivity_PubMed_PRAD$PubMed, Gleason_connectivity_PubMed_PRAD$GleasonScore, method = "kendall")

# Plot the correlation between Gleason score and number of PubMed publications:
pdf(file = "Plots/Figure_3a.pdf")
plot(data = Gleason_connectivity_PubMed_PRAD, PubMed ~ GleasonScore,
     xlab = "Gleason score", ylab = "Number of PubMed IDs", pch=19, main = "PRAD", xlim = c(-0.3, 0.4),
     col=ifelse(PubMed == max(PubMed)|  GleasonScore> 0.302 , "red", "black"))
with(data = Gleason_connectivity_PubMed_PRAD,
     text(PubMed ~GleasonScore, pos = 4, cex = 0.8,
          labels=ifelse(PubMed == max(PubMed)|  GleasonScore > 0.302, 
                        Gleason_connectivity_PubMed_PRAD$gemename, "")))
dev.off()
#################################################################################################

## 3. WGCNA Analysis for colon adenocarcinoma (COAD) data set

#Remove uninformative data:
COAD = COAD_DATA
row.names(COAD) = COAD_DATA[,1]; colnames(COAD) = COAD_DATA[1,]
COAD = COAD[-c(1,2),-1]

#Convert data frame to numeric matrix:
COAD1 = as.matrix(sapply(COAD, as.numeric))  
colnames(COAD1) = colnames(COAD)
row.names(COAD1) = row.names(COAD)
class(COAD1)
is.numeric(COAD1)

#Transpose matrix and filter for 10,000 most variant genes:
COADdata = t(COAD1[order(apply(COAD1,1,mad), decreasing = T)[1:10000],])

#Check if COADdata have many missing values:
gsg = goodSamplesGenes(COADdata, verbose = 3)
gsg$allOK

#The command returns "TRUE", so all genes have passed the cuts.
#################################################################################################

# Re-cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers:
sampleTree = hclust(dist(COADdata), method = "average")
# Plot the sample tree:
sizeGrWindow(12,9)

pdf(file = "Plots/Sample clustering to detect outliers_COAD.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()
# It appears there is 1 outlier.
#################################################################################################

#Outlier removal: 

# Plot a line to show the cut

sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf(file = "Plots/Sample clustering to detect outliers_outlier removal_COAD.pdf")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 1500000, col = "red");

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 1500000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
COADdata1 = COADdata[keepSamples, ]
nGenes = ncol(COADdata1)
nSamples = nrow(COADdata1)

#Now the COADdata1 matrix is ready for WGCNA analysis.
#################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(COADdata1, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf(file = "Plots/Soft threshold power for COAD.pdf")
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
dev.off()
#################################################################################################

# Based on the scale-free topology graph, the soft-thresholding power of 6 was chosen.

# Constructing the gene network and identifying modules:
net_COAD = blockwiseModules(COADdata1, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "COADTOM",
                       verbose = 3)

# To see how many modules were identified and what the module sizes are, one can use table(net$colors).
table(net_COAD$colors)
#################################################################################################

# Now we can visualize the modules.

# Convert labels to colors for plotting
mergedColors_COAD = labels2colors(net_COAD$colors)

# Plot the dendrogram and the module colors underneath
pdf("Plots/Cluster dendrogram for COAD.pdf")
plotDendroAndColors(net_COAD$dendrograms[[1]], mergedColors_COAD[net_COAD$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Cluster dendrogram for COAD")
dev.off()

# We now save the module assignment and module eigengene information necessary for subsequent analysis:
moduleLabels_COAD = net_COAD$colors
moduleColors_COAD = labels2colors(net_COAD$colors)
MEs_COAD = net_COAD$MEs;
geneTree_COAD = net_COAD$dendrograms[[1]];
save(MEs_COAD, moduleLabels_COAD, moduleColors_COAD, geneTree_COAD,
     file = "Data Outputs/COADnetwork_modulecolor_and_label.RData")

# Define numbers of genes and samples
nGenes = ncol(COADdata1)
nSamples = nrow(COADdata1)

# Recalculate MEs with color labels
MEs0_COAD = moduleEigengenes(COADdata1, moduleColors_COAD)$eigengenes
MEs_COAD = orderMEs(MEs0_COAD)

#Calculate module membership to identify important genes. 
# names (colors) of the modules
modNames_COAD = substring(names(MEs_COAD), 3)

geneModuleMembership_COAD = as.data.frame(cor(COADdata1, MEs_COAD, use = "p"));
MMPvalue_COAD = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_COAD), nSamples));

names(geneModuleMembership_COAD) = paste("MM", modNames_COAD, sep="");
names(MMPvalue_COAD) = paste("p.MM", modNames_COAD, sep="");
#################################################################################################

#We now create a data frame holding the following information for all genes: gene names, 
#module color,and module membership and p-values in all modules: 

# Create the starting data frame
geneInfoCOAD = data.frame(Genename = colnames(COADdata1),
                          moduleColor = moduleColors_COAD,
                          geneModuleMembership_COAD,
                          MMPvalue_COAD)

#Order the genes in the geneInfo variable by module color:
geneOrder_COAD = order(geneInfoCOAD$moduleColor)
geneInfoCOAD_1 = geneInfoCOAD[geneOrder_COAD, ]

save(geneInfoCOAD_1, file = "Data Outputs/COAD_geneMM.Rfile")
#################################################################################################


# Now we calculate scaled connectivity of genes in the COAD network:

#Create a TOM matrix:
tom_COAD = TOMsimilarityFromExpr(COADdata1)
fun = fundamentalNetworkConcepts(tom_COAD, GS = NULL) #Take > 30 minutes to compute

# We want the scaled connectivity k = connectivity/max(connectivity), which is an indication of hub gene significance.
connectivity_COAD = as.data.frame(fun$ScaledConnectivity)
row.names(connectivity_COAD) = colnames(COADdata1) 

#To visualize correlation between scaled connectivity of genes and number of publications, we merge the two variables into one dataframe. 

#Separate the row name of the connectivity_COAD dataframe into two columns (official gene symbol and Entrez ID) before merging with PubMed info:
connectivity_COAD = tibble::rownames_to_column(connectivity_COAD, "gene") #To make the rowname column into a new column.

#Separate the "gene" column:
connectivity_COAD = separate(connectivity_COAD, 1, into = c("genename", "entrez"), sep = "([|])", remove = FALSE)

#Reformat the connectivity_COAD dataframe:
row.names(connectivity_COAD) = connectivity_COAD$gene
connectivity_COAD = connectivity_COAD[,-1]

#Now we can merge genesig_COAD and PubMed:
connectivity_PubMed_COAD = merge(PubMed, connectivity_COAD, by = "entrez")
connectivity_PubMed_COAD = connectivity_PubMed_COAD[, -4]
colnames(connectivity_PubMed_COAD) = c("entrez", "genename", "PubMed", "scaledconnectivity")

#Save the data:
save(connectivity_PubMed_COAD, file = "Data Outputs/Scaled connectivity and PubMed ID for COAD.Rfile")

#Now we calculate correlation between scaled connectivity and PubMed publications for PRAD genes:

#Calcualte Kendall correlation between scaled connectivity and PubMed number:
corAndPvalue(connectivity_PubMed_COAD$PubMed, connectivity_PubMed_COAD$scaledconnectivity, method = "kendall")

#Plot a scatterplot to show the correlation between scaled connectivity and publications:

plot(data = connectivity_PubMed_COAD, PubMed ~ scaledconnectivity, xlim = c(0, 1.19),
     xlab = "Scaled connectivity of genes", ylab = "Number of PubMed IDs", pch=19, main = "COAD",
     col=ifelse(PubMed > 7000 | scaledconnectivity == max(scaledconnectivity), "red", "black"))
with(data = connectivity_PubMed_COAD,
     text(PubMed ~ scaledconnectivity, pos = 4, cex = 0.80,
          labels=ifelse(PubMed > 7000 | scaledconnectivity == max(scaledconnectivity), 
                        as.character(connectivity_PubMed_COAD$genename), "")))

# We also correlate scaled connectivity with aggressiveness score:

#Load aggressiveness score data:
aggressiveness_COAD = as.data.frame(COAD.aggressiveness[, c(3,7)])

#Add aggressiveness score to connectivity_PubMed_COAD dataframe:
aggressiveness_connectivity_PubMed_COAD = merge(connectivity_PubMed_COAD, aggressiveness_COAD, by = "genename")
dim(aggressiveness_connectivity_PubMed_COAD)
#Save the data:
save(aggressiveness_connectivity_PubMed_COAD, file = "Data Outputs/Scaled connectivity, PubMed IDs, and aggressivenss score for COAD.Rdata")

#Perform correlation between scaled connectivity and aggressiveness score:
corAndPvalue(aggressiveness_connectivity_PubMed_COAD$scaledconnectivity, aggressiveness_connectivity_PubMed_COAD$aggressiveness, method = "kendall")

# Create a figure for the correlation between scaled connectivity and aggressiveness score:

plot(data = aggressiveness_connectivity_PubMed_COAD, aggressiveness ~ scaledconnectivity, xlim = c(0, 1.19),
     xlab = "Scaled connectivity of genes", ylab = "Aggressiveness score", pch=19, main = "COAD",
     col=ifelse(aggressiveness > 11|scaledconnectivity == max(scaledconnectivity), "red", "black"))
with(data = aggressiveness_connectivity_PubMed_COAD, text(aggressiveness ~ scaledconnectivity, pos = 4, cex = 0.80,
                                           labels=ifelse(aggressiveness > 11|scaledconnectivity == max(scaledconnectivity), 
                                                        as.character(aggressiveness_connectivity_PubMed_COAD$genename), "")))

#Calculate correlation between aggressiveness score and number publications:

#Kendall correlation between aggressiveness score and number of PubMed publications:
corAndPvalue(aggressiveness_connectivity_PubMed_COAD$PubMed, aggressiveness_connectivity_PubMed_COAD$aggressiveness,method = "kendall")

# Create the corrrelation figure for aggressiveness score vs. number of PubMed publications:
plot(data = aggressiveness_connectivity_PubMed_COAD, PubMed ~ aggressiveness,
     xlab = "Aggressiveness score", ylab = "Number of PubMed IDs", pch=19, main = "COAD", xlim = c(-15, 15),
     col=ifelse(aggressiveness > 11|PubMed == max(PubMed), "red", "black"))
with(data = aggressiveness_connectivity_PubMed_COAD, text(PubMed ~ aggressiveness, pos = 4, cex = 0.80,
                                                   labels=ifelse(aggressiveness > 11|PubMed == max(PubMed), 
                                                                 as.character(aggressiveness_connectivity_PubMed_COAD$genename), "")))
#################################################################################################

## 4. WGCNA Analysis for glioma (GBMLGG) data set

#Remove uninformative data:
GBMLGG = GBMLGG_DATA
row.names(GBMLGG) = GBMLGG_DATA[,1]; colnames(GBMLGG) = GBMLGG_DATA[1,]
GBMLGG = GBMLGG[-c(1,2),-1]

#Convert data frame to numeric matrix:
GBMLGG1 = as.matrix(sapply(GBMLGG, as.numeric))  
colnames(GBMLGG1) = colnames(GBMLGG)
row.names(GBMLGG1) = row.names(GBMLGG)
class(GBMLGG1)
is.numeric(GBMLGG1)

#Transpose matrix and filter for 10,000 most variant genes:
GBMLGGdata = t(GBMLGG1[order(apply(GBMLGG1,1,mad), decreasing = T)[1:10000],])

#Check if GBMLGGdata have many missing values:
gsg = goodSamplesGenes(GBMLGGdata, verbose = 3)
gsg$allOK

#The command returns "TRUE", so all genes have passed the cuts.
#################################################################################################

# Re-cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers:
sampleTree = hclust(dist(GBMLGGdata), method = "average")
# Plot the sample tree:
sizeGrWindow(12,9)
pdf(file = "Plots/Sample clustering to detect outliers_GBMLGG.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
#################################################################################################

#It appears there is 1 outlier. Outlier removal: 

# Plot a line to show the cut
sizeGrWindow(12,9)
pdf(file = "Plots/Sample clustering to detect outliers_GBMLGG_outlier removal.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 2500000, col = "red")

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2500000, minSize = 10)
table(clust)
# Clust 1 contains the samples we want to keep
keepSamples = (clust==1)
GBMLGGdata1 = GBMLGGdata[keepSamples, ]

# Now the GBMLGGdata1 matrix is ready for WGCNA analysis.
#################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(GBMLGGdata1, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "Plots/Soft-thresholding pwoer for GBMLGG.pdf")

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

dev.off()

# Based on the scale-free topology graph, the soft-thresholding power of 6 was chosen.
#################################################################################################

# Constructing the gene network and identifying modules:
net_GBMLGG = blockwiseModules(GBMLGGdata1, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GBMLGGTOM",
                       verbose = 3)

# To see how many modules were identified and what the module sizes are, one can use table(net$colors).
table(net_GBMLGG$colors)

#Now we can visualize the modules.

# Convert labels to colors for plotting
mergedColors_GBMLGG = labels2colors(net_GBMLGG$colors)

# Plot the dendrogram and the module colors underneath
pdf(file = "Plots/Cluster dendrogram for GBMLGG.pdf")
plotDendroAndColors(net_GBMLGG$dendrograms[[1]], mergedColors_GBMLGG[net_GBMLGG$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Cluster dendrogram for GBMLGG")
dev.off()
#################################################################################################

# We now save the module assignment and module eigengene information necessary for subsequent analysis:
moduleLabels_GBMLGG = net_GBMLGG$colors
moduleColors_GBMLGG = labels2colors(net_GBMLGG$colors)
MEs_GBMLGG = net_GBMLGG$MEs;
geneTree_GBMLGG = net_GBMLGG$dendrograms[[1]];
save(MEs_GBMLGG, moduleLabels_GBMLGG, moduleColors_GBMLGG, geneTree_GBMLGG,
     file = "Data Outputs/GBMLGG_network_modulecolor_and_label.RData")

# Define numbers of genes and samples
nGenes = ncol(GBMLGGdata1)
nSamples = nrow(GBMLGGdata1)

# Recalculate MEs with color labels
MEs0_GBMLGG = moduleEigengenes(GBMLGGdata1, moduleColors_GBMLGG)$eigengenes
MEs_GBMLGG = orderMEs(MEs0_GBMLGG)

#Calculate module membership to identify important genes. 
# names (colors) of the modules
modNames_GBMLGG = substring(names(MEs_GBMLGG), 3)

geneModuleMembership_GBMLGG = as.data.frame(cor(GBMLGGdata1, MEs_GBMLGG, use = "p"));
MMPvalue_GBMLGG = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_GBMLGG), nSamples));

names(geneModuleMembership_GBMLGG) = paste("MM", modNames_GBMLGG, sep="");
names(MMPvalue_GBMLGG) = paste("p.MM", modNames_GBMLGG, sep="");

#We now create a data frame holding the following information for all genes: gene names, 
#module color,and module membership and p-values in all modules: 

# Create the starting data frame
geneInfoGBMLGG = data.frame(Genename = colnames(GBMLGGdata1),
                          moduleColor = moduleColors_GBMLGG,
                          geneModuleMembership_GBMLGG,
                          MMPvalue_GBMLGG)

#Order the genes in the geneInfo variable by module color:
geneOrder_GBMLGG = order(geneInfoGBMLGG$moduleColor)
geneInfoGBMLGG_1 = geneInfoGBMLGG[geneOrder_GBMLGG, ]

# Save the data frame into a text-format spreadsheet:
save(geneInfoGBMLGG_1, file = "Data Outputs/GBMLGG_geneMM.Rfile")
#################################################################################################

# Now we calculate scaled connectivity of genes in the GBMLGG network:

#Create a TOM matrix:
tom_GBMLGG = TOMsimilarityFromExpr(GBMLGGdata1)
fun = fundamentalNetworkConcepts(tom_GBMLGG, GS = NULL) #Take > 30 minutes to compute

# We want the scaled connectivity k = connectivity/max(connectivity), which is an indication of hub gene significance.
connectivity_GBMLGG = as.data.frame(fun$ScaledConnectivity)
row.names(connectivity_GBMLGG) = colnames(GBMLGGdata1) 

#To visualize correlation between scaled connectivity of genes and number of publications, we merge the two variables into one dataframe. 

#Separate the row name of the connectivity_COAD dataframe into two columns (official gene symbol and Entrez ID) before merging with PubMed info:
connectivity_GBMLGG = tibble::rownames_to_column(connectivity_GBMLGG, "gene") #To make the rowname column into a new column.

#Separate the "gene" column:
connectivity_GBMLGG = separate(connectivity_GBMLGG, 1, into = c("genename", "entrez"), sep = "([|])", remove = FALSE)

#Reformat the connectivity_GBMLGG dataframe:
row.names(connectivity_GBMLGG) = connectivity_GBMLGG$gene
connectivity_GBMLGG = connectivity_GBMLGG[,-1]

#Now we can merge genesig_COAD and PubMed:
connectivity_PubMed_GBMLGG = merge(PubMed, connectivity_GBMLGG, by = "entrez")
connectivity_PubMed_GBMLGG = connectivity_PubMed_GBMLGG[, -4]
colnames(connectivity_PubMed_GBMLGG) = c("entrez", "genename", "PubMed", "scaledconnectivity")

#Save the data:
save(connectivity_PubMed_GBMLGG, file = "Data Outputs/Scaled connectivity and PubMed ID for GBMLGG.Rfile")
#################################################################################################

#Now we calculate correlation between scaled connectivity and PubMed publications for GBMLGG genes:

#Calcualte Kendall correlation between scaled connectivity and PubMed number:
corAndPvalue(connectivity_PubMed_GBMLGG$PubMed, connectivity_PubMed_GBMLGG$scaledconnectivity, method = "kendall")

#Plot a scatterplot to show the correlation between scaled connectivity and publications:
pdf(file = "Plots/Figure_3a.pdf")
plot(data = connectivity_PubMed_GBMLGG, PubMed ~ scaledconnectivity, xlim = c(0, 1.19),
     xlab = "Scaled connectivity of genes", ylab = "Number of PubMed IDs", pch=19, main = "GBMLGG",
     col=ifelse(PubMed > 7000 | scaledconnectivity == max(scaledconnectivity), "red", "black"))
with(data = connectivity_PubMed_GBMLGG,
     text(PubMed ~ scaledconnectivity, pos = 4, cex = 0.80,
          labels=ifelse(PubMed > 7000 | scaledconnectivity == max(scaledconnectivity), 
                        as.character(connectivity_PubMed_GBMLGG$genename), "")))
dev.off()
#################################################################################################

#We also correlate scaled connectivity with C index score, a metric for glioma severity:
colnames(GBMLGG.Cindex) = c("genename", "entrez", "Q", "C_index")

#Add Cindex score to connectivity_PubMed_GBMLGG dataframe:
Cindex_connectivity_PubMed_GBMLGG = merge(connectivity_PubMed_GBMLGG, GBMLGG.Cindex, by = "entrez")
Cindex_connectivity_PubMed_GBMLGG = Cindex_connectivity_PubMed_GBMLGG[, -c(2)]
dim(Cindex_connectivity_PubMed_GBMLGG)
class(Cindex_connectivity_PubMed_GBMLGG)
colnames(Cindex_connectivity_PubMed_GBMLGG)[colnames(Cindex_connectivity_PubMed_GBMLGG) == "genename.y"] = "gemename"

#Save the data:
save(Cindex_connectivity_PubMed_GBMLGG, file = "Data Outputs/Scaled connectivity, PubMed number, and Cindex score for GBMLGG.RData")
#################################################################################################

# Calculate Kendall correlation between Cindex score and number of PubMed publications:
corAndPvalue(Cindex_connectivity_PubMed_GBMLGG$PubMed, Cindex_connectivity_PubMed_GBMLGG$C_index, method = "kendall")

# Plot the correlation between Cindex score and number of PubMed publications:
# Create the corrrelation figure for Cindex score vs. number of PubMed publications:
pdf(file = "Plots/Figure_3c.pdf")

plot(data = Cindex_connectivity_PubMed_GBMLGG, PubMed ~ C_index, xlab = "C-index score", 
     ylab = "Number of PubMed IDs", pch=19, main = "GBMLGG",
     col=ifelse(C_index > 0.83|PubMed == max(PubMed), "red", "black"))
     with(data = Cindex_connectivity_PubMed_GBMLGG, text(PubMed ~ C_index, pos = 4, cex = 0.80,
    labels=ifelse(C_index > 0.83|PubMed == max(PubMed), 
                  as.character(Cindex_connectivity_PubMed_GBMLGG$genename),"")))
          
dev.off()

#################################################################################################

## 5. Examine the distribution of underannotated genes in modules:

# 5a. Examine the distribution of underannotated genes in PRAD:
Modulemembership_PRAD = geneInfoPRAD_1[,c(1,2)]
colnames(Modulemembership_PRAD) = c("genename", "modulecolor")

#Separate the "genename" column:
Modulemembership_PRAD = separate(Modulemembership_PRAD, 1, into = c("genename", "entrez"), sep = "([|])", remove = TRUE)

#Merge the module color with PubMed ID:
module_PubMed_PRAD = merge(Modulemembership_PRAD, PubMed, by = "entrez")
module_PubMed_PRAD = module_PubMed_PRAD[,-4]
module_PubMed_PRAD = as.data.frame(module_PubMed_PRAD)

#Extract only underannoatation (PubMed <51) genes into a new data set:
underannotated_PRAD = module_PubMed_PRAD[which(module_PubMed_PRAD$PubMed < 51),]

#Count the number of all genes and underannotated genes in each module:
allannotation_bymodule = aggregate(module_PubMed_PRAD$PubMed ~ module_PubMed_PRAD$modulecolor, data = module_PubMed_PRAD,
                                  length)
colnames(allannotation_bymodule) = c("modulecolor", "PubMed")
underannotation_bymodule = aggregate(underannotated_PRAD$PubMed ~ underannotated_PRAD$modulecolor, data = underannotated_PRAD,
                                  length)
colnames(underannotation_bymodule) = c("modulecolor", "PubMed")
annotation_bymodule = merge(allannotation_bymodule, underannotation_bymodule, by = "modulecolor")          
colnames(annotation_bymodule) = c("Module.color", "Total.genes", "Functionally.enigmatic.genes")
annotation_bymodule = mutate(annotation_bymodule, Percent.functionally.enigmatic.genes = Functionally.enigmatic.genes*100/Total.genes)
annotation_bymodule = arrange(annotation_bymodule, desc(Percent.functionally.enigmatic.genes))
head(annotation_bymodule)
save(annotation_bymodule, file = "Data Outputs/Annotation by module PRAD.Rfile")
#################################################################################################

# 5b. Examine the distribution of underannotated genes in COAD:
Modulemembership_COAD = geneInfoCOAD_1[,c(1,2)]
colnames(Modulemembership_COAD) = c("genename", "modulecolor")

#Separate the "genename" column:
Modulemembership_COAD = separate(Modulemembership_COAD, 1, into = c("genename", "entrez"), sep = "([|])", remove = TRUE)

#Merge the module color with PubMed ID:
module_PubMed_COAD = merge(Modulemembership_COAD, PubMed, by = "entrez")
module_PubMed_COAD = module_PubMed_COAD[,-4]
module_PubMed_COAD = as.data.frame(module_PubMed_COAD)

#Extract only underannoatation (PubMed <51) genes into a new data set:
underannotated_COAD = module_PubMed_COAD[which(module_PubMed_COAD$PubMed < 51),]

#Count the number of all genes and underannotated genes in each module:
allannotation_bymodule = aggregate(module_PubMed_COAD$PubMed ~ module_PubMed_COAD$modulecolor, data = module_PubMed_COAD,
                                   length)
colnames(allannotation_bymodule) = c("modulecolor", "PubMed")
underannotation_bymodule = aggregate(underannotated_COAD$PubMed ~ underannotated_COAD$modulecolor, data = underannotated_COAD,
                                     length)
colnames(underannotation_bymodule) = c("modulecolor", "PubMed")
annotation_bymodule = merge(allannotation_bymodule, underannotation_bymodule, by = "modulecolor")          
colnames(annotation_bymodule) = c("Module.color", "Total.genes", "Functionally.enigmatic.genes")
annotation_bymodule = mutate(annotation_bymodule, Percent.functionally.enigmatic.genes = Functionally.enigmatic.genes*100/Total.genes)
annotation_bymodule = arrange(annotation_bymodule, desc(Percent.functionally.enigmatic.genes))
head(annotation_bymodule)
save(annotation_bymodule, file = "Data Outputs/Annotation by module COAD.Rfile")
#################################################################################################

# 5c. Examine the distribution of underannotated genes in GBMLGG:
Modulemembership_GBMLGG = geneInfoGBMLGG_1[,c(1,2)]
colnames(Modulemembership_GBMLGG) = c("genename", "modulecolor")

#Separate the "genename" column:
Modulemembership_GBMLGG = separate(Modulemembership_GBMLGG, 1, into = c("genename", "entrez"), sep = "([|])", remove = TRUE)

#Merge the module color with PubMed ID:
module_PubMed_GBMLGG = merge(Modulemembership_GBMLGG, PubMed, by = "entrez")
module_PubMed_GBMLGG = module_PubMed_GBMLGG[,-4]
module_PubMed_GBMLGG = as.data.frame(module_PubMed_GBMLGG)

#Extract only underannoatation (PubMed <51) genes into a new data set:
underannotated_GBMLGG = module_PubMed_GBMLGG[which(module_PubMed_GBMLGG$PubMed < 51),]

#Count the number of all genes and underannotated genes in each module:
allannotation_bymodule = aggregate(module_PubMed_GBMLGG$PubMed ~ module_PubMed_GBMLGG$modulecolor, data = module_PubMed_GBMLGG,
                                   length)
colnames(allannotation_bymodule) = c("modulecolor", "PubMed")
underannotation_bymodule = aggregate(underannotated_GBMLGG$PubMed ~ underannotated_GBMLGG$modulecolor, data = underannotated_GBMLGG,
                                     length)
colnames(underannotation_bymodule) = c("modulecolor", "PubMed")
annotation_bymodule = merge(allannotation_bymodule, underannotation_bymodule, by = "modulecolor")          
colnames(annotation_bymodule) = c("Module.color", "Total.genes", "Functionally.enigmatic.genes")
annotation_bymodule = mutate(annotation_bymodule, Percent.functionally.enigmatic.genes = Functionally.enigmatic.genes*100/Total.genes)
annotation_bymodule = arrange(annotation_bymodule, desc(Percent.functionally.enigmatic.genes))
head(annotation_bymodule)
save(annotation_bymodule, file = "Data Outputs/Annotation by module GBMLGG.Rfile")
#################################################################################################



## 6. Extract modules for network visualization in Cytoscape:

# 6a. Extract GBMLGG "saddlebrown" module:

# Select modules
modules = c("saddlebrown")

# Select module probes
Genes = colnames(GBMLGGdata1)
inModule = is.finite(match(moduleColors_GBMLGG, modules))
modGenes =Genes[inModule]

# Select the corresponding topological overlap
modTOM = tom_GBMLGG[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("GBMLGG_edges_", paste(modules, collapse="_"), ".txt", sep=""),
                               nodeFile = paste("GBMLGG_nodes_", paste(modules, collapse="_"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.06,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors_GBMLGG[inModule])

#################################################################################################
# 6b. Extract COAD "cyan" module containing APOL6:

# Select modules
modules = c("cyan")

# Select module probes
Genes = colnames(COADGdata1)
inModule = is.finite(match(moduleColors_COAD, modules))
modGenes =Genes[inModule]

# Select the corresponding Topological Overlap
modTOM = tom_COAD[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("COAD_edges_", paste(modules, collapse="_"), ".txt", sep=""),
                               nodeFile = paste("COAD_nodes_", paste(modules, collapse="_"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.06,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors_COAD[inModule])
#################################################################################################

# c. Extract COAD "blue" module containing C6orf48:

# Select modules
modules = c("blue")

# Select module probes
Genes = colnames(COADGdata1)
inModule = is.finite(match(moduleColors_COAD, modules))
modGenes =Genes[inModule]

# Select the corresponding Topological Overlap
modTOM = tom_COAD[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("COAD_edges_", paste(modules, collapse="_"), ".txt", sep=""),
                               nodeFile = paste("COAD_nodes_", paste(modules, collapse="_"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors_COAD[inModule])
