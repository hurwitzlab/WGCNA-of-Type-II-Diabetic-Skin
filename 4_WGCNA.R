library(matrixStats)
library(WGCNA)

sampleID = c(rownames(control_pheno), rownames(type2_pheno))
datExpr = tpm_clean[, colnames(tpm_clean) %in% sampleID]
datExpr = data.frame(t(datExpr))
pheno = subset(pheno, rownames(pheno) %in% sampleID)

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize = 9000)
sizeGrWindow(6.8, 3)
par(mfrow = c(1,2), mar = c(5, 5, 4, 2), ps = 10, cex = 1)
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("A. Scale independence"), cex.main = 1);
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.80,col="red")
#abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("B. Mean connectivity"), cex.main = 1)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#power estimate = 7

net_OneNetwork = blockwiseModules(datExpr, power = 7,
                                  TOMType = "unsigned", minModuleSize = 30,
                                  reassignThreshold = 0, mergeCutHeight = 0.20,
                                  numericLabels = TRUE, pamRespectsDendro = FALSE,
                                  saveTOMs = TRUE,
                                  saveTOMFileBase = "OneNetworkTOM",
                                  verbose = 3,
                                  maxBlockSize = 9000)
# open a graphics window
sizeGrWindow(12, 10)
# Convert labels to colors for plotting
mergedColors = labels2colors(net_OneNetwork$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_OneNetwork$dendrograms[[1]], mergedColors[net_OneNetwork$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.colorLabels = 1,
                    main = "A. Cluster Dendrogram",
                    cex.main = 1.2)
table(net_OneNetwork$colors) #39 modules

moduleLabels = net_OneNetwork$colors
moduleColors = labels2colors(net_OneNetwork$colors)
MEs = net_OneNetwork$MEs;
geneTree = net_OneNetwork$dendrograms[[1]];
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, pheno[1:5], use = "p") #whether modules are significant or not
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(4, 7, 2, 1), ps = 10, cex = 1) #c(bottom, left, top, right)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(pheno)[1:5],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("B. Module-Trait Relationships"),
               cex.main = 1.2)

