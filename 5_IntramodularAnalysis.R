T2DM = as.data.frame(pheno$T2DM)
# names (colors) of the modules
modNames = substring(names(MEs), 3) #module name: from 'MEyellow' to 'yellow'
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, T2DM, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(T2DM), sep="")
names(GSPvalue) = paste("p.GS.", names(T2DM), sep="")

module5 = "lightgreen" #p = 0.02
column5 = match(module5, modNames) #15
module5Genes = moduleColors==module5 #145

module6 = "magenta" #p = 0.05
column6 = match(module6, modNames) #4
module6Genes = moduleColors==module6 #213
sizeGrWindow(4.4, 2.2)
par(mfrow = c(1,2), mar = c(5, 5, 4, 2), ps = 10, cex = 1)

verboseScatterplot(abs(geneModuleMembership[module5Genes, column5]),
                   abs(geneTraitSignificance[module5Genes, 1]),
                   xlab = paste("Module Membership in", module5, "Module"),
                   ylab = "Gene Significance for T2DM",
                   main = paste("A. GS vs. MM in lightgreen Module\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module5,
                   pch = 20,
                   abline = TRUE,
                   abline.lty = "dashed")

verboseScatterplot(abs(geneModuleMembership[module6Genes, column6]),
                   abs(geneTraitSignificance[module6Genes, 1]),
                   xlab = paste("Module Membership in", module6, "Module"),
                   ylab = "Gene Significance for T2DM",
                   main = paste("B. GS vs. MM in magenta Module\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module6,
                   pch = 20,
                   abline = TRUE,
                   abline.lty = "dashed")


names(datExpr) #all geneIDs
names(datExpr)[moduleColors=="lightgreen"]
geneInfo0 = data.frame(gene_ID = annot$Description[ENSG2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
modOrder = order(-abs(cor(MEs, MHT2D, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pheno.MHT2D))
geneInfo = geneInfo0[geneOrder, ]
library('org.Hs.eg.db')
hs <- org.Hs.eg.db #homo sapiens gene database
library(dplyr)
detach("package:dplyr", unload=TRUE)
gene_symbol <- as.character(unique(geneInfo[, "gene_ID"])) #duplicates: SNORA63, 6043, keep ENSG00000200320.1
geneInfo = subset(geneInfo, rownames(geneInfo) != "ENSG00000200418.1") #remove row ENSG00000200418.1
ccols <- c("ENTREZID","SYMBOL")
gene_symbol2ENTREZID = select(hs, keys = gene_symbol,columns = ccols,keytype = "SYMBOL")
#gene HBD mapped to two ENTREZID, double check with ENSG_ID in GTEX ENSG, HBD = 3045, then remove the other duplicates
gene_symbol2ENTREZID = subset(gene_symbol2ENTREZID[duplicated(gene_symbol2ENTREZID$SYMBOL) == "FALSE",])
geneInfo = cbind(geneInfo, gene_symbol2ENTREZID)
geneInfo[, "SYMBOL"] = NULL
geneInfo = geneInfo[, c(1, 85, 2:84)] #reorder the columns
lightgreenGenes = geneInfo[geneInfo$moduleColor == "lightgreen", ]
magentaGenes = geneInfo[geneInfo$moduleColor == "magenta", ]
write.csv(lightgreenGenes, file = "./lightgreenGenes.csv")
write.csv(magentaGenes, file = "./magentaGenes.csv")
