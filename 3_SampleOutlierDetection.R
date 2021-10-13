#Sample outlier detection: standardized connectivity (Z.K method)
#datX: the rows correpsond to samples and the column correspond to network nodes (genes)
#if we use sample network for outlier detection, then the nodes should be samples and make sure nodes are columns
sample_control = tpm_clean[, colnames(tpm_clean) %in% colnames(tpm_control)]
SampleControlEuclid = adjacency(sample_control, type = "distance")
k_control = as.numeric(apply(SampleControlEuclid, 2, sum)) - 1
Z.k_control = scale(k_control)
Z.k_controlTab = cbind(rownames(SampleControlEuclid), Z.k_control)
row.names(Z.k_controlTab) = Z.k_controlTab[, 1]
Z.k_controlTab = Z.k_controlTab[, -1]
Z.k_controlTab = as.data.frame(Z.k_controlTab)

thresholdZ.k = -2.5 #often -2.5 or -2
# the color vector indicates outlyingness (red)
#for control samples
outlierColor=ifelse(Z.k_control<thresholdZ.k,"red","black")
# calculate the cluster tree using flashClust or hclust
SampleControlTree = hclust(as.dist(1-SampleControlEuclid), method = 'average')
control_pheno = subset(pheno, pheno$MHT2D == 0)
# Convert traits to a color representation: where red indicates high values
#TraitControlColors = data.frame(numbers2colors(control_pheno, signed = FALSE))
#dimnames(TraitControlColors)[[2]] = paste(names(control_pheno), 'C', sep = '')
datControlColors = data.frame(outlierC=outlierColor)

par(mfrow = c(2,1), mar = c(4, 4, 4, 2), mgp = c(1.5,0.5,0), ps = 10, cex = 1)
plot(Z.k_controlTab$Z.k_controlTab, xaxt = "n", xlab = "Sample", ylab = "Z.K", 
     col = ifelse(Z.k_control<thresholdZ.k,"red","black"),
     pch = 20,
     main = "A. Sample Network Outlier Detection in non-T2DM Group")
text(x=which(as.numeric(Z.k_controlTab$Z.k_controlTab) < thresholdZ.k) + 4,
     y=as.numeric(Z.k_controlTab[x,])-0.15,
     labels = rownames(Z.k_controlTab)[x],
     cex= 0.9,
     font = 0.5)

abline(h=-2.5,col="black", lty = 2)

#remove 3 outliers
removeControlSample = Z.k_control<thresholdZ.k | is.na(Z.k_control) #before '|' is true, then true
sample_control = sample_control[, removeControlSample == 'FALSE']
removeControlSampleID = subset(rownames(SampleControlEuclid), removeControlSample == 'TRUE')
#[1] "GTEX.11WQK.0526.SM.5EQLD" "GTEX.POMQ.2226.SM.2I5FK"  "GTEX.RUSQ.1526.SM.2TF62"
removeControlSamplePheno = subset(control_pheno, rownames(control_pheno) %in% removeControlSampleID)
#removeControlSampleTraits = subset(controlTraits, removeControlSample == 'TRUE')
removeControlSamplePheno$Z.k = Z.k_controlTab$Z.k_controlTab[match(rownames(removeControlSamplePheno), rownames(Z.k_controlTab))]
control_pheno = subset(control_pheno, ! (rownames(control_pheno) %in% removeControlSampleID))

#for type2 samples
sample_type2 = tpm_clean[, colnames(tpm_clean) %in% colnames(tpm_type2)]
SampleType2Euclid = adjacency(sample_type2, type = "distance")
k_type2 = as.numeric(apply(SampleType2Euclid, 2, sum)) - 1
Z.k_type2 = scale(k_type2)
Z.k_type2Tab = cbind(rownames(SampleType2Euclid), Z.k_type2)
row.names(Z.k_type2Tab) = Z.k_type2Tab[, 1]
Z.k_type2Tab = Z.k_type2Tab[, -1]
Z.k_type2Tab = as.data.frame(Z.k_type2Tab)

outlierColor=ifelse(Z.k_type2<thresholdZ.k,"red","black")
# calculate the cluster tree using flashClust or hclust
SampleType2Tree = hclust(as.dist(1-SampleType2Euclid), method = 'average')
type2_pheno = subset(pheno, pheno$MHT2D == 1)
datType2Colors = data.frame(outlierC=outlierColor)

plot(Z.k_type2Tab$Z.k_type2Tab, xaxt = "n", xlab = "Sample", ylab = "Z.K",
     col = ifelse(Z.k_type2<thresholdZ.k,"red","black"),
     pch = 20,
     xgap.axis = 0,
     main = "B. Sample Network Outlier Detection in T2DM Group"
     )
text(x=which(as.numeric(Z.k_type2Tab$Z.k_type2Tab) < thresholdZ.k) + 4,
     y=as.numeric(Z.k_type2Tab[x,])-0.15,
     labels = rownames(Z.k_type2Tab)[x],
     cex= 0.9,
     font = 0.5)

abline(h=-2.5,col="black", lty = 2)

#remove 2 outliers
removeType2Sample = Z.k_type2<thresholdZ.k | is.na(Z.k_type2) #before '|' is true, then true
sample_type2 = sample_type2[, removeType2Sample == 'FALSE']
removeType2SampleID = subset(rownames(SampleType2Euclid), removeType2Sample == 'TRUE')
removeType2SamplePheno = subset(type2_pheno, rownames(type2_pheno) %in% removeType2SampleID)
#removeType2SampleTraits = subset(type2Traits, removeType2Sample == 'TRUE')
removeType2SamplePheno$Z.k = Z.k_type2Tab$Z.k_type2Tab[match(rownames(removeType2SamplePheno), rownames(Z.k_type2Tab))]
type2_pheno = subset(type2_pheno, ! (rownames(type2_pheno) %in% removeType2SampleID))
#done

