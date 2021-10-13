library(limma)
library(sva)
library(WGCNA)
library(ggplot2)

#read raw counts table
counts_control <- read.delim("./FinalGeneReadfileType2isZero.txt", row.names = "Name")
ENSG_geneID = counts_control[, 'Description', drop = FALSE] ##'drop' keep the row name
counts_type2 <- read.delim("./FinalGeneReadfileType2isOne.txt", row.names = "Name")
#counts_all <- read.delim("./FinalGeneReadfile.txt", row.names = "Name")
final_75_150 = read.delim("./final_75_150.txt")
counts_control = counts_control[, gsub('GTEX[.]([^.]+)[.].*', '\\1', colnames(counts_control)) %in% final_75_150$SubID]
counts_type2$Description = NULL
counts_all = cbind(counts_control, counts_type2)

#read tpm value
tpm_control <- read.delim("./FinalGenetpmfileType2isZero.txt", row.names = "Name")
tpm_control = tpm_control[, gsub('GTEX[.]([^.]+)[.].*', '\\1', colnames(tpm_control)) %in% final_75_150$SubID]
tpm_type2 <- read.delim("./FinalGenetpmfileType2isOne.txt", row.names = "Name")
tpm_type2$Description = NULL
tpm_all = cbind(tpm_control, tpm_type2)

#remove low-expressed genes
counts_control_LowExpressed = subset(counts_control, rowQuantileC(counts_control, 0.1) < 10)
counts_type2_LowExpressed = subset(counts_type2, rowQuantileC(counts_type2, 0.1) < 10)
LowExpressedGenesID = intersect(rownames(counts_control_LowExpressed), rownames(counts_type2_LowExpressed))
RemainingGenesID = subset(rownames(counts_all),  ! (rownames(counts_all) %in% LowExpressedGenesID))
counts_all_remaining = subset(counts_all, rownames(counts_all) %in% RemainingGenesID)
counts_control_remaining = subset(counts_control, rownames(counts_control) %in% RemainingGenesID)
counts_type2_remaining = subset(counts_type2, rownames(counts_type2) %in% RemainingGenesID)
tpm_control_remaining = subset(tpm_control, rownames(tpm_control) %in% RemainingGenesID)
tpm_type2_remaining = subset(tpm_type2, rownames(tpm_type2) %in% RemainingGenesID)
tpm_all_remaining = subset(tpm_all, rownames(tpm_all) %in% RemainingGenesID)

#pick top varying genes based on MAD calculated from tpm matrix
#top MAD across all samples: 8812
tpm_all_remaining$MAD = apply(tpm_all_remaining, 1, mad)
tpm_all_TopVaryingGenesID = rownames(subset(tpm_all_remaining, tpm_all_remaining$MAD >= quantile(tpm_all_remaining$MAD, 0.5)))
tpm_network = subset(tpm_all_remaining, rownames(tpm_all_remaining) %in% tpm_all_TopVaryingGenesID) #min MAD = 7.44
tpm_network$MAD = NULL

#log transformation of tpm value
edata = log2(tpm_network + 1)

#sva: two models
#null model: includes terms for all of the adjustment variables (age, gender, BMI, ischemic time), but not the variables of interest (e.g. T2DM)
#full model: includes both the adjustment variables and the variables of interest
#read variables (metadata)
subject_trait <- read.delim("./Final331MHT1DMHT2Dfilter15AttributeData_ModifyID.txt")
subject_trait = subject_trait[, c(2, 4, 5, 6, 12, 15)]
rownames(subject_trait) = subject_trait[, 1]
subject_trait = subject_trait[, -1]
subject_trait$SubID = gsub('GTEX[.]([^.]+)[.].*', '\\1', rownames(subject_trait))
subject_trait = subset(subject_trait, subject_trait$SubID %in% final_75_150$SubID)
GRU = read.delim("./phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")

sub_ischemic_time = GRU[, c("SUBJID", "TRISCHD")]
row.names(sub_ischemic_time) = sub_ischemic_time$SUBJID
sub_ischemic_time$SUBJID = NULL
sub_ischemic_time = subset(sub_ischemic_time, rownames(sub_ischemic_time) %in% final_75_150$SUBJID)

sample_trait <- read.delim("./GTEx_v7_Annotations_SampleAttributesDS.txt")
sample_trait$SAMPID = gsub('\\-', '.', sample_trait$SAMPID)
sample_trait = subset(sample_trait, sample_trait$SAMPID %in% rownames(subject_trait))
row.names(sample_trait) = sample_trait$SAMPID
sample_trait$SAMPID = NULL
pheno = cbind(subject_trait, sub_ischemic_time$TRISCHD)
colnames(pheno)[1] = "Gender"
colnames(pheno)[2] = "Age"
colnames(pheno)[5] = "T2DM"
colnames(pheno)[7] = "Ischemic_Time"
sample_subject_IDs = pheno$SubID
pheno$SubID = NULL
#fill in mean values when NA
mean(pheno$Ischemic_Time, na.rm = TRUE) #528.8267
pheno$Ischemic_Time[is.na(pheno$SMTSISC)] = 528.8267
pheno = pheno[ , c(1,2,4,6,5)]

mod = model.matrix(~as.factor(T2DM) + as.factor(Gender) + as.numeric(BMI) + as.numeric(AGE) + as.numeric(Ischemic_Time), data = pheno)
mod0 = model.matrix(~as.factor(Gender) + as.numeric(BMI) + as.numeric(AGE) + as.numeric(Ischemic_Time), data = pheno)
svobj = sva(as.matrix(edata), mod, mod0, method="irw")
#Number of significant surrogate variables is:  23, pick top 3.

#batch correction
cleaningY = function(y, mod, svobj) {
  X=cbind(mod,svobj$sv[, 1:3])
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  P=ncol(mod)
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

tpm_clean = cleaningY(edata, mod = mod, svobj = svobj)
#WGCNA expression value matrix done.
