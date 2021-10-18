# WGCNA-of-Type-II-Diabetic-Skin
Scripts for paper: Network Analysis Reveals Dysregulated Functional Patterns in Type II Diabetic Skin

#1st step: sample selection
Use 'MatchIt' package to select T2DM and non-T2DM subjects matched for age, sex, and race.
Background factors can be added/adjusted under 'matchit' function.
m.out are the non-T2DM (control) samples whose backgroud is matching with T2DM samples.
Details can be found at: https://cran.r-project.org/web/packages/MatchIt/MatchIt.pdf

#2nd step: RNA-seq preprocessing
Remove low-expressed genes, here we use 'rwoQuantileC' to remove the gene whose raw count is less than 10 in more than 90% samples.
Pick out the genes exhibiting the most variation in expression levels across samples based on MAD (mean absolute deviation) value.
Log transformed TPM expression matrix of left genes.
Use 'sva' package to find out the surrogate variables (SVs): https://www.bioconductor.org/packages/devel/bioc/vignettes/sva/inst/doc/sva.pdf
Generally, build two models (mod0 and mod). mod0 is the null model that contains only known adjustment variables while mod is the full model that contains both adjustment variables and the variable of interest.
The output of 'sva' is SVs, and select top several SVs for batch correction.
The resulting matrix is batch corrected expression matrix.

#3rd step: sample outlier detection
Use standardized connectivity (Z.K method): https://link.springer.com/book/10.1007/978-1-4419-8819-5
In T2DM and non-T2DM group, calculate z.k values for each sample separately.
The threshold for z.k value is often -2 or -2.5.
The sample with z.k out of the threshold is considered outlier and will be removed from the clean expression matrix resulted from 3rd step.
The resulting expression matrix is ready for WGCNA.

#4th step: weghted gene co-expression network analysis
Use 'WGCNA' package to construct network, details can be found at https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
'pickSoftThreshold' is used to find out the best parameter for adjacency power function.
Once the parameter is determined, gene modules can be clustered with 'blockwiseModules'.
Each module has a module eigengen (ME) which represents that module.
The correlation between each ME and disease status (here is T2DM) is calculated and reflected in the heat map.
The modules that are significantly correlated with disease will be selected for intramodular analysis.

#5th step: intramodular analysis
Extract the gene modules of interest.
Calculate MM and GS to identify important genes (hubs).
Transform gene ID for enrichment analysis, here we use DAVID and g:Profiler (website tools).

#6th: betweenness centrality
Generate T2DM and non-T2DM network separately, same method as step 4.
Use 'igraph' package to calculate BC: https://igraph.org/
Note: As BC is defined as the number of shortest paths between every two other nodes in the module that pass through that node, here we use the inverse of topological overlap represents the length of an edge.
The result should be the BC values from two networks of each gene.
#done
