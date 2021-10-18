library(igraph)

#generate T2DM and non-T2DM networks separately
tpm_clean_control = tpm_clean[, colnames(tpm_clean) %in% colnames(sample_control)]
tpm_clean_type2 = tpm_clean[, colnames(tpm_clean) %in% colnames(sample_type2)]
datExprNonT2DM = data.frame(t(tpm_clean_control))
datExprT2DM = data.frame(t(tpm_clean_type2))
TOM_NonT2DM = TOMsimilarityFromExpr(datExprNonT2DM, power = 7)
TOM_T2DM = TOMsimilarityFromExpr(datExprT2DM, power = 7)
TOM = TOMsimilarityFromExpr(datExpr, power = 7)

inModule5 = (moduleColors == module5)
inModule6 = (moduleColors == module6)
TOM_NonT2DM_lightgreen = TOM_NonT2DM[inModule5, inModule5]
TOM_T2DM_lightgreen = TOM_T2DM[inModule5, inModule5]
TOM_lightgreen = TOM[inModule5, inModule5]
dimnames(TOM_NonT2DM_lightgreen) = list(intGenes5_symbol, intGenes5_symbol)
dimnames(TOM_T2DM_lightgreen) = list(intGenes5_symbol, intGenes5_symbol)
dimnames(TOM_lightgreen) = list(intGenes5_symbol, intGenes5_symbol)

TOM_NonT2DM_magenta = TOM_NonT2DM[inModule6, inModule6]
TOM_T2DM_magenta = TOM_T2DM[inModule6, inModule6]
TOM_magenta = TOM[inModule6, inModule6]
dimnames(TOM_NonT2DM_magenta) = list(intGenes6_symbol, intGenes6_symbol)
dimnames(TOM_T2DM_magenta) = list(intGenes6_symbol, intGenes6_symbol)
dimnames(TOM_magenta) = list(intGenes6_symbol, intGenes6_symbol)

TOM_T2DM_lightgreen = 1/TOM_T2DM_lightgreen
TOM_NonT2DM_lightgreen = 1/TOM_NonT2DM_lightgreen
TOM_lightgreen = 1/TOM_lightgreen
T2DM_lightgreen_graph = graph_from_adjacency_matrix(TOM_T2DM_lightgreen, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
NonT2DM_lightgreen_graph = graph_from_adjacency_matrix(TOM_NonT2DM_lightgreen, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
lightgreen_graph = graph_from_adjacency_matrix(TOM_lightgreen, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
T2DM_lightgreen_bet = betweenness(T2DM_lightgreen_graph, v = V(T2DM_lightgreen_graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)
NonT2DM_lightgreen_bet = betweenness(NonT2DM_lightgreen_graph, v = V(NonT2DM_lightgreen_graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)
lightgreen_bet = betweenness(lightgreen_graph, v = V(lightgreen_graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)

T2DM_lightgreen_bet = data.frame(keyName = names(T2DM_lightgreen_bet), value = T2DM_lightgreen_bet, row.names = NULL)
row.names(T2DM_lightgreen_bet) = T2DM_lightgreen_bet$keyName
T2DM_lightgreen_bet$keyName = NULL
colnames(T2DM_lightgreen_bet) = "T2DM_BC"

NonT2DM_lightgreen_bet = data.frame(keyName = names(NonT2DM_lightgreen_bet), value = NonT2DM_lightgreen_bet, row.names = NULL)
row.names(NonT2DM_lightgreen_bet) = NonT2DM_lightgreen_bet$keyName
NonT2DM_lightgreen_bet$keyName = NULL
colnames(NonT2DM_lightgreen_bet) = "NonT2DM_BC"

lightgreen_bet = data.frame(keyName = names(lightgreen_bet), value = lightgreen_bet, row.names = NULL)
row.names(lightgreen_bet) = lightgreen_bet$keyName
lightgreen_bet$keyName = NULL
colnames(lightgreen_bet) = "lightgreen_BC"

lightgreenBC = cbind(NonT2DM_lightgreen_bet, T2DM_lightgreen_bet, lightgreen_bet)
lightgreenBC = lightgreenBC[with(lightgreenBC, order(-lightgreenBC$NonT2DM_BC, -lightgreenBC$T2DM_BC)), ]
write.csv(lightgreenBC, file = "./lightgreenBC.csv")

TOM_T2DM_magenta = 1/TOM_T2DM_magenta
TOM_NonT2DM_magenta = 1/TOM_NonT2DM_magenta
TOM_magenta = 1/TOM_magenta
T2DM_magenta_graph = graph_from_adjacency_matrix(TOM_T2DM_magenta, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
NonT2DM_magenta_graph = graph_from_adjacency_matrix(TOM_NonT2DM_magenta, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
magenta_graph = graph_from_adjacency_matrix(TOM_magenta, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
T2DM_magenta_bet = betweenness(T2DM_magenta_graph, v = V(T2DM_magenta_graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)
NonT2DM_magenta_bet = betweenness(NonT2DM_magenta_graph, v = V(NonT2DM_magenta_graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)
magenta_bet = betweenness(magenta_graph, v = V(magenta_graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)

T2DM_magenta_bet = data.frame(keyName = names(T2DM_magenta_bet), value = T2DM_magenta_bet, row.names = NULL)
row.names(T2DM_magenta_bet) = T2DM_magenta_bet$keyName
T2DM_magenta_bet$keyName = NULL
colnames(T2DM_magenta_bet) = "T2DM_BC"

NonT2DM_magenta_bet = data.frame(keyName = names(NonT2DM_magenta_bet), value = NonT2DM_magenta_bet, row.names = NULL)
row.names(NonT2DM_magenta_bet) = NonT2DM_magenta_bet$keyName
NonT2DM_magenta_bet$keyName = NULL
colnames(NonT2DM_magenta_bet) = "NonT2DM_BC"

magenta_bet = data.frame(keyName = names(magenta_bet), value = magenta_bet, row.names = NULL)
row.names(magenta_bet) = magenta_bet$keyName
magenta_bet$keyName = NULL
colnames(magenta_bet) = "magenta_BC"

magentaBC = cbind(NonT2DM_magenta_bet, T2DM_magenta_bet, magenta_bet)
magentaBC = magentaBC[with(magentaBC, order(-magentaBC$NonT2DM_BC, -magentaBC$T2DM_BC)), ]
write.csv(magentaBC, file = "./magentaBC.csv")
