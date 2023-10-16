setwd("/Users/erinconnolly/Downloads/")
devtools::install_github("sqjin/CellChat")

library(CellChat)
library(patchwork)

Idents(Control_Eileen_scdata_postQC_MERGED_Tumor_IMMUNE) <-"celltype_broad2"

     
data.input <- GetAssayData(Control_Eileen_scdata_postQC_MERGED_Tumor_IMMUNE, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Control_Eileen_scdata_postQC_MERGED_Tumor_IMMUNE)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.human
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data


cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cell.type.fine")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
#set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#optional
cellchat <- projectData(cellchat, PPI.human)


df.net.control <- subsetCommunication(cellchat)
write.csv(df.net.control, file = "df.net.control.csv", quote = F)


cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- aggregateNet(cellchat)
netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, signaling = c("TGFb"), remove.isolate = FALSE, font.size = 20)
netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL, signaling = c("TGFb", "BMP", "GDF", "ACTIVIN"), lab.cex = 1)

cellchat <- computeNetSimilarity(cellchat, type = "functional",slot.name = "net",)
cellchat <- netEmbedding(cellchat, type = "functional",slot.name = "net")
?computeNetSimilarity()

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat@netP$pathways

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use[["interaction"]][["pathway_name"]]
pathways.show <- c("WNT") 

pathways.show <- c("BAFF")
pathways.show <- c("CXCL") 

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,5) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

#ORDER
levels(cellchat@idents)
vertex.receiver = seq(1,4)

netVisual_bubble(cellchat, sources.use = c(28), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(30,3,10,22,28,20,15,17,18,9,11,12), targets.use = c(6), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(10), targets.use = c(1:30), remove.isolate = FALSE)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2





# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = TRUE)

par(mfrow=c(1,1))
?netVisual_heatmap()
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = c("#2166ac", "#b2182b"))
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,9) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)


netVisual_individuaL(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",top=.01)

#CHORD
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", top=10)
#> Note: The first link end is drawn out of sector 'Inflam. FIB'

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(7) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

pathways.show.all <- cellchat@netP$pathways

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 12, targets.use = c(0:13), remove.isolate = FALSE)
#> Comparing communications on a single object



# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = c(0:13), targets.use = 12, signaling = c('BAFF'), remove.isolate = FALSE)
#> Comparing communications on a single object


netVisual_chord_gene(cellchat, sources.use = c(1,2), targets.use = 12, legend.pos.x = 15)

plotGeneExpression(cellchat, signaling = "TNFRSF17")
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

cellchat <- computeNetSimilarity(cellchat, type = "functional")
saveRDS(cellchat, file = "cellchat_BLOOD_HEALTHY.rds")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

plotGeneExpression(cellchat, signaling = "BAFF")


nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)


netAnalysis_river(cellchat, pattern = "outgoing")

netAnalysis_dot(cellchat, pattern = "incoming")
cellchat <- computeNetSimilarity(cellchat, type = "functional")



##########################################
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)

saveRDS(cellchat,'cellchat_CONTROL.RDS')
