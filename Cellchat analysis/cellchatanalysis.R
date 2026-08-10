
obj<-readRDS("mergedallcellsfetal.rds")
# =========================================================
# CellChat analysis: broad compartments
# epithelial, mesenchymal, immune, neuronal, RBC, endothelial
# =========================================================

library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)
library(future)

options(future.globals.maxSize = 64 * 1024^3)
plan("multisession", workers = 4)

# -----------------------------
# INPUT
# -----------------------------

# metadata column with your fine cell annotations
celltype_col <- "Cell_Type"  # change if needed

outdir <- "CellChat_broad_celltypes"
dir.create(outdir, showWarnings = FALSE)

# -----------------------------
# Make broad celltype labels
# -----------------------------
obj$cellchat_group <- case_when(
  obj@meta.data[[celltype_col]] %in% c(
    "Epithelial", "Enterocyte", "Stem Cells", "Goblet Cells", "EECs", "TA"
  ) ~ "Epithelial",
  
  obj@meta.data[[celltype_col]] %in% c(
    "Mesenchymal", "Fibroblast", "FibroblastBM", "Telocyte",
    "Telocyte_tip", "Smooth_muscle", "Pericyte"
  ) ~ "Mesenchymal",
  
  obj@meta.data[[celltype_col]] %in% c(
    "Immune cells", "Tcell", "Tcell_CD4", "Tcell_CD8", "Bcell",
    "Plasma_cell", "Macrophage", "cDC1", "cDC2", "Mast"
  ) ~ "Immune",
  
  obj@meta.data[[celltype_col]] %in% c(
    "Neuronal", "Neuron", "Enteric_neuron", "Glia", "Schwann"
  ) ~ "Neuronal",
  
  obj@meta.data[[celltype_col]] %in% c(
    "RBCs", "Erythrocyte", "Red_blood_cell"
  ) ~ "RBC",
  
  obj@meta.data[[celltype_col]] %in% c(
    "Endothelial", "EndoArt", "EndoVenous", "EndoLymph", "Endo"
  ) ~ "Endothelial",
  
  TRUE ~ NA_character_
)
obj<-JoinLayers(obj)
obj_sub <- subset(
  obj,
  subset = !is.na(cellchat_group)
)

obj_sub$cellchat_group <- factor(
  obj_sub$cellchat_group,
  levels = c(
    "Epithelial",
    "Mesenchymal",
    "Immune",
    "Neuronal",
    "RBC",
    "Endothelial"
  )
)

table(obj_sub$cellchat_group)

# Optional: downsample if very large
set.seed(123)
obj_sub <- subset(obj, downsample = 1000)

# -----------------------------
# Create CellChat object
# -----------------------------
saveRDS(obj_sub, "obj_sub.rds")
data.input <- GetAssayData(obj_sub, assay = "RNA", slot = "data")
meta <- obj_sub@meta.data
data.input <- GetAssayData(obj_sub, assay = "RNA", slot = "data")

meta <- obj_sub@meta.data
meta <- meta[colnames(data.input), , drop = FALSE]

identical(rownames(meta), colnames(data.input))
meta$cellchat_group <- droplevels(factor(meta$cellchat_group))

# remove any NA groups just in case
keep <- !is.na(meta$cellchat_group)

meta <- meta[keep, , drop = FALSE]
data.input <- data.input[, rownames(meta), drop = FALSE]

stopifnot(identical(colnames(data.input), rownames(meta)))
table(meta$cellchat_group)
cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cellchat_group"
)

# Human database
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# -----------------------------
# Run CellChat
# -----------------------------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)

saveRDS(cellchat, file.path(outdir, "cellchat_broad_groups.rds"))

# -----------------------------
# Export interaction tables
# -----------------------------
df.net <- subsetCommunication(cellchat)
write.csv(
  df.net,
  file.path(outdir, "all_ligand_receptor_interactions.csv"),
  row.names = FALSE
)

df.pathway <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(
  df.pathway,
  file.path(outdir, "pathway_level_interactions.csv"),
  row.names = FALSE
)

# -----------------------------
# Circle plots: number and strength
# -----------------------------
groupSize <- as.numeric(table(cellchat@idents))

pdf(file.path(outdir, "circle_number_of_interactions.pdf"), width = 7, height = 7)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
dev.off()

pdf(file.path(outdir, "circle_interaction_strength.pdf"), width = 7, height = 7)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength"
)
dev.off()

# -----------------------------
# Heatmaps
# -----------------------------
pdf(file.path(outdir, "heatmap_number_of_interactions.pdf"), width = 7, height = 6)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

pdf(file.path(outdir, "heatmap_interaction_strength.pdf"), width = 7, height = 6)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# -----------------------------
# Incoming / outgoing signaling roles
# -----------------------------
pdf(file.path(outdir, "signaling_role_scatter.pdf"), width = 8, height = 6)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

pdf(file.path(outdir, "signaling_role_heatmap_outgoing.pdf"), width = 8, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()

pdf(file.path(outdir, "signaling_role_heatmap_incoming.pdf"), width = 8, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

# -----------------------------
# Pathway-specific plots
# -----------------------------
pathways.show <- cellchat@netP$pathways
writeLines(pathways.show, file.path(outdir, "detected_pathways.txt"))

for (pw in pathways.show) {
  message("Plotting pathway: ", pw)
  
  pdf(file.path(outdir, paste0("pathway_", pw, "_circle.pdf")), width = 7, height = 7)
  print(
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      layout = "circle"
    )
  )
  dev.off()
  
  pdf(file.path(outdir, paste0("pathway_", pw, "_chord.pdf")), width = 8, height = 8)
  print(
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      layout = "chord"
    )
  )
  dev.off()
}

# -----------------------------
# Bubble plot: all interactions between broad groups
# -----------------------------
pdf(file.path(outdir, "bubble_all_interactions.pdf"), width = 12, height = 8)
netVisual_bubble(
  cellchat,
  sources.use = 1:6,
  targets.use = 1:6,
  remove.isolate = FALSE
)
dev.off()


plotdir <- file.path(outdir, "plots")
dir.create(plotdir, showWarnings = FALSE)

groups <- levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
pathways.show <- cellchat@netP$pathways

# =========================
# Global network plots
# =========================

pdf(file.path(plotdir, "01_circle_count.pdf"), 7, 7)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
dev.off()

pdf(file.path(plotdir, "02_circle_weight.pdf"), 7, 7)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction strength")
dev.off()

pdf(file.path(plotdir, "03_heatmap_count.pdf"), 7, 6)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

pdf(file.path(plotdir, "04_heatmap_weight.pdf"), 7, 6)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

pdf(file.path(plotdir, "05_bubble_all.pdf"), 12, 8)
netVisual_bubble(cellchat, sources.use = groups, targets.use = groups,
                 remove.isolate = FALSE)
dev.off()

pdf(file.path(plotdir, "06_signaling_role_scatter.pdf"), 8, 6)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

pdf(file.path(plotdir, "07_outgoing_role_heatmap.pdf"), 9, 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()

pdf(file.path(plotdir, "08_incoming_role_heatmap.pdf"), 9, 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

pdf(file.path(plotdir, "09_alluvial_outgoing.pdf"), 10, 7)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf(file.path(plotdir, "10_dot_outgoing.pdf"), 10, 7)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# =========================
# Sender-specific circle plots
# =========================

for (i in seq_along(groups)) {
  mat <- matrix(0, nrow = nrow(cellchat@net$count), ncol = ncol(cellchat@net$count),
                dimnames = dimnames(cellchat@net$count))
  mat[i, ] <- cellchat@net$count[i, ]
  
  pdf(file.path(plotdir, paste0("sender_", groups[i], "_count_circle.pdf")), 7, 7)
  netVisual_circle(mat, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste(groups[i], "outgoing interactions"))
  dev.off()
}

# =========================
# Receiver-specific circle plots
# =========================

for (i in seq_along(groups)) {
  mat <- matrix(0, nrow = nrow(cellchat@net$count), ncol = ncol(cellchat@net$count),
                dimnames = dimnames(cellchat@net$count))
  mat[, i] <- cellchat@net$count[, i]
  
  pdf(file.path(plotdir, paste0("receiver_", groups[i], "_count_circle.pdf")), 7, 7)
  netVisual_circle(mat, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste(groups[i], "incoming interactions"))
  dev.off()
}

# =========================
# Pathway-level plots
# =========================

writeLines(pathways.show, file.path(plotdir, "detected_pathways.txt"))

for (pw in pathways.show) {
  message("Plotting pathway: ", pw)
  
  safe_pw <- gsub("[^A-Za-z0-9_]", "_", pw)
  
  pdf(file.path(plotdir, paste0("pathway_", safe_pw, "_circle.pdf")), 7, 7)
  print(netVisual_aggregate(cellchat, signaling = pw, layout = "circle"))
  dev.off()
  
  pdf(file.path(plotdir, paste0("pathway_", safe_pw, "_chord.pdf")), 8, 8)
  print(netVisual_aggregate(cellchat, signaling = pw, layout = "chord"))
  dev.off()
  
  pdf(file.path(plotdir, paste0("pathway_", safe_pw, "_heatmap.pdf")), 7, 6)
  print(netVisual_heatmap(cellchat, signaling = pw))
  dev.off()
  
  pdf(file.path(plotdir, paste0("pathway_", safe_pw, "_bubble.pdf")), 10, 7)
  print(netVisual_bubble(cellchat, signaling = pw, remove.isolate = FALSE))
  dev.off()
}

# =========================
# Ligand-receptor pair plots
# =========================

df.net <- subsetCommunication(cellchat)
write.csv(df.net, file.path(plotdir, "all_ligand_receptor_interactions.csv"),
          row.names = FALSE)

top_lr <- df.net %>%
  arrange(desc(prob)) %>%
  distinct(interaction_name, .keep_all = TRUE) %>%
  head(50)

for (lr in top_lr$interaction_name) {
  safe_lr <- gsub("[^A-Za-z0-9_]", "_", lr)
  
  pdf(file.path(plotdir, paste0("LR_", safe_lr, "_bubble.pdf")), 10, 7)
  print(netVisual_bubble(cellchat, pairLR.use = lr, remove.isolate = FALSE))
  dev.off()
}

# =========================
# Contribution plots per pathway
# =========================

for (pw in pathways.show) {
  safe_pw <- gsub("[^A-Za-z0-9_]", "_", pw)
  
  pdf(file.path(plotdir, paste0("pathway_", safe_pw, "_contribution.pdf")), 9, 6)
  print(netAnalysis_contribution(cellchat, signaling = pw))
  dev.off()
}

# =========================
# Save summary tables
# =========================

count_mat <- cellchat@net$count
weight_mat <- cellchat@net$weight

write.csv(count_mat, file.path(plotdir, "interaction_count_matrix.csv"))
write.csv(weight_mat, file.path(plotdir, "interaction_weight_matrix.csv"))

outgoing_strength <- rowSums(weight_mat)
incoming_strength <- colSums(weight_mat)

summary_df <- data.frame(
  celltype = names(outgoing_strength),
  outgoing_strength = outgoing_strength,
  incoming_strength = incoming_strength,
  outgoing_count = rowSums(count_mat),
  incoming_count = colSums(count_mat)
)

write.csv(summary_df, file.path(plotdir, "cellchat_group_summary.csv"),
          row.names = FALSE)

saveRDS(cellchat, file.path(outdir, "cellchat_with_all_plots.rds"))
# -----------------------------
# Save object again
# -----------------------------
saveRDS(cellchat, file.path(outdir, "cellchat_broad_groups_final.rds"))

library(dplyr)
library(ggplot2)

plotdir <- file.path(outdir, "important_top_pathway_plots")
dir.create(plotdir, showWarnings = FALSE)

groups <- levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

# =========================================================
# 1. Rank pathways by total communication strength
# =========================================================

df.pathway <- subsetCommunication(cellchat, slot.name = "netP")

top.pathways <- df.pathway %>%
  group_by(pathway_name) %>%
  summarise(
    total_prob = sum(prob, na.rm = TRUE),
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_prob))

write.csv(
  top.pathways,
  file.path(plotdir, "ranked_pathways_by_total_strength.csv"),
  row.names = FALSE
)

topN <- min(20, nrow(top.pathways))
pathways.use <- top.pathways$pathway_name[1:topN]

# Barplot of top pathways
p <- ggplot(top.pathways[1:topN, ],
            aes(x = reorder(pathway_name, total_prob), y = total_prob)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    x = "Pathway",
    y = "Total communication probability",
    title = "Top CellChat signaling pathways"
  )

ggsave(
  file.path(plotdir, "top_pathways_ranked_barplot.pdf"),
  p,
  width = 7,
  height = 6
)

# =========================================================
# 2. Global important plots
# =========================================================

pdf(file.path(plotdir, "global_circle_number_interactions.pdf"), 7, 7)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
dev.off()

pdf(file.path(plotdir, "global_circle_interaction_strength.pdf"), 7, 7)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength"
)
dev.off()

pdf(file.path(plotdir, "global_heatmap_number_interactions.pdf"), 7, 6)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

pdf(file.path(plotdir, "global_heatmap_interaction_strength.pdf"), 7, 6)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

pdf(file.path(plotdir, "global_bubble_all_interactions.pdf"), 12, 8)
netVisual_bubble(
  cellchat,
  sources.use = groups,
  targets.use = groups,
  remove.isolate = FALSE
)
dev.off()

# =========================================================
# 3. Top pathway circle / chord / heatmap / bubble
# =========================================================

for (pw in pathways.use) {
  message("Plotting top pathway: ", pw)
  safe_pw <- gsub("[^A-Za-z0-9_]", "_", pw)
  
  pdf(file.path(plotdir, paste0("top_pathway_", safe_pw, "_circle.pdf")), 7, 7)
  print(
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      layout = "circle"
    )
  )
  dev.off()
  
  pdf(file.path(plotdir, paste0("top_pathway_", safe_pw, "_chord.pdf")), 8, 8)
  print(
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      layout = "chord"
    )
  )
  dev.off()
  
  pdf(file.path(plotdir, paste0("top_pathway_", safe_pw, "_heatmap.pdf")), 7, 6)
  print(
    netVisual_heatmap(
      cellchat,
      signaling = pw
    )
  )
  dev.off()
  
  pdf(file.path(plotdir, paste0("top_pathway_", safe_pw, "_bubble.pdf")), 10, 7)
  print(
    netVisual_bubble(
      cellchat,
      signaling = pw,
      remove.isolate = FALSE
    )
  )
  dev.off()
  
  pdf(file.path(plotdir, paste0("top_pathway_", safe_pw, "_contribution.pdf")), 8, 6)
  print(
    netAnalysis_contribution(
      cellchat,
      signaling = pw
    )
  )
  dev.off()
}

# =========================================================
# 4. Cell-type signaling role plots
# =========================================================

cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)

pdf(file.path(plotdir, "signaling_role_scatter.pdf"), 8, 6)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

pdf(file.path(plotdir, "outgoing_signaling_role_heatmap.pdf"), 9, 7)
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing"
)
dev.off()

pdf(file.path(plotdir, "incoming_signaling_role_heatmap.pdf"), 9, 7)
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming"
)
dev.off()

# =========================================================
# 5. Source-target specific important plots
# =========================================================

# Epithelial outgoing
pdf(file.path(plotdir, "epithelial_outgoing_bubble.pdf"), 10, 7)
netVisual_bubble(
  cellchat,
  sources.use = "Epithelial",
  targets.use = groups,
  remove.isolate = FALSE
)
dev.off()

# Epithelial incoming
pdf(file.path(plotdir, "epithelial_incoming_bubble.pdf"), 10, 7)
netVisual_bubble(
  cellchat,
  sources.use = groups,
  targets.use = "Epithelial",
  remove.isolate = FALSE
)
dev.off()

# Mesenchymal to epithelial
pdf(file.path(plotdir, "mesenchymal_to_epithelial_bubble.pdf"), 10, 7)
netVisual_bubble(
  cellchat,
  sources.use = "Mesenchymal",
  targets.use = "Epithelial",
  remove.isolate = FALSE
)
dev.off()

# Epithelial to mesenchymal
pdf(file.path(plotdir, "epithelial_to_mesenchymal_bubble.pdf"), 10, 7)
netVisual_bubble(
  cellchat,
  sources.use = "Epithelial",
  targets.use = "Mesenchymal",
  remove.isolate = FALSE
)
dev.off()

# Immune to epithelial
pdf(file.path(plotdir, "immune_to_epithelial_bubble.pdf"), 10, 7)
netVisual_bubble(
  cellchat,
  sources.use = "Immune",
  targets.use = "Epithelial",
  remove.isolate = FALSE
)
dev.off()

# Endothelial to epithelial
pdf(file.path(plotdir, "endothelial_to_epithelial_bubble.pdf"), 10, 7)
netVisual_bubble(
  cellchat,
  sources.use = "Endothelial",
  targets.use = "Epithelial",
  remove.isolate = FALSE
)
dev.off()

# =========================================================
# 6. Top ligand-receptor pairs overall
# =========================================================

df.lr <- subsetCommunication(cellchat)

top.lr <- df.lr %>%
  arrange(desc(prob)) %>%
  distinct(interaction_name, .keep_all = TRUE) %>%
  head(30)

write.csv(
  top.lr,
  file.path(plotdir, "top_30_ligand_receptor_pairs.csv"),
  row.names = FALSE
)

for (lr in top.lr$interaction_name) {
  safe_lr <- gsub("[^A-Za-z0-9_]", "_", lr)
  
  pdf(file.path(plotdir, paste0("top_LR_", safe_lr, "_bubble.pdf")), 10, 7)
  print(
    netVisual_bubble(
      cellchat,
      pairLR.use = lr,
      remove.isolate = FALSE
    )
  )
  dev.off()
}

# =========================================================
# 7. Save summary matrices
# =========================================================

write.csv(
  cellchat@net$count,
  file.path(plotdir, "cellchat_interaction_count_matrix.csv")
)

write.csv(
  cellchat@net$weight,
  file.path(plotdir, "cellchat_interaction_weight_matrix.csv")
)

summary.df <- data.frame(
  celltype = groups,
  outgoing_count = rowSums(cellchat@net$count),
  incoming_count = colSums(cellchat@net$count),
  outgoing_strength = rowSums(cellchat@net$weight),
  incoming_strength = colSums(cellchat@net$weight)
)

write.csv(
  summary.df,
  file.path(plotdir, "celltype_incoming_outgoing_summary.csv"),
  row.names = FALSE
)

top5.pathways <- top.pathways$pathway_name[1:5]

pdf(file.path(plotdir, "global_bubble_top5pathways.pdf"), 10, 8)

netVisual_bubble(
  cellchat,
  signaling = top5.pathways,
  remove.isolate = TRUE
)

dev.off()

pdf(file.path(plotdir, "epithelial_incoming_bubble.pdf"),
    width = 11,
    height = 9)

p <- netVisual_bubble(
  cellchat,
  sources.use = groups,
  targets.use = "Epithelial",
  remove.isolate = TRUE,
  angle.x = 45
) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(
      size = 6,          # reduce LR pair font
      face = "italic"
    ),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 25)
  ) +
  ggtitle("Incoming signaling to epithelial cells")

print(p)

dev.off()

pdf(file.path(plotdir, "epithelial_outgoing_bubble.pdf"),
    width = 11,
    height = 9)

p <- netVisual_bubble(
  cellchat,
  sources.use = "Epithelial",
  targets.use = groups,
  remove.isolate = TRUE,
  top = 30,            # show only top 30 interactions
  angle.x = 45
) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(
      size = 6,        # smaller LR labels
      face = "italic"
    ),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 25)
  ) +
  ggtitle("Outgoing signaling from epithelial cells")

print(p)

dev.off()