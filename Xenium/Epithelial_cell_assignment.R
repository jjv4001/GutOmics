
#Load object with only epithelial cells from different regions from both Patients 
epithelial<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/mergedepithelial.rds")

###############################################################################
# ASSIGN EPITHELIAL CELL TYPES BY MODULE SCORE
# For a Xenium / Seurat object called: epithelial
#
# Cell types:
# - Stemcell
# - TAcells
# - Enterocytes
# - Gobletcells
# - EECs
#
# Output metadata columns:
# - Stemcell_Score1, TAcells_Score1, ...
# - celltype_score
# - max_score
# - second_score
# - score_margin
# - celltype_score_filtered
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

## ---------------------------
## 1. Set assay and inspect genes
## ---------------------------
DefaultAssay(epithelial) <- DefaultAssay(epithelial)  # keep current assay
cat("Default assay:", DefaultAssay(epithelial), "\n")
cat("Number of genes:", nrow(epithelial), "\n")
cat("Number of cells:", ncol(epithelial), "\n")

## ---------------------------
## 2. Define marker sets
## ---------------------------
Stemcells_markers   <- c("LGR5", "PLAGL1", "MECOM")
TAcells_markers     <- c("UBE2C","MKI67")
Enterocytes_markers <- c("SLC26A3", "RARA")
Gobletcells_markers <- c("CLCA1", "IRF6", "TFCP2")
EECs_markers        <- c("CHGA", "ATF3", "ONECUT2", "FOXJ3")

marker_list <- list(
  Stemcell    = Stemcells_markers,
  TAcells     = TAcells_markers,
  Enterocytes = Enterocytes_markers,
  Gobletcells = Gobletcells_markers,
  EECs        = EECs_markers
)

## ---------------------------
## 3. Keep only genes present in object
## ---------------------------
marker_list_present <- lapply(marker_list, function(x) intersect(x, rownames(epithelial)))

cat("\nMarkers found in object:\n")
for (nm in names(marker_list_present)) {
  cat(nm, ":", paste(marker_list_present[[nm]], collapse = ", "), "\n")
}

## remove empty signatures
marker_list_present <- marker_list_present[sapply(marker_list_present, length) > 0]

if (length(marker_list_present) == 0) {
  stop("None of the marker genes were found in the epithelial object.")
}

## ---------------------------
## 4. Remove any old score columns from previous runs
## ---------------------------
old_score_cols <- grep("_Score[0-9]+$|celltype_score|max_score|second_score|score_margin|celltype_score_filtered",
                       colnames(epithelial@meta.data), value = TRUE)

if (length(old_score_cols) > 0) {
  epithelial@meta.data <- epithelial@meta.data[, setdiff(colnames(epithelial@meta.data), old_score_cols), drop = FALSE]
}

## ---------------------------
## 5. Run AddModuleScore for each cell type
## Xenium panels are small, so keep ctrl low
## ---------------------------
for (nm in names(marker_list_present)) {
  epithelial <- AddModuleScore(
    object = epithelial,
    features = list(marker_list_present[[nm]]),
    name = paste0(nm, "_Score"),
    ctrl = 5,
    search = FALSE
  )
}

## ---------------------------
## 6. Collect score columns automatically
## ---------------------------
score_cols <- paste0(names(marker_list_present), "_Score1")
score_cols <- intersect(score_cols, colnames(epithelial@meta.data))

if (length(score_cols) == 0) {
  stop("No module score columns were created.")
}

score_mat <- epithelial@meta.data[, score_cols, drop = FALSE]

## ---------------------------
## 7. Assign highest-scoring cell type
## ---------------------------
max_idx <- max.col(score_mat, ties.method = "first")

epithelial$celltype_score <- gsub(
  "_Score1$",
  "",
  colnames(score_mat)[max_idx]
)

## ---------------------------
## 8. Add confidence metrics
## ---------------------------
epithelial$max_score <- apply(score_mat, 1, max, na.rm = TRUE)

get_second_highest <- function(x) {
  x <- sort(as.numeric(x), decreasing = TRUE)
  if (length(x) < 2) return(NA_real_)
  x[2]
}

epithelial$second_score <- apply(score_mat, 1, get_second_highest)
epithelial$score_margin <- epithelial$max_score - epithelial$second_score

## ---------------------------
## 9. Filter ambiguous / weak assignments
##
## Adjust thresholds as needed:
## - max_score < 0.05  --> too weak
## - score_margin < 0.02 --> too ambiguous
## ---------------------------
epithelial$celltype_score_filtered <- ifelse(
  epithelial$max_score < 0.05 | epithelial$score_margin < 0.02,
  "Unassigned",
  epithelial$celltype_score
)

epithelial$celltype_score_filtered <- factor(
  epithelial$celltype_score_filtered,
  levels = c("Stemcell", "TAcells", "Enterocytes", "Gobletcells", "EECs", "Unassigned")
)

## ---------------------------
## 10. Summaries
## ---------------------------
cat("\nRaw assigned cell types:\n")
print(table(epithelial$celltype_score, useNA = "ifany"))

cat("\nFiltered assigned cell types:\n")
print(table(epithelial$celltype_score_filtered, useNA = "ifany"))

## ---------------------------
## 11. Inspect average scores by assigned class
## ---------------------------
avg_scores <- epithelial@meta.data %>%
  group_by(celltype_score_filtered) %>%
  summarise(across(all_of(score_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")

print(avg_scores)

## ---------------------------
## 12. Plot counts of assigned cell types
## ---------------------------
ggplot(epithelial@meta.data, aes(x = celltype_score_filtered)) +
  geom_bar() +
  theme_classic() +
  labs(x = "Assigned cell type", y = "Number of cells") +
  coord_flip()

## ---------------------------
## 13. Violin plots of module scores
## ---------------------------
VlnPlot(
  epithelial,
  features = score_cols,
  group.by = "celltype_score_filtered",
  pt.size = 0
)

## ---------------------------
## 14. UMAP / dimensional reduction plot if available
## ---------------------------
if ("umap" %in% names(epithelial@reductions)) {
  print(
    DimPlot(
      epithelial,
      reduction = "umap",
      group.by = "celltype_score_filtered",
      label = TRUE
    )
  )
}

## ---------------------------
## 15. Spatial plots for Xenium if image exists
## ---------------------------
## Uncomment if this is a Xenium object with image data
# ImageDimPlot(
#   epithelial,
#   group.by = "celltype_score_filtered",
#   size = 0.6
# )

## ---------------------------
## 16. Feature sanity checks
## ---------------------------
present_check_genes <- intersect(
  c("LGR5", "PLAGL1", "MECOM", "UBE2C", "SLC26A3", "RARA",
    "CLCA1", "IRF6", "TFCP2", "CHGA", "ATF3", "ONECUT2", "FOXJ3"),
  rownames(epithelial)
)

if (length(present_check_genes) > 0) {
  print(
    DotPlot(
      epithelial,
      features = present_check_genes,
      group.by = "celltype_score_filtered"
    ) + RotatedAxis()
  )
}

new_levels <- c("Unknown","EECs","Goblet cells","Enterocytes","TA cells","Stem cells")
Idents(epithelial) <- factor(Idents(epithelial), levels = new_levels)
DotPlot(epithelial, features=features, dot.scale=20)+RotatedAxis()+scale_color_viridis_c(option = "A", direction=-1)
new_levels <- c("Stem cells","TA cells","Enterocytes","Goblet cells","EECs","Unknown")
Idents(epithelial) <- factor(Idents(epithelial), levels = new_levels)

#save Seurat analysis
saveRDS(epithelial, "/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/newgutassignedmarker.rds")

###############################################################################
# OPTIONAL: save metadata table
###############################################################################
# write.csv(epithelial@meta.data, "epithelial_celltype_assignments.csv", row.names = TRUE)


