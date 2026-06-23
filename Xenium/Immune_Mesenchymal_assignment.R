#Increase memory limit 
options(future.globals.maxSize = 64 * 1024^3) 
#merge objects from same tissue type but different Patients into same object and save object
foregut169<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/foregut169.rds")
foregut168<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/foregut168.rds")
foregut<-merge(foregut169, y=foregut168)
saveRDS(foregut, "foregut.rds")

midgut169<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/midgut169.rds")
midgut168<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/midgut168.rds")
midgut<-merge(midgut169, y=midgut168)
saveRDS(midgut, "midgut.rds")

colon169<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/colon169.rds")
colon168<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/colon168.rds")
colon<-merge(colon169, y=colon168)
saveRDS(colon, "colon.rds")

#read one tissue type and assign identities
foregut<-readRDS("foregut.rds")
foregut$Identity<-Idents(foregut)

###############################################################################
# ASSIGN MESENCHYMAL AND IMMUNE CELL TYPES BY MODULE SCORE
# For a Xenium / Seurat object called: foregut
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

## ---------------------------
## 1. Set assay and inspect genes
## ---------------------------
DefaultAssay(foregut) <- DefaultAssay(foregut)  # keep current assay
cat("Default assay:", DefaultAssay(foregut), "\n")
cat("Number of genes:", nrow(foregut), "\n")
cat("Number of cells:", ncol(foregut), "\n")

## ---------------------------
## 2. Define marker sets
## ---------------------------
Fibroblast_markers   <- c("ACTA2","THY1","VCAN","TNC","HPLN1")
FibroblastBM_markers   <- c("ACTA2","THY1","VCAN","TNC","MMP14","VIM","FN1","PECAM1","PDPN")
EndoArt_markers   <- c("PECAM1")
EndoLymph_markers   <- c("PECAM1","PDPN")
EndoVenous_markers   <- c("PECAM1","SLC2A1")
Telocyte_markers  <- c("PDGFRA","GAP43")
Telocyte_tip_markers   <- c("PDGFRA","TNR","TNC","GAP43","HAPLN1")
Macrophage_markers   <- c("PTPRC","CD14","CD163","MRC1")
cDC2_markers   <- c("PTPRC","CD14","CD163","MRC1","CD1C","HLA-DR","ITGAX","CD4")
cDC1_markers   <- c("PTPRC","CD14","CD163","MRC1","HLA-DR","ITGAX","CD4")
FDC_markers   <- c("LAMP3","CD8A","PDCD1","CD69","CD19")
Tcell_CD4_markers   <- c("CD4","CD8A","PDCD1","CD69","CD19")
Tcell_CD8_markers  <- c("CD8A","CD69")
Bcell_markers  <- c("ITGAX","CD69","CD19")
Plasma_cell_markers   <- c("CD19","SDC1")
EpiGLUT_markers   <- c("SLC2A1")
EpiNDRG1_markers   <- c("NDRG1")
EpiCA9_markers   <- c("CA9","POSTN","MKI67")


marker_list <- list(
  Fibroblast    = Fibroblast_markers,
  FibroblastBM     = FibroblastBM_markers ,
  EndoArt = EndoArt_markers,
  EndoLymph = EndoLymph_markers,
  EndoVenous = EndoVenous_markers,
  Telocyte = Telocyte_markers ,
  Telocyte_tip=Telocyte_tip_markers,
  Macrophage    = Macrophage_markers,
  cDC2     = cDC2_markers,
  cDC1 = cDC1_markers,
  FDC = FDC_markers,
  Tcell_CD4 = Tcell_CD4_markers,
  Tcell_CD8 = Tcell_CD8_markers ,
  Bcell= Bcell_markers,
  Plasma_cell=Plasma_cell_markers,
  EpiGLUT    = EpiGLUT_markers ,
  EpiNDRG1    = EpiNDRG1_markers  ,
  EpiCA9 = EpiCA9_markers
)

## ---------------------------
## 3. Keep only genes present in object
## ---------------------------
marker_list_present <- lapply(marker_list, function(x) intersect(x, rownames(foregut)))

cat("\nMarkers found in object:\n")
for (nm in names(marker_list_present)) {
  cat(nm, ":", paste(marker_list_present[[nm]], collapse = ", "), "\n")
}

## remove empty signatures
marker_list_present <- marker_list_present[sapply(marker_list_present, length) > 0]

if (length(marker_list_present) == 0) {
  stop("None of the marker genes were found in the object.")
}

## ---------------------------
## 4. Remove any old score columns from previous runs
## ---------------------------
old_score_cols <- grep("_Score[0-9]+$|celltype_score|max_score|second_score|score_margin|celltype_score_filtered",
                       colnames(foregut@meta.data), value = TRUE)

if (length(old_score_cols) > 0) {
  foregut@meta.data <- foregut@meta.data[, setdiff(colnames(foregut@meta.data), old_score_cols), drop = FALSE]
}

## ---------------------------
## 5. Run AddModuleScore for each cell type
## Xenium panels are small, so keep ctrl low
## ---------------------------
for (nm in names(marker_list_present)) {
  foregut <- AddModuleScore(
    object = foregut,
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
score_cols <- intersect(score_cols, colnames(foregut@meta.data))

if (length(score_cols) == 0) {
  stop("No module score columns were created.")
}

score_mat <- foregut@meta.data[, score_cols, drop = FALSE]

## ---------------------------
## 7. Assign highest-scoring cell type
## ---------------------------
max_idx <- max.col(score_mat, ties.method = "first")

foregut$celltype_score <- gsub(
  "_Score1$",
  "",
  colnames(score_mat)[max_idx]
)

## ---------------------------
## 8. Add confidence metrics
## ---------------------------
foregut$max_score <- apply(score_mat, 1, max, na.rm = TRUE)

get_second_highest <- function(x) {
  x <- sort(as.numeric(x), decreasing = TRUE)
  if (length(x) < 2) return(NA_real_)
  x[2]
}

foregut$second_score <- apply(score_mat, 1, get_second_highest)
foregut$score_margin <- foregut$max_score - foregut$second_score

## ---------------------------
## 9. Filter ambiguous / weak assignments
##
## Adjust thresholds as needed:
## - max_score < 0.05  --> too weak
## - score_margin < 0.02 --> too ambiguous
## ---------------------------
foregut$celltype_score_filtered <- ifelse(
  foregut$max_score < 0.05 | foregut$score_margin < 0.02,
  "Unassigned",
  foregut$celltype_score
)

#Save seurat analysis for mesenchymal cell assigned object
saveRDS(foregut, "foregutimmunemesenchymeassigned.rds")

## ---------------------------
## 10. Summaries
## ---------------------------
cat("\nRaw assigned cell types:\n")
print(table(foregut$celltype_score, useNA = "ifany"))

cat("\nFiltered assigned cell types:\n")
print(table(foregut$celltype_score_filtered, useNA = "ifany"))

