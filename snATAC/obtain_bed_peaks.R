#Load required packages
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(Matrix)

#Load ATAC assay
DefaultAssay(epithelial) <- "ATAC"
celltypes <- unique(epithelial$celltype)

# Get full peak set with names
all_peaks <- granges(epithelial)
names(all_peaks) <- rownames(epithelial)

# Initialize list to store accessible GRanges per group
all_accessible <- list()

for (ct in celltypes) {
  message("Processing: ", ct)
  
  # Subset cells for this cell type
  sub_obj <- subset(epithelial, subset = celltype == ct)
  
  # Get peak accessibility matrix
  peak_matrix <- GetAssayData(sub_obj, layer = "counts")  # Use slot = "counts" for Seurat v4
  
  # Get peaks accessible in ≥1 cell
  accessible_peaks <- rownames(peak_matrix)[Matrix::rowSums(peak_matrix) > 0]
  
  # Match GRanges
  peaks_ct <- all_peaks[accessible_peaks]
  
  # Save for union later
  all_accessible[[ct]] <- peaks_ct
  
  # Export per cell type BED
  export(peaks_ct, con = paste0("peaks_", gsub(" ", "_", ct), ".bed"), format = "bed")
}

# Step 1: Combine all GRanges objects into a single GRanges
combined_gr <- do.call(c, unname(all_accessible))

# Step 2: Reduce overlapping/adjacent regions
library(GenomicRanges)
union_peaks <- reduce(combined_gr)

#export bed files
export(union_peaks, con = "union_peaks.bed", format = "bed")
