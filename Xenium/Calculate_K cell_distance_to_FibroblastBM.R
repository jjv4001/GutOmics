#Load the required packages 
library(Seurat)
library(FNN)
library(dplyr)
library(ggplot2)

#Load RDS object for foregut with immune and mesenchymal cells assigned 
foregut1<-readRDS("foregutimmunemesenchymeassigned.rds")

#Load eec object with all eecs assigned for foregut
eec<-readRDS("newguteec.rds")

#Subset only the eecs from foregut tissue 
eec<-subset(eec, tissue=="foregut")

# Specify identity column
foregut1$Identity <- as.character(foregut1$celltype_score_filtered)
eec$Identity<-eec$gut_endocrine_call

# =========================================================
# 0. Set metadata columns
# =========================================================

foregut1$Identity <- as.character(foregut1$celltype_score_filtered)

# choose cell type labels for foregut
eec$Identity <- as.character(eec$gut_endocrine_call)

# =========================================================
# 1. Function to get coordinates safely
# =========================================================

get_coords <- function(obj) {
  coords_list <- lapply(Images(obj), function(img) {
    tmp <- GetTissueCoordinates(obj, image = img)
    tmp$image <- img
    
    if (!"cell" %in% colnames(tmp)) {
      tmp$cell <- rownames(tmp)
    }
    
    tmp
  })
  
  coords <- bind_rows(coords_list) %>%
    filter(!is.na(cell), !is.na(x), !is.na(y)) %>%
    distinct(cell, .keep_all = TRUE)
  
  rownames(coords) <- coords$cell
  coords
}

# =========================================================
# 2. Get coordinates
# =========================================================

coords_eec  <- get_coords(eec)
coords_foregut1 <- get_coords(foregut1)

# =========================================================
# 3. Build query matrix from EEC cells
# =========================================================

query_cells <- intersect(colnames(eec), rownames(coords_eec))

coords_query <- coords_eec[query_cells, ]
coords_query <- coords_query[colnames(eec)[colnames(eec) %in% query_cells], ]

query_mat <- as.matrix(coords_query[, c("x", "y")])
storage.mode(query_mat) <- "numeric"

query_mat <- query_mat[complete.cases(query_mat), , drop = FALSE]

# =========================================================
# 4. Build reference matrix from FibroblastBM cells in foregut1
# =========================================================

ref_cells <- intersect(colnames(foregut1), rownames(coords_foregut1))

coords_ref <- coords_foregut1[ref_cells, ]
coords_ref <- coords_ref[colnames(foregut1)[colnames(foregut1) %in% ref_cells], ]

ref_mat <- as.matrix(coords_ref[, c("x", "y")])
storage.mode(ref_mat) <- "numeric"

ref_mat <- ref_mat[complete.cases(ref_mat), , drop = FALSE]

fibroBM_cells <- colnames(foregut1)[
  foregut1$Identity == "FibroblastBM"
]

fibroBM_cells <- intersect(fibroBM_cells, rownames(ref_mat))

fibroBM_mat <- ref_mat[fibroBM_cells, , drop = FALSE]

# check
length(fibroBM_cells)
dim(fibroBM_mat)

# =========================================================
# 5. Calculate distance from every EEC cell to nearest FibroblastBM
# =========================================================

nn <- get.knnx(
  data = fibroBM_mat,
  query = query_mat,
  k = 1
)

eec$dist_to_foregut1_FibroblastBM <- NA
eec$dist_to_foregut1_FibroblastBM[rownames(query_mat)] <- nn$nn.dist[, 1]

summary(eec$dist_to_foregut1_FibroblastBM)

# =========================================================
# 6. Create dataframe for plotting/summarizing
# =========================================================

df_dist <- data.frame(
  cell = colnames(eec),
  celltype = eec$Identity,
  dist_to_FibroblastBM = eec$dist_to_foregut1_FibroblastBM
) %>%
  filter(
    !is.na(celltype),
    celltype != "",
    !is.na(dist_to_FibroblastBM)
  )

table(df_dist$celltype)

# =========================================================
# 6. Create dataframe with distance + GIP expression
# =========================================================

gip_expr <- FetchData(eec, vars = "GIP")

df_combined1 <- data.frame(
  cell = colnames(eec),
  celltype = eec$Identity,
  dist_to_FibroblastBM = eec$dist_to_foregut1_FibroblastBM,
  GIP = gip_expr[colnames(eec), "GIP"]
) %>%
  filter(
    celltype == "K",
    !is.na(dist_to_FibroblastBM),
    !is.na(GIP)
  )

# Check
head(df_combined1)
summary(df_combined1$GIP)

# =========================================================
# 7. Summary by cell type
# =========================================================

summary_dist <- df_dist %>%
  group_by(celltype) %>%
  summarise(
    n_cells = n(),
    mean_dist = mean(dist_to_FibroblastBM),
    median_dist = median(dist_to_FibroblastBM),
    q25 = quantile(dist_to_FibroblastBM, 0.25),
    q75 = quantile(dist_to_FibroblastBM, 0.75),
    .groups = "drop"
  ) %>%
  arrange(median_dist)

summary_dist

# =========================================================
# 9. Histogram by cell type
# =========================================================
df_combined1<-subset(df_combined1, celltype=="K")
ggplot(df_combined1, aes(x = dist_to_FibroblastBM)) +
  geom_histogram(
    binwidth = 2,
    fill = "#4C8DCB",
    color = "white",
    linewidth = 0.2
  )  +
  theme_classic(base_size = 15) +
  labs(
    x = "Distance to nearest foregut1 FibroblastBM (µm)",
    y = "Cell count"
  )