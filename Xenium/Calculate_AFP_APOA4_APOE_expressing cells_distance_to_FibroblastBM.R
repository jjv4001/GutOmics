#Load required packages 
library(Seurat)
library(FNN)
library(dplyr)
library(ggplot2)
library(tidyr)

#Load foregut object with immune and mesenchymal cells specified 
foregutIM<-readRDS("foregutimmunemesenchymeassigned.rds")

#Load foregut object and assign categories and relevant identities 
foregut<-readRDS("foregut.rds")
foregut$category<-Idents(foregut)
foregutIM$category<-foregut$category
Idents(foregutIM)<-foregutIM$celltype_score_filtered

# =========================================================
# 1. Remove zoom image
# =========================================================

foregutIM@images[["zoom"]] <- NULL

# =========================================================
# 2. Get coordinates
# =========================================================

coords_list <- lapply(Images(foregutIM), function(img) {
  tmp <- GetTissueCoordinates(foregutIM, image = img)
  tmp$image <- img
  
  if (!"cell" %in% colnames(tmp)) {
    tmp$cell <- rownames(tmp)
  }
  
  tmp
})

coords <- bind_rows(coords_list) %>%
  filter(!is.na(cell), !is.na(x), !is.na(y))

rownames(coords) <- coords$cell

# Add identity metadata
foregutIM$Identity <- as.character(Idents(foregutIM))

# Keep only cells with coordinates
common_cells <- intersect(colnames(foregutIM), rownames(coords))

coords <- coords[common_cells, ]
foregut_sub <- subset(foregutIM, cells = common_cells)

coords <- coords[colnames(foregut_sub), ]

coords_mat <- as.matrix(coords[, c("x", "y")])
storage.mode(coords_mat) <- "numeric"

# =========================================================
# 3. Define mesenchyme/base cells
# =========================================================

mes_cells <- WhichCells(
  foregut_sub,
  idents = "FibroblastBM"
)

mes_coords_mat <- coords_mat[mes_cells, , drop = FALSE]

# =========================================================
# 4. Calculate distance to nearest mesenchyme/base
# =========================================================

nn <- get.knnx(
  data = mes_coords_mat,
  query = coords_mat,
  k = 1
)

foregut_sub$dist_to_base <- nn$nn.dist[, 1]

summary(foregut_sub$dist_to_base)

# =========================================================
# 5. Restrict to epithelial cells only
# =========================================================

epi_cells <- WhichCells(
  foregut_sub,
  expression = category == "Epithelial"
)

foregut_epi <- subset(
  foregut_sub,
  cells = epi_cells
)


# =========================================================
# Cells in top quartile of ALL 3 genes simultaneously
# AFP + APOA4 + APOE
# =========================================================
#get gene expression matrix for cells expressing AFP, APOA4, and APOE
expr_df <- FetchData(
  foregut_epi,
  vars = c("AFP", "APOA4", "APOE")
)

#make dataframe with the relevant information
df_expr <- data.frame(
  cell = colnames(foregut_epi),
  dist_to_FibroblastBM = foregut_epi$dist_to_base,
  AFP = expr_df$AFP,
  APOA4 = expr_df$APOA4,
  APOE = expr_df$APOE
)%>%
  filter(!is.na(dist_to_FibroblastBM))

# =========================================================
# Calculate cutoffs
# =========================================================

afp_cutoff <- quantile(
  df_expr$AFP[df_expr$AFP > 0],
  0.75,
  na.rm = TRUE
)

apoa4_cutoff <- quantile(
  df_expr$APOA4[df_expr$APOA4 > 0],
  0.75,
  na.rm = TRUE
)

apoe_cutoff <- quantile(
  df_expr$APOE[df_expr$APOE > 0],
  0.75,
  na.rm = TRUE
)

# =========================================================
# Keep cells high for ALL 3 genes
# =========================================================

df_top3 <- df_expr %>%
  filter(
    AFP >= afp_cutoff,
    APOA4 >= apoa4_cutoff,
    APOE >= apoe_cutoff,
    !is.na(dist_to_FibroblastBM)
  )

# =========================================================
# Summary statistics
# =========================================================

summary_top3 <- df_top3 %>%
  summarise(
    n_cells = n(),
    
    mean_dist = mean(dist_to_FibroblastBM, na.rm = TRUE),
    
    median_dist = median(dist_to_FibroblastBM, na.rm = TRUE),
    
    min_dist = min(dist_to_FibroblastBM, na.rm = TRUE),
    
    q25 = quantile(dist_to_FibroblastBM, 0.25, na.rm = TRUE),
    
    q75 = quantile(dist_to_FibroblastBM, 0.75, na.rm = TRUE),
    
    max_dist = max(dist_to_FibroblastBM, na.rm = TRUE),
    
    sd_dist = sd(dist_to_FibroblastBM, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), round, 2))

summary_top3

#Plot histogram for top quartile cells expressing all 3 genes and the distance from mesenchyme

ggplot(df_top3, aes(x = dist_to_FibroblastBM)) +
  geom_histogram(
    binwidth = 5,
    fill = "darkblue",
    color = "white"
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Distance to nearest FibroblastBM",
    y = "Cell count",
    title = "Top quartile AFP/APOA4/APOE cells"
  )