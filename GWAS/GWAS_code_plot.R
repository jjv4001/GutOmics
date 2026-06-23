############################################################
## 1. Load libraries
############################################################
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(dendsort)

sort_hclust <- function(...) {
  as.hclust(dendsort(as.dendrogram(...)))
}

############################################################
## 2. Read LDSC summary data
############################################################
results <- read_tsv(
  "res.summary.tsv",
  col_names = c(
    "Cluster",
    "Coefficient",
    "Coefficient_std_error",
    "Coefficient_P_value",
    "Trait",
    "Set"
  )
)

############################################################
## 3. Define metabolic GWAS traits (FINAL)
############################################################
traits_metabolic <- c(
  ## BMI
  "Locke.Nature.2015.BMI",
  "UKB.Neale_v2.2018.Body_mass_index",
  "Yengo.biorxiv.2018.BMI",
  
  ## Diabetes
  "Mahajan.NatGenet.2018.Type_2_Diabetes",
  "Aylward.2018.biorxiv.Type_1_Diabetes",
  "Chiou.Nature.2021.Type_1_diabetes",
  
  ## Glycemic traits
  "Manning.NatGenet.2012.Fasting_Glucose",
  "Manning.NatGenet.2012.Fasting_Insulin",
  
  ## Liver
  "Namjou.BMCMed.2019.Nonalcoholic_fatty_liver_disease"
)

############################################################
## 4. Keep ALL clusters
############################################################
cls <- results %>%
  filter(Trait %in% traits_metabolic) %>%
  distinct(Cluster) %>%
  pull()

############################################################
## 5. Compute FDR per trait (for stars)
############################################################
cluster_traits <- results %>%
  filter(
    Trait %in% traits_metabolic,
    Cluster %in% cls
  ) %>%
  group_by(Trait) %>%
  mutate(FDR = p.adjust(Coefficient_P_value, method = "fdr")) %>%
  ungroup()

############################################################
## 6. Build enrichment matrix (Z = coef / SE)
############################################################
mat <- results %>%
  filter(
    Trait %in% traits_metabolic,
    Cluster %in% cls
  ) %>%
  mutate(enrich = Coefficient / Coefficient_std_error) %>%
  select(Cluster, Trait, enrich) %>%
  pivot_wider(
    names_from  = Trait,
    values_from = enrich
  ) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

############################################################
## 7. Build FDR matrix
############################################################
mat_p <- cluster_traits %>%
  filter(
    Trait %in% traits_metabolic,
    Cluster %in% cls
  ) %>%
  select(Cluster, Trait, FDR) %>%
  pivot_wider(
    names_from  = Trait,
    values_from = FDR
  ) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

############################################################
## 8. Clean column names
############################################################
colnames(mat)   <- gsub(".*[.]", "", colnames(mat))
colnames(mat_p) <- gsub(".*[.]", "", colnames(mat_p))

############################################################
## 9. Force cluster (row) order
############################################################
desired_cluster_order <- c(
  "SCs",
  "Enterocytes",
  "Goblet",
  "EECs"
)

desired_cluster_order <- desired_cluster_order[
  desired_cluster_order %in% rownames(mat)
]

mat   <- mat[desired_cluster_order, , drop = FALSE]
mat_p <- mat_p[desired_cluster_order, , drop = FALSE]

############################################################
## 10. Create annotation matrix (FDR stars)
############################################################
anno <- matrix(
  "",
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = list(rownames(mat), colnames(mat))
)

anno[mat_p < 0.1]   <- "-"
anno[mat_p < 0.05]  <- "*"
anno[mat_p < 0.01]  <- "**"
anno[mat_p < 0.001] <- "***"

############################################################
## 11. Column clustering only (rows fixed)
############################################################
hclust_method <- "ward.D2"

hclust_rows <- FALSE

hclust_cols <- if (ncol(mat) > 1) {
  sort_hclust(hclust(dist(t(mat)), method = hclust_method))
} else {
  FALSE
}

############################################################
## 12. Heatmap parameters
############################################################
upper <- 4
lower <- -2
n_cols <- (upper - lower) * 10

############################################################
## 13. Plot heatmap
############################################################
options(repr.plot.width = 30, repr.plot.height = 3.8)

pheatmap(
  mat,
  color = colorRampPalette(
    rev(brewer.pal(n = 7, name = "RdBu"))
  )(n_cols),
  breaks = seq(lower, upper, by = 0.1),
  cluster_rows = hclust_rows,
  cluster_cols = hclust_cols,
  fontsize_row = 10,
  fontsize_col = 10,
  display_numbers = anno,
  fontsize_number = 6,
  labels_col = gsub("_", " ", colnames(mat)),
  labels_row = rownames(mat),
  angle_col = 45,
  border_color = NA,
  main = "GWAS on snATAC-seq clusters"
)
