#Load required packages 
library(Seurat)
library(FNN)
library(dplyr)
library(ggplot2)
library(tidyr)
#Load RDS object with immune and mesenchymal cells assigned 
foregut<-readRDS("foregutimmunemesenchymeassigned.rds")
#Load RDS object for foregut and assign identities for both objects
tissue<-readRDS("foregut.rds")
tissue$Identity<-Idents(tissue)
foregut$category<-tissue$Identity
Idents(foregut)<-foregut$celltype_score_filtered

# =========================================================
# 1. Remove zoom image
# =========================================================

foregut@images[["zoom"]] <- NULL

# =========================================================
# 2. Get coordinates
# =========================================================

coords_list <- lapply(Images(foregut), function(img) {
  tmp <- GetTissueCoordinates(foregut, image = img)
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
foregut$Identity <- as.character(Idents(foregut))

# Keep only cells with coordinates
common_cells <- intersect(colnames(foregut), rownames(coords))

coords <- coords[common_cells, ]
foregut_sub <- subset(foregut, cells = common_cells)

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
# 6. Optional: remove extreme outlier distances
# =========================================================

max_dist <- quantile(
  foregut_epi$dist_to_base,
  probs = 0.95,
  na.rm = TRUE
)

foregut_epi <- subset(
  foregut_epi,
  subset = dist_to_base <= max_dist
)

# =========================================================
# 7. Plot gene expression with increasing distance from base
# =========================================================

genes_use <- c("APOA4","DMBT1","RBP2","CCK","GCG","GIP","MLN","TPH1","CLCA1","SH2D6","TRPM5")
genes_use <- c("TRAC","TRBC1","TRBC2")
genes_use <- c("CA9","NDRG1","SLC2A1")

bin_size <- 10

expr <- FetchData(
  foregut_epi,
  vars = genes_use
)

df <- cbind(
  expr,
  dist = foregut_epi$dist_to_base
) %>%
  as.data.frame() %>%
  filter(!is.na(dist))

df_binned <- df %>%
  mutate(
    dist_bin = floor(dist / bin_size) * bin_size
  ) %>%
  group_by(dist_bin) %>%
  summarise(
    across(all_of(genes_use), ~ mean(.x, na.rm = TRUE)),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(dist_bin)

df_long <- df_binned %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene",
    values_to = "expression"
  )

gene_colors <- c(
  "APOA4" = "darkblue",
  "DMBT1" = "darkblue",
  "RBP2" = "darkblue",
  "CCK"="darkmagenta",
  "GCG"="darkmagenta",
  "GIP"="darkmagenta",
  "MLN"="darkmagenta",
  "TPH1"="darkmagenta",
  "CLCA1"="red",
  "SH2D6"="orange",
  "TRPM5"="orange"
)
p <- ggplot(df_long, aes(x = dist_bin, y = expression, color = gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(size = n_cells)) +
  facet_wrap(~gene, scales = "fixed") +
  scale_color_manual(values = gene_colors) +
  scale_x_continuous(
    limits = c(0, 150),
    breaks = seq(0, 150, by = bin_size),
    expand = expansion(mult = c(0.05, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    x = "Distance from nearest mesenchyme / crypt base",
    y = "Average expression in epithelial cells",
    title = "Epithelial gene expression vs distance from crypt/base",
    size = "Cells/bin"
  )

print(p)

ggsave(
  "foregut_gene_expression_distance.pdf",
  plot = p,
  width = 20,
  height = 15,
  device = cairo_pdf   # nicer fonts
)


p <- ggplot(df_long, aes(x = dist_bin, y = expression, color = gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(size = n_cells)) +
  scale_x_continuous(
    limits = c(0, 150),
    breaks = seq(0, 150, by = bin_size)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    x = "Distance from nearest mesenchyme / crypt base",
    y = "Average expression in epithelial cells",
    size = "Cells/bin"
  )

print(p)


genes_use <- c(
  "APOA4","DMBT1","RBP2",
  "CCK","GCG","GIP","MLN","TPH1",
  "CLCA1","SH2D6","TRPM5"
)

bin_size <- 10

expr <- FetchData(foregut_epi, vars = genes_use)

df <- cbind(
  expr,
  dist = foregut_epi$dist_to_base
) %>%
  as.data.frame() %>%
  filter(!is.na(dist))

df_binned <- df %>%
  mutate(dist_bin = floor(dist / bin_size) * bin_size) %>%
  group_by(dist_bin) %>%
  summarise(
    across(all_of(genes_use), ~ mean(.x, na.rm = TRUE)),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(dist_bin)

df_long <- df_binned %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(
    gene = factor(gene, levels = genes_use)  # 🔑 FORCE ORDER
  )

gene_colors <- c(
  "APOA4" = "darkblue",
  "DMBT1" = "darkblue",
  "RBP2" = "darkblue",
  "CCK"="darkmagenta",
  "GCG"="darkmagenta",
  "GIP"="darkmagenta",
  "MLN"="darkmagenta",
  "TPH1"="darkmagenta",
  "CLCA1"="red",
  "SH2D6"="orange",
  "TRPM5"="orange"
)

p <- ggplot(df_long, aes(x = dist_bin, y = expression, color = gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(size = n_cells)) +
  facet_wrap(~gene, ncol = 5, scales = "fixed") +   # 🔑 layout
  scale_color_manual(values = gene_colors) +
  scale_x_continuous(
    limits = c(0, 150),                              # 🔑 match your figure
    breaks = seq(0, 150, by = bin_size),
    expand = expansion(mult = c(0.05, 0.02))
  ) +
  scale_y_continuous(limits = c(0, 5)) +            # optional: match scale
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Distance from mesenchyme (µm)",
    y = NULL,
    title = "Jejunum",
    size = "Cells/bin"
  )

print(p)

ggsave(
  "colon_gene_expression_distance.pdf",
  plot = p,
  width = 12,
  height = 8,
  device = cairo_pdf
)