##############################################################
# Load libraries
##############################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

##############################################################
# 1. File paths and replicate mapping (Load duplicates in bulk RNA-seq)
##############################################################

dir <- "/athena/chenlab/scratch/jjv4001/gutrevision/KcellGSEA"

files  <- c("DMSO1.csv","DMSO2.csv","cAMP1.csv","cAMP2.csv","CREB1.csv","CREB2.csv")
paths  <- file.path(dir, files)

# Explicitly define conditions for each sample
conditions <- c(
  "DMSO", "DMSO",
  "cAMP", "cAMP",
  "CREB", "CREB"
)

##############################################################
# 2. Define Marker sets for each lineage
##############################################################


K_markers <- c("GIP","GIPR","FFAR1","FFAR4","RFX6","FOXA2","ISL1","FEV","GCK","STXBP1","CPE")

marker_lists <- list(
  K = K_markers
)

##############################################################
# 3. Module score calculator (Function for module score)
##############################################################

compute_module_score <- function(df, genes) {
  g <- intersect(genes, df$Gene)
  if (length(g) < 3) return(NA_real_)
  df_sub <- df[df$Gene %in% g, ]
  mean(log2(df_sub$avg_expr + 0.1), na.rm = TRUE)
}


##############################################################
# 4. Compute module scores (Compute module score)
##############################################################

all_scores <- list()

for (i in seq_along(paths)) {
  df <- read.csv(paths[i], stringsAsFactors = FALSE)
  
  lineage_scores <- sapply(marker_lists, compute_module_score, df = df)
  
  all_scores[[i]] <- data.frame(
    replicate = files[i],
    lineage   = names(lineage_scores),
    score     = as.numeric(lineage_scores),
    stringsAsFactors = FALSE
  )
}

scores_df <- bind_rows(all_scores)

# Add correct condition assignment
scores_df$condition <- conditions[ match(scores_df$replicate, files) ]

##############################################################
# 5. Bar plot: each replicate shown separately
##############################################################

scores_df$replicate <- factor(scores_df$replicate, levels = files)

ggplot(scores_df, aes(x = replicate, y = score, fill = condition)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c("cAMP" = "#F8766D", "DMSO" = "#00B0F6", "CREB" = "darkgreen")) +
  theme_bw(base_size = 16) +
  ylab("K-cell module score (mean log2(avg_expr + 0.1))") +
  xlab("Replicate") +
  ggtitle("K-cell module score per replicate") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

##############################################################
# 6. Heatmap of replicate-level module scores
##############################################################

heat_df <- scores_df %>%
  select(lineage, replicate, score) %>%
  pivot_wider(names_from = replicate, values_from = score) %>%
  as.data.frame()

rownames(heat_df) <- heat_df$lineage
heat_mat <- as.matrix(heat_df[ , -1])

pheatmap(
  heat_mat,
  color = colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Replicate-level K-cell module scores"
)
