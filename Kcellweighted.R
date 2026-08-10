##############################################################
# Load packages
##############################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

##############################################################
# 1. File paths and replicate mapping
##############################################################

dir <- "/athena/chenlab/scratch/jjv4001/gutrevision/KcellGSEA"

files  <- c("DMSO1.csv","DMSO2.csv","cAMP1.csv","cAMP2.csv")
paths  <- file.path(dir, files)

# Explicitly define conditions
conditions <- c(
  "DMSO", "DMSO",
  "cAMP", "cAMP"
  
)

##############################################################
# 2. Marker sets for each lineage (K-cell module)
##############################################################

K_markers <- c("GIP","FFAR4","GCK","RFX6","FOXA2",
               "ISL1","TCF7L2","STXBP1","CPE")

# Assign weights to each gene
K_weights <- c(
  GIP = 3.0,
  FFAR1 = 2.0,
  FFAR4 = 2.0,
  GCK = 2.0,
  RFX6 = 2.0,
  ISL1 = 2.0,
  TCF7L2 = 2.0,
  FOXA2=1.0,
  STXBP1=1.0,
  CPE=1.0
)

##############################################################
# 3. Weighted module score calculator
##############################################################

compute_module_score <- function(df, genes, weights) {
  
  # Genes present in CSV AND in weight vector
  g <- intersect(genes, df$Gene)
  g <- intersect(g, names(weights))   # ensures weights match
  
  if (length(g) < 1) return(NA_real_)
  
  df_sub <- df[df$Gene %in% g, ]
  
  # Reorder weights to match df_sub gene order
  w <- weights[df_sub$Gene]
  
  # Weighted average of log2 expression
  score <- sum(w * log2(df_sub$avg_expr + 0.1), na.rm = TRUE) / sum(w)
  
  return(score)
}

##############################################################
# 4. Compute weighted module scores (one value per replicate)
##############################################################

all_scores <- list()

for (i in seq_along(paths)) {
  df <- read.csv(paths[i], stringsAsFactors = FALSE)
  
  score_value <- compute_module_score(
    df = df,
    genes = K_markers,
    weights = K_weights
  )
  
  all_scores[[i]] <- data.frame(
    replicate = files[i],
    lineage = "K",
    score = score_value,
    stringsAsFactors = FALSE
  )
}

scores_df <- bind_rows(all_scores)

# Add condition column
scores_df$condition <- conditions[ match(scores_df$replicate, files) ]

##############################################################
# 5. Bar plot: weighted K-cell module score per replicate
##############################################################

scores_df$replicate <- factor(scores_df$replicate, levels = files)

ggplot(scores_df, aes(x = replicate, y = score, fill = condition)) +
  geom_col(width = 0.4) +
  scale_fill_manual(values = c("cAMP" = "#F8766D",
                               "DMSO" = "#00B0F6"
                              )) +
  theme_bw(base_size = 16) +
  ylab("Weighted K-cell module score") +
  xlab("Replicate") +
  ggtitle("Weighted K-cell module score per replicate") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

##############################################################
# 6. Heatmap of weighted module scores
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
  main = "Weighted K-cell module scores per replicate"
)
