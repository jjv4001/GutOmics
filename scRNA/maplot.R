library(ggplot2)
library(ggrepel)  # for text labels that don't overlap

# Assuming your data frame is called `df` and contains:
# - avg_log2FC: log2 fold change
# - avg_expr: average expression (log10 scale)
# - padj: adjusted p-value
# - gene: gene symbol

# Create a column to color points by significance
markers_filtered$significance <- ifelse(markers_filtered$p_val_adj < 0.05, "padj < 0.05", "padj > 0.05")

# Optional: Select genes to label (e.g., top N)
top_genes <- markers_filtered[markers_filtered$p_val_adj  < 0.05 & (markers_filtered$avg_log2FC) > 0.5 & markers_filtered$logLI>0.2, ]  # adjust as needed

# Plot
ggplot(markers_filtered, aes(x = logLI, y = avg_log2FC)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1) +
  scale_color_manual(values = c("padj < 0.05" = "red", "padj > 0.05" = "grey")) +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = top_genes,
                  aes(label = Gene),
                  size = 3, max.overlaps = 50) +
  theme_minimal() +
  labs(x = "Log₁₀(Average Expression)",
       y = "Log₂(Fold Change)",
       color = "Significance")


