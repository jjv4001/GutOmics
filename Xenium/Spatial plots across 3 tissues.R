#Load required packages 
library(Seurat)
library(ggplot2)
library(patchwork)
library(future)

#increase memory limit 
plan(sequential)
options(future.globals.maxSize = 64 * 1024^3)

# =========================================================
# OUTPUT FOLDER
# =========================================================
#Create output directory 
outdir <- "spatial_gene_panels"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Specify genes for plotting 
genes_use <- features$rownames.features.   # your 422 genes

# =========================================================
# CROP COORDINATES
# =========================================================
duo_x <- c(1500, 3500)
duo_y <- c(1800, 4300)

mid_x <- c(1000, 3000)
mid_y <- c(4100, 6300)

col_x <- c(3500, 6000)
col_y <- c(2500, 5000)

#Specify colors 
duo_col <- "blue"
mid_col <- "blue"
col_col <- "red"

# =========================================================
# HELPER FUNCTION
# =========================================================
make_spatial_panel <- function(obj, source_fov, xlim, ylim,
                               gene, title_text, high_col,
                               reverse_y = TRUE) {
  
  cropped_fov <- Crop(
    obj[[source_fov]],
    x = xlim,
    y = ylim,
    coords = "plot"
  )
  
  obj[["temp_zoom"]] <- cropped_fov
  DefaultBoundary(obj[["temp_zoom"]]) <- "segmentation"
  
  p <- ImageFeaturePlot(
    obj,
    fov = "temp_zoom",
    features = gene,
    cols = c("grey85", high_col),
    dark.background = FALSE,
    border.color = NA,
    axes = FALSE,
    coord.fixed = FALSE,
    size = 0.4
  ) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(size = 22, hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
  
  if (reverse_y) {
    p <- p + scale_y_reverse()
  }
  
  return(p)
}

# =========================================================
# LOOP THROUGH GENES
# =========================================================
for (gene_use in genes_use) {
  
  message("Plotting: ", gene_use)
  
  tryCatch({
    
    p_duo <- make_spatial_panel(
      obj = epithelial,
      source_fov = "foregut.fov.1691",
      xlim = duo_x,
      ylim = duo_y,
      gene = gene_use,
      title_text = "Duodenum",
      high_col = duo_col
    )
    
    p_mid <- make_spatial_panel(
      obj = epithelial,
      source_fov = "midgut.fov.1691",
      xlim = mid_x,
      ylim = mid_y,
      gene = gene_use,
      title_text = "Jejunum",
      high_col = mid_col
    )
    
    p_col <- make_spatial_panel(
      obj = epithelial,
      source_fov = "colon.fov.1691",
      xlim = col_x,
      ylim = col_y,
      gene = gene_use,
      title_text = "Colon",
      high_col = col_col
    )
    
    gene_label <- ggplot() +
      annotate(
        "text",
        x = 1, y = 1,
        label = gene_use,
        angle = 90,
        size = 8,
        fontface = "bold.italic",
        colour = "black"
      ) +
      theme_void() +
      xlim(0, 2) +
      ylim(0, 2)
    
    final_plot <- gene_label + p_duo + p_mid + p_col +
      plot_layout(widths = c(0.12, 1, 1, 1))
    
    safe_gene <- gsub("[^A-Za-z0-9_.-]", "_", gene_use)
    
    ggsave(
      filename = file.path(outdir, paste0(safe_gene, "_duodenum_midgut_colon_spatial.png")),
      plot = final_plot,
      width = 10,
      height = 3.5,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(outdir, paste0(safe_gene, "_duodenum_midgut_colon_spatial.pdf")),
      plot = final_plot,
      width = 13,
      height = 5
    )
    
  }, error = function(e) {
    message("Skipping ", gene_use, " due to error: ", e$message)
  })
  
  gc()
}