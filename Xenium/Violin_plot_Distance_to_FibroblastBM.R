#Load RDS object with epithelial cells
epithelial<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/newgutassignedmarker.rds")
#subset only the foregut tissue 
foregut<-subset(epithelial, tissue=="foregut")
#Load immune and mesenchymal cell assignments 
foregutIM<-readRDS("foregutimmunemesenchymeassigned.rds")
colnames(foregutIM@meta.data)
colnames(foregut@meta.data)
coords <- GetTissueCoordinates(foregutIM)

#Get coordinates from foregut object with immune and mesenchymal cell assignments 
foregutIM$x <- coords[,1]
foregutIM$y <- coords[,2]
coords1 <- GetTissueCoordinates(foregut)

#Get coordinates from foregut object 
foregut$x <- coords1[,1]
foregut$y <- coords1[,2]

#find the coordinates for fibroblast BM cells 
foregutIM$Identity<-Idents(foregutIM)
mesenchyme_cells <- WhichCells(foregut, expression = Identity == "FibroblastBM")
mes_coords <- foregut@meta.data[mesenchyme_cells, c("x","y")]

#find the coordinates for foregut epithelial cell types
foregut_coords <- foregut@meta.data[, c("x","y")]
foregut_types <- foregut$celltype_score_filtered   # or your subtype column

#Load required package
library(FNN)

#find the nearest neighbor distance from each epithelial cell type in foregut to fibroblast BM
nn <- get.knnx(
  data = as.matrix(mes_coords),
  query = as.matrix(foregut_coords),
  k = 1
)

foregut$dist_to_mes <- nn$nn.dist[,1]
df <- data.frame(
  dist = foregut$dist_to_mes,
  subtype = foregut$celltype_score_filtered
)

#Load ggplot package and order cell types, plot violin plot showing distance of epithelial cell types to Fibroblast BM
library(ggplot2)

subtype_order <- c(
  # Epithelial
  "EpiCA9", "EpiGLUT", "EpiNDRG1",
  
  # Immune
  "Bcell", "cDC1", "cDC2", "FDC",
  "Macrophage", "Plasma_cell",
  "Tcell_CD4", "Tcell_CD8",
  
  # Mesenchymal
  "Fibroblast", "FibroblastBM",
  "EndoArt", "EndoLymph", "EndoVenous",
  "Telocyte", "Telocyte_tip",
  
  # Other
  "Unassigned"
)

df$subtype <- factor(df$subtype, levels = subtype_order)

p <- ggplot(df, aes(x = subtype, y = dist, fill = subtype)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "black",
    linewidth = 0.25,
    alpha = 0.9,
    adjust = 0.7,        # 🔑 LESS smoothing → more true shape
    bw = "nrd0"          # 🔑 stable bandwidth (like Seurat defaults)
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.8,
    color = "black"
  ) +
  scale_fill_hue() +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Distance to nearest mesenchyme (µm)"
  )

p +
  geom_vline(xintercept = c(3.5, 11.5), linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 2, y = max(df$dist)*1.05, label = "Epithelial", size = 5) +
  annotate("text", x = 7, y = max(df$dist)*1.05, label = "Immune", size = 5) +
  annotate("text", x = 15, y = max(df$dist)*1.05, label = "Mesenchyme", size = 5)