
library(Seurat)

genes <- c("EPCAM", "CDX2","COL1A1","COL6A1","TUBB2B", 
           "ELAVL4","HBA1","HBB","PECAM1","GSN",
           "CCL3","CD74")
# markers

id <- "PMID33406409"
message(paste("sample:", id))
# PMID of the dataset

dirs <- dir(path="./", pattern="run_count")
seurat <- lapply(dirs, function(dir){
  message(dir)
  dir <- paste("./", dir, "/outs/filtered_feature_bc_matrix/", sep="")
  seurat <- Read10X(dir)
  seurat <- CreateSeuratObject(counts=seurat, project=gsub("/outs/filtered_feature_bc_matrix/", "", gsub("./run_count_", "", dir)))
  seurat})
seurat <- merge(seurat[[1]], y=seurat[2:length(seurat)], project=id, add.cell.ids=gsub("run_count_", "", dirs))
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern="^MT-")
# data preparation

seurat <- subset(seurat, subset=nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 40)
# data filtering

seurat <- NormalizeData(seurat, normalization.method="LogNormalize", scale.factor=10000)
seurat <- FindVariableFeatures(seurat, selection.method="vst", nfeatures=400)
seurat <- ScaleData(seurat, features=rownames(seurat))
# data normalization

seurat <- RunPCA(seurat , features=VariableFeatures(object=seurat))
seurat <- RunUMAP(seurat, dims=1:10, n.components=3)
seurat@misc$umap3d <- seurat@reductions$umap
seurat <- RunUMAP(seurat, dims=1:10, n.components=2)
seurat@misc$umap2d <- seurat@reductions$umap
# dimension reduction

seurat <- FindNeighbors(seurat, dims=1:10)
pdf(paste(id, "_umap_clusters_iter_1.pdf", sep=""), width=5, height=5)
for (resolution in seq(0, 1, 0.05)){
  g <- FindClusters(seurat, resolution=resolution)
  g <- DimPlot(g, reduction="umap", label=TRUE, raster=TRUE)
  print(g)}
dev.off()


pdf(paste(id, "_umap_QC_iter_1.pdf", sep=""), width=9, height=8)
FeaturePlot(seurat, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), raster=TRUE)
dev.off()

pdf(paste(id, "_umap_markers_iter_1.pdf", sep=""), width=16, height=48)
FeaturePlot(seurat, features=genes, raster=TRUE)
dev.off()

pdf(paste(id, "_umap_group_iter_1.pdf", sep=""))
DimPlot(seurat, reduction="umap", group.by="orig.ident", label=TRUE, raster=TRUE)
dev.off()
# visualization

saveRDS(seurat, file=paste(id, "_iter_1_seurat.rds", sep=""))
# save Seurat analysis

####################################################
##########          iteration #1          ##########
####################################################

merged.gut<-readRDS("/home/fs01/jjv4001/PMID33406409_iter_1_seurat.rds")
merged.gut <- FindClusters(merged.gut, resolution=0.05)
cells <- structure(c("Mesenchymal", "Epithelial", "Neuronal", "Endothelial", "Immune cells", "Epithelial","Mesenchymal","7","Mesenchymal","RBCs","Mesenchymal","Mesenchymal","RBCs","RBCs","Mesenchymal"), names=levels(merged.gut))
merged.gut <- RenameIdents(merged.gut, cells)
merged.gut$Cell_Type <- as.character(merged.gut@active.ident)
pdf(paste("PMID33406409_umap_cluster_iter_6.pdf", sep=""))
DimPlot(merged.gut, reduction="umap", label=TRUE)
dev.off()
# clustering analysis
merged.gut<-subset(merged.gut, idents="Epithelial")
saveRDS(merged.gut, file=paste("PMID33406409_integrated_iter_2_seurat.rds", sep=""))
# save Seurat analysis

####################################################
##########          iteration #2          ##########
####################################################