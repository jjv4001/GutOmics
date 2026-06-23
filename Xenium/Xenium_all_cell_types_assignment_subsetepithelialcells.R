#Load required packages
library(Seurat)
#increase memory limit
options(future.globals.maxSize = 64 * 1024^3)  # 10 GiB
#Load Xenium directories for Patient 169
foregut1691 <- LoadXenium("output-XETG00621__0079276__foregut_3__20251219__212612", fov = "foregut-fov-1691", segmentations = "cell")
foregut1692 <- LoadXenium("output-XETG00621__0079276__foregut_4__20251219__212612", fov = "foregut-fov-1692", segmentations = "cell")
midgut1691 <- LoadXenium("output-XETG00621__0079276__midgut1__20251219__212612", fov = "midgut-fov-1691", segmentations = "cell")
midgut1692 <- LoadXenium("output-XETG00621__0079276__midgut2_hindgut__20251219__212612", fov = "midgut-fov-1692", segmentations = "cell")
colon1691 <- LoadXenium("output-XETG00621__0079276__colon_3__20251219__212612", fov = "colon-fov-1691", segmentations = "cell")
colon1692 <- LoadXenium("output-XETG00621__0079276__colon_4__20251219__212612", fov = "colon-fov-1692", segmentations = "cell")
#Assign relevant metadata labels
foregut1691$tissue<-rep("foregut")
foregut1691$patient<-rep("169")
foregut1691$view<-rep("fov1")
foregut1692$tissue<-rep("foregut")
foregut1692$patient<-rep("169")
foregut1692$view<-rep("fov2")
midgut1691$tissue<-rep("midgut")
midgut1691$patient<-rep("169")
midgut1691$view<-rep("fov1")
midgut1692$tissue<-rep("midgut")
midgut1692$patient<-rep("169")
midgut1692$view<-rep("fov2")
colon1691$tissue<-rep("colon")
colon1691$patient<-rep("169")
colon1691$view<-rep("fov1")
colon1692$tissue<-rep("colon")
colon1692$patient<-rep("169")
colon1692$view<-rep("fov2")
#merge tissues from different FOV but from same patients and same region together
foregut169<-merge(x=foregut1691, y=foregut1692)
midgut169<-merge(x=midgut1691, y=midgut1692)
colon169<-merge(x=colon1691, y=colon1692)
#remove cells with 0 counts
foregut169 <- subset(foregut169, subset = nCount_Xenium > 0)
midgut169 <- subset(midgut169, subset = nCount_Xenium > 0)
colon169 <- subset(colon169, subset = nCount_Xenium > 0)
#Join Layers, normalize, scale, run non-linear dimensional reduction and run clustering for foregut in Patient 169
foregut169<-JoinLayers(foregut169)
foregut169<-NormalizeData(foregut169)
foregut169 <- SCTransform(foregut169, assay = "Xenium")
foregut169 <- RunPCA(foregut169, npcs = 30, features = rownames(foregut169))
foregut169 <- RunUMAP(foregut169, dims = 1:30)
foregut169 <- FindNeighbors(foregut169, reduction = "pca", dims = 1:30)
foregut169 <- FindClusters(foregut169, resolution = 0.3)
foregut169.markers <- FindAllMarkers(foregut169, only.pos = TRUE)
saveRDS(foregut169,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut169.rds")
write.csv(foregut169.markers,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut169markers.csv")
#Join Layers, normalize, scale, run non-linear dimensional reduction and run clustering for midgut in Patient 169
midgut169<-JoinLayers(midgut169)
midgut169<-NormalizeData(midgut169)
midgut169 <- SCTransform(midgut169, assay = "Xenium")
midgut169 <- RunPCA(midgut169, npcs = 30, features = rownames(midgut169))
midgut169 <- RunUMAP(midgut169, dims = 1:30)
midgut169 <- FindNeighbors(midgut169, reduction = "pca", dims = 1:30)
midgut169 <- FindClusters(midgut169, resolution = 0.3)
midgut169.markers <- FindAllMarkers(midgut169, only.pos = TRUE)
saveRDS(midgut169,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut169.rds")
write.csv(midgut169.markers,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut169markers.csv")
#Join Layers, normalize, scale, run non-linear dimensional reduction and run clustering for colon in Patient 169
colon169<-JoinLayers(colon169)
colon169<-NormalizeData(colon169)
colon169 <- SCTransform(colon169, assay = "Xenium")
colon169 <- RunPCA(colon169, npcs = 30, features = rownames(colon169))
colon169 <- RunUMAP(colon169, dims = 1:30)
colon169 <- FindNeighbors(colon169, reduction = "pca", dims = 1:30)
colon169 <- FindClusters(colon169, resolution = 0.3)
colon169.markers <- FindAllMarkers(colon169, only.pos = TRUE)
saveRDS(colon169,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon169.rds")
write.csv(colon169.markers,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon169markers.csv")
#Load Xenium directories for Patient 168
foregut1681 <- LoadXenium("output-XETG00621__0079106__foregut_1__20251219__212612", fov = "foregut-fov-1681", segmentations = "cell")
foregut1682 <- LoadXenium("output-XETG00621__0079106__foregut_2__20251219__212612", fov = "foregut-fov-1682", segmentations = "cell")
midgut1681 <- LoadXenium("output-XETG00621__0079106__midgut_1__20251219__212612", fov = "midgut-fov-1681", segmentations = "cell")
midgut1682 <- LoadXenium("output-XETG00621__0079106__midgut_2__20251219__212612", fov = "midgut-fov-1682", segmentations = "cell")
colon1681 <- LoadXenium("output-XETG00621__0079106__colon_1__20251219__212612", fov = "colon-fov-1681", segmentations = "cell")
colon1682 <- LoadXenium("output-XETG00621__0079106__colon_2__20251219__212612", fov = "colon-fov-1682", segmentations = "cell")
#Assign relevant metadata labels
foregut1681$tissue<-rep("foregut")
foregut1681$patient<-rep("168")
foregut1681$view<-rep("fov1")
foregut1682$tissue<-rep("foregut")
foregut1682$patient<-rep("168")
foregut1682$view<-rep("fov2")
midgut1681$tissue<-rep("midgut")
midgut1681$patient<-rep("168")
midgut1681$view<-rep("fov1")
midgut1682$tissue<-rep("midgut")
midgut1682$patient<-rep("168")
midgut1682$view<-rep("fov2")
colon1681$tissue<-rep("colon")
colon1681$patient<-rep("168")
colon1681$view<-rep("fov1")
colon1682$tissue<-rep("colon")
colon1682$patient<-rep("168")
colon1682$view<-rep("fov2")
#merge tissues from different FOV but from same patients and same region together
foregut168<-merge(x=foregut1681, y=foregut1682)
midgut168<-merge(x=midgut1681, y=midgut1682)
colon168<-merge(x=colon1681, y=colon1682)
#remove cells with 0 counts
foregut168 <- subset(foregut168, subset = nCount_Xenium > 0)
midgut168 <- subset(midgut168, subset = nCount_Xenium > 0)
colon168 <- subset(colon168, subset = nCount_Xenium > 0)
#Join Layers, normalize, scale, run non-linear dimensional reduction and run clustering for foregut in Patient 168
foregut168<-JoinLayers(foregut168)
foregut168<-NormalizeData(foregut168)
foregut168 <- SCTransform(foregut168, assay = "Xenium")
foregut168 <- RunPCA(foregut168, npcs = 30, features = rownames(foregut168))
foregut168 <- RunUMAP(foregut168, dims = 1:30)
foregut168 <- FindNeighbors(foregut168, reduction = "pca", dims = 1:30)
foregut168 <- FindClusters(foregut168, resolution = 0.3)
foregut168.markers <- FindAllMarkers(foregut168, only.pos = TRUE)
saveRDS(foregut168,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut168.rds")
write.csv(foregut168.markers,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut168markers.csv")

#Join Layers, normalize, scale, run non-linear dimensional reduction and run clustering for midgut in Patient 168
midgut168<-JoinLayers(midgut168)
midgut168<-NormalizeData(midgut168)
midgut168 <- SCTransform(midgut168, assay = "Xenium")
midgut168 <- RunPCA(midgut168, npcs = 30, features = rownames(midgut168))
midgut168 <- RunUMAP(midgut168, dims = 1:30)
midgut168 <- FindNeighbors(midgut168, reduction = "pca", dims = 1:30)
midgut168 <- FindClusters(midgut168, resolution = 0.3)
midgut168.markers <- FindAllMarkers(midgut168, only.pos = TRUE)
saveRDS(midgut168,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut168.rds")
write.csv(midgut168.markers,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut168markers.csv")

#Join Layers, normalize, scale, run non-linear dimensional reduction and run clustering for colon in Patient 168
colon168<-JoinLayers(colon168)
colon168<-NormalizeData(colon168)
colon168 <- SCTransform(colon168, assay = "Xenium")
colon168 <- RunPCA(colon168, npcs = 30, features = rownames(colon168))
colon168 <- RunUMAP(colon168, dims = 1:30)
colon168 <- FindNeighbors(colon168, reduction = "pca", dims = 1:30)
colon168 <- FindClusters(colon168, resolution = 0.3)
colon168.markers <- FindAllMarkers(colon168, only.pos = TRUE)
saveRDS(colon168,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon168.rds")


#Assign clusters in foregut in Patient 169
foregut169<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut169.rds")
foregut169 <- FindClusters(foregut169, resolution = 2.0)
new.cluster.ids <- c("Endothelial", "Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Immune cells","Immune cells","Mesenchymal","Mesenchymal","Endothelial","Epithelial",
                     "Neuronal","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal","Mesenchymal","Endothelial","Epithelial",
                     "Epithelial","Mesenchymal","Neuronal","Mesenchymal","Endothelial","Epithelial","Epithelial","Immune cells","Epithelial","Mesenchymal",
                     "Mesenchymal","Mesenchymal","Endothelial","Epithelial","Mesenchymal","Immune cells","Epithelial","Mesenchymal","Mesenchymal","Immune cells",
                     "Mesenchymal","Neuronal","Endothelial","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal",
                     "Mesenchymal","Immune cells","Epithelial","Mesenchymal")
names(new.cluster.ids) <- levels(foregut169)
foregut169 <- RenameIdents(foregut169, new.cluster.ids)
DimPlot(foregut169, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(foregut169,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut169.rds")


#Assign clusters in midgut in Patient 169
midgut169<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut169.rds")
midgut169 <- FindClusters(midgut169, resolution = 2.0)
new.cluster.ids <- c("Immune cells","Mesenchymal","Endothelial","Mesenchymal","Immune cells","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Neuronal",
                     "Mesenchymal","Endothelial","Mesenchymal","Epithelial","Mesenchymal","Endothelial","Immune cells","Epithelial","Mesenchymal","Immune cells",
                     "Epithelial","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Immune cells","Endothelial","Endothelial",
                     "Endothelial","Epithelial","Endothelial","Mesenchymal","Endothelial","Epithelial","Immune cells","Epithelial","Mesenchymal","Mesenchymal",
                     "Mesenchymal","Epithelial","Mesenchymal","Immune cells","Neuronal","Epithelial","Mesenchymal","Mesenchymal","Immune cells","Epithelial",
                     "Immune cells","Neuronal","Epithelial","Immune cells")
names(new.cluster.ids) <- levels(midgut169)
midgut169 <- RenameIdents(midgut169, new.cluster.ids)
DimPlot(midgut169, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(midgut169,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut169.rds")

#Assign clusters in colon in Patient 169
colon169<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon169.rds")
colon169 <- FindClusters(colon169, resolution = 2.0)
new.cluster.ids <- c("Mesenchymal","Mesenchymal","Endothelial","Mesenchymal","Immune cells","Mesenchymal","Epithelial","Mesenchymal","Epithelial","Epithelial","Mesenchymal",
                     "Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Endothelial","Epithelial","Immune cells",
                     "Immune cells","Immune cells","Mesenchymal","Neuronal","Endothelial","Neuronal","Neuronal","Epithelial","Mesenchymal","Endothelial",
                     "Endothelial","Endothelial","Mesenchymal","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Epithelial","Mesenchymal","Epithelial",
                     "Endothelial","Endothelial","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Mesenchymal",
                     "Epithelial","Epithelial","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal")
names(new.cluster.ids) <- levels(colon169)
colon169 <- RenameIdents(colon169, new.cluster.ids)
DimPlot(colon169, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(colon169,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon169.rds")

#Assign clusters in foregut in Patient 168
foregut168<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut168.rds")
foregut168 <- FindClusters(foregut168, resolution = 2.0)
new.cluster.ids <- c("Mesenchymal","Mesenchymal","Mesenchymal","Endothelial","Mesenchymal","Immune cells","Endothelial","Endothelial","Immune cells","Epithelial","Epithelial",
                     "Mesenchymal","Epithelial","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Epithelial","Neuronal","Immune cells","Mesenchymal",
                     "Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Immune cells","Mesenchymal",
                     "Epithelial","Mesenchymal","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Endothelial","Epithelial",
                     "Neuronal","Endothelial","Epithelial","Immune cells","Epithelial","Mesenchymal","Epithelial","Epithelial","Endothelial","Endothelial",
                     "Mesenchymal")
names(new.cluster.ids) <- levels(foregut168)
foregut168 <- RenameIdents(foregut168, new.cluster.ids)
DimPlot(foregut168, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(foregut168,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut168.rds")

#Assign clusters in midgut in Patient 168
midgut168<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut168.rds")
midgut168 <- FindClusters(midgut168, resolution = 2.0)
new.cluster.ids <- c("Mesenchymal","Endothelial","Mesenchymal","Neuronal","Mesenchymal","Mesenchymal","Mesenchymal","Immune cells","Endothelial","Mesenchymal","Immune",
                     "Epithelial","Epithelial","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Epithelial","Epithelial","Mesenchymal","Epithelial",
                     "Endothelial","Endothelial","Mesenchymal","Endothelial","Epithelial","Epithelial","Mesenchymal","Mesenchymal","Endothelial","Neuronal",
                     "Epithelial","Immune cells","Epithelial","Mesenchymal","Epithelial","Epithelial","Epithelial","Immune cells","Mesenchymal","Mesenchymal",
                     "Epithelial","Immune cells","Epithelial","Mesenchymal","Epithelial","Mesenchymal","Immune cells","Immune cells","Immune cells")
names(new.cluster.ids) <- levels(midgut168)
midgut168 <- RenameIdents(midgut168, new.cluster.ids)
DimPlot(midgut168, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(midgut168,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut168.rds")

#Assign clusters in colon in Patient 168
colon168<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon168.rds")
colon168 <- FindClusters(colon168, resolution = 2.0)
new.cluster.ids <- c("Mesenchymal","Mesenchymal","Mesenchymal","Endothelial","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Immune cells","Endothelial",
                     "Neuronal","Immune cells","Mesenchymal","Mesenchymal","Endothelial","Epithelial","Mesenchymal","Endothelial","Epithelial","Immune cells",
                     "Epithelial","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal","Neuronal","Epithelial","Neuronal","Mesenchymal",
                     "Endothelial","Mesenchymal","Mesenchymal","Immune cells","Epithelial","Immune cells","Epithelial","Mesenchymal","Neuronal","Mesenchymal",
                     "Immune cells","Epithelial","Mesenchymal","Epithelial","Immune cells","Neuronal")
names(new.cluster.ids) <- levels(colon168)
colon168 <- RenameIdents(colon168, new.cluster.ids)
DimPlot(colon168, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(colon168,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon168.rds")

#Subset only epithelial cells in foregut in Patient 169 and normalize, scale, and cluster 
foregut169<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut169.rds")
epiforegut169<-subset(foregut169, idents=="Epithelial")
epiforegut169<-NormalizeData(epiforegut169)
epiforegut169<-SCTransform((epiforegut169), assay = "Xenium")
epiforegut169 <- RunPCA(epiforegut169, npcs = 30, features = rownames(epiforegut169))
epiforegut169 <- RunUMAP(epiforegut169, dims = 1:30)
epiforegut169 <- FindNeighbors(epiforegut169, reduction = "pca", dims = 1:30)
epiforegut169 <- FindClusters(epiforegut169, resolution = 2.0)
saveRDS(epiforegut169 ,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/epiforegut169.rds")

#Subset only epithelial cells in midgut in Patient 169 and normalize, scale, and cluster 
midgut169<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut169.rds")
epimidgut169<-subset(foregut169, idents=="Epithelial")
epimidgut169<-NormalizeData(epimidgut169)
epimidgut169<-SCTransform((epimidgut169), assay = "Xenium")
epimidgut169 <- RunPCA(epimidgut169, npcs = 30, features = rownames(epimidgut169))
epimidgut169 <- RunUMAP(epimidgut169, dims = 1:30)
epimidgut169 <- FindNeighbors(epimidgut169, reduction = "pca", dims = 1:30)
epimidgut169 <- FindClusters(epimidgut169, resolution = 2.0)
saveRDS(epimidgut169 ,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/epimidgut169.rds")

#Subset only epithelial cells in colon in Patient 169 and normalize, scale, and cluster 
colon169<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon169.rds")
epicolon169<-subset(foregut169, idents=="Epithelial")
epicolon169<-NormalizeData(epicolon169)
epicolon169<-SCTransform((epicolon169), assay = "Xenium")
epicolon169 <- RunPCA(epicolon169, npcs = 30, features = rownames(epicolon169))
epicolon169 <- RunUMAP(epicolon169, dims = 1:30)
epicolon169 <- FindNeighbors(epicolon169, reduction = "pca", dims = 1:30)
epicolon169 <- FindClusters(epicolon169, resolution = 2.0)
saveRDS(epicolon169 ,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/epicolon169.rds")

#Subset only epithelial cells in foregut in Patient 168 and normalize, scale, and cluster 
foregut168<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/foregut168.rds")
epiforegut168<-subset(foregut169, idents=="Epithelial")
epiforegut168<-NormalizeData(epiforegut168)
epiforegut168<-SCTransform((epiforegut168), assay = "Xenium")
epiforegut168 <- RunPCA(epiforegut168, npcs = 30, features = rownames(epiforegut168))
epiforegut168 <- RunUMAP(epiforegut168, dims = 1:30)
epiforegut168 <- FindNeighbors(epiforegut168, reduction = "pca", dims = 1:30)
epiforegut168 <- FindClusters(epiforegut168, resolution = 2.0)
saveRDS(epiforegut168 ,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/epiforegut168.rds")

#Subset only epithelial cells in midgut in Patient 168 and normalize, scale, and cluster 
midgut168<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/midgut168.rds")
epimidgut168<-subset(foregut169, idents=="Epithelial")
epimidgut168<-NormalizeData(epimidgut168)
epimidgut168<-SCTransform((epimidgut168), assay = "Xenium")
epimidgut168 <- RunPCA(epimidgut168, npcs = 30, features = rownames(epimidgut168))
epimidgut168 <- RunUMAP(epimidgut168, dims = 1:30)
epimidgut168 <- FindNeighbors(epimidgut168, reduction = "pca", dims = 1:30)
epimidgut168 <- FindClusters(epimidgut168, resolution = 2.0)
saveRDS(epimidgut168 ,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/epimidgut168.rds")

#Subset only epithelial cells in colon in Patient 168 and normalize, scale, and cluster 
colon168<-readRDS("/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/colon168.rds")
epicolon168<-subset(foregut169, idents=="Epithelial")
epicolon168<-NormalizeData(epicolon168)
epicolon168<-SCTransform((epicolon168), assay = "Xenium")
epicolon168 <- RunPCA(epicolon168, npcs = 30, features = rownames(epicolon168))
epicolon168 <- RunUMAP(epicolon168, dims = 1:30)
epicolon168 <- FindNeighbors(epicolon168, reduction = "pca", dims = 1:30)
epicolon168 <- FindClusters(epicolon168, resolution = 2.0)
saveRDS(epicolon168 ,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/epicolon168.rds")

#merge all epithelial cells
merged<-merge(epiforegut169, y=c(epimidgut169, epicolon169, epiforegut168, epimidgut168, epicolon168))

#save epithelial cells
saveRDS(merged,"/athena/chenlab/scratch/jjv4001/Xeniumfetalintestine/mergedepithelial.rds")
