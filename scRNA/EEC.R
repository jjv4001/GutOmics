
#Load Seurat
library(Seurat)
#Load the Seurat dataset
merged.gut<-readRDS("/home/fs01/jjv4001/integrated_iter_2_seurat.rds")
#Subset the enteroendocrine cells
merged.gut <- subset(merged.gut, idents="EECs")
# normalize data
merged.gut <- NormalizeData(merged.gut, normalization.method="LogNormalize", scale.factor=10000)
# identification of highly variable features 
merged.gut <- FindVariableFeatures(merged.gut, selection.method="vst", nfeatures=400)
# scale data
merged.gut <- ScaleData(merged.gut, features=rownames(merged.gut))
# run linear dimensional reduction
merged.gut <- RunPCA(merged.gut , features=VariableFeatures(object=merged.gut))
# run non-linear dimensional reduction
merged.gut <- RunUMAP(merged.gut, dims=1:10, n.components=3)
merged.gut@misc$umap3d <- merged.gut@reductions$umap
merged.gut <- RunUMAP(merged.gut, dims=1:10, reduction="pca")
# cluster cells
merged.gut <- FindNeighbors(merged.gut, dims=1:10)
# find cells with high GIP expression (K cells)
Kcell<-WhichCells(merged.gut, expression=GIP >4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% Kcell, "K cell","gene-neg")
write.csv(merged.gut@meta.data, "Kcell.csv")
# find cells with high GHRL or MLN expression (M/X/A cells)
MXAcell<-WhichCells(merged.gut, expression=GHRL >4 & MLN>4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% MXAcell, "M/X/A cell", "gene-neg")
write.csv(merged.gut@meta.data, "MXAcell.csv")
# find cells with high CCK expression (I cells)
Icell<-WhichCells(merged.gut, expression=CCK >4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% Icell, "I cell", "gene-neg")
write.csv(merged.gut@meta.data, "Icell.csv")
# find cells with high SST expression (D cells)
Dcell<-WhichCells(merged.gut, expression=SST >4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% Dcell, "D cell", "gene-neg")
write.csv(merged.gut@meta.data, "Dcell.csv")
# find cells with high GCG/PYY expression (L cells)
Lcell<-WhichCells(merged.gut, expression=GCG >4 & PYY>4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% Lcell, "L cell", "gene-neg")
write.csv(merged.gut@meta.data, "Lcell.csv")
# find cells with high NTS expression (N cells)
Ncell<-WhichCells(merged.gut, expression=NTS >4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% Ncell, "N cell", "gene-neg")
write.csv(merged.gut@meta.data, "Ncell.csv")
# find cells with high TPH1 expression (EC cells)
ECcell<-WhichCells(merged.gut, expression=TPH1 >4)
merged.gut$gene <- ifelse(colnames(merged.gut) %in% ECcell, "EC cell", "gene-neg")
write.csv(merged.gut@meta.data, "ECcell.csv")
#assign cell identities
Idents(merged.gut)<-merged.gut$gene
#save Seurat analysis
saveRDS(merged.gut, file=paste("EEChormoneassignment", sep=""))
