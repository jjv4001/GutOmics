#Load Seurat object and relevant packages
library(Seurat)
library(readxl)
#Load age-assigned metadata
metadata<-read_xlsx("/athena/chenlab/scratch/jjv4001/gutrevision/metadata.xlsx")
#Load epithelial data corresponding to both adult and fetal samples
all<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/assignedintegrated.RDS")
#Assign metadata to the dataset
all$SampleType<-metadata$SampleType
#subset fetal samples
fetal<-subset(all, SampleType=="fetal")
#assign Ages of fetal samples as identities 
Idents(fetal)<-fetal$Age
#subset age-matched samples ranging between 11 and 17 weeks post-conceptioon
fetalcorrected<-subset(fetal, idents=c("101 days","11.2 weeks","11.9 weeks","12 weeks","12.2 weeks","122 days","12weeks","13weeks","15 weeks","15weeks","16 weeks","16weeks",
                                       "17 weeks","17weeks","80 days"))
#rename cell identities
new.cluster.ids <- c("12.2 weeks", "11.2 weeks", "11.9 weeks", "12 weeks", "13 weeks", "17 weeks",
                     "12 weeks", "16 weeks", "15 weeks","17.4 weeks","11.4 weeks","14.4 weeks","17 weeks","15 weeks","16 weeks")
names(new.cluster.ids) <- levels(fetalcorrected)
fetalcorrected <- RenameIdents(fetalcorrected, new.cluster.ids)
new_levels <- c("11.2 weeks", "11.4 weeks", "11.9 weeks", "12 weeks","12.2 weeks","13 weeks","14.4 weeks","15 weeks","16 weeks","17 weeks","17.4 weeks")
Idents(fetalcorrected) <- factor(Idents(fetalcorrected), levels = new_levels)
#save subsetted age-matched fetal dataset
saveRDS(fetalcorrected, "/athena/chenlab/scratch/jjv4001/gutrevision/fetalgutcorrected.rds")

#visualize elbowplot of fetal samples
ElbowPlot(fetalcorrected)
#cluster fetal epithelial cells
fetalcorrected <- FindNeighbors(fetalcorrected, dims = 1:15)
fetalcorrected <- FindClusters(fetalcorrected, resolution = 3.5)
#choose features and visualize feature expression
features<-c("LGR5","MYC","RGMB","TOP2A","UBE2C","FABP1","ANPEP","CLCA1","MUC2","CHGA","NEUROD1","LYZ","PRSS2","LYZ1")
DotPlot(fetalcorrected, features=features)+RotatedAxis()

#assign fetal epithelial cells
new.cluster.ids <- c("Enterocytes","Enterocytes", "Enterocytes", "Enterocytes","Enterocytes", "Enterocytes", 
                     "Stem cells", "Enterocytes","Enterocytes","Enterocytes","TA cells",
                     "Enterocytes","Enterocytes","Enterocytes","Enterocytes","Enterocytes",
                     "Enterocytes","Stem cells","Enterocytes","Enterocytes","Enterocytes",
                     "Enterocytes","Enterocytes","Enterocytes","Enterocytes","Enterocytes",
                     "Enterocytes","Enterocytes","Enterocytes","Enteroendocrine cells","Enterocytes",
                     "Stem cells","Stem cells","Goblet cells","TA cells","Stem cells",
                     "Enterocytes","Enterocytes","Paneth cells","Enterocytes","TA cells",
                     "Enterocytes","TA cells","Enterocytes","Stem cells","TA cells",
                     "TA cells","Stem cells","Enterocytes","Enterocytes","Stem cells",
                     "Enterocytes","Enterocytes","Enteroendocrine cells","Stem cells","Enterocytes",
                     "Enterocytes","Enteroendocrine cells","Enteroendocrine cells","Enterocytes","Enterocytes",
                     "Enterocytes","Enterocytes","Enterocytes","Goblet cells","Enterocytes",
                     "Goblet cells","Enterocytes","Enterocytes")

names(new.cluster.ids) <- levels(fetalcorrected)
fetalcorrected <- RenameIdents(fetalcorrected, new.cluster.ids)
#change levels of fetal epithelial cells
new_levels <- c("Paneth cells", "Enteroendocrine cells", "Goblet cells", "Enterocytes","TA cells","Stem cells")
Idents(fetalcorrected) <- factor(Idents(fetalcorrected), levels = new_levels)
#Load viridis package
library(viridis)
#Visualize dotplot with markers
DotPlot(fetalcorrected, features=features, dot.scale=20)+RotatedAxis()+scale_color_viridis_c(option = "A")
#Assign Paneth cells in the large intestine to enterocytes
fetalcorrected$identity<-Idents(fetalcorrected)
fetalcorrected$identity[fetalcorrected$identity == "Paneth cells" & fetalcorrected$region == "large intestine"] <- "Enterocytes"
Idents(fetalcorrected)<-fetalcorrected$identity
#Visualize dotplot with markers
DotPlot(fetalcorrected, features=features, dot.scale=20)+RotatedAxis()+scale_color_viridis_c(option = "A")
#Visualize umap of epithelial cells
DimPlot(fetalcorrected, reduction = "umap.integrated", group.by="Age") +
  scale_color_viridis_d(option = "B") +
  theme_bw(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )+labs(x = "UMAP_1", y = "UMAP_2")
#Visualize dotplot with markers
DotPlot(fetalcorrected, features=features, dot.scale=20)+RotatedAxis()+scale_color_viridis_c(option = "A") 

#subset enteroendocrine cells from fetal epithelial cells
merged.gut<-subset(fetalcorrected, idents="Enteroendocrine cells")

# Create a a new gene column with "gene-neg" as entry
merged.gut$gene <- "gene-neg"

# Assign different enteroendocrine cell subtypes
Kcell <- WhichCells(merged.gut, expression = GIP > 1)
merged.gut$gene[Kcell] <- "K cell"

MXAcell <- WhichCells(merged.gut, expression = GHRL > 1 & MLN > 1)
merged.gut$gene[MXAcell] <- "M/X/A cell"

Icell <- WhichCells(merged.gut, expression = CCK > 1)
merged.gut$gene[Icell] <- "I cell"

Dcell <- WhichCells(merged.gut, expression = SST > 1)
merged.gut$gene[Dcell] <- "D cell"

Lcell <- WhichCells(merged.gut, expression = GCG > 1 & PYY > 1)
merged.gut$gene[Lcell] <- "L cell"

Ncell <- WhichCells(merged.gut, expression = NTS > 1)
merged.gut$gene[Ncell] <- "N cell"

ECcell <- WhichCells(merged.gut, expression = TPH1 > 1)
merged.gut$gene[ECcell] <- "EC cell"

#Assign eec cell sub-type identities
Idents(merged.gut) <- merged.gut$gene

#Save fetal eec metadata 
write.csv(merged.gut@meta.data, "AllCellTypes.csv")



#subset adult samples
adult<-subset(all, SampleType=="adult")
#normalize data 
adult <- NormalizeData(adult)
#find highly variable features 
adult <- FindVariableFeatures(adult, selection.method = "vst", nfeatures = 2000)
#scale data
adult <- ScaleData(adult)
#run linear dimensional reduction
adult <- RunPCA(adult, features = VariableFeatures(object = adult))
#cluster cells
adult <- FindNeighbors(adult, dims = 1:15)
adult <- FindClusters(adult, resolution = 1.0)
#choose features and visualize feature expression
features<-c("LGR5","MYC","RGMB","TOP2A","UBE2C","ALPI","FABP1","CLCA1","MUC2","SPINK4","CHGA","NEUROD1","INSM1","LYZ","DEFA5","DEFA6")
#Visualize dotplot with markers
DotPlot(adult, features=features, dot.scale=20)+RotatedAxis()+scale_color_viridis_c(option = "A", direction=-1)
#assign adult epithelial cells
new.cluster.ids <- c("Enterocytes","Stem cells","Enterocytes","Stem cells","Enterocytes","Enterocytes",
                     "Enterocytes","Goblet cells","Enterocytes","Stem cells","Enterocytes",
                     "Enterocytes","Paneth cells","Enterocytes","Enterocytes","Enterocytes",
                     "Goblet cells","TA cells","TA cells","Stem cells","Goblet cells",
                     "Enterocytes","Enterocytes","Enterocytes","Enterocytes","Enterocytes",
                     "Goblet cells","Enterocytes","Enterocytes","Enterocytes","Goblet cells",
                     "Enterocytes","Enterocytes","TA cells","Goblet cells","Stem cells",
                     "Enterocytes","EECs","Goblet cells","Enterocytes","Enterocytes",
                     "Enterocytes","Goblet cells","Enterocytes","Enterocytes","Goblet cells",
                     "Enterocytes","Paneth cells","Goblet cells","Enterocytes")
names(new.cluster.ids) <- levels(adult)
adult <- RenameIdents(adult, new.cluster.ids)
new_levels <- c("Stem cells","TA cells","Enterocytes","Goblet cells","EECs","Paneth cells")
Idents(adult) <- factor(Idents(adult), levels = new_levels)
new_levels <- c("Paneth cells","EECs","Goblet cells","Enterocytes","TA cells","Stem cells")
Idents(adult) <- factor(Idents(adult), levels = new_levels)
new_levels <- c("25 to 30 years","45 to 50 years","50 to 55 years","60 to 65 years","65 to 70 years", "70 to 75 years")
Idents(adult) <- factor(Idents(adult), levels = new_levels)

#Visualize umap of epithelial cells
DimPlot(adult, reduction = "umap.integrated") +
  scale_color_viridis_d(option = "B") +
  theme_bw(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )+labs(x = "UMAP_1", y = "UMAP_2")


