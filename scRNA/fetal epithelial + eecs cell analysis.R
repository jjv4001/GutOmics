#Load Seurat object and relevant packages
library(Seurat)
library(readxl)
#Load age-assigned metadata
metadata<-read_xlsx("/athena/chenlab/scratch/jjv4001/gutrevision/metadata.xlsx")
#Load epithelial data corresponding to both adult and fetal samples
all<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/assignedintegrated.rds")
#Assign metadata to the dataset
all$SampleType<-metadata$SampleType
#subset fetal samples
fetal<-subset(all, SampleType=="fetal")
#assign Ages of fetal samples as identities 
Idents(fetal)<-fetal$Age
#subset age-matched samples ranging between 11 and 17 weeks post-conceptioon
fetalcorrected<-subset(fetal, idents=c("101 days","11.2 weeks","11.9 weeks","12 weeks","12.2 weeks","122 days","12weeks","13weeks","15 weeks","15weeks","16 weeks","16weeks",
                                     "17 weeks","17weeks","80 days"))
#rename cell identities and assigned identities 
new.cluster.ids <- c("12.2 weeks", "11.2 weeks", "11.9 weeks", "12 weeks", "13 weeks", "17 weeks",
                     "12 weeks", "16 weeks", "15 weeks","17.4 weeks","11.4 weeks","14.4 weeks","17 weeks","15 weeks","16 weeks")
names(new.cluster.ids) <- levels(fetalcorrected)
fetalcorrected <- RenameIdents(fetalcorrected, new.cluster.ids)
new_levels <- c("11.2 weeks", "11.4 weeks", "11.9 weeks", "12 weeks","12.2 weeks","13 weeks","14.4 weeks","15 weeks","16 weeks","17 weeks","17.4 weeks")
Idents(fetalcorrected) <- factor(Idents(fetalcorrected), levels = new_levels)
#save subsetted age-matched fetal dataset
saveRDS(fetalcorrected, "/athena/chenlab/scratch/jjv4001/gutrevision/fetalgutcorrected.rds")
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

#Assignment of EECs according to hormone expression
# 1. Initialize as "gene-neg"
fetalcorrected$gene <- "gene-neg"

# 2. Assign eec cell types
Kcell <- WhichCells(fetalcorrected, expression = GIP > 1)
fetalcorrected$gene[Kcell] <- "K cell"

MXAcell <- WhichCells(fetalcorrected, expression = GHRL > 1 & MLN > 1)
fetalcorrected$gene[MXAcell] <- "M/X/A cell"

Icell <- WhichCells(fetalcorrected, expression = CCK > 1)
fetalcorrected$gene[Icell] <- "I cell"

Dcell <- WhichCells(fetalcorrected, expression = SST > 1)
fetalcorrected$gene[Dcell] <- "D cell"

Lcell <- WhichCells(fetalcorrected, expression = GCG > 1 & PYY > 1)
fetalcorrected$gene[Lcell] <- "L cell"

Ncell <- WhichCells(fetalcorrected, expression = NTS > 1)
fetalcorrected$gene[Ncell] <- "N cell"

ECcell <- WhichCells(fetalcorrected, expression = TPH1 > 1)
fetalcorrected$gene[ECcell] <- "EC cell"

#Set identity class
Idents(fetalcorrected) <- fetalcorrected$gene

# 4. (Optional) Save final metadata
write.csv(fetalcorrected@meta.data, "AllCellTypes.csv")
