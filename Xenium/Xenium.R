Week20<-merge(Reg1, Reg4)
Week17<-merge(Reg6, Reg9)
merge2<-merge(Week17, Week20)


library(Seurat)
merge2 <- SCTransform(merge2, assay = "Xenium")
merge2 <- RunPCA(merge2, npcs = 30, features = rownames(merge2))
merge2 <- RunUMAP(merge2, dims = 1:30, reduction="pca")
merge2 <- FindNeighbors(merge2, dims = 1:30)
merge2 <- RunHarmony(object = merge2, group.by.vars = 'sample')
merge2 <- RunUMAP(merge2, dims = 1:30, reduction="harmony")
merge2 <- FindNeighbors(merge2, dims = 1:30)
merge2 <- FindClusters(merge2, resolution = 4.15)


new.cluster.ids <- c("Enterocytes","Enterocytes", "2","Enterocytes", "4", "5",
                     "6", "7", "8", "9", "10",
                     "Enterocytes","12", "13", "14", "15", "16",
                     "17", "18", "19", "20", "21", "22",
                     "23", "24", "25","26", "27", "28",
                     "29", "30", "31", "32", "33", "34", "35", "Enterocytes",
                     "37", "38", "39", "40", "41","42","43",
                     "Enterocytes","45","46","47","48","49","50","51",
                     "52","53","54","55","56","57","58","59","60",
                     "61","62","63","64","65","66","67","68","69","70",
                     "71","72","73","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)


new.cluster.ids <- c("Enterocytes", "2", "4", "5",
                     "6", "7", "8", "9", "10",
                     "12", "13", "14", "15", "16",
                     "17", "Stem cells", "Stem cells", "20", "21", "22",
                     "23", "24", "25","26", "27", "28",
                     "29", "EECs", "31", "32", "33", "Stem cells", "35",
                     "37", "38", "39", "40", "41","42","43",
                     "45","46","47","48","49","50","51",
                     "EECs","53","EECs","55","56","57","58","59","60",
                     "61","62","63","64","65","66","67","68","69","70",
                     "71","72","73","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "2", "4", "5",
                     "6", "7", "8", "9", "10",
                     "12", "13", "14", "15", "16",
                     "17", "Stem cells", "20", "21", "22",
                     "23", "Stem cells", "25","Enterocytes", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "37", "38", "39", "40", "41","42","43",
                     "45","46","47","48","49","50","Enterocytes",
                     "53","55","56","57","58","EECs","60",
                     "61","62","63","64","65","66","67","68","69","70",
                     "Stem cells","72","Stem cells","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "2", "4", "5",
                     "6", "7", "8", "9", "10",
                     "12", "Goblet cells", "14", "Goblet cells", "16",
                     "17", "Stem cells", "20", "21", "22",
                     "23", "25", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "37", "38", "39", "Goblet cells", "41","42","43",
                     "45","46","47","48","49","50",
                     "53","55","56","57","58","Goblet cells",
                     "61","62","63","64","65","66","67","68","69","70",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "2", "4", "5",
                     "6", "7", "8", "9", "10",
                     "12", "Goblet cells", "14", "16",
                     "17", "Stem cells", "20", "21", "22",
                     "23", "25", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "37", "38", "39", "41","42","43",
                     "45","Goblet cells","Goblet cells","48","49","50",
                     "53","55","Goblet cells","57","58",
                     "61","62","63","64","65","66","67","68","69","70",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "2", "4", "5",
                     "6", "7", "8", "9", "10",
                     "12", "Goblet cells", "14", "16",
                     "17", "Stem cells", "20", "21", "22",
                     "23", "25", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "37", "38", "39", "41","42","43",
                     "45","Enterocytes","Enterocytes","TA",
                     "53","55","57","58",
                     "61","62","63","64","65","66","67","68","69","EECs",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "2", "4", "5",
                     "6", "7", "8", "9", "10",
                     "12", "Goblet cells", "14", "16",
                     "17", "Stem cells", "20", "21", "22",
                     "23", "25", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "Enterocytes", "38", "39", "41","42","43",
                     "Enterocytes","TA",
                     "53","55","57","Enterocytes",
                     "61","62","63","64","65","66","67","68","69",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "2", "Enterocytes", "5",
                     "6", "7", "8", "Stem cells", "10",
                     "12", "Goblet cells", "Stem cells", "Stem cells",
                     "17", "Stem cells", "20", "21", "22",
                     "23", "25", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "38", "39", "41","42","43",
                     "TA",
                     "53","55","57",
                     "61","62","63","64","65","66","67","68","69",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "Enterocytes", "5",
                     "6", "Enterocytes", "8", "Stem cells", "Enterocytes",
                     "Enterocytes", "Goblet cells",
                     "17", "20", "21", "22",
                     "23", "25", "27", "28",
                     "29", "EECs", "31", "32", "33", "35",
                     "38", "39", "41","42","43",
                     "TA",
                     "53","55","57",
                     "61","62","63","64","65","66","67","68","69",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)


new.cluster.ids <- c("Enterocytes", "5",
                     "6", "8", "Stem cells",
                     "Goblet cells",
                     "17", "20", "21", "22",
                     "23", "25", "27", "Enterocytes",
                     "Enterocytes", "EECs", "31", "32", "33", "35",
                     "38", "39", "41","42","43","TA",
                     "53","55","57",
                     "61","Enterocytes","63","64","Enterocytes","66","67","68","69",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "Enterocytes",
                     "6", "8", "Stem cells",
                     "Goblet cells",
                     "17", "20", "Enterocytes", "22",
                     "23", "25", "27",
                     "EECs", "31", "32", "33", "35",
                     "Stem cells", "39", "41","42","43","TA",
                     "53","55","57",
                     "61","63","64","66","67","68","69",
                     "72","74","75")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes",
                     "6", "8", "Stem cells",
                     "Goblet cells",
                     "17", "20", "22",
                     "23", "25", "27",
                     "EECs", "31", "32", "33", "35",
                     "39", "41","42","43","TA",
                     "53","55","57",
                     "61","Enterocytes","Goblet cells","66","67","68","69",
                     "72","74","EECs")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes",
                     "6", "8", "Stem cells",
                     "Goblet cells",
                     "17", "20", "22",
                     "23", "25", "27",
                     "EECs", "31", "32", "33", "35",
                     "39", "41","42","Enterocytes","TA",
                     "53","55","57",
                     "61","66","67","Enterocytes","69",
                     "72","Enterocytes")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", 
                     "6", "8", "Stem cells",
                     "Goblet cells",
                     "17", "20", "22",
                     "23", "25", "Enterocytes",
                     "EECs", "31", "32", "Goblet cells", "Goblet cells",
                     "39", "41","Enterocytes","TA",
                     "53","55","57",
                     "61","66","67","69",
                     "72")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", 
                     "Enterocytes", "8", "Stem cells",
                     "Goblet cells",
                     "17", "20", "22",
                     "23", "25",
                     "EECs", "31", "32",
                     "Enterocytes", "41","TA",
                     "53","55","57",
                     "61","66","67","69",
                     "72")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", 
                     "Enterocytes", "Stem cells",
                     "Goblet cells",
                     "Enterocytes", "20", "22",
                     "23", "25",
                     "EECs", "31", "32",
                     "41","TA",
                     "53","55","57",
                     "Enterocytes","66","67","69",
                     "72")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", 
                     "Stem cells",
                     "Goblet cells",
                     "Enterocytes", "22",
                     "23", "Enterocytes",
                     "EECs", "31", "32",
                     "41","TA",
                     "53","55","Enterocytes",
                     "66","Enterocytes","69",
                     "72")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", 
                     "Stem cells",
                     "Goblet cells",
                     "Enterocytes",
                     "Enterocytes",
                     "EECs", "31", "32",
                     "Stem cells","TA",
                     "53","55",
                     "66","69",
                     "72")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", 
                     "Stem cells",
                     "Goblet cells",
                     "EECs", "Goblet cells", "Goblet cells","TA",
                     "Goblet cells","Enterocytes",
                     "Enterocytes","Enterocytes",
                     "Enterocytes")
names(new.cluster.ids) <- levels(merge2)
merge2 <- RenameIdents(merge2, new.cluster.ids)
features<-c("CHGA","RAP1GAP2","RIMBP2",
            "RBP2","APOA1","APOC3","KHK",
            "LGR5","ASCL2","RGMB","PREP",
            "CLCA1","MUC2","TCEA3","ATOH1",
            "NUSAP1","TOP2A","UBE2C")
DotPlot(
  merge2,
  features=features)+RotatedAxis()
epithelial<-merge2
cropped.coords <- Crop(merge2[["fov1"]], x = c(2000, 3500), y = c(2000, 2750), coords = "plot")
cropped.coords1 <- Crop(epithelial[["fov4"]], x = c(2150, 2750), y = c(1000, 1800), coords = "plot")
merge2[["zoom"]] <- cropped.coords
epithelial[["zoom"]] <- cropped.coords1
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(merge2[["zoom"]]) <- "segmentation"
DefaultBoundary(epithelial[["zoom"]]) <- "segmentation"
options(future.globals.maxSize = 8000 * 1024^2)
custom_colors <- c("#FFCC99", "#FF99FF", "#99CCFF", "#99FFCC", "#FF6666")
ImageDimPlot(merge2, fov = "zoom", axes = TRUE, border.color = "NA", border.size = 0.1, cols = custom_colors,dark.background = FALSE,coord.fixed = FALSE, nmols = 10000, flip_xy = TRUE)
ImageFeaturePlot(merge2,fov ="zoom",features = c("RBP2"), size = 2, cols = c("blue", "red"), min.cutoff = c("1"), max.cutoff=c("4"), dark.background = FALSE)

features<-c("APOA4","APOB","AADAC","AFP","CCL25","CPVL","CPS1","ALDOB","ADCY2","CYP2W1","CD68","APOE","AIG1","ARID3A","ALDH1A1","ACE2","CDKN1C","CHST9","CTSH","DAB1","ADGRG7","ATRN","ATP11A","BDH2")

p1<-ImageFeaturePlot(merge2,fov ="zoom",features = c("AFP"),max.cutoff =2,  size = 2, cols = c("grey", "#ff6699"), dark.background = FALSE, border.color = NA)

p2<-ImageFeaturePlot(epithelial,fov ="zoom",features = c("AFP"),max.cutoff =2,  size = 2, cols = c("grey", "#99ccff"), dark.background = FALSE, border.color = NA)

pdf(paste(".pdf", sep=""), width=32, height=24)
for (gene in features) {
  p1<-ImageFeaturePlot(merge2,fov ="zoom",features = features,max.cutoff =2,  size = 2, cols = c("grey", "#ff6699"), dark.background = FALSE, border.color = NA)+ ggtitle("Duodenum")
  p2<-ImageFeaturePlot(epithelial,fov ="zoom",features = features,max.cutoff =2,  size = 2, cols = c("grey", "#99ccff"), dark.background = FALSE, border.color = NA)+ ggtitle("Colon")
  combined_plot<-p1+p2+plot_layout(ncol=2)+ plot_annotation(title = "AADAC",
                                                            theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold.italic", color = "black")))
  ggsave(
    filename = paste0(gene,".pdf"),  # Create a unique file for each gene
    plot = combined_plot,
    width = 15,
    height = 7
  )
}

saveRDS(merge2, "/Users/chenlab/scratch/jjv4001/16445-selected/xeniumepithelialfinal.rds")

EECs<-subset(merge2, idents="EECs")

EECs <- SCTransform(EECs, assay = "Xenium")
EECs <- RunPCA(EECs, npcs = 30, features = rownames(EECs))
EECs <- RunUMAP(EECs, dims = 1:30, reduction="pca")
EECs <- FindNeighbors(EECs, dims = 1:30)
EECs <- FindClusters(EECs, resolution = 1.0)
features<-c("GIP","PYY","TPH1","SST","CCK")
DotPlot(EECs, features=features)

new.cluster.ids <- c("EC cells", "EEC Progenitors","EC cells","EC cells","EC cells",
                     "K cells","EC cells","I cells","EC cells","D cells","D cells",
                     "EC cells", "EC cells","EC cells","K cells","EC cells","L cells")

names(new.cluster.ids) <- levels(EECs)
EECs <- RenameIdents(EECs, new.cluster.ids)
EECs$tissue <- EECs$sample
EECs$tissue[EECs$tissue == "Reg1"] <- "Duodenum"
EECs$tissue[EECs$tissue == "Reg4"] <- "Colon"
EECs$tissue[EECs$tissue == "Reg6"] <- "Duodenum"
EECs$tissue[EECs$tissue == "Reg9"] <- "Colon"
EECs$identity<-EECs@active.ident
EECs$celltype<-paste0(EECs$identity, EECs$tissue)

new.cluster.ids <- c("D cells", "EC cells","I cells","K cells","EEC Progenitors",
                     "L cells","EC cells","D cells","EEC Progenitors","L cells","EC cells",
                     "L cells")
names(new.cluster.ids) <- levels(EECs)
EECs <- RenameIdents(EECs, new.cluster.ids)

saveRDS(merge2, "/Users/chenlab/scratch/jjv4001/16445-selected/eecfinal.rds")