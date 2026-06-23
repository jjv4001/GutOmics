#Load required packages 
library(Seurat)
library(dplyr)
library(FNN)
library(ggplot2)
library(dplyr)
library(purrr)

#Load Seurat object with epithelial cells assigned 
newgutassigned<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Xeniumfetalintestine/newgutassignedmarker.rds")
newgutassigned$Identity<-Idents(newgutassigned)
#Subset enteroendocrine cell population from epithelial cells 
newguteec<-subset(newgutassigned, Identity=="EECs")

#Pull expression of enteroendocrine-specific genes (use "data" slot if normalized)
genes <- c("GIP","PYY","CCK","SST","TPH1")
expr <- FetchData(newguteec, vars = genes)

#Set Thresholds: example uses >0 (detected) AND above gene-specific percentile
thr <- sapply(genes, function(g) quantile(expr[[g]], 0.95, na.rm = TRUE))
is_pos <- sweep(expr, 2, thr, `>`)  # TRUE if above 95th percentile

#Assign labels/identity to EECs with priority / ambiguity handling based on enteroendocrine-specific gene expression
label <- rep("Early endocrine progenitors", nrow(expr))
label[is_pos[,"TPH1"]] <- "Enterochromaffin"
label[is_pos[,"SST"]]  <- "D"
label[is_pos[,"CCK"]]  <- "I"
label[is_pos[,"PYY"]]  <- "L"
label[is_pos[,"GIP"]]  <- "K"

#Flag multi-positive and assign multi-positive cells as late endocrine prohenitors 
n_pos <- rowSums(is_pos)
label[n_pos >= 2] <- "Late Endocrine Progenitors"

#Assign eec identity according to gut_endocrine_call
newguteec$gut_endocrine_call <- label
table(newguteec$gut_endocrine_call)
Idents(newguteec)<-newguteec$gut_endocrine_call
newguteec$group1<-paste0(newguteec$gut_endocrine_call,newguteec$tissue)

#Rename cluster labels 
new.cluster.ids <- c("Enterochromaffin cells", "Late Endocrine Progenitors", "Early Endocrine Progenitors", "K cells", "I cells", "D cells",
                     "L cells", "Early endocrine progenitors", "Late Endocrine Progenitors","L cells","K cells","D cells","I cells","Enterochromaffin cells","Early Endocrine Progenitors",
                     "Enterochromaffin cells","L cells","D cells","Late Endocrine Progenitors","Late Endocrine Progenitors","Late Endocrine Progenitors")
names(new.cluster.ids) <- levels(newguteec)
newguteec <- RenameIdents(newguteec, new.cluster.ids)

#Rename cluster labels 
new.cluster.ids <- c("Enterochromaffin cells", "Late Endocrine Progenitors", "Early Endocrine Progenitors", "K cells", "I cells", "D cells",
                     "L cells", "Early Endocrine Progenitors")
names(new.cluster.ids) <- levels(newguteec)
newguteec <- RenameIdents(newguteec, new.cluster.ids)

#Save seurat analysis for EEC assignment
saveRDS(newguteec, "newguteecassigned.rds")


