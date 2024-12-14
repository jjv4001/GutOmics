# GutOmics

## scRNA

seurat33209721.R, seurat33406409.R, seurat34497389.R, seurat33278341.R, seurat31753849.R, and seurat35176508.R provide codes for initial processing and Seurat analysis of scRNA data from fetal and adult human intestinal samples. 
Cell types from different germ layers are assigned and epithelial cells which are of interest are subset. 
dataintegration.R integrates epithelial cells from multiple scRNA datasets and assigns different epithelial cell type identities. 
EEC.R is used to subset enteroendocrine cells from the epithelial cell population and assign enteroendocrine subtypes. 

##snATAC

ATAC.R generates a unified set of peaks across all fetal intestinal samples from different regions including the foregut (duodenum), midgut (jejunum), hindgut (ileum), and colon using the reduce function, 
creates Seurat objects with chromatin assays, assigns cells originating from different germ layers, and identifies different epithelial cell types based on open chromatin accessibility peaks of characteristic 
genetic markers.

##Xenium

Xenium.R assigns epithelial cell types in intestinal villi and crypt stuctures from different regions of interest of the fetal intestine including the duodenum and colon from 17 and 20 week donors using Seurat. 
Moreover, enteroendocrine cells are also subset and enteroendocrine cell subtypes are classified.

