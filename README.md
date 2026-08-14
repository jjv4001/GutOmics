# GutOmics

## Bulk RNA-seq
Kcellweighted.R provides codes for calculating weighted gene module score for relevant genes of interest from bulk RNA-seq data for hPSC-derived intestinal organoids treated with different small molecules.

## GWAS 
run_cell_type_individual.sh, step1.sh, step2.sh, summary.sh provide codes for GWAS analysis of the epithelial clusters detected in the snATAC-seq data to datasets from previous literature to understand their genetic correlation to metabolic disorders. plot_ldsc_results.ipynb and GWAS_code_plot.R provide codes for plotting the GWAS results. 

## scRNA

seurat33209721.R, seurat33406409.R, seurat34497389.R, seurat33278341.R, seurat31753849.R, and seurat35176508.R provide codes for initial processing and Seurat analysis of scRNA data from fetal and adult human intestinal samples. 
Cell types from different germ layers are assigned and epithelial cells which are of interest are subset. 
dataintegration.R integrates epithelial cells from multiple scRNA datasets and assigns different epithelial cell type identities. 
EEC.R is used to subset enteroendocrine cells from the epithelial cell population and assign enteroendocrine subtypes. fetal + adult epithelial cell analysis.R provides the updated code for subsetting samples corressponding to fetal developmental timepoints of interest and reanalysis as well as separate analysis of the scRNA-seq data corresponding to adult samples. fetal epithelial + eecs cell analysis.R provides the updated code for the assignment of epithelial and enteroendocrine cell types in the fetal samples corresponding to developmental timepoints of interest (11-17 weeks PCW).  

## snATAC

ATAC.R generates a unified set of peaks across all fetal intestinal samples from different regions including the foregut (duodenum), midgut (jejunum), hindgut (ileum), and colon using the reduce function, 
creates Seurat objects with chromatin assays, assigns cells originating from different germ layers, and identifies different epithelial cell types based on open chromatin accessibility peaks of characteristic 
genetic markers. ciceroanalysis-allcells.R, cicero and motif enrichment analysis.R, cicero and motif enrichment analysis-epithelial cells.R and cicero analysis-epithelial cells-v1.R provide codes for the cicero and motif enrichment analysis of all cell types as well as epithelial cell types. obtain_bed_peaks.R provides codes for obtaining the bed file for peaks from snATAC-seq data.

## Xenium

Xenium.R assigns epithelial cell types in intestinal villi and crypt stuctures from different regions of interest of the fetal intestine including the duodenum and colon from 17 and 20 week donors using Seurat. 
Moreover, enteroendocrine cells are also subset and enteroendocrine cell subtypes are classified. Xenium_all_cell_types_assignment_subsetepithelialcells.R provides the updated code for assigning all cell types in the new dataset corresponding to duodenum, jejunum, and colon samples from 18 and 20 PCW donors. Epithelial_cell_assignment.R specifically assigns the epithelial cell types detected in the duodenum, jejunum, and colon samples from 18 and 20 PCW donors. EEC assignment.R provides codes for the assignment of enteroendocrine cells in the duodenum, jejunum, and colon samples from 18 and 20 PCW donors. Immune_Mesenchymal_assignment.R specifically assigns cell types in the different tissue sections based on marker expression and assignment of cell types in Harnik Y. et al., Nature. (2024);632(8027):1101-1109. Gene_Expression_plot_to_FibroblastBM.R, Gene_Expression_plot_to_FibroblastBM_v1.R, and Violin_plot_Distance_to_FibroblastBM.R provide code for generating plots for the visualization of gene expression with increasing distance from the mesenchyme or Fibroblast BM. Spatial plots across 3 tissues.R generates spatial plots of all of the 3 different tissues: the duodenum, the jejunum, and the colon and the expression of different genes of interest across these tissue sections. Calculate_K cell_distance_to_FibroblastBM.R and Calculate_AFP_APOA4_APOE_expressing cells_distance_to_FibroblastBM.R provide codes for the calculation of the distance of the K cells to the fibroblast BM cells/mesenchyme and the distance of the top quartile of cells expressing AFP/APOA4/APOE to the fibroblast BM cells/mesenchyme.

## CellChat analysis

CellChat_analysis.R performs cell–cell communication analysis across broad fetal intestinal cell compartments, including epithelial, mesenchymal, immune, neuronal, red blood cell, and endothelial populations. The script extracts normalized RNA expression data from the Seurat object, and uses the human CellChat ligand–receptor interaction database (CellChatDB.human) to identify overexpressed signaling genes and ligand–receptor interactions. It calculates communication probabilities between cell populations, aggregates ligand–receptor interactions at the signaling pathway level, and computes overall interaction numbers and communication strengths between cell types. The script further evaluates incoming and outgoing signaling roles and network centrality, ranks signaling pathways and ligand–receptor pairs based on communication strength

