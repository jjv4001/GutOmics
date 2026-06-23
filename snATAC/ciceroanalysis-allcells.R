#Load required packages
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(monocle3)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cicero)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)
library(txdbmaker)
#read integrated object with all cell types in ATAC-seq
integrated<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/integratedgutcategory.rds")
#Calculate gene activity and add to RNA assay, normalize and scale
gene.activities <- GeneActivity(integrated)
integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)
DefaultAssay(integrated)<-"RNA"

# Get Expression matrix (genes x cells)
peak_counts <- GetAssayData(integrated, assay = "ATAC", slot = "counts")

# Get cell metadata
cell_metadata <- integrated@meta.data

# Get feature metadata — must include 'gene_short_name'
peak_names <- rownames(peak_counts)
peak_metadata <- data.frame(gene_short_name = peak_names, row.names = peak_names)

common <- intersect(colnames(peak_counts), rownames(cell_metadata))
peak_counts <- peak_counts[, common]
cell_metadata <- cell_metadata[common, , drop = FALSE]
identical(colnames(peak_counts), rownames(cell_metadata))
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = peak_metadata)

#Create new cds dataset from all cell atac-seq dataset
cds <- newCellDataSet(
  as(peak_counts, "dgCMatrix"),  # ensure it's a sparse matrix
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

umap_coords <- Embeddings(integrated, "umap")
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome.df <- data.frame(V1 = names(genome), V2 = as.numeric(genome))
edb <- EnsDb.Hsapiens.v86
genes <- genes(edb)
transcripts <- transcripts(edb)


# Proceed with cicero
cicero_cds <- make_cicero_cds(cds, reduced_coordinates = umap_coords)
conns <- run_cicero(cicero_cds, genomic_coords = genome.df)
#Find cis-co-accessible networks (CCANs)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(integrated) <- links
edb <- EnsDb.Hsapiens.v86
genes <- genes(edb)
#save co-accessible links and cis-co-accessible networks (CCANs)
saveRDS(conns, "/athena/chenlab/scratch/jjv4001/gutrevision/connsallcells.rds")
saveRDS(ccans, "/athena/chenlab/scratch/jjv4001/gutrevision/ccansallcells.rds")

# Download the GTF file for human (GRCh38)
temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz", temp)

# Read it using rtracklayer
library(rtracklayer)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# Add standardized columns
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name


#### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

#annotate the cds with gene names
cds <- annotate_cds_by_site(cds, gene_annotation_sub)


# Create a Seurat assay from gene activity matrix
unnorm_ga <- build_gene_activity_matrix(cds, conns)
gene_activity_assay <- CreateAssayObject(counts = unnorm_ga)

# Add to Seurat object
integrated[["GeneActivity"]] <- gene_activity_assay

# Normalize the gene activity data and scale the data
integrated <- NormalizeData(integrated, assay = "GeneActivity")
integrated <- ScaleData(integrated, assay = "GeneActivity")

#save ATAC-seq data for all cells with cicero links
saveRDS(integrated, "/athena/chenlab/scratch/jjv4001/gutrevision/integratedciceroatacallcells.rds")
