#Load relevant packages 
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

#Load fetal epithelial ATAC dataset 
epithelial<-readRDS("/athena/chenlab/scratch/jjv4001/gutrevision/Epithelialnew.rds")
DefaultAssay(epithelial)<-"ATAC"
#Idents(epithelial)<-epithelial$predicted.id

#Find differential accessible regions of all epithelial cell types and save differential accessible regions
da_peaks <- FindAllMarkers(
  object = epithelial,
  assay = "ATAC",
  only.pos = TRUE,               # show only enriched peaks per cluster
  min.pct = 0.1,                 # minimum fraction of cells accessible in each group
  test.use = 'LR',               # logistic regression recommended for sparse data
  latent.vars = 'nCount_ATAC'  # optional: account for sequencing depth
)
write.csv(da_peaks, "/athena/chenlab/scratch/jjv4001/gutrevision/olddata/da_peaks.csv")


da_peaks <- FindMarkers(
  epithelial,
  ident.1="Goblet cells",
  ident.2="Other",           # show only enriched peaks per cluster
  min.pct = 0.05,                 # minimum fraction of cells accessible in each group
  test.use = 'LR',               # logistic regression recommended for sparse data
  latent.vars = 'nCount_ATAC'  # optional: account for sequencing depth
)


# Expression matrix (genes x cells)
peak_counts <- GetAssayData(epithelial, assay = "ATAC", slot = "counts")

# Cell metadata
cell_metadata <- epithelial@meta.data

# Feature metadata — must include 'gene_short_name'
peak_names <- rownames(peak_counts)
peak_metadata <- data.frame(gene_short_name = peak_names, row.names = peak_names)

common <- intersect(colnames(peak_counts), rownames(cell_metadata))
peak_counts <- peak_counts[, common]
cell_metadata <- cell_metadata[common, , drop = FALSE]
identical(colnames(peak_counts), rownames(cell_metadata))
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = peak_metadata)

#Create new cds dataset from epithelial atac-seq dataset
cds <- newCellDataSet(
  as(peak_counts, "dgCMatrix"),  # ensure it's a sparse matrix
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

umap_coords <- Embeddings(epithelial, "umap")
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
Links(epithelial) <- links
edb <- EnsDb.Hsapiens.v86
genes <- genes(edb)
#save co-accessible links and cis-co-accessible networks (CCANs)
saveRDS(conns, "/athena/chenlab/scratch/jjv4001/gutrevision/olddata/conns.rds")
saveRDS(ccans, "/athena/chenlab/scratch/jjv4001/gutrevision/olddata/ccans.rds")

# Download the GTF file for human (GRCh38)
temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz", temp)

# Read it using rtracklayer
library(rtracklayer)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)
seqlevelsStyle(gene_anno) <- "UCSC"
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

# Create a Seurat assay from gene activity matrix using Cicero data
fData(cds)$gene <- fData(cds)$gene_symbol
unnorm_ga <- build_gene_activity_matrix(cds, conns)
gene_activity_assay <- CreateAssayObject(counts = unnorm_ga)

# Add to Seurat object
epithelial[["GeneActivity"]] <- gene_activity_assay

# Normalize and scale the gene activity data
epithelial <- NormalizeData(epithelial, assay = "GeneActivity")
epithelial <- ScaleData(epithelial, assay = "GeneActivity")

#save epithelial ATAC data with links
saveRDS(epithelial, "/athena/chenlab/scratch/jjv4001/gutrevision/olddata/epithelialciceroatac.rds")

#Create gene activity matrix using fragment counts intersecting promoter and gene body
gene.activities <- GeneActivity(epithelial)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
epithelial[['RNA']] <- CreateAssayObject(counts = gene.activities)
epithelial <- NormalizeData(
  object = epithelial,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(epithelial$nCount_RNA)
)
#save epithelial ATAC data with links
saveRDS(epithelial, "/athena/chenlab/scratch/jjv4001/gutrevision/olddata/epithelialciceroatac.rds")


#load related packages
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(dplyr)
library(EnsDb.Hsapiens.v86)

# Combine and get unique peak names
all_peaks <- unique(c(conns$Peak1, conns$Peak2))

# Separate into chrom, start, end
peak_tbl <- data.frame(peak = all_peaks) %>%
  mutate(original_peak = peak) %>%  # keep original peak string
  separate(peak, into = c("chr", "start", "end"), sep = "-", convert = TRUE, extra = "merge", fill = "right")

peak_tbl_clean <- peak_tbl %>% 
  filter(!is.na(start), !is.na(end))

gr_peaks <- GRanges(
  seqnames = peak_tbl_clean$chr,
  ranges = IRanges(start = peak_tbl_clean$start, end = peak_tbl_clean$end)
)
names(gr_peaks) <- as.character(peak_tbl_clean$original_peak)
genes <- genes(EnsDb.Hsapiens.v86)
genes_clean <- keepStandardChromosomes(genes, pruning.mode = "coarse")

# Extend to cover promoter region: 2kb upstream, 500bp downstream
promoters_extended <- trim(promoters(genes_clean, upstream = 2000, downstream = 500))
seqlevelsStyle(promoters_extended) <- "UCSC"

# Overlap with promoter-extended gene annotations
hits <- findOverlaps(gr_peaks, promoters_extended)

# Create a peak-to-gene mapping
peak_gene_links <- data.frame(
  peak = names(gr_peaks)[queryHits(hits)],
  gene_id = promoters_extended$gene_id[subjectHits(hits)]
)

# Optional: Add gene symbols
library(org.Hs.eg.db)
peak_gene_links$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = peak_gene_links$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Extract fData as a data.frame
library(Biobase)
fdata_df <- Biobase::fData(cds)
#fdata_df <- as.data.frame(fData(cds))

# If rownames are peak names, add as a column for join
fdata_df$peak <- rownames(fdata_df)

# Join with peak_gene_links to add gene info
library(dplyr)
fdata_df <- fdata_df %>%
  left_join(peak_gene_links %>% select(peak, gene_symbol), by = "peak")

# Add back gene column named 'gene' as Cicero expects
fdata_df$gene <- fdata_df$gene_symbol

# Update the fData of cds
fData(cds) <- fdata_df

# If DA peaks are stored as rownames
da_peak_names <- rownames(da.markers)

# Filter links to DA peaks
linked_da_genes <- subset(peak_gene_links, peak %in% da_peak_names)

# Unique genes
unique_da_genes <- unique(linked_da_genes$gene_symbol)

write.csv(linked_da_genes, "/athena/chenlab/scratch/jjv4001/gutrevision/olddata/DA_Peak_Gene_Links_SI.csv", row.names = FALSE)
