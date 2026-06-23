

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
saveRDS(conns, "/athena/chenlab/scratch/jjv4001/gutrevision/conns.rds")
saveRDS(ccans, "/athena/chenlab/scratch/jjv4001/gutrevision/ccans.rds")
# Extend gene ranges by 2kb upstream to include promoters
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(epithelial) <- links
# Add gene annotations as standardized columns
gene_anno <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene_anno) <- "UCSC"
# Convert to data frame for manipulation
gene_anno_df <- as.data.frame(gene_anno)
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
#### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos_df <- as.data.frame(pos)
pos_df <- pos_df[order(pos_df$start), ]
# remove all but the first exons per transcript
pos_df <- pos_df[!duplicated(pos_df$transcript),] 
# make a 1 base pair marker of the TSS
pos_df$end <- pos_df$start + 1 
neg <- subset(gene_anno, strand == "-")
neg_df <- as.data.frame(neg)
neg_df <- neg_df[order(neg_df$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg_df <- neg_df[!duplicated(neg_df$transcript),] 
neg_df$start <- neg_df$end - 1

gene_annotation_sub <- bind_rows(pos_df, neg_df) %>%
  select(chromosome, start, end, symbol)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

#annotate the cds with gene names
cds <- annotate_cds_by_site(cds, gene_annotation_sub)

# View structure
fData(cds)$gene <- fData(cds)$gene_symbol
unnorm_ga <- build_gene_activity_matrix(cds, conns)
gene_activity_assay <- CreateAssayObject(counts = unnorm_ga)

# Add geneactivity to Seurat object
epithelial[["GeneActivity"]] <- gene_activity_assay

# Normalize and scale the gene activity data
epithelial <- NormalizeData(epithelial, assay = "GeneActivity")
epithelial <- ScaleData(epithelial, assay = "GeneActivity")

#save Seurat analysis
saveRDS(epithelial, "/athena/chenlab/scratch/jjv4001/gutrevision/epithelialciceroatac.rds")

# Assuming conns has Peak1 and Peak2 in "chr-start-end" format

# Convert peak names to GRanges
peaks_to_gr <- function(peaks) {
  gr_list <- lapply(peaks, function(p) {
    parts <- unlist(strsplit(p, "-"))
    GRanges(seqnames=parts[1], ranges=IRanges(as.numeric(parts[2]), as.numeric(parts[3])))
  })
  do.call(c, gr_list)
}

# Make GRanges for peaks in conns
peak1_gr <- peaks_to_gr(conns$Peak1)
peak2_gr <- peaks_to_gr(conns$Peak2)

# Combine peaks GRanges (unique)
all_peaks_gr <- unique(c(peak1_gr, peak2_gr))

genes_clean <- genes[seqnames(genes) != "KI270750.1"]
extended_genes <- trim(promoters(genes_clean, upstream = 2000, downstream = 500))

# Find overlaps between peaks and extended genes
hits <- findOverlaps(all_peaks_gr, extended_genes)

# Create a data.frame mapping peaks to genes
peak_gene_links <- data.frame(
  peak = names(all_peaks_gr)[queryHits(hits)],
  gene_id = extended_genes$gene[subjectHits(hits)]
)

# Map gene IDs (assuming ENSEMBL IDs) to gene symbols using org.Hs.eg.db
peak_gene_links$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = peak_gene_links$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Now get unique gene symbols linked to peaks
cicero_genes <- unique(peak_gene_links$gene_symbol)

# Assuming da_peaks$gene (or your DA peaks gene column) exists, find cicero linked DA genes
intersect_genes <- intersect(da_peaks$gene, cicero_genes)
write.csv(intersect_genes, "/athena/chenlab/scratch/jjv4001/gutrevision/DA_Cicero_Linked_Genes.csv", row.names = FALSE)
          



