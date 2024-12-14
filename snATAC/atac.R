library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
plan("multisession", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

elevenwkcolon <- read.table(file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_colon/filtered_peak_bc_matrix/peaks.bed",
                            col.names = c("chr", "start", "end"))

elevenwkforegut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_foregut/filtered_peak_bc_matrix/peaks.bed",
                               col.names = c("chr", "start", "end")
)

elevenwkmidgut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_midgut/filtered_peak_bc_matrix/peaks.bed",
                              col.names = c("chr", "start", "end")
)

eighteenwkcolon <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_colon/filtered_peak_bc_matrix/peaks.bed",
                               col.names = c("chr", "start", "end")
)
eighteenwkhindgut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_hindgut/filtered_peak_bc_matrix/peaks.bed",
                                 col.names = c("chr", "start", "end")
)
eighteenwkmiddgut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_midgut/filtered_peak_bc_matrix/peaks.bed",
                                 col.names = c("chr", "start", "end")
)
fifteenwkforegut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/F1_3834_foregut/filtered_peak_bc_matrix/peaks.bed",
                                col.names = c("chr", "start", "end")
)
fifteenwkmidgut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/M1_3834_midgut/filtered_peak_bc_matrix/peaks.bed",
                               col.names = c("chr", "start", "end")
)
fifteenwkhindgut <- read.table( file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3834_hindgut/filtered_peak_bc_matrix/peaks.bed",
                                col.names = c("chr", "start", "end")
)

gr.eighteenwkcolon <- makeGRangesFromDataFrame(eighteenwkcolon)
gr.eighteenwkhindgut <- makeGRangesFromDataFrame(eighteenwkhindgut)
gr.eighteenwkmidgut <- makeGRangesFromDataFrame(eighteenwkmiddgut)
gr.elevenwkcolon <- makeGRangesFromDataFrame(elevenwkcolon)
gr.elevenwkforegut <- makeGRangesFromDataFrame(elevenwkforegut)
gr.elevenwkmidgut <- makeGRangesFromDataFrame(elevenwkmidgut)
gr.fifteenwkforegut <- makeGRangesFromDataFrame(fifteenwkforegut)
gr.fifteenwkhindgut <- makeGRangesFromDataFrame(fifteenwkhindgut)
gr.fifteenwkmidgut <- makeGRangesFromDataFrame(fifteenwkmidgut)

combined.peaks <- reduce(x = c(gr.eighteenwkcolon, gr.eighteenwkhindgut, gr.eighteenwkmidgut, gr.elevenwkcolon, gr.elevenwkforegut, gr.elevenwkmidgut, gr.fifteenwkforegut, gr.fifteenwkhindgut, gr.fifteenwkmidgut))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

md.eighteenwkcolon <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_colon/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.eighteenwkhindgut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_hindgut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.eighteenwkmidgut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_midgut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.elevenwkcolon <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_colon/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.elevenwkforegut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_foregut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.elevenwkmidgut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_midgut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.fifteenwkforegut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/F1_3834_foregut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.fifteenwkhindgut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3834_hindgut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.fifteenwkmidgut <- read.table(
  file = "/athena/chenlab/scratch/jjv4001/gut/source-selected/M1_3834_midgut/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 


md.eighteenwkcolon <- md.eighteenwkcolon[md.eighteenwkcolon$passed_filters > 500, ]
md.eighteenwkhindgut <- md.eighteenwkhindgut[md.eighteenwkhindgut$passed_filters > 500, ]
md.eighteenwkmidgut <- md.eighteenwkmidgut[md.eighteenwkmidgut$passed_filters > 500, ]
md.elevenwkcolon <- md.elevenwkcolon[md.elevenwkcolon$passed_filters > 500, ]
md.elevenwkforegut <- md.elevenwkforegut[md.elevenwkforegut$passed_filters > 500, ]
md.elevenwkmidgut <- md.elevenwkmidgut[md.elevenwkmidgut$passed_filters > 500, ]
md.fifteenwkhindgut <- md.fifteenwkhindgut[md.fifteenwkhindgut$passed_filters > 500, ]
md.fifteenwkforegut <- md.fifteenwkforegut[md.fifteenwkforegut$passed_filters > 500, ]
md.fifteenwkmidgut <- md.fifteenwkmidgut[md.fifteenwkmidgut$passed_filters > 500, ]

frags.eighteenwkcolon <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_colon/fragments.tsv.gz",
  cells = rownames(md.eighteenwkcolon)
)


frags.eighteenwkhindgut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_hindgut/fragments.tsv.gz",
  cells = rownames(md.eighteenwkhindgut)
)


frags.eighteenwkmidgut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3824_midgut/fragments.tsv.gz",
  cells = rownames(md.eighteenwkmidgut)
)



frags.elevenwkcolon <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_colon/fragments.tsv.gz",
  cells = rownames(md.elevenwkcolon)
)


frags.elevenwkforegut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_foregut/fragments.tsv.gz",
  cells = rownames(md.elevenwkforegut)
)


frags.elevenwkmidgut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3767_midgut/fragments.tsv.gz",
  cells = rownames(md.elevenwkmidgut)
)

frags.fifteenwkforegut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/F1_3834_foregut/fragments.tsv.gz",
  cells = rownames(md.fifteenwkforegut)
)


frags.fifteenwkhindgut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/3834_hindgut/fragments.tsv.gz",
  cells = rownames(md.fifteenwkhindgut)
)


frags.fifteenwkmidgut <- CreateFragmentObject(
  path = "/athena/chenlab/scratch/jjv4001/gut/source-selected/M1_3834_midgut/fragments.tsv.gz",
  cells = rownames(md.fifteenwkmidgut)
)

colon3824.counts <- FeatureMatrix(
  fragments = frags.eighteenwkcolon,
  features = combined.peaks,
  cells = rownames(md.eighteenwkcolon)
)

hindgut3824.counts <- FeatureMatrix(
  fragments = frags.eighteenwkhindgut,
  features = combined.peaks,
  cells = rownames(md.eighteenwkhindgut)
)

midgut3824.counts <- FeatureMatrix(
  fragments = frags.eighteenwkmidgut,
  features = combined.peaks,
  cells = rownames(md.eighteenwkmidgut)
)

midgut3834.counts <- FeatureMatrix(
  fragments = frags.fifteenwkmidgut,
  features = combined.peaks,
  cells = rownames(md.fifteenwkmidgut)
)

hindgut3834.counts <- FeatureMatrix(
  fragments = frags.fifteenwkhindgut,
  features = combined.peaks,
  cells = rownames(md.fifteenwkhindgut)
)

foregut3834.counts <- FeatureMatrix(
  fragments = frags.fifteenwkforegut,
  features = combined.peaks,
  cells = rownames(md.fifteenwkforegut)
)

midgut3767.counts <- FeatureMatrix(
  fragments = frags.elevenwkmidgut,
  features = combined.peaks,
  cells = rownames(md.elevenwkmidgut)
)

colon3767.counts <- FeatureMatrix(
  fragments = frags.elevenwkcolon,
  features = combined.peaks,
  cells = rownames(md.elevenwkcolon)
)

foregut3767.counts <- FeatureMatrix(
  fragments = frags.elevenwkforegut,
  features = combined.peaks,
  cells = rownames(md.elevenwkforegut)
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

colon3824_assay <- CreateChromatinAssay(colon3824.counts, fragments = frags.eighteenwkcolon)
colon3824 <- CreateSeuratObject(colon3824_assay, assay = "ATAC")
midgut3824_assay <- CreateChromatinAssay(midgut3824.counts, fragments = frags.eighteenwkmidgut)
midgut3824 <- CreateSeuratObject(midgut3824_assay, assay = "ATAC")
hindgut3824_assay <- CreateChromatinAssay(hindgut3824.counts, fragments = frags.eighteenwkhindgut)
hindgut3824 <- CreateSeuratObject(hindgut3824_assay, assay = "ATAC")
colon3767_assay <- CreateChromatinAssay(colon3767.counts, fragments = frags.elevenwkcolon)
colon3767 <- CreateSeuratObject(colon3767_assay, assay = "ATAC")
midgut3767_assay <- CreateChromatinAssay(midgut3767.counts, fragments = frags.elevenwkmidgut)
midgut3767 <- CreateSeuratObject(midgut3767_assay, assay = "ATAC")
foregut3767_assay <- CreateChromatinAssay(foregut3767.counts, fragments = frags.elevenwkforegut)
foregut3767 <- CreateSeuratObject(foregut3767_assay, assay = "ATAC")
foregut3834_assay <- CreateChromatinAssay(foregut3834.counts, fragments = frags.fifteenwkforegut)
foregut3834 <- CreateSeuratObject(foregut3834_assay, assay = "ATAC")
midgut3834_assay <- CreateChromatinAssay(midgut3834.counts, fragments = frags.fifteenwkmidgut)
midgut3834 <- CreateSeuratObject(midgut3834_assay, assay = "ATAC")
hindgut3834_assay <- CreateChromatinAssay(hindgut3834.counts, fragments = frags.fifteenwkhindgut)
hindgut3834 <- CreateSeuratObject(hindgut3834_assay, assay = "ATAC")

genome(annotations) <- "hg38"

Annotation(colon3824) <- annotations

colon3824 <- NucleosomeSignal(object = colon3824)
colon3824 <- TSSEnrichment(object = colon3824, fast = FALSE)
colon3824$peak_region_fragments<-md.eighteenwkcolon$peak_region_fragments
colon3824$passed_filters<-md.eighteenwkcolon$passed_filters
colon3824$blacklist_region_fragments<-md.eighteenwkcolon$blacklist_region_fragments
colon3824$pct_reads_in_peaks <- colon3824$peak_region_fragments / colon3824$passed_filters * 100
colon3824$blacklist_ratio <- colon3824$blacklist_region_fragments / colon3824$peak_region_fragments
VlnPlot(
  object = colon3824,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

colon3824 <- subset(
  x = colon3824,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 10000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(colon3824,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/colon3824.rds")

Annotation(colon3767) <- annotations

colon3767 <- NucleosomeSignal(object = colon3767)
colon3767 <- TSSEnrichment(object = colon3767, fast = FALSE)
colon3767$peak_region_fragments<-md.elevenwkcolon$peak_region_fragments
colon3767$passed_filters<-md.elevenwkcolon$passed_filters
colon3767$blacklist_region_fragments<-md.elevenwkcolon$blacklist_region_fragments
colon3767$pct_reads_in_peaks <- colon3767$peak_region_fragments / colon3767$passed_filters * 100
colon3767$blacklist_ratio <- colon3767$blacklist_region_fragments / colon3767$peak_region_fragments
VlnPlot(
  object = colon3767,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

colon3767 <- subset(
  x = colon3767,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 15000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(colon3767,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/colon3767.rds")

Annotation(midgut3767) <- annotations
midgut3767 <- NucleosomeSignal(object = midgut3767)
midgut3767 <- TSSEnrichment(object = midgut3767, fast = FALSE)
midgut3767$peak_region_fragments<-md.elevenwkmidgut$peak_region_fragments
midgut3767$passed_filters<-md.elevenwkmidgut$passed_filters
midgut3767$blacklist_region_fragments<-md.elevenwkmidgut$blacklist_region_fragments
midgut3767$pct_reads_in_peaks <- midgut3767$peak_region_fragments / midgut3767$passed_filters * 100
midgut3767$blacklist_ratio <- midgut3767$blacklist_region_fragments / midgut3767$peak_region_fragments
VlnPlot(
  object = midgut3767,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

midgut3767 <- subset(
  x = midgut3767,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 15000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(midgut3767,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/midgut3767.rds")

Annotation(foregut3767) <- annotations
foregut3767 <- NucleosomeSignal(object = foregut3767)
foregut3767 <- TSSEnrichment(object = foregut3767, fast = FALSE)
foregut3767$peak_region_fragments<-md.elevenwkforegut$peak_region_fragments
foregut3767$passed_filters<-md.elevenwkforegut$passed_filters
foregut3767$blacklist_region_fragments<-md.elevenwkforegut$blacklist_region_fragments
foregut3767$pct_reads_in_peaks <- foregut3767$peak_region_fragments / foregut3767$passed_filters * 100
foregut3767$blacklist_ratio <- foregut3767$blacklist_region_fragments / foregut3767$peak_region_fragments
VlnPlot(
  object = foregut3767,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

foregut3767 <- subset(
  x = foregut3767,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 13000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(foregut3767,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/foregut3767.rds")

Annotation(midgut3824) <- annotations
midgut3824 <- NucleosomeSignal(object = midgut3824)
midgut3824 <- TSSEnrichment(object = midgut3824, fast = FALSE)
midgut3824$peak_region_fragments<-md.eighteenwkmidgut$peak_region_fragments
midgut3824$passed_filters<-md.eighteenwkmidgut$passed_filters
midgut3824$blacklist_region_fragments<-md.eighteenwkmidgut$blacklist_region_fragments
midgut3824$pct_reads_in_peaks <- midgut3824$peak_region_fragments / midgut3824$passed_filters * 100
midgut3824$blacklist_ratio <- midgut3824$blacklist_region_fragments / midgut3824$peak_region_fragments
VlnPlot(
  object = midgut3824,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

midgut3824 <- subset(
  x = midgut3824,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 15000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(midgut3824,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/midgut3824.rds")

Annotation(hindgut3824) <- annotations
hindgut3824 <- NucleosomeSignal(object = hindgut3824)
hindgut3824 <- TSSEnrichment(object = hindgut3824, fast = FALSE)
hindgut3824$peak_region_fragments<-md.eighteenwkhindgut$peak_region_fragments
hindgut3824$passed_filters<-md.eighteenwkhindgut$passed_filters
hindgut3824$blacklist_region_fragments<-md.eighteenwkhindgut$blacklist_region_fragments
hindgut3824$pct_reads_in_peaks <- hindgut3824$peak_region_fragments / hindgut3824$passed_filters * 100
hindgut3824$blacklist_ratio <- hindgut3824$blacklist_region_fragments / hindgut3824$peak_region_fragments
VlnPlot(
  object = hindgut3824,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

hindgut3824 <- subset(
  x = hindgut3824,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 5000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1
)

saveRDS(hindgut3824,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/hindgut3824.rds")

Annotation(hindgut3834) <- annotations
hindgut3834 <- NucleosomeSignal(object = hindgut3834)
hindgut3834 <- TSSEnrichment(object = hindgut3834, fast = FALSE)
hindgut3834$peak_region_fragments<-md.fifteenwkhindgut$peak_region_fragments
hindgut3834$passed_filters<-md.fifteenwkhindgut$passed_filters
hindgut3834$blacklist_region_fragments<-md.fifteenwkhindgut$blacklist_region_fragments
hindgut3834$pct_reads_in_peaks <- hindgut3834$peak_region_fragments / hindgut3834$passed_filters * 100
hindgut3834$blacklist_ratio <- hindgut3834$blacklist_region_fragments / hindgut3834$peak_region_fragments
VlnPlot(
  object = hindgut3834,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

hindgut3834 <- subset(
  x = hindgut3834,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 12500 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(hindgut3834,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/hindgut3834.rds")

Annotation(midgut3834) <- annotations
midgut3834 <- NucleosomeSignal(object = midgut3834)
midgut3834 <- TSSEnrichment(object = midgut3834, fast = FALSE)
midgut3834$peak_region_fragments<-md.fifteenwkmidgut$peak_region_fragments
midgut3834$passed_filters<-md.fifteenwkmidgut$passed_filters
midgut3834$blacklist_region_fragments<-md.fifteenwkmidgut$blacklist_region_fragments
midgut3834$pct_reads_in_peaks <- midgut3834$peak_region_fragments / midgut3834$passed_filters * 100
midgut3834$blacklist_ratio <- midgut3834$blacklist_region_fragments / midgut3834$peak_region_fragments
VlnPlot(
  object = midgut3834,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

midgut3834 <- subset(
  x = midgut3834,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 10000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(midgut3834,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/midgut3834.rds")

Annotation(foregut3834) <- annotations
foregut3834 <- NucleosomeSignal(object = foregut3834)
foregut3834 <- TSSEnrichment(object = foregut3834, fast = FALSE)
foregut3834$peak_region_fragments<-md.fifteenwkforegut$peak_region_fragments
foregut3834$passed_filters<-md.fifteenwkforegut$passed_filters
foregut3834$blacklist_region_fragments<-md.fifteenwkforegut$blacklist_region_fragments
foregut3834$pct_reads_in_peaks <- foregut3834$peak_region_fragments / foregut3834$passed_filters * 100
foregut3834$blacklist_ratio <- foregut3834$blacklist_region_fragments / foregut3834$peak_region_fragments
VlnPlot(
  object = foregut3834,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

foregut3834 <- subset(
  x = foregut3834,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 7500 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

saveRDS(foregut3834,"/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/foregut3834.rds")

library(Seurat)
library(Signac)

foregut3767<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/foregut3767.rds")
midgut3767<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/midgut3767.rds")
colon3767<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/colon3767.rds")
midgut3834<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/midgut3834.rds")
hindgut3834<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/hindgut3834.rds")
foregut3834<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/foregut3834.rds")
midgut3824<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/midgut3824.rds")
hindgut3824<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/hindgut3824.rds")
colon3824<-readRDS("/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/colon3824.rds")


colon3824$dataset <- 'colon3824'
hindgut3824$dataset <- 'hindgut3824'
midgut3824$dataset <- 'midgut3824'
colon3767$dataset <- 'colon3767'
midgut3767$dataset <- 'midgut3767'
foregut3767$dataset <- 'foregut3767'
hindgut3834$dataset <- 'hindgut3834'
foregut3834$dataset <- 'foregut3834'
midgut3834$dataset <- 'midgut3834'


combined<-merge(
  x=midgut3834, 
  y=list(hindgut3824, midgut3824, colon3767, midgut3767, foregut3767, hindgut3834, foregut3834, colon3824),
  add.cell.ids=c("colon3824","hindgut3824","midgut3824","colon3767","midgut3767","foregut3767","hindgut3834","foregut3834","midgut3834")
  )

combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)

library(harmony)
hm.integrated <- RunHarmony(object = combined, group.by.vars = 'dataset', reduction = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1, raster=FALSE)
hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'lsi', dims = 2:30)
hm.integrated <- FindClusters(object = hm.integrated, verbose = FALSE, algorithm = 3, resolution=0.5)

new.cluster.ids <- c("Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial",
                     "Epithelial", "Epithelial", "8","Epithelial","Epithelial","11","Epithelial","13","14","15","Epithelial","17","18","19")
names(new.cluster.ids) <- levels(hm.integrated)
hm.integrated <- RenameIdents(hm.integrated, new.cluster.ids)

new.cluster.ids <- c("Epithelial", "Mesenchymal","11","Mesenchymal","Mesenchymal","15","17","18","19")
names(new.cluster.ids) <- levels(hm.integrated)
hm.integrated <- RenameIdents(hm.integrated, new.cluster.ids)

new.cluster.ids <- c("Epithelial", "Mesenchymal","Neuronal","15","Neuronal","18","19")
names(new.cluster.ids) <- levels(hm.integrated)
hm.integrated <- RenameIdents(hm.integrated, new.cluster.ids)

new.cluster.ids <- c("Epithelial", "Mesenchymal","Neuronal","15","Endothelial","Immune cells")
names(new.cluster.ids) <- levels(hm.integrated)
hm.integrated <- RenameIdents(hm.integrated, new.cluster.ids)

new.cluster.ids <- c("Epithelial", "Mesenchymal","Neuronal","RBCs","Endothelial","Immune cells")
names(new.cluster.ids) <- levels(hm.integrated)
hm.integrated <- RenameIdents(hm.integrated, new.cluster.ids)

saveRDS(hm.integrated, "/athena/chenlab/scratch/jjv4001/gut/newgutanalysis/hmintegratedcategory.rds")

Week11<- FindNeighbors(object = Week11, reduction = 'lsi', dims = 2:30)
Week11<- FindClusters(object = Week11, verbose = FALSE, algorithm = 3, resolution=1)

Epithelial<-subset(gut, idents='Epithelial')
Epithelial <- FindTopFeatures(Epithelial, min.cutoff = 'q0')
Epithelial <- RunSVD(Epithelial)
Epithelial<- RunUMAP(object = Epithelial, reduction = 'lsi', dims = 2:30)
Epithelial<- FindNeighbors(object = Epithelial, reduction = 'lsi', dims = 2:30)
Epithelial<- FindClusters(object = Epithelial, verbose = FALSE, algorithm = 3, resolution=7)
DimPlot(Epithelial, reduction="umap", label=TRUE)

new.cluster.ids <- c("0", "1", "Enterocytes", "Enterocytes", "Enterocytes", "Enterocytes","6", "7", "8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6", "Enterocytes", "Enterocytes","Enterocytes","10","Enterocytes","Enterocytes","Enterocytes","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","10","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","EECs","32","33","34","35","36","37","38","EECs","40","41","42","43","44","45","Enterocytes","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","Enterocytes","63","64","65","66")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","10","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","EECs","32","33","Enterocytes","35","36","37","38","Enterocytes","41","42","43","44","45","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","63","64","65","66")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","10","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","EECs","32","33","35","36","37","38","Enterocytes","Enterocytes","43","44","45","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","63","64","65","66")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Enterocytes","Stem cells","Enterocytes","16","17","18","Enterocytes","20","21","22","23","24","25","26","27","28","29","30","EECs","32","33","35","36","37","38","43","44","45","47","48","49","50","51","52","53","54","55","56","57","58","59","Enterocytes","61","63","64","65","Enterocytes")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","17","18","20","21","Enterocytes","23","24","25","26","27","28","29","30","EECs","32","33","35","36","37","38","43","44","45","47","48","49","50","51","52","53","54","55","56","57","58","59","61","63","Enterocytes","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","17","Enterocytes","20","21","23","24","25","26","27","28","29","30","EECs","Goblet cells","33","35","36","37","38","43","44","45","47","Enterocytes","Goblet cells","50","51","Goblet cells","53","54","55","56","57","58","59","61","63","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","17","20","21","23","24","25","26","Enterocytes","28","29","Enterocytes","EECs","Goblet cells","33","35","36","37","38","43","44","45","47","Enterocytes","51","53","54","55","56","57","58","59","61","63","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","Enterocytes","20","21","23","24","25","26","28","29","EECs","Goblet cells","Enterocytes","35","36","37","38","43","44","45","47","51","53","Enterocytes","55","56","57","58","59","61","63","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","20","Enterocytes","23","24","25","26","28","29","EECs","Goblet cells","35","36","37","38","43","44","45","47","51","53","55","56","57","58","Enterocytes","61","63","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","20","Stem cells","Stem cells","25","26","28","29","EECs","Goblet cells","35","Stem cells","37","38","43","44","45","47","51","53","55","56","Stem cells","58","61","63","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","20","25","26","28","29","EECs","Goblet cells","35","37","38","43","44","45","47","51","53","55","56","Enterocytes","61","63","65")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("0", "1","Enterocytes","6","Stem cells","16","20","25","26","28","29","EECs","Goblet cells","35","37","38","43","44","45","47","51","53","55","56","Enterocytes","Enterocytes","Enterocytes")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("Enterocytes", "Enterocytes","Enterocytes","6","Stem cells","16","20","25","26","28","29","EECs","Goblet cells","35","37","38","43","44","45","47","51","53","55","56")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)

new.cluster.ids <- c("Enterocytes","Enterocytes","Stem cells","Enterocytes","Enterocytes","Enterocytes","Enterocytes","Stem cells","Stem cells","EECs","Goblet cells","Stem cells","Enterocytes","Enterocytes","Stem cells","Enterocytes","Enterocytes","Enterocytes","Stem cells","Enterocytes","Enterocytes","Enterocytes")
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)



