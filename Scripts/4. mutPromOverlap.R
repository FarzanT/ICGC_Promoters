# mutPromOverlap.R

library(BiocInstaller)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome)
library(stringi)
library(data.table)

eSet <- readRDS("Data/ExpressionSet.rds")

# ==== Extract promoter info ====
# Load the promoters
load("Promoters-Encode.rsav")
# Extract the important columns
proms <- promoters[, c(1,2,3,4,6)]
# Split the 'name' column
temp = str_replace_all(unlist(proms[,4]), ".*::.*::(.*)::(.*)\\..*",
                       replacement = "\\1::\\2")
tempGene = str_split(temp, "::", simplify = TRUE)[,1]
tempENSG = str_split(temp, "::", simplify = TRUE)[,2]
proms[,"6"] <- tempGene
proms[,"7"] <- tempENSG
# Rearrange columns
proms <- proms[, c(1,2,3,6,7,5)]
# Assign colmun names
colnames(proms) <- c("chromosome", "start", "end", "hgnc_symbol",
                     "ensembl_gene_id", "strand")
# Prepend 'chr' to chromosome names
proms$chromosome <- paste0("chr", proms$chromosome)
proms
# Save 
saveRDS(proms, "Data/promoters.rds")
# TODO: create a lookup function to query
# expression data based on features (genes) and samples (donors)

# ==== Create a GenomicRanges object for these promoters ====
GRproms <- GenomicRanges::makeGRangesFromDataFrame(
    df = proms,
    keep.extra.columns = TRUE,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end",
    strand.field = "strand")
# Save
saveRDS(GRproms, "Data/GRproms.rds")

# ==== Find overlaps of mutations in promoters ====
# Load the updated mutation scores files (test)
muts.chr22 <- readRDS("Data/Updated_mutation_scores/Chr22.rds")
muts.chr22
# Create a GRanges object
gr.chr22 <- makeGRangesFromDataFrame(
    df = muts.chr22,
    keep.extra.columns = TRUE,
    seqnames.field = "chromosome",
    start.field = "position",
    end.field = "position",
    ignore.strand = TRUE)

ov = subsetByOverlaps(gr.chr22, GRproms)



# Load the hg19 genome
BSgenome::available.genomes()
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
# Load annotations for hg19
library(AnnotationHub)
AH <- AnnotationHub()
# UCSC
q <- query(AH, c("homo", "genes"))
q
# Ensembl
q2 <- query(AH, c("homo", "ensg"))
q2
ensemblGenes <- q2[[1]]
seqlevels(ensemblGenes)
ensg <- keepStandardChromosomes(ensemblGenes, pruning.mode = "coarse")
ensg
# Create a view
Views(Hsapiens$chr22, 19278657, width = 1)
Views(Hsapiens$chr22, 19275639, width = 1)
Views(Hsapiens$chr22, 40290257, width = 1)
# Confirmed, mid's are base conversions on the hg19

# ==== Subset transcripts that have mutations in their promoters ====
# As indicated by the mid in the mut_score files:

chr22 = Views(Hsapiens$chr22, 1, width = 3 * 10 ^ 8)
GenomicFeatures::promoters(ensg)
