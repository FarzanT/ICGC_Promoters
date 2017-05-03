# Useful packages
install.packages(c("readr", "data.table", "igraph"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("DESeq", "limma", "GenomicRanges", "IRanges", "BSGenome",
           "Biostrings", "AnnotationHub", "GenomicFeatures",
           "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19",
           "Biobase", "hgu95av2.db", "airway", "GEOquery", "biomaRt",
           "minfiData", "oligo"))

library(data.table)
library(readr)
library(Biobase)
# Read the expression data
row.names(temp)
eSetTSV <- data.table::fread("joint_fpkm_uq.tsv")
# ENSG IDs are in the first column
(eSetTSV)[1:10, 1]
class(eSetTSV)
# Append row names
row.names(eSetTSV) <- unlist(eSetTSV[, 1])
# Check column names
colnames(eSetTSV) # Correct

# Create ExpressionSet:
# Must convert to a matrix first
?`expression-class`
eSetMatrix <- as.matrix(eSetTSV)


# featureData is a matrix where each row contains information regarding each
# feature



eSet <- ExpressionSet(assayData = eSetMatrix, )
class(eSet)

# Now we can use Biobase functions:
# sampleNames() returns patient IDs in this case
sampleNames(eSet)[1:5]
ncol(eSet)

# featureNames returns the feature IDs, which may need to be translated
featureNames(eSet)[1:5]
nrow(eSet)
max(as.numeric(featureNames(eSet)))


# Load the updated mutation scores files
mut_scores <- readRDS("SummerProject/Summer_2017/Updated_mutation_scores/Chr22.rds")
mut_scores

# Load the lookup table
lookup <- read_tsv("SummerProject/Summer_2017/pcawg_sample_sheet.v1.4.2016-09-14.tsv")
# The elements that are mapped are the columns
lookup[1,]
colnames(lookup)



# TODO: Create featureData, experimentData etc. and attach to the expession set
