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
library(fastmatch)
library(stringr)

# Read the expression data
eSetTSV <- data.table::fread("joint_fpkm_uq.tsv")

# ==== Assay Data Preparation ====
# ENSG IDs are in the first column
(eSetTSV)[1:10, 1]
class(eSetTSV)
# Append row names
# Must drop decimals in the ENSG IDs
?str_replace
rowNames = stringr::str_replace_all(string = unlist(eSetTSV[, 1]),
                                    pattern = "(ENSG\\d*)\\..*",
                                    replacement = "\\1")
row.names(eSetTSV) <- rowNames
# Check column names
colnames(eSetTSV) # Correct

# Create ExpressionSet:
# Must convert to a matrix first
?ExpressionSet
eSetMatrix <- as.matrix(eSetTSV)
dimnames(eSetMatrix) <- list(row.names(eSetTSV), colnames(eSetTSV))

# ==== Feature Data Preparation ====

# featureData is an AnnotatedDataFram where each row contains information
# regarding each feature, e.g. Gene in this case. Row names must match.
# Create a table where the ENSG IDs are mapped to ENST IDs and HGNC symbols
# Can append other information later on
# Match column names with the lookup table, and retrieve other information 

# Read the (updated) gene lookup table
geneLookup <- data.table::fread("LookupTables/IDs_convert_Table_updated.txt")
# Retrieve Gene IDs
ENSGIDs = row.names(eSetMatrix)
# Be sure to replace NA with 0, so that further subsetting works properly
tableMatchIdx = fastmatch::fmatch(x = ENSGIDs,
                                  table = geneLookup$`Ensembl Gene ID`,
                                  nomatch = 0)
myFeatureData <- Biobase::AnnotatedDataFrame(data = geneLookup[tableMatchIdx])
# Recall that row names must match
Biobase::featureNames(myFeatureData) <- geneLookup[tableMatchIdx]$`Ensembl Gene ID`

# ==== Phenotype Data Preparation ====
# TODO: What to do with unmatchable sample IDs
# Phenotype data contains sample information, e.g. tumor type for each patient
# This data can be retrieved from the PCAWG file
pcawg <- data.table::fread("LookupTables/pcawg_sample_sheet.v1.4.2016-09-14.tsv")

# We can match using the aliquot_id:
# Keep that NA rows since not all aliquot_ids can be matched (currently)
sampleIdx = fastmatch::fmatch(x = colnames(eSetMatrix), table = pcawg$aliquot_id)

# From the man page:
# The number of rows in phenoData must match the number of columns in
# assayData. Row names of phenoData must match column names of the
# matrix / matricies in assayData.
myPhenoData <- Biobase::AnnotatedDataFrame(data = pcawg[sampleIdx])

# Check
nrow(myPhenoData) == ncol(eSetMatrix)
Biobase::featureNames(myPhenoData) <- pcawg[sampleIdx]$aliquot_id
  
  
  
  
  
  
eSet <- ExpressionSet(assayData = eSetMatrix, featureData = myFeatureData, phenoData = )
class(eSet)

# Now we can test Biobase functions:
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
