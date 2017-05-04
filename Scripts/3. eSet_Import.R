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
                                    pattern = "(ENSG.*)\\..*",
                                    replacement = "\\1")
row.names(eSetTSV) <- rowNames
# Check column names
colnames(eSetTSV) # Correct

# Must convert assayData to a matrix
eSetMatrix <- as.matrix(eSetTSV)
# Discard the "feature" column
eSetMatrix <- eSetMatrix[,-1]
ncol(eSetMatrix)
# Discard the name as well
dimnames(eSetMatrix) <- list(row.names(eSetTSV), colnames(eSetTSV)[-1])

### ====  Update ID conversion table ====
# Read the (updated) gene lookup table
geneLookup <- data.table::fread("LookupTables/IDs_convert_Table_updated.txt")
# Retrieve Gene IDs
ENSGIDs = row.names(eSetMatrix)

# See which ENSG IDs are not matched
# Be sure to replace NA with 0, so that further subsetting works properly
tableMatchIdx = fastmatch::fmatch(x = ENSGIDs,
                                  table = geneLookup$`Ensembl Gene ID`,
                                  nomatch = 0)
unmatched = ENSGIDs[which(tableMatchIdx == 0)]
# These must be added to the GeneIDTranslation file
# Check random ENSG to ensure fmatch worked correctly
any(geneLookup$`Ensembl Gene ID` == "ENSG00000221062")
head(unmatched)
# Query biomaRt
library(biomaRt)
mart = useMart("ensembl", "hsapiens_gene_ensembl")
query <- getBM(attributes =c("ensembl_gene_id", "ensembl_transcript_id",
                       "hgnc_symbol"),
filters = "ensembl_gene_id", values = unmatched, mart = mart)
# A quick online search indicates that these ENSG IDs are now deprecated, and 
# thus an archived version of the ENSEMBL database has to be queried
# Use version 54
mart = useMart(host = 'may2009.archive.ensembl.org', 
              biomart ='ENSEMBL_MART_ENSEMBL', 
              dataset ='hsapiens_gene_ensembl')
query <- getBM(attributes =c("ensembl_gene_id", "ensembl_transcript_id",
                       "ensembl_peptide_id", "hgnc_symbol"),
         filters = "ensembl_gene_id", values = unmatched, mart = mart)

# Number of ENSG IDs left
length(unmatched) - nrow(query)
head(query)
class(query)

colnames(query) <- colnames(geneLookup)
# Check
before = nrow(geneLookup)
geneLookup <- rbindlist(list(geneLookup, query))
after = nrow(geneLookup)
after == before + nrow(query)

# Again see which ENSG IDs were not matched
tableMatchIdx = fastmatch::fmatch(x = ENSGIDs,
                                  table = geneLookup$`Ensembl Gene ID`,
                                  nomatch = 0)
unmatched = ENSGIDs[which(tableMatchIdx == 0)]
length(unmatched)
head(unmatched)
# Use version 67 from May 2012
mart = useMart(host = 'may2012.archive.ensembl.org', 
               biomart = 'ENSEMBL_MART_ENSEMBL', 
               dataset = 'hsapiens_gene_ensembl')
query <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                             "ensembl_peptide_id", "hgnc_symbol"),
               filters = "ensembl_gene_id", values = unmatched, mart = mart)
length(unmatched) - nrow(query)
# Update
colnames(query) <- colnames(geneLookup)
# Check
before = nrow(geneLookup)
geneLookup <- rbindlist(list(geneLookup, query))
after = nrow(geneLookup)
after == before + nrow(query)

# And again
tableMatchIdx = fastmatch::fmatch(x = ENSGIDs,
                                  table = geneLookup$`Ensembl Gene ID`,
                                  nomatch = 0)
unmatched = ENSGIDs[which(tableMatchIdx == 0)]
length(unmatched)
head(unmatched)
# Use version 74 from Dec 2013
mart = useMart(host = 'dec2013.archive.ensembl.org', 
               biomart = 'ENSEMBL_MART_ENSEMBL', 
               dataset = 'hsapiens_gene_ensembl')
query <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                             "ensembl_peptide_id", "hgnc_symbol"),
               filters = "ensembl_gene_id", values = unmatched, mart = mart)
length(unmatched) - nrow(query)
# Update
colnames(query) <- colnames(geneLookup)
# Check
before = nrow(geneLookup)
geneLookup <- rbindlist(list(geneLookup, query))
after = nrow(geneLookup)
after == before + nrow(query)

# One last time
tableMatchIdx = fastmatch::fmatch(x = ENSGIDs,
                                  table = geneLookup$`Ensembl Gene ID`,
                                  nomatch = 0)
unmatched = ENSGIDs[which(tableMatchIdx == 0)]
length(unmatched)
head(unmatched)
# ENSGR IDs, for genes in Pseudoautosomal Regions on ChrY
# Let featureData contain info for the ENSG version
# i.e. create rows that have ENSGR as gene IDs, but have other gene related data
# for the ENSG gene IDs
nonENSGR = stringr::str_replace_all(unmatched, "(.*)R(.*)", replacement = "\\10\\2")
# Check lookup table for existence of nonENSGR versions
tempIdx = fastmatch::fmatch(nonENSGR, geneLookup$`Ensembl Gene ID`)
# They all match, no need to query biomaRt
new = data.table("Ensembl Gene ID" = unmatched, geneLookup[tempIdx, 2:4])
geneLookup <- data.table::rbindlist(list(geneLookup, new))
# Check e.g. ENSGR0000264819
a = geneLookup[fmatch("ENSG00000264819", geneLookup$`Ensembl Gene ID`), 2]
b = geneLookup[fmatch("ENSGR0000264819", geneLookup$`Ensembl Gene ID`), 2]
a == b

# ==== Feature Data Preparation ====

# featureData is an AnnotatedDataFram where each row contains information
# regarding each feature, e.g. Gene in this case. Row names must match.
# Create a table where the ENSG IDs are mapped to ENST IDs and HGNC symbols
# Can append other information later on
# Match column names with the lookup table, and retrieve other information 

# Check
tableMatchIdx = fastmatch::fmatch(x = ENSGIDs,
                                  table = geneLookup$`Ensembl Gene ID`,
                                  nomatch = 0)
sum(tableMatchIdx != 0) == nrow(eSetMatrix)

myFeatureData <- Biobase::AnnotatedDataFrame(data = geneLookup[tableMatchIdx])
# Recall that row names must match
Biobase::featureNames(myFeatureData) <- geneLookup[tableMatchIdx]$`Ensembl Gene ID`
myFeatureData

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
# matrix / matricies in assayData. i.e. pcawg[sampleIdx]
prePhenoData = data.table::data.table(pcawg[sampleIdx,])

# Update: Can also append other sample data e.g. file names from the
# Samples_info.rsav file.
load("Samples_info.rsav")
# View(vcfs)
class(vcfs)
# Convert to data.table
vcfs <- data.table::data.table(vcfs)

# We can use "sample ID" (which is icgc_donor_id in prePhenoData) to attach
# related data to data on PCAWG table.
# Note: Must convert NA values in prePhenoData$icgc_donor_id to "" so that they
# are not matched
prePhenoData[is.na(prePhenoData$icgc_donor_id), "icgc_donor_id"] <- ""
vcfsIdx = fastmatch::fmatch(x = vcfs$Sample_ID,
                            table = prePhenoData$icgc_donor_id, nomatch = 0)

# Attach other columns in VCFS to their respective samples in the prePhenoData
# Note: sampleIdx != vcfsIdx
# Add VCFS data (except for "icgc_donor_id")
prePhenoData[vcfsIdx,
             c("File", "vcf_file", "NormalBAM", "TumorBAM")] <- vcfs[, -3]

# Check
temp = prePhenoData[1:10, "icgc_donor_id"]
tempIdx = fmatch(x = unlist(temp), table = vcfs$Sample_ID, nomatch = NA)
prePhenoData[1:10, c("icgc_donor_id", "File")]$File == vcfs$File[tempIdx]

myPhenoData <- Biobase::AnnotatedDataFrame(data = prePhenoData)
# Check row and column numbers
nrow(myPhenoData) == ncol(eSetMatrix)

# Assign the column names of eSetMatrix to row names
Biobase::featureNames(myPhenoData) <- colnames(eSetMatrix)

# Check
sampleNames(myPhenoData)
tempIdx = fmatch(pData(myPhenoData)[1:10, "aliquot_id"], pcawg$aliquot_id)
pcawg$donor_unique_id[tempIdx] == pData(myPhenoData)$`donor_unique_id`[1:10]
myPhenoData


# ==== Create Expression Set ====
?ExpressionSet
# Check
length(featureNames(myFeatureData)) == length(row.names(eSetMatrix))
length(featureNames(myPhenoData)) == length(colnames(eSetMatrix))

all(row.names(eSetMatrix) == featureNames(myFeatureData))
eSet <- ExpressionSet(assayData = eSetMatrix, featureData = myFeatureData,
                      phenoData = myPhenoData)
class(eSet)

# Now we can test Biobase functions:
# sampleNames() returns patient IDs in this case
sampleNames(eSet)[1:5]
ncol(eSet)

# featureNames returns the feature IDs, which may need to be translated
featureNames(eSet)[1:5]
nrow(eSet)

# Save eSet
saveRDS(eSet, "Data/ExpressionSet.rds")


# Now we can load the ExpressionSet and modify it directly
eSet <- readRDS("Data/ExpressionSet.rds")

temp = pData(eSet)[1,]

# TODO: Add cancer type to each sample, create a lookup function to query
# expression data based on features (genes) and samples (donors)


# Load the updated mutation scores files
mut_scores <- readRDS("Data/Updated_mutation_scores/Chr22.rds")
mut_scores

