# Merge_Mutation_Info.R

library(readr)
library(data.table)
install.packages("qdapTools")
library(qdapTools)
library(fastmatch)

# Load the lookup table
lookupTable <- read.csv("pcawg_sample_sheet.v1.4.2016-09-14.tsv", sep = "\t")
View(lookupTable)

# PWM index and name conversion table
load("pwms_2016_01.rsav")
head(pwms)
View(pwms)
# Tumors are separated by type in the meta_tumor_cohorts
# Mutation scores (by chromosome) are in mutation_scores

load("Promoters-Encode.rsav")
head(promoters)

# P values 
load("PValues_new.rsav")
View(pValues)

# Loop over the files in the meta_tumor_cohorts
metaTumor <- list.files(path = "meta_tumor_cohorts/", pattern = ".*whitelist.txt", full.names = T)
allTumors <- vector("list", length(metaTumor))
for (i in 1:length(metaTumor)) {
  allTumors[[i]] <- read.delim(metaTumor[i])
}
head(allTumors)


nonNormalized <- read_tsv("tophat_star_fpkm_uq.v2_sample_gl.tsv")
head(nonNormalized)

load("Samples_info.rsav")
head(vcfs)

# Submitter donor ID in lookup table (take sample ID) + 
# allTumors (meta_tumor_cohort) + 
head(lookupTable)
head(allTumors)

promoters

load("ICGC_All_SNVs_Final.rsav")
head(all_mut)


# Append new column for cancer type to the lookup table
any(allTumors[[1]][1,] == lookupTable$submitter_donor_id)
X <- F
for (i in 1:14) {
  if (any(unlist(allTumors[[i]]) == "096b4f32-10c1-4737-a0dd-cae04c54ee33")) {
    x <- TRUE
  }
  
}
X
print(lookupTable$aliquot_id[1])
(lookupTable$aliquot_id[[1]])
class(lookupTable$aliquot_id[[1]])


# ==== Append wildtype and mutant p values and name of pwm to each mutation_score file ====

# Create a named vector as a dictionary to translate PWM index to its name
PWMdict <- vector("character", length = length(pwms))
for (i in 1:length(pwms)) {
  names(PWMdict)[i] <- i
  PWMdict[i] <- names(pwms)[i]
}
PWMdict

# Create the same type of dictionary for column names in pValues for faster
# access
pValCols <- new.env(size = ncol(pValues))
for (i in 2:ncol(pValues)) {
  pValCols[[colnames(pValues)[i]]] <- i - 1
}
ls(pValCols)
# Now append the names of PWMs, p values and chromosome number to each row of 
# mutation_scores files:
# Create a new directory for the updated mutation_score files
dir.create("Data")
dir.create("Data/Updated_mutation_scores")
# Loop over the files in mutation_scores directory
chrScores <- list.files(path = "mutation_scores/", full.names = T)
genomeScores <- vector("list", 24)

for (i in 1:length(genomeScores)) {
  # Extract chromosome number
  chromNum = gsub(".*Chr(.*)\\.rds", "\\1", chrScores[i])
  # Read chromosome mutation_score file
  tempChrom = readRDS(chrScores[i])
  # Prepend chromosome number to 'mid'
  tempChrom$mid <- paste0(chromNum, tempChrom$mid)
  # Convert character columns to numeric for easier operation
  tempChrom$pwm_index <- as.numeric(tempChrom$pwm_index)
  tempChrom$wt.scores_p <- as.numeric(tempChrom$wt.scores_p)
  tempChrom$mt.scores_p <- as.numeric(tempChrom$mt.scores_p)
  # Append pwm_id for each row based on the pwm_index
  tempChrom$pwm_id <- PWMdict[tempChrom$pwm_index]
  # Get the p values for each PWM from the pValues file
  
  # Find the closest score to wt.scores_p for each row, and use the pwm_id to
  # append wt_pval and mt_pval:
  
  # Append P value for each PWM (this is faster):
  PWMids = unique(tempChrom$pwm_id)
  # Run in parallel
  result <- mclapply(PWMids, function(i) {
    # Get a subset of this chromosome for a specific PWM
    tempSubset = tempChrom[tempChrom$pwm_id == i]
    # Get p values associated with this PWM
    tempPVals = data.frame(seq(1, 101, by = 1), pValues[ , get(i)])
    # Convert scores to indices in the p value table
    wtIdx = tempSubset$wt.scores_p * 100
    mtIdx = tempSubset$mt.scores_p * 100
    # Lookup p values corresponding scores
    wt_Pvals = lookup(wtIdx, tempPVals, key.reassign = NULL)
    mt_Pvals = lookup(mtIdx, tempPVals, key.reassign = NULL)
    # Bind results
    tempSubset = cbind(tempSubset, wt_Pvals, mt_Pvals)
    return(tempSubset)
  }, mc.cores = 4)
  
  # Bind all results (list of data.tables)
  final <- rbindlist(result)
  # Save to the new folder that was created above
  saveRDS(object = final, file = paste0("Data/Updated_mutation_scores/Chr", chromNum, ".rds"))
  
}

# Update May 8th, 2017
# ==== Add chromosome number, base conversion and position as separate columns ====
allFiles = list.files("Data/Updated_mutation_scores/", full.names = TRUE)
for (file in allFiles[1:22]) {
  # Load the file
  tempFile = readRDS(file)
  # Get the mid's, separate chromosome, bases and positions by a '.'
  tempString = stringr::str_replace_all(string = tempFile$mid,
                                        pattern = "(\\d{1,2})(\\w)(.*)(\\w)",
                                        replacement = "\\1\\.\\2\\4\\.\\3")
  # Split by '.'
  splits = stringr::str_split(string = tempString, pattern = "\\.", 3, simplify = TRUE)
  # Assign as new columns
  tempFile$`chromosome` <- splits[, 1]
  tempFile$`mutation` <- splits[, 2]
  tempFile$`position` <- splits[, 3]
  # Rearrange columns
  tempFile <- tempFile[,c("mid", "chromosome", "mutation", "position",
                          "wt.scores_p", "mt.scores_p", "diff", "pwm_id",
                          "pwm_index", "wt_Pvals", "mt_Pvals")]
  saveRDS(tempFile, file = file)
}

# Do the same for X and Y
for (file in allFiles[23:24]) {
  tempFile = readRDS(file)
  tempString = stringr::str_replace_all(string = tempFile$mid,
                                        pattern = "(\\w)(\\w)(.*)(\\w)",
                                        replacement = "\\1\\.\\2\\4\\.\\3")
  splits = stringr::str_split(string = tempString, pattern = "\\.", 3, simplify = TRUE)
  tempFile$`chromosome` <- splits[, 1]
  tempFile$`mutation` <- splits[, 2]
  tempFile$`position` <- splits[, 3]
  tempFile <- tempFile[,c("mid", "chromosome", "mutation", "position",
                          "wt.scores_p", "mt.scores_p", "diff", "pwm_id",
                          "pwm_index", "wt_Pvals", "mt_Pvals")]
  saveRDS(tempFile, file = file)
}

# ==== Prepend 'chr' to chromosome names ====
allFiles = list.files("Data/Updated_mutation_scores/", full.names = TRUE)
for (file in allFiles) {
  tempFile <-  readRDS(file)
  tempFile$chromosome <- paste0("chr", tempFile$chromosome)
  saveRDS(tempFile, file)
}

test = readRDS(file)
test

# Update
# ==== Separate bases in the mutation ====
allFiles <- list.files("Data/Updated_mutation_scores/", full.names = TRUE, ignore.case = T)
for (file in allFiles) {
  tempChr <- readRDS(file)
  tempMut = stringr::str_split(string = tempChr$mutation,
                               pattern = "", n = 2, simplify = TRUE)
  tempChr$`from` <- tempMut[,1]
  tempChr$`to` <- tempMut[,2]
  tempChr <- tempChr[,c(1,2,12,13,4,5,6,7,8,9,10,11)]
  saveRDS(tempChr, file)
}


# ==== Attach strand information (+ or -) based on base conversion ====
library(BiocInstaller)
library(AnnotationHub)

AH <- AnnotationHub()
q <- query(AH, c("hg19", "gene", "gc"))
q

# Load the hg19 genome
# BSgenome::available.genomes()
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# library(BSgenome.Hsapiens.UCSC.hg19)
# # Load annotations for hg19
# library(AnnotationHub)
# AH <- AnnotationHub()
# # UCSC
# q <- query(AH, c("homo", "genes"))
# q
# # Ensembl
# q2 <- query(AH, c("homo", "ensg"))
# q2
# ensemblGenes <- q2[[1]]
# seqlevels(ensemblGenes)
# ensg <- keepStandardChromosomes(ensemblGenes, pruning.mode = "coarse")
# ensg
# # Create a view
# a = Views(Hsapiens$chr22, 19278657, width = 1)
# Views(Hsapiens$chr22, 19275639, width = 1)
# Views(Hsapiens$chr22, 40290257, width = 1)
## Confirmed, mid's are base conversions on the hg19

# Read and update mutation_score files with strand information
allFiles <- list.files("Data/Updated_mutation_scores/", full.names = TRUE)
for (file in allFiles) {
  tempChr <- readRDS(file)
  # Create a GRanges object
  tempGR <- makeGRangesFromDataFrame(
    df = tempChr,
    keep.extra.columns = TRUE,
    seqnames.field = "chromosome",
    start.field = "position",
    end.field = "position")
  tempGR$from
  tempGR
  # Hsapiens has the "+" strand
  vi = Views(Hsapiens, tempGR)
  if (!(all(DNAStringSet(vi) == tempGR$from))) {
    X <- 100
  }
}
X  
# Confirmed, all mutations correspond to the positive strand

# ==== Extract promoter info from Promoters.bed file ====
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

# ==== Attach ENSG/ENST IDs ====
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)
library(data.table)
library(plyr)
# Read promoter data
GRproms <- readRDS("Data/GRproms.rds")
# # Disjoin the promoters
# reducedProms <- GenomicRanges::reduce(GRproms, with.revmap = TRUE)
# revmap <- mcols(reducedProms)$revmap
# # Retrieve the corresponding IDs
# ENSG <- IRanges::extractList(GRproms$ensembl_gene_id, revmap)
# HGNC <- IRanges::extractList(GRproms$hgnc_symbol, revmap)
# Must assign which gene's promoter each mutation belongs to.
dir.create("Data/PromoterOverlaps")
allFiles <- list.files("Data/Updated_mutation_scores/", full.names = TRUE)
for (file in allFiles) {
  tempChr <- readRDS(file)
  # Order the mutations by position (to correspond to the GRanges), add strand
  data.table::setorder(tempChr,-position)
  tempChr$strand <- "+"
  # Create a GRanges object
  tempGR <- GenomicRanges::makeGRangesFromDataFrame(
    df = tempChr,
    start.field = "position",
    end.field = "position",
    strand.field = "strand",
    seqnames.field = "chromosome",
    keep.extra.columns = TRUE
  )
  # Subset those that fall inside a promoter region, as defined in GRproms
  # Ignore strand, as a mutation on one will lead to a mutation on the other
  ov = IRanges::subsetByOverlaps(tempGR, GRproms, ignore.strand = TRUE)
  ovIdx = GenomicRanges::findOverlaps(tempGR, GRproms, ignore.strand = TRUE)
  ovIdx = as.data.table(ovIdx)[, by = list(queryHits, subjectHits)]
  ovIdx$`HGNC` = GRproms$hgnc_symbol[ovIdx$subjectHits]
  ovIdx$`ENSG` = GRproms$ensembl_gene_id[ovIdx$subjectHits]
  # Aggregate rows by queryHits (i.e. all HGNC and ENSG IDs for a hit can be 
  # captured)
  merged = aggregate(ovIdx, by = list(ovIdx$queryHits), c)
  
  # These indices are rows in the GRproms object
  
  # Attach IDs
  ov$`hgnc_symbol` <- merged$HGNC
  ov$`ensemble_gene_id` <- merged$ENSG
  
  # Save
  saveRDS(tempChr, file)
  saveRDS(ov, paste0("Data/PromoterOverlaps/", basename(file)))
}

# [END]
