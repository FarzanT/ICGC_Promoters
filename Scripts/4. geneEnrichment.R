# geneEnrichment.R

library(Biobase)
library(limma)
library(oligo)
library(clusterProfiler)
library(igraph)

# Load eSet
eSet <- readRDS("Data/ExpressionSet.rds")
eSet
featureData(eSet)
# Load a chromosome's mutation data
gr.chrY <- readRDS("Data/PromoterOverlaps/ChrY.rds")
gr.chrY

Biobase::exprs(eSet)[1,]

# Type and number of tumors
table(eSet$tumor_type)
which(is.na(eSet$tumor_type))

eSet$tumor_type

anyDuplicated(sampleNames(eSet))

# Subset the data for Kidney-RCC cancers, separating normal tissue from tumor
curData <- eSet[, eSet$tumor_type == "Kidney-RCC"]
varLabels(curData)
# Discard NAs for now
curData <- curData[,!(is.na(curData$dcc_specimen_type))]

design <- model.matrix(~ curData$dcc_specimen_type)
head(design)

fit <- limma::lmFit(curData, design)
fit <- limma::eBayes(fit)
topTable(fit)
View(fit)


# Find the mean expression of a feature (gene) in a certain group
gene = featureNames(eSet)[1]
?tapply
typeMean = tapply(exprs(curData)[gene,], curData$tumor_type, mean)

boxplot(curData)

# ==== clusterProfiler ====
clusterProfiler::enricher(featureNames(eSet), featureData(eSet)$HGNC)

# ==== Create a lookup function ====
# Must retrieve samples based on tumor type
varLabels(eSet)
(pData(eSet)[1,])


#'enrichSamples
#'
#'Find and enrich samples based on tumor type, specimen type,
enrichSamples <- function(myTumorType, mySpecimenType = "contrast", myESet) {
  if(missing(myTumorType)) {
    stop("Must provide a valid tumor type")
  } else {
    if (!(myTumorType %in% unique(pData(eSet)$tumor_type))) {
      stop(cat("Tumor type is invalid, available tumor types are:",
                      unique(pData(eSet)$tumor_type), sep = "\n"))
    }
  }
  # if (missing(mySpecimenType)) {
  #   stop("Must provide a valid specimen type")
  # } else {
  #   if (!(mySpecimenType %in% unique(pData(eSet)$dcc_specimen_type))) {
  #     stop(cat("Specimen type is invalid, available specimen types are:",
  #              unique(pData(eSet)$dcc_specimen_type), sep = "\n"))
  #   }
  # }
  if (missing(myESet)) {
    stop("Must provide a valid expression set")
  }
  # Subset the expression set based on given tumor_type
  curData <- myESet[ , pData(myESet)$tumor_type == myTumorType]
  # Subset based on dcc_specimen_type
  if (mySpecimenType != "contrast") {
    curData <- curData[ , pData(myESet)$dcc_specimen_type == mySpecimenType]
  } else {
    # Find compare normal and tumor samples (constrast)
      
  #   curData <- curData[ , pData(myESet)$dcc_specimen_type %in% c("Normal - solid tissue",
  #                                                                "Normal - tissue adjacent to primary")]
  }
  
}
