# geneEnrichment.R

library(Biobase)
library(limma)

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

# Subset the data for different Lymphoid cancers
curData <- eSet[, eSet$tumor_type %in% c("Lymph-BNHL", "Lymph-CLL")]

design <- model.matrix(~ curData$tumor_type)
head(design)
class(exprs(eSet)[1,])
fit <- limma::lmFit(curData, design)
fit <- limma::eBayes(fit)
topTable(fit)


# Find the mean expression of a feature (gene) in a certain group
gene = featureNames(eSet)[1]
typeMean = tapply(exprs(curData)[gene,], curData$tumor_type, mean)

