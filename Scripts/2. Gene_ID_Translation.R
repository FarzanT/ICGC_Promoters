# Gene_ID_Translation.R

# ==== Translate between ENSG, ENST and HGNC IDs ====
# Packages (install if needed):
library(fastmatch)
library(readr)
library(biomaRt)
library(stringr)

# Load ID translation data
table = read_tsv("LookupTables/IDs_convert_Table.txt")
table
colnames(table)

# Create a dictionary where ENSG and ENST can be converted to HGNC symbols
listMarts()
# ENSEMBL_MART_ENSEMBL
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")
listDatasets(mart)

# No.41, 'hsapiens_gene_ensembl'
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")
grep(".*transcript.*", listAttributes(mart)$name, value = T)
# No.99, "ensembl_transcript_id"

ensgNames <- table[, 1]
enstNames <- table[, 2]

HGNCfromG <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                   filters = "ensembl_gene_id", values = ensgNames, mart = mart)
HGNCfromT <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id"),
                   filters = "ensembl_gene_id", values = ensgNames, mart = mart)

# Ordering has changed
ensgNames[1:10,] == HGNCfromG$ensembl_gene_id[1:10]


# Find the rows where there is HGNC symbol, get rid of failed queries
any(is.na(table$`HGNC symbol`))
sum(is.na(table$`HGNC symbol`))
idx <- which(is.na(table$`HGNC symbol`))
# Get the ENSG and ENSTs
ENS <- table[idx, 1:2]
# Find these ENS IDs in the retrieved data from biomaRt to see if any HGNC 
# symbols were found
ENSGidx <- fastmatch::fmatch(x = unlist(ENS[,1]), table = HGNCfromG$ensembl_gene_id)
ENSTidx <- fastmatch::fmatch(x = unlist(ENS[,2]), table = HGNCfromT$ensembl_transcript_id)

# Delete null ENSG rows
HGNCfromG <- HGNCfromG[-which(is.na(HGNCfromG[ENSGidx,])), ]
HGNCfromG <- HGNCfromG[-which(HGNCfromG[ENSGidx,] == ""), ]

# Delete null ENST rows
HGNCfromT <- HGNCfromT[-which(is.na(HGNCfromT[ENSTidx,])), ]
HGNCfromT <- HGNCfromT[-which(HGNCfromT[ENSTidx,] == ""), ]

# ==== Update HGNC symbols by ENST ====
ENSTmatchIdx = fmatch(x = HGNCfromT$ensembl_transcript_id,
                      table = table$`Ensembl Transcript ID`)
# If the HGNC slot (in the table) is empty, update it
table[which(is.na(table[,4]))[1],]
# Get ENST IDs from rows with missing HGNC symbols
missingIdx = which(is.na(table[,4]))
# Confirm
table[missingIdx, ][1:10, ]
enstOfMissing = table[missingIdx, 2]
# Check HGNCfromT table to see if HGNC symbol was retrieved
tempIdx = fmatch(unlist(enstOfMissing), HGNCfromT$ensembl_transcript_id)
temp = HGNCfromT[tempIdx, ]
temp = temp[!(temp$hgnc_symbol == ""), ]
temp = temp[!(is.na(temp$hgnc_symbol)), ]
temp = temp[!(is.na(temp$ensembl_transcript_id)), ]
temp

# Now find ENSTs from temp in the table
updateIdx = fmatch(x = temp$ensembl_transcript_id,
                   table = table$`Ensembl Transcript ID`)
# Update
table$`HGNC symbol`[updateIdx] <- temp$hgnc_symbol
table


# ==== Update HGNC symbols by ENSG ====
ENSGmatchIdx = fmatch(x = HGNCfromG$ensembl_gene_id,
                      table = table$`Ensembl Gene ID`)
# If the HGNC slot (in the table) is empty, update it
table[which(is.na(table[,4]))[1],]
# Get ENST IDs from rows with missing HGNC symbols
missingIdx = which(is.na(table[,4]))
# Confirm
table[missingIdx, ][1:10, ]
ensgOfMissing = table[missingIdx, 1]
# Check HGNCfromT table to see if HGNC symbol was retrieved
tempIdx = fmatch(unlist(ensgOfMissing), HGNCfromG$ensembl_gene_id)
temp = HGNCfromG[tempIdx, ]
temp = temp[!(temp$hgnc_symbol == ""), ]
temp = temp[!(is.na(temp$hgnc_symbol)), ]
temp = temp[!(is.na(temp$ensembl_gene_id)), ]
temp

# Now find ENSGs from temp in the table
updateIdx = fmatch(x = temp$ensembl_gene_id, table = table$`Ensembl Gene ID`)
# Update
table$`HGNC symbol`[updateIdx] <- temp$hgnc_symbol
table


# Reload the original table for comparison
origTable = read_tsv("LookupTables/IDs_convert_Table.txt")
sum(is.na(origTable$`HGNC symbol`)) - sum(is.na(table$`HGNC symbol`))
# 1699 Change (hopefully in a good way)

# Save the table
write_delim(table, "LookupTables/IDs_convert_Table_updated.txt", col_names = T)

# [END]
