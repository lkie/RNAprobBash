options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

library(RNAprobR)

# Read inputs
merged <- args[1]
unique_barcodes <- args[2]
output <- args[3]

k2n_calc( merged, unique_barcodes, output )

q("no")
