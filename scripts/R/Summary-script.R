library("GenomicTools")
library("data.table")

samplesheet <- args[1]
output.file <- args[2]

cat("Use the samplesheet :", samplesheet, "\n")
cat("Use the output file :", output.file, "\n")

