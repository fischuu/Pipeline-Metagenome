# This is a development script to summarise the data it will be later generalised and added to the pipeline

projFolder <- "/scratch/project_2006940/MASTER"
coas <- 2

library("GenomicTools")
library("data.table")

filename <- file.path(projFolder, "PRODIGAL", paste0("final.contigs_group_", coas), paste0("final.contigs_group_",coas,".prodigal_plain.gtf"))

prodGFF <- importGFF(filename, level="CDS")

# Now import the FeatureCounts quantification

fc.files <- list.files(file.path(projFolder, "QUANTIFICATION", "PRODIGAL_FC", paste0("final.contigs_group_", coas)), pattern="*.txt$")

# This is now for me to speed up things, remove it later for production
fc.files <- fc.files[grep("FEC", fc.files)]

tmp <- importFeatureCounts(file.path(projFolder, "QUANTIFICATION", "PRODIGAL_FC", paste0("final.contigs_group_", coas), fc.files[1]))
tmpExp <- tmp$expValues
colnames(tmpExp) <- gsub("X.scratch.project_2006940.MASTER.BAM.final.contigs_coas2.LUKE_", "", colnames(tmpExp))
colnames(tmpExp) <- gsub("_mega.bam", "", colnames(tmpExp))

fcExp <- matrix(0, ncol=length(fc.files)+2, nrow=nrow(tmpExp))
fcExp <- as.data.table(fcExp)
rownames(fcExp) <- tmpExp[,1]
colnames(fcExp)[1] <- "Geneid"
colnames(fcExp)[2] <- "Eggnog_id"
colnames(fcExp)[3] <- colnames(tmpExp)[2]
fcExp[,1] <- tmp$geneInfo$Geneid
fcExp[,2] <- paste0(tmp$geneInfo$Chr,"_",sapply(strsplit(tmp$geneInfo$Geneid, "_"),"[",2))
fcExp[,3] <- tmpExp[,2]

for(i in 2:length(fc.files)){
  tmp <- importFeatureCounts(file.path(projFolder, "QUANTIFICATION", "PRODIGAL_FC", paste0("final.contigs_group_", coas), fc.files[i]))
  tmpExp <- tmp$expValues
  colnames(tmpExp) <- gsub("X.scratch.project_2006940.MASTER.BAM.final.contigs_coas2.LUKE_", "", colnames(tmpExp))
  colnames(tmpExp) <- gsub("_mega.bam", "", colnames(tmpExp))

  colnames(fcExp)[i+2] <- colnames(tmpExp)[2]
  fcExp[,i+2] <- tmpExp[,2]
  
  
    cat("Sample",i,"/", length(fc.files),"processed -",date(),"\n")
}

fwrite(fcExp, file.path(projFolder,"fcExp_coas2.csv"))
