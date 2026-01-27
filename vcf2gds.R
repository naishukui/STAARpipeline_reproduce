library(SeqVarTools)
vcf_directory <- "/rsrch6/home/biostatistics/nkui/vcfnew2/"
gds_directory <- "/rsrch6/home/biostatistics/nkui/ukbiobank/output6/agds/"
i <- as.numeric(commandArgs(TRUE)[1])
vcffile <- paste0(vcf_directory, "chr", i, ".vcf")
gdsfile <- paste0(gds_directory, "chr", i, ".gds")
seqVCF2GDS(vcffile, gdsfile, fmt.import="GT", storage.option="LZMA_RA")

