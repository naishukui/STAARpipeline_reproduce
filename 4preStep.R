###########################################################
# Pre-step for running STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/21/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
###########################################################
#           User Input
###########################################################
## file directory of aGDS file (genotype and annotation data)
dir.geno <- "/rsrch6/home/biostatistics/nkui/ukbiobank_gds/ukb_panc/"
agds_file_name_1 <- "ukb_imp_chr"
agds_file_name_2 <- "_v3_qc.gds"
QC_label <- "annotation/filter"
output_path <- "/rsrch6/home/biostatistics/nkui/ukbiobank/output12/"



###########################################################
#           Main Function
###########################################################
#### aGDS directory
agds_dir <- paste0(dir.geno,agds_file_name_1,seq(1,22),agds_file_name_2)
save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

#### Annotation name catalog (alternatively, can skip this part by providing Annotation_name_catalog.csv with the same information)
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
          "aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info",
         "/genecode_comprehensive_exonic_category","/metasvm_pred",
         "/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))

#### Number of jobs for each chromosome
jobs_num <- matrix(rep(0,66),nrow=22)

#run this instead
#modify filter
for(chr in 1:22)
{
  print(chr)
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path,readonly=FALSE)
  position <- as.numeric(seqGetData(genofile, "position"))
  # filter <- read.gdsn(index.gdsn I 12(genofile, "annotation/oldfilter"))   14ss34erw
  filter <- read.gdsn(index.gdsn(genofile, "annotation/filter"))

  new<-cbind.data.frame(filter,position)
  new$newfilter[position!=0]<-'PASS'


  newfilter=as.factor(new$newfilter)

  #filter<- index.gdsn(genofile, "annotation/filter")
  #delete.gdsn(filter,force=TRUE)
  #rename.gdsn(filter,"oldfilter")
  Anno <- index.gdsn(genofile, "annotation")
  add.gdsn(Anno, "filter", val=newfilter, compress="LZMA_ra", replace=TRUE,closezip=TRUE)

  seqClose(genofile)
}

for(chr in 1:22)
{
  print(chr)
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  position <- as.numeric(seqGetData(genofile, "position"))
  filter <- seqGetData(genofile, QC_label)
  SNVlist <- filter == "PASS"

  jobs_num[chr,1] <- chr
  jobs_num[chr,2] <- min(position[SNVlist],na.rm = TRUE)
  jobs_num[chr,3] <- max(position[SNVlist],na.rm = TRUE)

  seqClose(genofile)
}



# Individual Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/10e6))
# Sliding Window Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/5e6))
# Dynamic Window Analysis (SCANG-STAAR)
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/1.5e6))

colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num","sliding_window_num","scang_num")
jobs_num <- as.data.frame(jobs_num)

save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))




#-----------------------previous code :do not use--------------------------------------------------------#


for(chr in 1:22)
{
  print(chr)
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  position <- as.numeric(seqGetData(genofile, "position"))
  filter <- seqGetData(genofile, QC_label)
  table(filter)
  seqClose(genofile)
}

#previous code:do nto use
####previous code before fix lift over####

for(chr in 1:22)
{
  print(chr)
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path,readonly=FALSE)

  position <- as.numeric(read.gdsn(index.gdsn(genofile, "position")))
  filter <- seqGetData(genofile, QC_label)


  new<-cbind.data.frame(filter,position)
  nezw$newfilter[filter %in% c('VQSRTrancheSNP90.00to95.00', 'VQSRTrancheSNP95.00to98.00', 'VQSRTrancheSNP98.00to99.00',
                               'VQSRTrancheSNP99.00to99.10', 'VQSRTrancheSNP99.10to99.20', 'VQSRTrancheSNP99.20to99.30',
                               'VQSRTrancheSNP99.30to99.40', 'VQSRTrancheSNP99.40to99.50', 'VQSRTrancheSNP99.50to99.60',
                               'VQSRTrancheSNP99.60to99.70', 'VQSRTrancheSNP99.70to99.80', 'VQSRTrancheSNP99.80to99.90',
                               'VQSRTrancheSNP99.90to99.91', 'VQSRTrancheSNP99.91to99.92', 'VQSRTrancheSNP99.92to99.93',
                               'VQSRTrancheSNP99.93to99.94', 'VQSRTrancheSNP99.94to99.95', 'VQSRTrancheSNP99.95to99.96',
                               'VQSRTrancheSNP99.96to99.97', 'VQSRTrancheSNP99.97to99.98', 'VQSRTrancheSNP99.98to99.99',
                               'VQSRTrancheSNP99.99to100.00+', 'VQSRTrancheSNP99.99to100.00', 'PASS') & position!=0]<-'PASS'

  newfilter=as.factor(new$newfilter)

  SNVlist <- newfilter == "PASS"



  jobs_num[chr,1] <- chr
  jobs_num[chr,2] <- min(position[SNVlist],na.rm = TRUE)
  jobs_num[chr,3] <- max(position[SNVlist],na.rm = TRUE)
  jobs_num[chr,4] <- max(position[SNVlist],na.rm = TRUE)
  filter<- index.gdsn(genofile, "annotation/filter")
  rename.gdsn(filter,"oldfilter")
  Anno <- index.gdsn(genofile, "annotation")
  add.gdsn(Anno, "filter", val=newfilter, compress="LZMA_ra", replace=TRUE,closezip=TRUE)

  seqClose(genofile)
}
