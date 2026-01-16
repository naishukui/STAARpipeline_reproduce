# STAARpipeline Protocal
This is a tutorial for (1) automatically functionally annotating the variants of whole-genome/whole-exome sequencing (WGS/WES) studies and integrating the functional annotations with the genotype data using **FAVORannotator** and (2) performing association analysis of WGS/WES studies, summarizing and visualization results using **STAARpipeline** and **STAARpipelineSummary**. The software prerequisites, dependencies and installation can be found in <a href="https://github.com/xihaoli/STAARpipeline">**STAARpipeline**</a> and <a href="https://github.com/xihaoli/STAARpipelineSummary">**STAARpipelineSummary**</a> packages.

**FAVORannotator**, **STAARpipeline** and **STAARpipelineSummary** are implemented as a collection of apps for cloud-based platforms. Please see the following apps

**favorannotator** (<a href="https://platform.sb.biodatacatalyst.nhlbi.nih.gov/public/apps/admin/sbg-public-data/favorannotator-1-0-0">**BDC-Seven Bridges**</a>, <a href="https://github.com/xihaoli/favorannotator-rap">**RAP-DNAnexus**</a>) <br>
**staarpipeline** (<a href="https://platform.sb.biodatacatalyst.nhlbi.nih.gov/public/apps/admin/sbg-public-data/staarpipeline-0-9-6">**BDC-Seven Bridges**</a>, <a href="https://github.com/xihaoli/staarpipeline-rap">**RAP-DNAnexus**</a>) <br>
**staarpipelinesummary_varset** (<a href="https://platform.sb.biodatacatalyst.nhlbi.nih.gov/public/apps/admin/sbg-public-data/staarpipelinesummary-varset-0-9-6">**BDC-Seven Bridges**</a>, <a href="https://github.com/xihaoli/staarpipelinesummary_varset-rap">**RAP-DNAnexus**</a>) <br>
**staarpipelinesummary_indvar** (<a href="https://platform.sb.biodatacatalyst.nhlbi.nih.gov/public/apps/admin/sbg-public-data/staarpipelinesummary-indvar-0-9-6">**BDC-Seven Bridges**</a>, <a href="https://github.com/xihaoli/staarpipelinesummary_indvar-rap">**RAP-DNAnexus**</a>)

that run on the <a href="https://biodatacatalyst.nhlbi.nih.gov">**NIH/NHLBI BioData Catalyst (BDC) ecosystem**</a> and the <a href="https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform">**UK Biobank Research Analysis Platform (RAP)**</a> for more details (<a href="tinyurl.com/staarpipelineapps">**user manual and tutorial**</a>).
## Pre-step of association analysis using STAARpipeline 
### Generate Genomic Data Structure (GDS) file
R/Bioconductor package **SeqArray** provides functions to convert the genotype data (in VCF/BCF/PLINK BED/SNPRelate format) to SeqArray GDS format. For more details on usage, please see the R/Bioconductor package <a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">**SeqArray**</a> [<a href="https://bioconductor.org/packages/release/bioc/manuals/SeqArray/man/SeqArray.pdf">manual</a>]. A wrapper for the `seqVCF2GDS`/`seqBCF2GDS` function in the SeqArray package can be found <a href="convertVCF2GDS.R">**here**</a> (**Credit: Michael R. Brown and Jennifer A. Brody**).

R package **gds2bgen** provides functions to convert the genotype data (in BGEN format) to SeqArray GDS format. For more details on usage, please see the R package <a href="https://github.com/zhengxwen/gds2bgen">**gds2bgen**</a>. An example for the `seqBGEN2GDS` function in the gds2bgen package can be found <a href="https://github.com/zhengxwen/gds2bgen#examples">**here**</a> (**Credit: Xiuwen Zheng**).

Note 1: As a file integrity check, it is expected that variant in the GDS file can be **uniquely identified** based on its **CHR-POS-REF-ALT** combination. That is, there shouldn't be two variants in the GDS file with identical CHR-POS-REF-ALT records. It is also expected that the physical positions of variants in the GDS file (of each chromosome) should be sorted in **ascending order**.

Note 2: After the GDS file is generated, there is supposed to be a channel in the GDS file (default is `annotation/filter`) where all variants passing the quality control (QC) should be labeled as `"PASS"`. If there is no such channel for a given post-QC GDS file (where all variants in the GDS file are pass variants), one can create a new channel in the GDS file by setting the value of all variants as `"PASS"`. An example script can be found <a href="Add_QC_label.R">**here**</a>. Then, in all scripts of STAARpipeline, `QC_label <- "annotation/filter"` should be updated to `QC_label <- "annotation/info/QC_label"`.

### Generate annotated GDS (aGDS) file using FAVORannotator
#### Prerequisites:
**FAVORannotator** (CSV version 1.0.0) depends on the **xsv software** and the **FAVOR database** in CSV format. Please install the <a href="https://github.com/BurntSushi/xsv">**xsv software**</a> and download the **FAVOR essential database CSV files** from <a href="http://favor.genohub.org">**FAVOR website**</a> (under the "FAVORannotator" tab's top panel, 31.2 GB for chr1 CSV) or <a href="https://doi.org/10.7910/DVN/1VGTJI">**Harvard Dataverse**</a> before using **FAVORannotator** (CSV version 1.0.0).
#### Step 0: Install xsv
The following steps are for the widely used operating system (Ubuntu) on a virtual machine.

1. Install Rust and Cargo:
 - ```$ curl https://sh.rustup.rs -sSf | sh```
2. Source the environment: 
 - ```$ source $HOME/.cargo/env``` 
3. Install xsv using Cargo:
 - ```$ cargo install xsv```
#### Step 1: Generate the variants list to be annotated
##### Script: <a href="FAVORannotator_csv/Varinfo_gds.R">**Varinfo_gds.R**</a>
##### Input: GDS files of each chromosome and the FAVOR database information <a href="FAVORannotator_csv/FAVORdatabase_chrsplit.csv">**FAVORdatabase_chrsplit.csv**</a>. For more details, please see the R script.
##### Output: CSV files of the variants list. For each chromosome, the number of CSV files is listed in <a href="FAVORannotator_csv/FAVORdatabase_chrsplit.csv">**FAVORdatabase_chrsplit.csv**</a>.
Note: The physical positions of variants in the GDS file (of each chromosome) should be sorted in ascending order.

#### Step 2: Annotate the variants using the FAVOR database through xsv software
##### Script: <a href="FAVORannotator_csv/Annotate.R">**Annotate.R**</a>
##### Input: CSV files of the variants list to be annotated, the FAVOR database information <a href="FAVORannotator_csv/FAVORdatabase_chrsplit.csv">**FAVORdatabase_chrsplit.csv**</a>,
the FAVOR database, and the directory xsv software. For more details, please see the R script.
##### Output: CSV files of the annotated variants list. 
* `Anno_chrXX.csv`: a CSV file containing annotated variants list of chromosome XX. <br>
* `Anno_chrXX_STAARpipeline.csv`: a CSV file containing the variants list with annotations required for STAARpipeline of chromosome XX. 
The annotations in this file is a subset of `Anno_chrXX.csv`. <br>

#### Step 3: Generate the annotated GDS (aGDS) file
##### Script: <a href="FAVORannotator_csv/gds2agds.R">**gds2agds.R**</a>
##### Input: GDS files and the CSV files of annotated variants list (`Anno_chrXX.csv` or `Anno_chrXX_STAARpipeline.csv`). For more details, please see the R script.
##### Output: aGDS files including both the genotype and annotation information.
Note: FAVORannotator also supports the database in SQL format. Please see the <a href="https://github.com/zhouhufeng/FAVORannotator">**FAVORannotator** tutorial</a> for detailed usage of **FAVORannotator** (SQL version).

### Generate sparse Genetic Relatedness Matrix (GRM)
R package **FastSparseGRM** provides functions and a pipeline to efficiently calculate genetic principal components (PCs) and the ancestry-adjusted sparse genetic relatedness matrix (GRM). It accounts for population heterogeneity using genetic PCs which are automatically calculated as part of the pipeline. The genetic PCs can be used as fixed effect covariates to account for the population stratification and the sparse GRM can be used to model the random effects to account for the sample relatedness in a mixed effects phenotype-genotype association testing model implemented in STAARpipeline. For more details on usage, please see the R package <a href="https://github.com/rounakdey/FastSparseGRM">**FastSparseGRM**</a> and <a href="https://doi.org/10.21203/rs.3.rs-5343361/v1">manuscript</a>.

## Association analysis using STAARpipeline
### Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
#### Script: <a href="Association_Analysis_PreStep.r">**Association_Analysis_PreStep.r**</a>
#### Input: aGDS files of all 22 chromosomes. For more details, please see the R script.
#### Output: `agds_dir.Rdata`, `Annotation_name_catalog.Rdata`, `jobs_num.Rdata`.
* `agds_dir.Rdata`: a vector containing directory of GDS/aGDS files of all chromosomes. <br>
* `Annotation_name_catalog.Rdata`: a data frame containing the annotation name and the corresponding channel name in the aGDS file. Alternatively, one can skip this part in the R script by providing `Annotation_name_catalog.csv` with the same information. An example of `Annotation_name_catalog.csv` can be found <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/Annotation_name_catalog.csv">here</a>. <br>
* `jobs_num.Rdata`: a data frame containing the number of jobs for association analysis, including individual analysis, sliding window analysis and dynamic window analysis (SCANG-STAAR).

### Step 1: Fit STAAR null model
#### Script: <a href="STAARpipeline_Null_Model.r">**STAARpipeline_Null_Model.r**</a> or <a href="STAARpipeline_Null_Model_GENESIS.r">**STAARpipeline_Null_Model_GENESIS.r**</a> or <a href="STAARpipeline_Null_Model_Multi.r">**STAARpipeline_Null_Model_Multi.r**</a>
* `STAARpipeline_Null_Model.r` fits the STAAR null model using the STAARpipeline package. <br>
* `STAARpipeline_Null_Model_GENESIS.r` fits the null model using the GENESIS package and convert it to the STAAR null model using the STAARpipeline package. <br>
* `STAARpipeline_Null_Model_Multi.r` fits the MultiSTAAR null model using the STAARpipeline package. <br>
#### Input: Phenotype data and (sparse) genetic relatedness matrix. For more details, please see the R scripts.
#### Output: a Rdata file of the STAAR or MultiSTAAR null model.
Note: Once the <a href="https://github.com/xihaoli/STAAR">STAAR</a> or <a href="https://github.com/xihaoli/MultiSTAAR">MultiSTAAR</a> null model is fit, all the remaining steps of STAARpipeline and STAARpipelineSummary share the same scripts (the information of single-trait or multi-trait analysis being considered is automatically retrieved from the null model object).

### Step 2: Individual (single-variant) analysis
#### Script: <a href="STAARpipeline_Individual_Analysis.r">**STAARpipeline_Individual_Analysis.r**</a>
Perform single-variant analysis for common and low-frequency variants across the genome using the STAARpipeline package. 
#### Input: aGDS files and the STAAR or MultiSTAAR null model. For more details, please see the R script.
#### Output: Rdata files with the user-defined names.
The number of output files is the summation of the column "individual_analysis_num" for the object in `jobs_num.Rdata`.

### Step 3.1: Gene-centric coding analysis
#### Script: <a href="STAARpipeline_Gene_Centric_Coding.r">**STAARpipeline_Gene_Centric_Coding.r**</a> and <a href="STAARpipeline_Gene_Centric_Coding_Long_Masks.r">**STAARpipeline_Gene_Centric_Coding_Long_Masks.r**</a>
Perform gene-centric analysis for coding rare variants using the STAARpipeline package. The gene-centric coding analysis provides five functional categories to aggregate coding rare variants of each protein-coding gene: (1) putative loss of function (stop gain, stop loss, and splice) RVs, (2) missense RVs, (3) disruptive missense RVs, (4) putative loss of function and disruptive missense RVs, and (5) synonymous RVs. <br>
* `STAARpipeline_Gene_Centric_Coding.r` performs gene-centric coding analysis for all protein-coding genes across the genome. There are 379 jobs using this script. <br>
* `STAARpipeline_Gene_Centric_Coding_Long_Masks.r` performs gene-centric coding analysis for some specific long masks, and might require larger memory compared to `STAARpipeline_Gene_Centric_Coding.r`. There are 2 jobs using this script.
#### Input: aGDS files and the STAAR or MultiSTAAR null model. For more details, please see the R scripts.
#### Output: 381 Rdata files with the user-defined names.
