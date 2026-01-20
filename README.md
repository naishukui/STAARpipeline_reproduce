# STAARpipeline Protocol

This repository contains the implementation of the **STAAR (Variant Set Test for Association using Annotation information)** protocol. The pipeline manages data conversion, variant standardization, functional annotation, and gene-centric association analysis.

## 📝 Background & Reproducibility
The protocol provided here is designed to reproduce the results presented in the following study:
> **"Large Impact of Genetic Data Processing Steps on Stability and Reproducibility of Set-Based Analyses in Genome-Wide Association Studies"** > [bioRxiv (2025)](https://www.biorxiv.org/content/10.1101/2025.07.21.665850v1)

While based on the original STAAR framework, this version includes **enhanced features** designed to improve data stability and consistency across different genomic sources.
---

## 🛠 Pipeline Workflow

| Step | Script | Description |
| :--- | :--- | :--- |
| **Step 0** | `0vcf2gds.R` | Convert VCF to GDS format. |
| **Step 1** | `1Varinfo_gds.R` | Extract variant info and standardize REF/ALT. |
| **Step 2** | `2Annotate.R` | Annotate variants via FAVOR database. |
| **Step 3** | `3gds2agds.R` | Integrate annotations into GDS (create aGDS). |
| **Step 4** | `4preStep.R` | Calculate job scaling based on variant density. |
| **Step 5** | `5pipeline_coding_combine.R` | Perform core STAAR association analysis. |
| **Step 6** | `summary_gene_centric_coding.R` | Result summarization and Manhattan plotting. |

---

## 📖 Step-by-Step Instructions

### Step 0: VCF to GDS Conversion
`0vcf2gds.R`  
STAAR requires data in the **GDS (Genomic Data Structure)** format. This script handles the conversion from VCF.
> **Note:** If your data is already in GDS format, you may skip this step.

### Step 1: Variant Standardization
`1Varinfo_gds.R`  
Extracts variant metadata (CHR-POS-REF-ALT) into a CSV. 
* We added a custom algorithm to standardize the format of **REF** and **ALT** alleles and update the **POS** accordingly. 
* This ensures consistency across disparate data sources and prevents mismatches during annotation.

### Step 2: Annotation
`2Annotate.R`  
Annotates the variant information from Step 1 using the **FAVOR** database. The output is saved as a CSV file containing functional scores and categories.

### Step 3: Integrate Annotations (aGDS)
`3gds2agds.R`  
Appends the annotation data back into the original GDS file. This "Annotated GDS" (aGDS) is the primary input for the STAAR statistical tests.

### Step 4: Job Preparation
`4preStep.R`  
Analyzes the total variant count in the aGDS file to calculate the required number of parallel jobs for high-performance computing (HPC) efficiency.

### Step 5: STAAR Association Analysis
`5pipeline_coding_combine.R` | `5pipeline_coding_combine_long.R`  
The core analytical engine of the pipeline.
* **Selection:** Use `_long.R` for high-complexity jobs requiring extended processing time.
* **Configurations:** Users can adjust Minor Allele Frequency (MAF) cutoffs, functional categories, and coding categories.
* **New Feature:** We have added a function to combine different coding categories; refer to `coding_combine.R` for specific logic.
This is the core analytical step. To reproduce the specific results from the paper, you must adjust the script parameters according to the settings below.

| Setting | Raw Data Type | `rare_maf_cutoff` | `variant_type` | Coding Categories (in `coding_combine.R`) |
| :--- | :--- | :--- | :--- | :--- |
| **Setting 1** | Splitted Multi-Allelic | `0.005` | `"variant"` | Original + Non-frameshift INDELs |
| **Setting 2** | Splitted Multi-Allelic | `0.005` | `"SNV"` | Original + Non-frameshift INDELs |
| **Setting 3** | Splitted Multi-Allelic | `0.01` | `"variant"` | Original + Non-frameshift INDELs |
| **Setting 4** | Splitted Multi-Allelic | `0.005` | `"variant"` | Original |
| **Setting 5** | Splitted Multi-Allelic | `0.01` | `"SNV"` | Original |
| **Setting 6** | Merged Multi-Allelic* | `0.01` | `"SNV"` | Original |

> **Note on Setting 6:** This requires the raw data to be processed such that multiple alternative alleles are merged into a single composite alternative allele before running the pipeline.

#### 🛠 How to Modify:
1.  **MAF & Variant Type:** Open `5pipeline_coding_combine.R` and update the `rare_maf_cutoff` and `variant_type` variables at the end of the script.
2.  **Coding Categories:** Open `coding_combine.R` and adjust the logical filters. For example, for **Setting 1**, ensure the logic follows:
 ```r
  coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
  (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")|
  (GENCODE.EXONIC.Category=="nonframeshift deletion")|(GENCODE.EXONIC.Category=="nonframeshift insertion")|
  (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|
  (GENCODE.EXONIC.Category=="nonsynonymous SNV")#|(GENCODE.EXONIC.Category=="synonymous SNV")
```

for **Setting 6**, delete the two categories for non-frameshift INDELS the logic:
 ```r
  coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
  (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")|
  #(GENCODE.EXONIC.Category=="nonframeshift deletion")|(GENCODE.EXONIC.Category=="nonframeshift insertion")|
  (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|
  (GENCODE.EXONIC.Category=="nonsynonymous SNV")#|(GENCODE.EXONIC.Category=="synonymous SNV")
```
3.  **Execution:** Use `5pipeline_coding_combine_long.R` for any jobs that exceed standard cluster walltime limits.
   
### Step 6: Summary and Visualization
`summary_gene_centric_coding.R`  
Generates the final summary of the gene-centric analysis for coding regions. This script automatically produces a **Manhattan plot** of gene-based p-values for visual interpretation.

---

