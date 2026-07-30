# MHDetect

**MHDetect** is an R package designed for detecting and classifying structural variants, particularly indels (insertions, deletions, and multi-nucleotide variants, MNVs). The algorithm uses microhomology-mediated end joining (MMEJ) to classify and analyze the indels, which is particularly useful in cancer genomics and studies involving large-scale genomic datasets such as those from The Cancer Genome Atlas (TCGA).

<div align="center">
  <img width="4968" height="4290" alt="image" src="https://github.com/user-attachments/assets/4fd08e72-1c34-48f4-b690-226760d78f88" />
</div>

## Features

- Detects and classifies indels (deletions, insertions, and MNVs) from VCF files.
- Identifies and analyzes microhomology-mediated deletions and insertions.
- Simulates double-strand breaks (DSBs) for MNVs and evaluates their repair mechanism via MMEJ.
- Provides a comprehensive output with results in separate data frames for deletions, insertions, and MNVs.
- Compatible with human genomic references such as `BSgenome.Hsapiens.UCSC.hg19`.

## Algorithm Overview

### Input Data
- **VCF File**: The input to the algorithm is a VCF file, which is read into R using the `readVcf` function from the **VariantAnnotation** package.
- **Parameters**:
  - `k`: Specifies the length (in nucleotides) of the sequence before and after the indel.
  - `N`: Defines the minimum number of matching nucleotides in the sequence after the indel for it to be classified as MMEJ-dependent.
  - `genome`: Indicates the genome used for analysis (e.g., `BSgenome.Hsapiens.UCSC.hg19`).
  - `Interval`: Defines the region before and after the MNV to simulate the Template Switching mechanism.

### Steps of the Algorithm

1. **VCF Data Processing**: 
   - Selects sequences marked as "PASS" in the VCF file and classifies them into deletions, insertions, single-nucleotide variants (SNVs), and multi-nucleotide variants (MNVs).
   - These sequences are saved in separate data frames: `DEL`, `INS`, `SNV`, and `MNV`.

2. **Deletions and Insertions**:
   - For deletions and insertions, the algorithm compares the indel sequences with the reference genome.
   - Sequences before and after the indel are extracted based on the parameter `k` (default set to 25).
   - Indels are classified as microhomology-dependent if the number of matching nucleotides exceeds or equals the value `N`.

3. **MNVs (Multi-Nucleotide Variants)**:
   - The algorithm simulates double-strand breaks (DSBs) within a defined range around the MNV and analyzes the repair process using microhomologies.
   - It evaluates the possible repair outcomes and calculates the positions of deleted and added DNA fragments.
   - The results are compiled in a summary table for further analysis.

### Output

The algorithm outputs comprehensive data for deletions, insertions, and MNVs, including their classification as microhomology-mediated and simulation results for MNVs. The result is a **DataFrame** containing:
- Data on deletions and insertions.
- Data on MNVs and their similarity to simulation results.
- Results of the repair simulations.

## Installation

You can install `MHDetect` directly from GitHub. To do so, you will need the `devtools` package.

### Prerequisites

Before installing the package, ensure that you have R (version 4.0 or higher) and the `devtools` package installed:

```r
install.packages("devtools")
devtools::install_github("darkoss566/MHDetect")
install.packages(c("GenomicRanges", "Biostrings", "BSgenome", "VariantAnnotation"))
# Load the necessary libraries
library(VariantAnnotation)
library(BSgenome)
library(GenomicRanges)
library(MHDetect)
library(BSgenome.Hsapiens.UCSC.hg38)

# Specify the path to your VCF file
vcf_file <- "path/to/your/file.vcf.gz"

# Read the VCF data using the VariantAnnotation package
vcf_data <- readVcf(vcf_file, genome = "hg38")

# Run the MHDetect algorithm
result <- MHDetect(vcf_data, k = 25, N = 2, genome = BSgenome.Hsapiens.UCSC.hg38, Interval = 25)

# View the results
head(result)
```


# compute_LoF()

The compute_LoF() function evaluates gene inactivation through a modular framework, identifying loss-of-function (LoF) events by integrating data across single-nucleotide variants (SNVs), short indels, structural variants (SVs), and copy number variations (CNVs).

<div align="center">
  <img width="3138" height="1936" alt="image" src="https://github.com/user-attachments/assets/8a7988ba-cdac-4e53-bdf4-27ac3bd31252" />
</div>


Genome-Agnostic Input
While the original methodology mapped variants using the GRCh38 reference genome (Ensembl release 115, chromosomes 1-22, X, and Y), the function itself is genome-agnostic. Users can utilize any reference genome or custom gene annotation. The algorithm simply requires a pre-constructed input table (genes × samples × mutations) formatted with the required variables, allowing it to seamlessly process your custom variant mappings.

## How it works
The algorithm evaluates the provided input table and assigns a binary score (1 = inactivated, 0 = retained function) based on the following rules:

- SNVs / Indels (LoF_IMPACT): The gene is marked as LoF if the variant's Ensembl VEP impact score is classified as HIGH or MODERATE.

- Copy Number Variations (LoF_CopyNum): Significant loss is identified when the gene dosage (TotalCopyNumMin) falls strictly below the rounded sample ploidy level.

- Structural Variants (LoF_SV):

- Deletions (<DEL>): Marked as LoF if the deletion overlaps any portion of the gene.

- Other SVs (<DUP>, <INV>, <TRA>): Marked as LoF if they disrupt the structural integrity of the gene (i.e., one breakpoint is located within the gene boundaries and the other externally).

Final Status (LoF): A maximum-score integration is applied. A single positive criterion across any of the data streams is sufficient to classify the gene as at least partially inactivated.

### Example Usage
Below is a practical example demonstrating how compute_LoF() processes different variant scenarios from a user-supplied table:

```r
R
library(dplyr)

# 1. Create an example dataset with the required variables
example_genes <- tibble::tibble(
  gene_id         = c("ENSG_1", "ENSG_2", "ENSG_3", "ENSG_4", "ENSG_5"),
  gene_name       = c("OR4F5", "AL627309", "GENE3", "GENE4", "GENE5"),
  chr             = c("1", "1", "1", "2", "2"),
  start           = c(69000, 134000, 367000, 621000, 738000),
  end             = c(70000, 139000, 368000, 622000, 739000),
  sample_id       = rep("DO46325", 5),
  impact          = c(NA, "HIGH", NA, NA, "MODERATE"),
  TotalCopyNumMin = c(1, 2, 0, 1, 2),
  ploidy          = rep(2.12, 5),
  sv_start        = c(NA, NA, NA, 621500, 738500),
  sv_end          = c(NA, NA, NA, 623000, 740000),
  sv_type         = c(NA, NA, NA, "<DEL>", "<DUP>")
)

# How the algorithm evaluates each row:
# Row 1: No variants -> LoF = 0
# Row 2: SNV impact is "HIGH" -> LoF_IMPACT = 1, Final LoF = 1
# Row 3: TotalCopyNumMin is 0 (below rounded ploidy 2) -> LoF_CopyNum = 1, Final LoF = 1
# Row 4: <DEL> starts in the gene, ends outside -> LoF_SV = 1, Final LoF = 1
# Row 5: <DUP> starts in the gene, ends outside -> LoF_SV = 1, Final LoF = 1

# 2. Execute the function
res <- compute_LoF(example_genes)

# 3. View the calculated LoF scores
res %>% 
  select(gene_name, sample_id, LoF_IMPACT, LoF_CopyNum, LoF_SV, LoF)

#> # A tibble: 5 × 6
#>   gene_name sample_id LoF_IMPACT LoF_CopyNum LoF_SV   LoF
#>   <chr>     <chr>          <dbl>       <dbl>  <dbl> <dbl>
#> 1 OR4F5     DO46325            0           0      0     0
#> 2 AL627309  DO46325            1           0      0     1
#> 3 GENE3     DO46325            0           1      0     1
#> 4 GENE4     DO46325            0           0      1     1
#> 5 GENE5     DO46325            1           0      1     1
```
