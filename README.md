# MHDetect

**MHDetect** is an R package designed for detecting and classifying structural variants, particularly indels (insertions, deletions, and multi-nucleotide variants, MNVs). The algorithm uses microhomology-mediated DNA repair (MMEJ) to classify and analyze the indels, which is particularly useful in cancer genomics and studies involving large-scale genomic datasets such as those from The Cancer Genome Atlas (TCGA).

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
library(MHDetect)
library(BSgenome.Hsapiens.UCSC.hg19)

# Specify the path to your VCF file
vcf_file <- "path/to/your/file.vcf.gz"

# Read the VCF data using the VariantAnnotation package
vcf_data <- readVcf(vcf_file, genome = "hg19")

# Run the MHDetect algorithm
result <- MHDetect(vcf_data, k = 25, N = 2, genome = BSgenome.Hsapiens.UCSC.hg19, Interval = 25)

# View the results
head(result)
```

