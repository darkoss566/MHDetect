# library(testthat)
# library(MHDetect)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(VariantAnnotation)

VCF_test <- system.file("extdata", "EXAMPLE_fixed.vcf.gz", package = "MHDetect")
mh_indels_expected_output <- system.file("extdata", "MH_indels_test_data.tsv", package = "MHDetect")
MNV_MH_dependent_expected_data <- system.file("extdata", "MNV_MH_dependent_test_data.tsv", package = "MHDetect")
MNV_variant_caller_detected_expected_data <- system.file("extdata", "MNV_variant_caller_detected_test_data.tsv", package = "MHDetect")
Template_switch_expected_data <- system.file("extdata", "Template_switch_test_data.tsv", package = "MHDetect")

test_that("MHDetect works", {

  VCF<-readVcf(VCF_test)
  result <- MHDetect(VCF, k = 25, N = 2, genome = BSgenome.Hsapiens.UCSC.hg38, Interval = 25)

  mh_indels_expected <- read.delim(mh_indels_expected_output, header = TRUE, stringsAsFactors = FALSE)
  MNV_MH_dependent_expected <- read.delim(MNV_MH_dependent_expected_data, header = TRUE, stringsAsFactors = FALSE)
  MNV_variant_caller_detected_expected <- read.delim(MNV_variant_caller_detected_expected_data, header = TRUE, stringsAsFactors = FALSE)
  Template_switch_expected <- read.delim(Template_switch_expected_data, header = TRUE, stringsAsFactors = FALSE)

  result$class$seqnames <- as.character(result$class$seqnames)
  result$class$repeat_count <- as.integer(result$class$repeat_count)
  result$class$matching_nucleotides_count <- as.integer(result$class$matching_nucleotides_count)

  mh_indels_expected$seqnames <- as.character(mh_indels_expected$seqnames)
  mh_indels_expected$repeat_count <- as.integer(mh_indels_expected$repeat_count)
  mh_indels_expected$matching_nucleotides_count <- as.integer(mh_indels_expected$matching_nucleotides_count)

  result$MNV_MH_dependent$seqnames <- as.character(result$MNV_MH_dependent$seqnames)
  MNV_MH_dependent_expected$seqnames <- as.character(MNV_MH_dependent_expected$seqnames)
  tolerance <- 1e-6
  row.names(result$MNV_MH_dependent) <- NULL
  row.names(MNV_MH_dependent_expected) <- NULL

  result$MNV_variant_caller_detected$seqnames <- as.character(result$MNV_variant_caller_detected$seqnames)
  MNV_variant_caller_detected_expected$seqnames <- as.character(MNV_variant_caller_detected_expected$seqnames)

  actual_cleaned <- as.data.frame(result$Template_switch)
  attr(actual_cleaned, "groups") <- NULL
  cols_to_convert <- sapply(actual_cleaned, is.double)
  actual_cleaned[ , cols_to_convert] <- lapply(actual_cleaned[ , cols_to_convert], as.integer)
  Template_switch_expected <- as.data.frame(Template_switch_expected)

  expect_identical(result$class, mh_indels_expected)
  expect_identical(result$MNV_MH_dependent, MNV_MH_dependent_expected,tolerance = tolerance)
  expect_identical(result$MNV_variant_caller_detected, MNV_variant_caller_detected_expected)
  expect_identical(actual_cleaned, Template_switch_expected)

})
