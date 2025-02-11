# library(testthat)
# library(MHDetect)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(VariantAnnotation)

VCF_test <- system.file("extdata", "EXAMPLE_fixed.vcf.gz", package = "MHDetect")
mh_indels_expected_output <- system.file("extdata", "MH_indels_test_data.rds", package = "MHDetect")
MNV_MH_dependent_expected_data <- system.file("extdata", "MNV_MH_dependent_test_data.rds", package = "MHDetect")
MNV_variant_caller_detected_expected_data <- system.file("extdata", "MNV_variant_caller_detected_test_data.rds", package = "MHDetect")
Template_switch_expected_data <- system.file("extdata", "Template_switch_test_data.rds", package = "MHDetect")

test_that("MHDetect works", {

  VCF<-readVcf(VCF_test)
  result <- MHDetect(VCF, k = 25, N = 2, genome = BSgenome.Hsapiens.UCSC.hg38, Interval = 25)

  mh_indels_expected <- readRDS(mh_indels_expected_output)
  MNV_MH_dependent_expected <- readRDS(MNV_MH_dependent_expected_data)
  MNV_variant_caller_detected_expected <- readRDS(MNV_variant_caller_detected_expected_data)
  Template_switch_expected <- readRDS(Template_switch_expected_data)

  expect_identical(result$class, mh_indels_expected)
  expect_identical(result$MNV_MH_dependent, MNV_MH_dependent_expected)
  expect_identical(result$MNV_variant_caller_detected, MNV_variant_caller_detected_expected)
  expect_identical(result$Template_switch, Template_switch_expected)

})
