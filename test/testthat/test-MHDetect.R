# tests/testthat/test-MHDetect.R
library(testthat)
library(MHDetect)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(stringr)

test_that("MHDetect działa poprawnie", {
  # Przygotuj dane wejściowe (przykładowy plik VCF)
  vcf <- "C:\\Users\\daria\\OneDrive - Politechnika Śląska\\MH\\MHDetect\\test\\testthat\\VCF.gz"
  vcf_data <- readVcf(vcf)

  # Ustaw parametry
  k <- 25
  N <- 2
  genome <- BSgenome.Hsapiens.UCSC.hg19
  Interval <- 25

  # Uruchom funkcję MHDetect
  result <- MHDetect(vcf_data, k = k, N = N, genome = genome, Interval = Interval)


})
