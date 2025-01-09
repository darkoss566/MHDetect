library(testthat)
library(MHDetect)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)

test_that("MHDetect działa poprawnie", {
  # Przygotuj dane wejściowe
  vcf <- "MHDetect/tests/testthat/EXAMPLE.gz"
  expect_true(file.exists(vcf)) # Upewnij się, że plik istnieje
  vcf_data <- readVcf(vcf)

  # Ustaw parametry
  k <- 25
  N <- 2
  genome <- BSgenome.Hsapiens.UCSC.hg19
  Interval <- 25

  # Uruchom funkcję MHDetect
  result <- MHDetect(vcf_data, k = k, N = N, genome = genome, Interval = Interval)

  # Sprawdź wyniki
  expect_type(result, "list")        # Sprawdź, czy wynik to lista
  expect_true(length(result) > 0)   # Sprawdź, czy lista nie jest pusta
})
