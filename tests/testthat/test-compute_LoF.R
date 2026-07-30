test_that("compute_LoF correctly identifies Loss of Function across all genes", {

  # Create a test table with your variables
  test_genes <- tibble::tibble(
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

  # What exactly we are testing in each row:
  # Row 1: No variants (like in your preview) -> LoF should be 0
  # Row 2: impact is "HIGH" -> LoF_IMPACT = 1, LoF = 1
  # Row 3: TotalCopyNumMin is 0 -> LoF_CopyNum = 1, LoF = 1
  # Row 4: <DEL> starts in the gene, ends outside -> LoF_DEL = 1, LoF = 1
  # Row 5: <DUP> starts in the gene, ends outside -> LoF_DUP_INV_TRA = 1, LoF = 1

  # Execute the function
  res <- compute_LoF(test_genes)

  # Assertions (checking results)
  expect_equal(res$LoF_IMPACT, c(0, 1, 0, 0, 1))
  expect_equal(res$LoF_CopyNum, c(0, 0, 1, 0, 0))
  expect_equal(res$LoF_SV, c(0, 0, 0, 1, 1))
  expect_equal(res$LoF, c(0, 1, 1, 1, 1))
})
