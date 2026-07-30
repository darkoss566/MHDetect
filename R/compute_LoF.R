#' Compute Ploidy
#'
#' @param ploidy Ploidy value
#' @export
rnd <- function(ploidy) {
  ifelse(ploidy <= 2, 2, floor(ploidy))
}

#' Compute Loss of Function gene
#'
#' @param df Data frame containing variants
#' @return A data frame with appended LoF columns
#' @export
#' @importFrom dplyr %>%
compute_LoF <- function(df) {
  df %>%
    dplyr::mutate(
      LoF_IMPACT = ifelse(impact %in% c("HIGH", "MODERATE"), 1, 0)
    ) %>%
    dplyr::mutate(
      LoF_CopyNum = ifelse(TotalCopyNumMin == 0, 1, 0)
    ) %>%
    dplyr::mutate(
      is_DEL = sv_type == "<DEL>",
      start_in_gene = sv_start >= start & sv_start <= end,
      end_in_gene   = sv_end   >= start & sv_end   <= end,
      LoF_DEL = ifelse(is_DEL & (start_in_gene | end_in_gene), 1, 0)
    ) %>%
    dplyr::mutate(
      is_DUP_INV_TRA = sv_type %in% c("<DUP>", "<INV>", "<TRA>"),
      LoF_DUP_INV_TRA = ifelse(is_DUP_INV_TRA & (start_in_gene != end_in_gene), 1, 0)
    ) %>%
    dplyr::mutate(
      LoF_DUP_INV_TRA = ifelse(is.na(LoF_DUP_INV_TRA), 0, LoF_DUP_INV_TRA),
      LoF_DEL = ifelse(is.na(LoF_DEL), 0, LoF_DEL),
      LoF_CopyNum = ifelse(is.na(LoF_CopyNum), 0, LoF_CopyNum),
      LoF_SV = pmax(LoF_DEL, LoF_DUP_INV_TRA, na.rm = TRUE),
      LoF_IMPACT = ifelse(is.na(LoF_IMPACT), 0, LoF_IMPACT),
      LoF = pmax(LoF_IMPACT, LoF_CopyNum, LoF_DEL, LoF_DUP_INV_TRA)
    ) %>%
    dplyr::select(-is_DEL, -is_DUP_INV_TRA, -LoF_DUP_INV_TRA,
                  -start_in_gene, -end_in_gene, -LoF_DEL)
}
