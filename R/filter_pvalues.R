# --------- filter_pvalues ----------

#' Filter GO enrichment results using a p-value cutoff
#'
#' \code{filter_pvalues} takes a data frame with the results of GO enrichment
#' analysis and performs row filtering based on provided \code{cutoff} value.
#' The data can be filtered either using the raw (\strong{pvalue} column) or the
#' adjusted (default) p-values (\strong{p_bonferroni} column).
#'
#' @param dat Data frame with the results of GO enrichment analysis.
#'
#' @param which_pvalue A character string specifying p-value for filtering.
#'  Must of one of two: "raw" or "adjusted"
#'
#' @param cutoff A numeric value specifying the cutoff. All rows with p-values
#'  that are less or equal to the cutoff will be retained.
#'
#' @return Data frame.
#'
#' @seealso \code{get_ontology}.
#'
#' @export

filter_pvalues <- function(dat, which_pvalue = "adjusted", cutoff){
  if (missing(cutoff)) {
    message(
      crayon::green("No filtering on p-values will be performed.")
    )
    return(dat)
  } else {
    if (cutoff > 1 | cutoff <= 0) {
      stop(
        crayon::red(
          paste0("ERROR: Invalid value for the p-value `cutoff` ",
                 "provided. The value should be more than zero and ",
                 "less than 1.")
        )
      )
    } else {
      pval_type <- verify_input(input_name = which_pvalue,
                                input_choices = c("raw", "adjusted"),
                                input_default = "adjusted")
      column_filter <- switch(pval_type, "raw" = "pvalue",
                              "adjusted" = "p_bonferroni")
      dat_filtered <- dat[dat[[column_filter]] <= cutoff, ]
      if (nrow(dat_filtered) == 0){
        message(
          crayon::yellow(
            paste0("ATTENTION: Using p-value cutoff on ",
                   crayon::underline(column_filter), " column and cutoff ",
                   crayon::underline(cutoff), " resulted in a data frame with ",
                   "0 rows. Returning original data frame with unfiltered ",
                   "p-values."))
        )
        return(dat)
      } else {
        return(dat_filtered)
      }
    }
  }
}
