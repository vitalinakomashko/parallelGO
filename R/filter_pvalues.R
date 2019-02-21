# --------- filter_pvalues ----------

#' Filter results based on a p-value cutoff

filter_pvalues <- function(dat, which_pvalue = c("raw", "adjusted"), cutoff){
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
    }
    pval_type <- verify_input(input_name = which_pvalue,
                              input_choices = c("raw", "adjusted"),
                              input_default = "adjusted")
  }
}
