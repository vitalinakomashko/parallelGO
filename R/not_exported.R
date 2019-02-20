# ---------verify_input----------
#' Verify validity of the provided function input.
#'
#' `verify_input` provides checks given the name of the input parameter,
#' possible choices and the default choice.  If the `input_name` is
#' not missing and is valid it is then returned; if the `input_name` is not
#' missing and not valid the execution is stopped; if the `input_name` is missing
#' the `input_default` is then returned. This is an internal function.
#'
#' @param input_name String containing the parameter name.
#' @param input_choices String containing the possible choices.
#' @param input_default String containing the default parameter.
#'
#' @seealso \code{\link{map_genes}}, \code{\link{run_parallel_go}}.
#'
#' @return String

verify_input <- function(input_name, input_choices, input_default){
  if (!missing(input_name)) {
    if (length(input_name) == 1) {
      if (input_name %in% input_choices) {
        return(input_name)
      } else {
        stop(
          stringr::str_wrap(
            crayon::red(
              paste0("Unexpected value for ", crayon::underline(input_name),
                     " was provided; please provide one of ",
                     crayon::underline(paste0(input_choices,
                                              collapse = " or ")),
                     ".")
            )
          )
        )
      }
    } else {
      stop(
        stringr::str_wrap(
          crayon::red(
            paste0("More than one value for ", crayon::underline(input_name),
                   " parameter was provided; please provide either ",
                   crayon::underline(paste0(input_choices, collapse = " or ")),
                   ".")
          )
        )
      )
    }
  } else {
    message(
      paste0("Using default setting for the parameter ",
             input_name, ": ",
             input_default)
    )
    return(input_default)
  }
}


# ---------return_error_result----------
#' Output a data frame for sets that generated an error or a warning
#'
#' `return_error_result` creates a one-row data frame with the same columns
#' as the regular output of get_ontology. The column \strong{count} is
#' represented by the number of genes, the column \strong{term} provides the
#' message, and the column \strong{set_label} provides the label of the gene set
#' for which a condition was encountered. This is an internal function.
#'
#' @seealso \code{\link{run_parallel_go}}, \code{\link{remove_errors}},
#'  \code{\link{get_ontology}}.
#'
#' @return Data frame with 1 row.

return_error_result <- function(error_warning, input_genes, input_label){
  df <- data.frame(goid = NA,
                   pvalue = NA,
                   odds_ratio = NA,
                   exp_count = NA,
                   count = length(input_genes),
                   size = NA,
                   term = conditionMessage(error_warning),
                   p_bonferroni = NA,
                   ontology = "error_warning",
                   set_label = input_label,
                   stringsAsFactors = FALSE)
  return(df)
}


# ---------remove_errors----------
#' Remove sets which generated errors or warnings by `get_all_ontologies`.
#'
#' `remove_errors` removes rows in the output of `get_all_ontologies` for sets
#' which generated errors or warnings. Also provides a helpful message listing
#' the labels of the removed sets. This is an internal function.
#'
#' @param dat Data frame, output of the parallelized run.
#'
#' @seealso \code{\link{run_parallel_go}}.
#'
#' @return Data frame.

remove_errors <- function(dat){
  if("error_warning" %in% dat$ontology) {
    k <- which(dat$ontology == "error_warning")
    message(
      stringr::str_wrap(
        crayon::yellow(
          paste0("GO enrichment generated errors or warnings for the following ",
                 "sets: ",
                 crayon::underline(paste0(as.character(dat$set_label[k]),
                                          collapse = ", ")),
                 ". These sets are excluded from the output.")
        )
      )
    )
    dat <- dat[-k, ]
  }
}
