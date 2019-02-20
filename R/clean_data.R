# --------- deduplicate_rows ----------
#' Remove duplicated rows.
#'
#' \code{deduplicate_rows} function removes duplicated rows prior to mapping
#' to ENTREZ gene identifiers.
#'
#' @param dat Data frame object with two columns: \strong{id} and set
#'  labels \strong{set_label}.
#'
#' @return Data frame.
#'
#' @export


deduplicate_rows <- function(dat){
  if (any(duplicated(dat))) {
    message(
      stringr::str_wrap(
        crayon::yellow("ATTENTION: Found duplicated rows in provided data,
                       removing duplicates.")
        )
      )
    dat <- unique(dat)
    message(
      stringr::str_wrap(
        crayon::green(
          paste0("Number of rows after removing duplicated rows is ",
                 crayon::underline(nrow(dat)), ".")
        )
      )
    )
  } else {
    message(
      stringr::str_wrap(
        crayon::green("No duplicate rows have been identified.")
        )
      )
  }
  return(dat)
}

# --------- remove_small_sets ----------
#' Remove gene sets where membership is smaller than a cutoff
#'
#' \code{remove_small_sets} function removes sets that have fewer genes
#' identifiers than \code{min_set_size}.
#'
#' @param dat Data frame object with two columns: ENTREZ gene idenitifiers
#'  \code{entrez} and set labels \code{set_label}.
#'
#' @param min_set_size Positive integer that defines the minimum number of genes
#' in a set. All sets that don't meet this criteria will be removed.
#'
#' @return Data frame object.
#'
#' @export

remove_small_sets <- function(dat, min_set_size){
  if (!missing(min_set_size)) {
    if (min_set_size <= 0){
      stop(
        stringr::str_wrap(
          crayon::red(
            paste0("You provided ", crayon::underline(min_set_size),
                   " as 'min_set_size'. Parameter 'min_set_size' ",
                   "needs to be a positive number.")
          )
        )
      )
    } else {
      message(
        stringr::str_wrap(
          crayon::green(
            paste0("Total number of sets in the dataset is ",
                   crayon::underline(length(unique(dat$set_label))), ".")
          )
        )
      )
      # calculate the number of genes in each set
      gene_count_by_set <- dat %>%
        dplyr::group_by(.data$set_label) %>%
        dplyr::summarise(n = n())
      # extract sets where membership is less than min_set_size
      low_count_set_list <- gene_count_by_set %>%
        dplyr::filter(.data$n <= min_set_size)
      # calculate how many sets need to be removed
      low_count_set_total <- low_count_set_list %>%
        dplyr::summarize(total = n())
      message(
        stringr::wrap(
          crayon::green(
            paste0("You provided ", crayon::underline(min_set_size),
                   " genes as 'min_set_size'. Found ",
                   crayon::underline(as.character(low_count_set_total)),
                   " sets with number of genes less than ",
                   crayon::underline(min_set_size), ".")
          )
        )
      )
      if (as.numeric(low_count_set_total) > 0) {
        dat <- dat[!(dat$set_label %in% low_count_set_list$set_label), ]
        if (nrow(dat) == 0) {
          stop(crayon::red("After removing clusters with fewer than
                           'min.set.size' no rows remained in the data."))
        }
        message(
          stringr::str_wrap(
            crayon::yellow(
              paste0("ATTENTION: Removed sets where the number of genes is
                     less than ", crayon::underline(min_set_size),
                     ". Number of rows after filtering is ",
                     crayon::underline(nrow(dat)),
                     ". Number of sets after filtering is ",
                     crayon::underline(length(unique(dat$set_label))), "."))
              )
              )
        } else {
          message(crayon::green("No sets will be removed."))
      }
    }
    } else {
      message(
        stringr::str_wrap(
          crayon::green("Parameter 'min_set_size' is not provided.
                        No gene sets will be removed.")
          )
        )
    }
  return(dat)
}
