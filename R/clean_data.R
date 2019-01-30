#' Basic data cleaning before gene mapping
#'
#' `clean_data` function removes duplicated rows, sets that have fewer genes
#' than `min_set_size`.
#'
#' @param dat Data frame object with two columns: genes (`id`) and set labels (`set.name`),
#'  output of the function \code{\link{read_data}} call.
#' @param min_set_size Positive integer. If provided will be used to remove
#'  all gene sets where then number of genes is fewer than \code{min.set.size}.
#'
#' @return Data frame object
#' @export


clean_data <- function(dat, min_set_size){
  # check the data for duplicated rows
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
        paste0("Number of rows after removing duplicated rows is ",
               nrow(dat), ".")
        )
      )
  }
  # remove modules
  if (!missing(min_set_size)) {
    if (min_set_size <= 0){
      stop(
        stringr::str_wrap(
          paste0("You provided ", min_set_size,
                 " as 'min_set_size'. Parameter 'min_set_size' should be a
                 positive number.")))
    } else {
      message(
        stringr::str_wrap(
          paste0("Total number of sets in the dataset is ",
                 length(unique(dat$set.)), ".")
          )
        )
      gene_count_by_set <- dat %>%
        dplyr::group_by(.data$set_label) %>%
        dplyr::summarise(n = n())
      low_count_set_list <- gene_count_by_set %>%
        dplyr::filter(.data$n <= min_set_size)
      low_count_set_total <- low_count_set_list %>%
        dplyr::summarize(total = n())
      message(
        paste0("You provided ", min_set_size,
               " genes as 'min_set_size'. \nFound ",
               as.character(low_count_set_total),
               " sets with number of genes less than ", min_set_size, "."))
      if (as.numeric(low_count_set_total) > 0) {
        dat <- dat[!(dat$set_label %in% low_count_set_list$set_label), ]
        if (nrow(dat) == 0) {
          stop("After removing clusters with fewer than 'min.set.size' no
               rows remained in the data.")
        }
        message(
          stringr::str_wrap(
            crayon::yellow(
              paste0("ATTENTION: Removed sets where the number of genes is
                      less than ", min_set_size, ". Number of rows after
                      filtering is ", nrow(dat), ". Number of sets after
                      filtering is ", length(unique(dat$set_label)), "."))
            )
          )
      } else {
        message("No sets will be removed.")
      }
    }
  } else {
    message(
      stringr::str_wrap("Parameter 'min_set_size' is not provided.
                        No gene sets will be removed."))
  }
  return(dat)
}
