# --------- read_data ----------

#' Read data with genes and gene set labels to a data frame
#'
#' `read_file` reads a file provided a path, column names and file field
#' delimiter and returns a data frame. For the file to be processed correctly
#' ensure that the gene identifiers are in the first columns and the set labels
#' are in the second column in the file.
#'
#' @param path The name of the fle which the data are to be read from.
#' @param col_names A logical value indicating whether the file contains the names
#'  of the variables as its first line.
#' @param delim The file field separator character.
#' @return A data frame with two columns, where the first column is 'id' (gene
#'  identifiers) and the second column is 'set.label' (gene set labels).


read_file <- function(path, col_names, delim){
  if (missing(path)) {
    stop("File path must be provided")
  } else {
    if (!file.exists(path)) {
      stop("The file doesn't exist, please provide path to an existing file")
    } else {
      message(
        crayon::yellow("ATTENTION: For the file to be processed correctly,",
                         "ensure that the gene identifiers are in the first",
                         "column and the gene set labels are in the second column.")
        )
      if (missing(col_names)) {
        stop("Please provide TRUE (present) or FALSE (absent) for the
             'col_names' parameter.")
      } else {
        col_names <- col_names
      }
      if (missing(delim)) {
        stop("Please provide a value for 'delim' parameter.")
      } else {
        delim <- delim
      }
      dat <- read.delim(file = path, header = col_names, sep = delim,
                        stringsAsFactors = FALSE,
                        colClasses = "character",
                        row.names = NULL)
      message(
        stringr::str_wrap(
          crayon::green(
            paste0("SUCCESS: Read the file '", path, "' with ",
                   nrow(dat), " and ", ncol(dat), " columns.\n")
          )
        )
      )
      if (ncol(dat) > 2) {
        message(
          stringr::str_wrap(
            crayon::yellow("ATTENTION: found more than two columns in the file,
                           retaining only the first two.")))
        dat <- dat[, 1:2]
      }
      colnames(dat) <- c("id", "set_label")
      return(dat)
    }
      }
  }

