#' Map gene identifiers to ENTREZ gene
#'
#' \code{map_genes} provides mapping to ENTREZ gene identifiers for Hugo or
#' Ensembl gene identifiers using either \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db} or \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db} Bioconductor annotation
#' packages. Function \code{\link[annotate:getSYMBOL]{annotate::lookUp()}}
#' extract the mappings. We use \code{ALIAS2EG} for Hugo - ENTREZ
#' mappings and \code{ENSEMBL} for Ensembl - ENTREZ gene mappings.
#'
#' @param dat Data frame with two columns: \strong{id} and \strong{set_label}.
#'
#' @param id A character string specifying the type of gene identifier in the
#'  data, must of be one of two "symbol" (default) or "ensembl".
#'
#' @param species A character string specifying the species, must be one of two
#'  "human" (default) or "mouse".
#'
#' @return Data frame with two columns: \strong{entrez} and \strong{set_label}.
#' The column with the original identifiers is removed.
#'
#' @examples
#'
#' \dontrun{
#' # load test data:
#' data("human_symbol")
#' # remove duplicate rows:
#' dat_clean <- deduplicate_rows(human_symbol)
#' # map symbol to ENTREZ gene id:
#' dat_mapped <- map_genes(dat_clean, id = "symbol", species = "human")
#' }
#'
#' @export

map_genes <- function(dat, id = "symbol", species = "human"){
  # verify 'species' parameter
  species <- verify_input(input_name = species,
                          input_choices = c("human", "mouse"),
                          input_default = "human")
  # verify 'id' parameter
  id <- verify_input(input_name = id, input_choices = c("symbol", "ensembl"),
                     input_default = "symbol")
  if (species == "mouse") {
    metadata_basename <- "org.Mm.eg"
  } else {
    metadata_basename <- "org.Hs.eg"
  }
  if (id == "symbol") {
    annotation_element <- "ALIAS2EG"
  } else {
    annotation_element <- "ENSEMBL"
    # also, remove versions in ENSEMBL ids in the input dat
    dat$id <- stringr::str_remove(dat$id, "\\..{1,}")
  }
  # this returns a list, where the name is the symbol/ensembl and the value is
  # the character vector of all ENTREZ ids this symbol maps to.
  xx <- annotate::lookUp(dat$id, metadata_basename, annotation_element,
                         load = TRUE)
  # each symbol/ensembl might map to more than one ENTREZ id.
  map_length <- sapply(xx, length)
  map_count <- data.frame(matches = unname(map_length), id = names(map_length))
  # summarize
  map_summary <- map_count %>%
    dplyr::group_by(.data$matches) %>%
    dplyr::summarise(count = n())
  message(
    stringr::str_wrap(
      crayon::green(
        paste0("Summary of the mapping results from the provided ",
               crayon::underline(id), " to ENTREZ gene identifiers:")
      )
    )
  )
  print(map_summary)
  message(
    stringr::str_wrap(
      crayon::green(
        paste0("If more than 1 matches are found to ENTREZ ids, only the ",
               "first match is used.")
      )
    )
  )
  # extract the first matches only
  first_match_entrez <- sapply(xx, "[", 1)
  # add ENTREZ ids to dat
  dat_with_match <- dplyr::mutate(dat, entrez = unname(first_match_entrez))
  # drop "id" column
  dat_with_match <- dat_with_match[, c("entrez", "set_label")]
  # remove rows with missing ENTREZ (where ENTREZ wasn't found)
  if(any(is.na(dat_with_match$entrez))) {
    k <- which(is.na(dat_with_match$entrez))
    dat_with_match <- dat_with_match[!is.na(dat_with_match$entrez), ]
    if (nrow(dat_with_match) == 0) {
      stop(
        stringr::str_wrap(
          crayon::red(
            paste0("ERROR: The number of rows after mapping identifiers to ",
                   "ENTREZ gene identifiers is 0. Please verify that you used ",
                   "correct input values for 'species' and 'id' parameters. ",
                   "You used ", crayon::underline(species),
                   " for the 'species' parameter and ",
                   crayon::underline(id), " for the 'id' parameter.")
          )
        )
      )
    } else if (nrow(dat_with_match) < nrow(dat) / 2) {
      warning(
        stringr::str_wrap(
          crayon::yellow(
            paste0("ATTENTION: The number of rows after mapping identifiers ",
                   "to ENTREZ gene identifiers is less than a half of the ",
                   " original identifiers provided.  Please verify that you ",
                   "used correct input values for 'species' or 'id' ",
                   "parameters. You used ", crayon::underline(species),
                   " for the 'species' parameter and ", crayon::underline(id),
                   " for the 'id' parameter.")
          )
        )
      )
    } else {
      message(
        stringr::str_wrap(
          crayon::green(
            paste0("SUCCESS: Performed mapping to ENTREZ gene identifiers. ",
                   crayon::underline(nrow(dat_with_match)),
                   " out of the original ",
                   crayon::underline(nrow(dat)), " identifiers were mapped.")
          )
        )
      )
    }
  } else {
    message(
      stringr::str_wrap(
        crayon::green(
          paste0("SUCCESS: Performed mapping to ENTREZ gene identifiers. ",
                 "All identifiers have been mapped")
        )
      )
    )
  }
  # remove duplicates
  if (any(duplicated(dat_with_match))) {
    dat_with_match <- unique(dat_with_match)
    message(
      stringr::str_wrap(
        crayon::yellow(
          paste0("ATTENTION: Removed duplicated rows after mapping to ENTREZ",
                 " gene identifiers. ",
                 crayon::underline(nrow(dat_with_match)),
                 " rows remained.")
        )
      )
    )
  }
  return(dat_with_match)
}


