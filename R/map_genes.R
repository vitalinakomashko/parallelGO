#' Map gene identifiers to ENTREZ gene
#'
#' \code{map_genes} provides mapping to ENTREZ gene identifiers for Hugo or Ensembl
#' gene identifiers using either \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db} or \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db} Bioconductor annotation
#' packages. Function \code{\link[AnnotationDb]{toTable}} is used to
#' extract the mappings. We use org.Xx.eg.db::org.Xx.egALIAS2EG for the mappings
#' between the common gene symbol idenditifers and org.Xx.eg.db::org.Xx.egENSEMBL for
#' the mappings between ENTREZ Gene identifiers and Exnsembl gene
#' accession numbers.
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
    if (id == "symbol") {
      xx <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egALIAS2EG)
    } else {
      xx <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL)
    }
  }
  if (species == "human") {
    if (id == "symbol") {
      xx <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egALIAS2EG)
    } else {
      xx <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
    }
  }
  # check if ensembl id has versions in it and remove them
  if (id == "ensembl") {
    dat$id <- stringr::str_remove(dat$id, "\\..{1,}")
  }
  colnames(xx) <- c("entrez", "mapping_id")
  # since there could be multiple matches, do alphabetic sort
  xx <- xx[order(xx$mapping_id), ]
  # match ids in the data to the mapping_id in the annotation data
  # frame extracting only the first match
  k <- match(dat$id, xx$mapping_id)
  dat_with_match <- dplyr::mutate(dat, index = k)
  # remove rows without any match if such are available
  if(any(is.na(k))){
    dat_with_match <- dat_with_match[!is.na(dat_with_match$index), ]
  }
  # add the matched entrez ID and select relevant columns
  dat_with_genes <- dplyr::bind_cols(dat_with_match,
                                     xx[dat_with_match$index, ]) %>%
    dplyr::select(.data$entrez, .data$set_label)
  if (nrow(dat_with_genes) == 0) {
    stop(
      stringr::str_wrap(
        crayon::red(
          paste0("ERROR: The number of rows after mapping identifiers to ENTREZ",
                 " gene identifiers is 0. Please verify that you used correct",
                 " input values for 'species' and 'id' parameters. You used ",
                 crayon::underline(species),
                 " for the 'species' parameter and ",
                 crayon::underline(id), " for the 'id' parameter.")
          )
        )
      )
  } else if (nrow(dat_with_genes) < nrow(dat) / 2) {
    warning(
      stringr::str_wrap(
        crayon::yellow(
          paste0("ATTENTION: The number of rows after mapping identifiers
                  to ENTREZ gene identifiers is significantly smaller than
                  in the original data provided. Please verify that you used
                  correct input values for 'species' or 'id' parameters.
                  You used ", crayon::underline(species), " for the 'species'
                  parameter and ", crayon::underline(id),
                 " for the 'id' parameter.")
          )
        )
      )
    return(dat_with_genes)
  } else {
    message(
      stringr::str_wrap(
        crayon::green(
          paste0("SUCCESS: Performed mapping to ENTREZ gene identifiers. ",
                 crayon::underline(nrow(dat_with_genes)),
                 " out of the original ",
                 crayon::underline(nrow(dat)), " identifiers were mapped.")
        )
      )
    )
    if (any(duplicated(dat_with_genes))) {
      dat_with_genes <- unique(dat_with_genes)
      message(
        stringr::str_wrap(
          crayon::yellow(
            paste0("ATTENTION: Removed duplicated rows after mapping to ENTREZ",
                   " gene identifiers. ",
                   crayon::underline(nrow(dat_with_genes)),
                   " rows remained.")
          )
        )
      )
    }
    return(dat_with_genes)
  }
}
