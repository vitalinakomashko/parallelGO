#' Map gene identifiers to Entrez Gene identifiers.
#'
#' @param dat Data frame with two columns: `id` (genes) and `set.label` (gene
#'  sets labels). The data frame is supposed to be returned by the call to the
#'  function \code{\link{clean_data}}.
#' @param id A character string specifying the type of gene identifier in the
#'  data, must of be one of two "hugo" (default) or "ensembl".
#' @param species A character string specifying the species, must be one of two
#'  "human" (default) or "mouse". For mapping to Entrez Gene identifiers we use
#'  'org.Hs.eg.db' (human) and 'org.Mm.eg.db' (mouse) Bioconductor annotation
#'  packages.
#' @return Data frame with an additional column `entrez` (column with the
#' original identifiers will be removed).
#' @export

map_genes <- function(dat, id = "hugo", species = "human"){
  # verify 'species' parameter
  species <- verify_input(input_name = species,
                          input_choices = c("human", "mouse"),
                          input_default = "human")
  # verify 'id' parameter
  id <- verify_input(input_name = id, input_choices = c("hugo", "ensembl"),
                     input_default = "hugo")
  if (species == "mouse") {
    if (id == "hugo") {
      xx <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egALIAS2EG)
    } else {
      xx <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL)
    }
  }
  if (species == "human") {
    if (id == "hugo") {
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
  dat_with_genes <- dplyr::bind_cols(dat_with_match, xx[dat_with_match$index, ]) %>%
    dplyr::select(.data$entrez, .data$set_label)
  if (nrow(dat_with_genes) == 0) {
    stop(
      stringr::str_wrap(
        crayon::red(
          paste0("The number of rows after mapping identifiers to Entrez Gene
                identifiers is 0. Please verify that you used correct input
                values for 'species' and 'id' parameters. You used ",
                crayon::underline(species), " for the 'species' parameter and ",
                crayon::underline(id), " for the 'id' parameter.")
          )
        )
      )
  } else if (nrow(dat_with_genes) < nrow(dat) / 2) {
    warning(
      stringr::str_wrap(
        crayon::yellow(
          paste0("ATTENTION: The number of rows after mapping identifiers
                  to Entrez Gene identifiers is significantly smaller than
                  in the original data provided. Please verify that you used
                  correct input values for 'species' or 'id' parameters.
                  You used ", crayon::underline(species), " for the 'species'
                  parameter and ", crayon::underline(id),
                 " for the 'id' parameter.")
          )
        )
      )
  } else {
    message(
      stringr::str_wrap(
        paste0("Performed mapping to Entrez Gene identifiers. ",
               nrow(dat_with_genes), " out of the original ",
             nrow(dat), " identifiers were mapped.")
        )
      )
    if (any(duplicated(dat_with_genes))) {
      dat_with_genes <- unique(dat_with_genes)
      message(
        stringr::str_wrap(
          paste0("Removed duplicated rows after mapping to Entrez Gene
                 identifiers. ",
                 nrow(dat_with_genes),
                 " rows remained.")
          )
        )
    }
    return(dat_with_genes)
  }
}
