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
  species <- verify_input(input.name = species,
                          input.choices = c("human", "mouse"),
                          input.default = "human")
  # verify 'id' parameter
  id <- verify_input(input.name = id, input.choices = c("hugo", "ensembl"),
                     input.default = "hugo")
  if (species == "mouse") {
    if (id == "hugo") {
      xx <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egALIAS2EG)
      colnames(xx) <- c("entrez", "maping_id")
    } else {
      xx <- AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL)
      colnames(xx) <- c("entrez", "maping_id")
    }
  }
  if (species == "human") {
    if (id == "hugo") {
      xx <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egALIAS2EG)
      colnames(xx) <- c("entrez", "maping_id")
    } else {
      xx <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
      colnames(xx) <- c("entrez", "mapping_id")
    }
  }
  dat_with_genes <- merge(dat, xx, by.x = "id", by.y = "mapping_id")
  dat_with_genes <- dat_with_genes[, c("entrez", "gene_set_label")]
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
