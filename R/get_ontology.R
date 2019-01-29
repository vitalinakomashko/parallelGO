#' Run analysis for a specified ontology
#'
#' @param gene_id set of unique Entrez Gene identifiers for which GO annotations
#' need to be calculated
#' @param background set of unique Entrez Gene identifiers
#' @param ontology a character string specifying the ontology, must be either "BP",
#' "CC" or "MF"
#' @param gene_set_label a string to add as a column to the final result. Might
#' be used as identification of the results in case there are many generated
#' using this function
#' @param annotation a string indicating annotation data package name
#' shouldn't be exported or maybe it should?


get_ontology <- function(gene_id, background, annotation,
                         ontology, gene_set_label){
  params <- new("GOHyperGParams", geneIds = gene_id,
                universeGeneIds = background,
                annotation = annotation,
                ontology = ontology,
                pvalueCutoff = 1,
                conditional = FALSE,
                testDirection = "over")
  hg_over_test <- GOstats::hyperGTest(params)
  hg_over_results <- GOstats::summary(hg_over_test)
  hg_over_results <- hg_over_results %>%
    dplyr::mutate(p_bonferroni = p.adjust(.data$Pvalue, method = "bonferroni"),
                            ontology = ontology)
  colnames(hg_over_results)[1] <- "GOID"
  if (!missing(gene_set_label)) {
    hg_over_results <- hg_over_results %>%
      dplyr::mutate(gene_set_label = gene_set_label)
  }
  return(hg_over_results)
}
