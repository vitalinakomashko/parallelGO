#' Run analysis for a specified ontology
#' 
#' @param gene.id set of unique Entrez Gene identifiers for which GO annotations
#' need to be calculated
#' @param background set of unique Entrez Gene identifiers
#' @param ontology a character string specifying the ontology, must be either "BP",
#' "CC" or "MF" 
#' @param module a string to add as a column to the final result. Might
#' be used as identification of the results in case there are many generated
#' using this function
#' @param annotation a string indicating annotation data package name
#' shouldn't be exported or maybe it should?


get_ontology <- function(gene.id, background, annotation, ontology, module){
  params <- new("GOHyperGParams", geneIds = gene.id,
                universeGeneIds = background, 
                annotation = annotation,
                ontology = ontology, 
                pvalueCutoff = 1, 
                conditional = FALSE, 
                testDirection = "over")
  hgOver <- GOstats::hyperGTest(params)
  hgOver.results <- GOstats::summary(hgOver)
  hgOver.results <- hgOver.results %>% 
    dplyr::mutate(P.Bonferroni = p.adjust(.data$Pvalue, method="bonferroni"),
                            Ontology = ontology)
  colnames(hgOver.results)[1] <- "GOID"
  if(!missing(module)){
    hgOver.results <- hgOver.results %>% dplyr::mutate(module = module)
  }
  return(hgOver.results)
}
