#' Run gene ontology given a list of ENTREZ GENE IDs and backgroun
#'
#' Returns BP, MF and CC results without p value filtering
#' @param dat data frame with 2 columns: gene_id and module. Module provides information
#' about group assignment, gene_id is the list of ENTREZ GENE IDs. Module is useful
#' for parallel computation
#' @param background vector of unique ENTREZ gene ids to use as a background
#' @param species a character string specifying the species, must be one of two "human"
#' (default) or "mouse". For GO analysis we use 'org.Hs.eg.db' for
#' human and 'org.Mm.eg.db' for mouse.


get_all_ontologies <- function(dat, background, species = c("human", "mouse"),
                               ontology){
  if (species == "human") {
    annotation <- "org.Hs.eg.db"
  } else {
    annotation <- "org.Mm.eg.db"
  }
  genes <- unique(dat$gene_id)
  module <- unique(dat$module)
  res.list <- lapply(ontology, function(x) get_ontology(gene.id = genes,
                                                        background = background,
                                                        annotation = annotation,
                                                        ontology = x,
                                                        module = module))
  res.df <- do.call("rbind", res.list)
  return(res.df)
}
