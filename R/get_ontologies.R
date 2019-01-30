# --------- get_ontology ----------

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

# --------- get_all_ontologies ----------

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


# --------- run_parallel_go ----------
#' Run GO analysis on parallel
#'
#' @param annotation string providing the name of the annotation package to use,
#'  Can we use anything else but the 'org.Hs.eg.db' forhuman and 'org.Mm.eg.db' for
#'  mouse? I haven't investigated
#'  @param ontology vector of ontologies for which to run the analysis. If not provided
#' then the analysis will be run for all ontologies: CC, BP and MF


run_parallel_go <- function(dat, background, species = c("human", "mouse"),
                            cores, ontology){
  if (!missing(cores)) {
    doParallel::registerDoParallel(cores = cores)
    workers <- doParallel::getDoParWorkers()
    message(
      stringr::str_wrap(
        paste0("Will run GO computation on ", workers, ".")
      )
    )
  } else {
    doParallel::registerDoParallel()
    workers <- doParallel::getDoParWorkers()
    message(
      stringr::str_wrap(
        paste0("Will run GO computation on ", workers, ".")
      )
    )
  }
  species <- verify_input(input_name = species,
                          input_choices = c("human", "mouse"),
                          input_default = "human")
  if (missing(background)) {
    stop("Please provide a vector with Entrez Gene identifiers to serve
         as a background")
  } else {
    background <- unique(background)
  }
  if (missing(ontology)) {
    ontologies <- c("BP", "CC", "MF")
  } else {
    if (!all(ontology %in% c("BP", "CC", "MF"))) {
      stop(
        stringr::str_wrap(
          paste0("Please provide valid values for the ontology parameter.
                 Possible valid values: CC, BP, MF. You provided: ",
                 crayon::underline(ontology), ".")
        )
        )
    } else {
      ontologies <- ontology
    }
  }
  iterated_df <- iterators::isplit(dat, dat$module)
  res <- foreach::foreach(a = iterated_df,
                          .combine = rbind,
                          .packages = c("GOstats", "dplyr"),
                          .verbose = TRUE) %dopar% {
                            get_all_ontologies(a, background = background,
                                               species = species,
                                               ontology = ontologies)
                          }
  return(res)
}

