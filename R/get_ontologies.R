# --------- get_ontology ----------

#' Run analysis for a specified ontology
#'
#' @param gene_id set of unique Entrez Gene identifiers for which GO annotations
#' need to be calculated
#' @param universe set of unique Entrez Gene identifiers
#' @param ontology a character string specifying the ontology, must be either "BP",
#' "CC" or "MF"
#' @param set_label a string to add as a column to the final result. Might
#' be used as identification of the results in case there are many generated
#' using this function
#' @param annotation a string indicating annotation data package name
#' shouldn't be exported or maybe it should?


get_ontology <- function(gene_id, universe, annotation,
                         ontology, set_label){
  params <- new("GOHyperGParams", geneIds = gene_id,
                universeGeneIds = universe,
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
  if (!missing(set_label)) {
    hg_over_results <- hg_over_results %>%
      dplyr::mutate(set_label = set_label)
  }
  return(hg_over_results)
}

# --------- get_all_ontologies ----------

#' Run gene ontology given a list of ENTREZ GENE IDs and backgroun
#'
#' Returns BP, MF and CC results without p value filtering, but with Bonferroni
#' adjustment.
#'
#' @param dat data frame with two columns: gene_id and module. Module provides information
#' about group assignment, gene_id is the list of ENTREZ GENE IDs. Module is useful
#' for parallel computation
#' @param universe vector of unique ENTREZ gene ids to use as a background
#' @param species a character string specifying the species, must be one of two "human"
#' (default) or "mouse". For GO analysis we use 'org.Hs.eg.db' for
#' human and 'org.Mm.eg.db' for mouse.
#' @param ontology


get_all_ontologies <- function(dat, universe, species = c("human", "mouse"),
                               ontology){
  if (species == "human") {
    annotation <- "org.Hs.eg.db"
  } else {
    annotation <- "org.Mm.eg.db"
  }
  genes <- unique(dat$entrez)
  set_label <- unique(dat$set_label)
  res_list <- lapply(ontology, function(x) get_ontology(gene_id = genes,
                                                        universe = universe,
                                                        annotation = annotation,
                                                        ontology = x,
                                                        module = set_label))
  res.df <- do.call("rbind", res_list)
  return(res.df)
}


# --------- run_parallel_go ----------
#' Run GO analysis on parallel
#'
#' This is a parallel implementation of the function to run GO analysis. The
#' function will split the supplied data into pieces based on the unique values
#' in the column "set_label" in the data and will send them to the "cores" for
#' calculations.
#'
#'
#' @param dat data frame with two columns: entrez and set_label. The column
#'  "entrez" should contain Entrez gene identifiers; the column "set_label"
#'  should contain labels for genes that will be analyzed as a group.
#' @param species a character string specifying the species, must be one of two "human"
#'  (default) or "mouse". For GO analysis we use 'org.Hs.eg.db' for
#'  human and 'org.Mm.eg.db' for mouse
#' @param ontology vector of ontologies for which to run the analysis. If not provided
#'  then the analysis will be run for all ontologies: cellular component (CC),
#'  biological process (BP) and molecular function (MF).
#' @param cores numeric value representing the number of cores to use.
#'


run_parallel_go <- function(dat, species = c("human", "mouse"),
                            universe,
                            ontology, cores){
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
  # currently universe is calculated using the full set of genes
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
                            get_all_ontologies(a, universe = universe,
                                               species = species,
                                               ontology = ontologies)
                          }
  return(res)
}

