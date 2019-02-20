# --------- get_ontology ----------

#' Run GO enrichment analysis for a specified ontology.
#'
#' \code{get_ontology} uses \code{\link[GOstats]{hyperGTest}} to calculate
#' hypergeometric statistics (for over-represenation of each GO term);
#' \code{\link[GOstats]{summary}} to return the summarized results with adjusted
#' p values using \code{\link[stats]{p.adjust}} function and method
#' "bonferroni".
#'
#' @param gene_id Character vector of unique ENTREZ gene identifiers of
#'  interest.
#'
#' @param universe Character vector of unique ENTREZ gene identifiers to be used
#'  as a universe.
#'
#' @param ontology Character string specifying the ontology, must be either
#'  "BP", "CC" or "MF".
#'
#' @param set_label Character string used to populate the column
#'  \strong{set_label}. This parameter is optional. It is used for
#'  identification of results when \code{get_ontology} is used in a loop.
#'
#' @param annotation Character string with the name of the annotation data
#'  package.
#'
#' @return Data frame with the following columns: \strong{goid},
#'  \strong{pvalue}, \strong{odds_ratio}, \strong{exp_count}, \strong{count},
#'  \strong{size}, \strong{term}, \strong{p_bonferroni}, \strong{ontology},
#'  \strong{set_label} (if provided).
#'
#' @seealso \code{\link[GOstats]{hyperGTest}},
#' \code{\link[GOstats]{summary,GOHyperGResult-method}}.
#'
#' @export


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
    dplyr::mutate(p_bonferroni = stats::p.adjust(.data$Pvalue, method = "bonferroni"),
                  ontology = ontology)
  # rename the columns
  colnames(hg_over_results)[1:7] <- c("goid", "pvalue", "odds_ratio",
                                      "exp_count", "count", "size", "term")
  if (!missing(set_label)) {
    hg_over_results <- hg_over_results %>%
      dplyr::mutate(set_label = set_label)
  }
  return(hg_over_results)
}

# --------- get_all_ontologies ----------

#' Perform GO enrichment analysis for a list of ontologies.
#'
#' \code{get_all_ontologies} loops through the list of ontologies and calls
#' \code{get_ontology} for each.
#'
#' @param gene_id Character vector of unique ENTREZ gene identifiers of
#'  interest.
#'
#' @param universe Character vector of unique ENTREZ gene identifiers to be used
#'  as a universe.
#'
#' @param species Character string specifying the species, must be one of two
#'  "human" (default) or "mouse". This will define which Bioconductor package
#'  is used as a source of GO annotations. For human we use
#'  \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db}, and for mouse \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db}.
#'
#' @param ontologies Character vector specifying ontologies of interest, such
#'  as "MF", "BP" and "CC".
#'
#' @param set_label Character string specifying the label characterizing the
#'  set of genes being analyzed. This parameter is optional.
#'
#' @return Data frame with the results of GO enrichment and an optional column
#'  populated with the provided "set_label".
#'
#' @seealso \code{get_ontology} to understand the output format.
#'
#' @export

get_all_ontologies <- function(gene_id, universe, species = c("human", "mouse"),
                               ontologies, set_label){
  if (species == "human") {
    annotation <- "org.Hs.eg.db"
  } else {
    annotation <- "org.Mm.eg.db"
  }
  res_list <- base::lapply(ontologies, function(x) get_ontology(gene_id = gene_id,
                                                        universe = universe,
                                                        annotation = annotation,
                                                        ontology = x,
                                                        set_label = set_label))
  res.df <- base::do.call("rbind", res_list)
  return(res.df)
}


# --------- run_parallel_go ----------
#' Run GO analysis on parallel.
#'
#' \code{run_parallel_go} splits provided data frame into a list of data frames
#' based on the values in the column \strong{set_label} and
#' using \code{\link[iterators]{isplit}} function. It then runs
#' \code{get_all_ontologies} in parallel using \pkg{foreach} by sending each
#' data frame onto a worker.
#'
#' @param dat Data frame with two columns: \strong{entrez} and
#'  \strong{set_label}. The column \strong{entrez} should contain ENTREZ gene
#'  identifiers; the column \strong{set_label} should contain labels for
#'  identifiers that will be analyzed as a group.
#'
#' @param species Character string specifying the species, must be one of two
#'  "human" (default) or "mouse". For GO analysis we use
#'  \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db}, and for mouse \href{https://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db}.
#'
#' @param universe Character vector of unique ENTREZ gene identifiers to be used
#'  as a universe.
#'
#' @param ontologies Character vector of ontologies for which to run the
#'  analysis. If not provided the analysis will be run for all ontologies:
#'  cellular component (CC), biological process (BP) and molecular function
#'  (MF).
#'
#' @param cores Integer value representing the number of cores to use. This
#'  parameter is optional. If not provided
#'  \code{\link[foreach]{getDoParWorkers}} function will be called to determine
#'  the number of workers.
#'
#' @return Data frame with the results of GO enrichment.
#'
#' @seealso \code{get_ontology} to understand the output format.
#'
#' @examples
#' \dontrun{
#' # load test data:
#' data("human_symbol")
#' remove duplicate rows:
#' dat_clean <- deduplicate_rows(human_symbol)
#' # map symbol to ENTREZ gene id:
#' dat_mapped <- map_genes(dat_clean, id = "symbol", species = "human")
#' # extract universe
#' universe <- unique(dat_mapped$entrez)
#' # run GO enrichment in parallel
#' res <- run_parallel_go(dat_mapped, species = "human", universe = universe)
#' }
#'
#' @export


run_parallel_go <- function(dat, species = c("human", "mouse"),
                            universe,
                            ontologies, cores){
  if (!missing(cores)) {
    doParallel::registerDoParallel(cores = cores)
    workers <- foreach::getDoParWorkers()
    message(
      stringr::str_wrap(
        crayon::green(
          paste0("Will run GO computation on ", workers, " cores.")
        )
      )
    )
  } else {
    doParallel::registerDoParallel()
    workers <- foreach::getDoParWorkers()
    message(
      stringr::str_wrap(
        crayon::green(
          paste0("Parameter 'cores' is not provided. Getting the number of ",
                 "available cores with foreach::getDoParWorkers(). GO enrichment ",
                 "will be run on ", crayon::underline(workers), " cores.")
        )
      )
    )
  }
  species <- verify_input(input_name = species,
                          input_choices = c("human", "mouse"),
                          input_default = "human")
  if (missing(ontologies)) {
    ontologies <- c("BP", "CC", "MF")
  } else {
    if (!all(ontologies %in% c("BP", "CC", "MF"))) {
      stop(
        stringr::str_wrap(
          crayon::red(
            paste0("ERROR: Please provide valid values for the ontology ",
                   "parameter. Possible valid values: CC, BP, MF. ",
                   "You provided: ", crayon::underline(ontologies), ".")
          )
        )
      )
    } else {
      ontologies <- ontologies
    }
  }
  iterated_df <- iterators::isplit(dat, as.factor(dat$set_label))
  # to prevent complaining that %dopar% is not found
  `%dopar%` <- foreach::`%dopar%`
  res <- foreach::foreach(a = iterated_df,
                          .combine = rbind) %dopar% {
                            tryCatch({
                              get_all_ontologies(a$value$entrez,
                                                 species = species,
                                                 universe = universe,
                                                 ontologies = ontologies,
                                                 set_label = a$key[[1]])},
                              warning = function(w){
                                return_error_result(w,
                                                    a$value$entrez,
                                                    a$key[[1]])
                              },
                              error = function(e){
                                return_error_result(e,
                                                    a$value$entrez,
                                                    a$key[[1]])
                              }
                            )
                          }
  # if any jobs created errors or warnings, remove these sets from the output
  res <- remove_errors(res)
  return(res)
}


