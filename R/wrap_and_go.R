#' Main function for running GO analysis in parallel.
#'
#' \code{wrap_and_go} serves as a wrapper function for reading the file with
#' gene identifiers and set labels, applying minor cleaning, mapping to ENTREZ
#' gene identifiers and calculating GO enrichment analysis in parallel.
#'
#' @param path Character string specifying the path to file. This parameter is
#'  mandatory. For the file to be processed correctly ensure that gene
#'  identifiers are in the first column and the gene set labels are in the
#'  second column of the file.
#'
#' @param col_names A logical value indicating presence of column names in the
#'  first line of the file. This parameter is mandatory.
#'
#' @param delim Character string indicating field separator character in the
#'  file. This parameter is mandatory.
#'
#' @param id A character string specifying the type of gene identifier in the
#'  data, must of be one of two "symbol" (default) or "ensembl". This parameter
#'  is mandatory.
#'
#' @param species A character string specifying the species, must be one of two
#'  "human" (default) or "mouse". This parameter is mandatory.
#'
#' @param min_set_size Positive integer that defines a minimum number of genes
#'  in a set. This parameter is optional. If provided all sets that don't meet
#'  the membership criteria will be removed before GO enrichment analysis.
#'
#' @param cores Positive integer specifying number of cores to be used for
#'  parallel computation. This parameter is optional.
#'
#' @param run_parallel Boolean indicating whether to run the execution in
#'  parallel. Default is TRUE. If FALSE parameter \code{cores} will be ignored.
#'
#' @param ontologies Optional character vector with ontologies. Valid values:
#'  CC, BP, MF. If not provided GO enrichment will be run for all categories.
#'
#'
#' @return Data frame with the results of GO enrichment analysis.
#'
#' @examples
#' \dontrun{
#' # Example 1: provide all possible parameters
#'
#'  out <- wrap_and_go("path_to_my_file", col_names = TRUE, delim = ",",
#'  id = "symbol", species = "human", min_set_size = 30, cores = 4)
#'
#' # Example 2: provide minimum necessary set of parameters (no sets will be
#' removed, number of cores will be determined by foreach::getDoParWorkers)
#'
#' out <- wrap_and_go("path_to_my_file", col_names = TRUE, delim = ",",
#' id = "symbol", species = "human")
#'
#' # Example 3: with ensembl gene ids instead of gene symbols:
#'
#' out <- wrap_and_go("path_to_my_file", col_names = TRUE, delim = ",",
#' id = "ensembl", species = "human")
#'
#' }
#'
#' @export




wrap_and_go <- function(path, col_names, delim, id, species,
                        min_set_size, cores, run_parallel = TRUE,
                        ontologies){
  # read file
  dat <- read_file(path = path, col_names = col_names, delim = delim)
  # remove duplicate rows
  dat_clean <- deduplicate_rows(dat)
  # map to ENTREZ gene identifiers
  dat_mapped <- map_genes(dat_clean, id = id, species = species)
  # extract the universe for hypergeometric test
  universe <- unique(dat_mapped$entrez)
  # remove sets with a small number of genes
  dat_large_sets <- remove_small_sets(dat_mapped, min_set_size = min_set_size)
  # run GO analysis in parallel
  res <- run_go(dat_large_sets, species = species, universe = universe,
                cores = cores, run_parallel = run_parallel,
                ontologies = ontologies)
  return(res)
}
