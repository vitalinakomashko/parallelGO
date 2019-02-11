#' Main function for running GO analysis
#'
#' @param path The name of the fle which the data are to be read from
#' @param col_names a logical value indicating whether the file contains the names
#' of the variables as its first line
#' @param delim The file field separator character.
#' @param min_set_size Positive integer. If provided will be used to remove
#'  all gene sets where then number of genes is fewer than \code{min.set.size}.
#' @param id A character string specifying the type of gene identifier in the
#'  data, must of be one of two "hugo" (default) or "ensembl".
#' @param species A character string specifying the species, must be one of two
#'  "human" (default) or "mouse". For mapping to Entrez Gene identifiers we use
#'  'org.Hs.eg.db' (human) and 'org.Mm.eg.db' (mouse) Bioconductor annotation
#'  packages.
#' @export




run_GO <- function(path, col_names, delim, min_set_size, species, id){
  # read the file
  dat <- read_file(path = path, col_names = col_names, delim = delim)
  # do basic cleaning
  dat_clean <- clean_data(dat)
  # map to Entrez Gene identifiers
  dat_mapped <- map_genes(dat_clean, id, species)
  universe <- unique(dat_mapped$entrez)
  dat_large_sets <- remove_small_sets(dat_mapped, min_set_size)
  # run GO
  res <- run_parallel_go(dat_large_sets)
}
