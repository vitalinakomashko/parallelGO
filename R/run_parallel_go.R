#' Run GO analysis on parallel
#'
#' @param annotation string providing the name of the annotation package to use,
#'  Can we use anything else but the 'org.Hs.eg.db' forhuman and 'org.Mm.eg.db' for
#'  mouse? I haven't investigated
#'  @param ontology vector of ontologies for which to run the analysis. If not provided
#' then the analysis will be run for all ontologies: CC, BP and MF


run_parallel_go <- function(dat, background, species = c("human", "mouse"), cores,
                            ontology){
  if(!missing(cores)){
    doParallel::registerDoParallel(cores = cores)
    workers <- getDoParWorkers()
    message(stringr::str_wrap(paste0("Will run GO computation on ", workers, ".")))
  }else{
    doParallel::registerDoParallel()
    workers <- getDoParWorkers()
    message(stringr::str_wrap(paste0("Will run GO computation on ", workers, ".")))
  }
  species <- verify_input(input.name = species, input.choices = c("human", "mouse"),
                          input.default = "human")
  if(missing(background)){
    stop("Please provide a vector with Entrez Gene identifiers to server as a background")
  }else{
    background <- unique(background)
  }
  if(missing(ontology)){
    ontologies <- c("BP", "CC", "MF")
  }else{
    if(!all(ontology %in% c("BP", "CC", "MF"))){
      stop(stringr::str_wrap(paste0("Please provide valid values for the ontology parameter.
                                    Possible valid values: CC, BP, MF. You provided: ",
                                    crayon::underline(ontology), ".")))
    }else{
      ontologies <- ontology
    }
  }
  df.iter <- iterators::isplit(dat, dat$module)
  res <- foreach::foreach(a = df.iter, .combine = rbind, .packages = c("GOstats", "dplyr"), .verbose = TRUE) %dopar% {
    get_all_ontologies(a, background = background, species = species, ontology = ontologies)
  }
  return(res)
}
