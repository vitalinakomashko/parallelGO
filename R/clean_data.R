#' Basic data cleaning before gene mapping
#' 
#' `clean_data` function removes duplicated rows, sets that have fewer genes
#' than `min.set.size`.
#'
#' @param dat Data frame object with two columns: genes (`id`) and set labels (`set.name`),
#'  output of the function \code{\link{read_data}} call.
#' @param min.set.size Positive integer. If provided will be used to remove
#'  all gene sets where then number of genes is fewer than \code{min.set.size}.
#'
#' @return Data frame object
#' @export


clean_data <- function(dat, min.set.size){
  # check the data for duplicated rows
  if(any(duplicated(dat))){
    message(stringr::str_wrap(crayon::yellow("ATTENTION: Found duplicated rows in 
                                             provided data, removing duplicates.")))
    dat <- unique(dat)
    message(stringr::str_wrap(paste0("Number of rows after removing
                                             duplicated rows is ", nrow(dat), ".")))
  }
  # remove modules 
  if(!missing(min.set.size)){
    if(min.set.size <= 0){
      stop(stringr::str_wrap(paste0("You provided ", min.set.size, 
                                    " as 'min.set.size'. Parameter 'min.set.size' 
                                    should be a positive number.")))
    }else{
      message(stringr::str_wrap(paste0("Total number of sets in the dataset is ", 
                                       length(unique(dat$set.)), ".")))
      gene.count.by.set <- dat %>% group_by(.data$set.label) %>% 
        summarise(n = n())
      low.count.set.list <- gene.count.by.set %>% 
        filter(.data$n <= min.set.size)
      low.count.set.total <- low.count.set.list %>% summarize(total = n())
      message(paste0("You provided ", min.set.size, 
                                       " genes as 'min.set.size'. \nFound ", 
                                       as.character(low.count.set.total), 
                                       " sets with number of genes less than ", 
                                       min.set.size, "."))
      if(as.numeric(low.count.set.total) > 0){
        dat <- dat[!(dat$set.label %in% low.count.set.list$set.label), ]
        if(nrow(dat) == 0){
          stop("After removing clusters with fewer than 'min.set.size' no rows remained in the data.")
        }
        message(stringr::str_wrap(crayon::yellow(paste0("ATTENTION: Removed sets where 
                                                      the number of genes is less than ", 
                                                        min.set.size, 
                                                        ". Number of rows after filtering is ", 
                                                        nrow(dat), ". Number of sets after filtering is ",
                                                        length(unique(dat$set.label)), "."))))
      }else{
        message("No sets will be removed.")
      }
    }
  }else{
    message(stringr::str_wrap("Parameter 'min.set.size' is not provided. No gene
                              sets will be removed."))
  }
  return(dat)
}