#' Supporting function for verifying function parameters
#'
#' Not to be used by the user. Not exported

verify_input <- function(input.name, input.choices, input.default){
  if(!missing(input.name)){
    if(length(input.name) == 1){
      if(input.name %in% input.choices){
        return(input.name)
      }else{
        stop(paste0("Unexpected value for ", input.name,
                    " was provided; please provide one of ",
                    paste0(input.choices, collapse = " or "), "."))
      }
    }else{
      stop(paste0("More than one value for ", input.name,
                  " parameter was provided; please provide either ",
                  paste0(input.choices, collapse = " or "), "."))
    }
  }else{
    message(paste0("Using default setting for the parameter ", input.name, ": ",  input.default))
    return(input.default)
  }
}
