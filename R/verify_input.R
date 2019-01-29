#' Supporting function for verifying function parameters
#'
#' Not to be used by the user. Not exported

verify_input <- function(input_name, input_choices, input_default){
  if (!missing(input_name)) {
    if (length(input_name) == 1) {
      if (input_name %in% input_choices) {
        return(input_name)
      } else {
        stop(
          paste0("Unexpected value for ", input_name,
                 " was provided; please provide one of ",
                 paste0(input_choices, collapse = " or "), ".")
          )
      }
    } else {
      stop(
        paste0("More than one value for ", input_name,
               " parameter was provided; please provide either ",
               paste0(input_choices, collapse = " or "), "."))
    }
  } else {
    message(
      paste0("Using default setting for the parameter ",
             input_name, ": ",
             input_default)
      )
    return(input_default)
  }
}
