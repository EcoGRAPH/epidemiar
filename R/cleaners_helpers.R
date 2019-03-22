#Data Cleaning and Helper function

## Data prep & cleaning functions
#interpolate NA data -- new case field -- cases_epidemiar & val_epidemiar
#' Interpolates missing epi data
#'
#' @return To fill in
#'
#'
epi_NA_interpolate <- function(epi_data, quo_casefield, quo_groupfield){
  epi_data %>%
    dplyr::group_by(!!quo_groupfield) %>%
    #confirm date sorting
    dplyr::arrange(obs_date) %>%
    #interpolate
    dplyr::mutate(cases_epidemiar = epidemiar::na_approx(!!quo_casefield)) %>%
    #finish by ungrouping
    dplyr::ungroup()
}

#' Interpolates missing env data
#'
env_NA_interpolate <- function(env_data, quo_obsfield, quo_valuefield, quo_groupfield){
  env_data %>%
    #first, mark which ones are observed versus (will be) interpolated
    dplyr::mutate(data_source = ifelse(!is.na(!!quo_valuefield), "Observed", "Interpolated")) %>%
    #two levels of group_by
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
    #confirm date sorting
    dplyr::arrange(obs_date) %>%
    #interpolate
    dplyr::mutate(val_epidemiar = !!quo_valuefield,
                  val_epidemiar = epidemiar::na_approx(val_epidemiar)) %>%
    #finish by ungrouping
    dplyr::ungroup()
}



#Helper functions
#' Create a named list
#' @export
#'
create_named_list <- function(...){
  list_to_name <- list(...)
  named_chr <- sapply(substitute(list(...)),deparse)[-1]
  given_names <- names(list_to_name)
  #if given no names, then just use names of original items
  if (is.null(given_names)) {
    names_to_use <- named_chr
  }
  #logical vector of which items did not have names given/to be assigned
  not_named_logical <- given_names == ""
  if (any(not_named_logical)) {
    #set up, names that were given
    names_to_use <- given_names
    #add names for things that had names originally
    names_to_use[not_named_logical] <- named_chr[not_named_logical]
  }
  #set the names
  stats::setNames(list_to_name, names_to_use)
}

#' Mode
#'
#' Calculate the mode of a set of values, for numeric or character/factor data. In ties, returns the first tied value.
#'
#' @param x A vector.
#' @param na.rm Logical indicating whether \code{NA} values should be excluded.
#'
#' @return A vector of length 1 and the same class as \code{x}.
#' @export
#'
#' @examples
#' Mode(c(1,1,2,3))
#' Mode(c(1,2,2,3))
#' Mode(c(1,1,3,3))
#' Mode(c(3,3,1,1))
#' Mode(c(1,NA,NA))
#' Mode(c(1,NA,NA), na.rm = TRUE)
#'
Mode <- function(x, na.rm = FALSE) {
  if (!is.vector(x)) stop("x is not a vector, but is class ",
                          paste(class(x), collapse = ", "))
  ux <- unique(x)
  if (na.rm) ux <- stats::na.exclude(ux)
  ux[which.max(tabulate(match(x, ux)))]
}

