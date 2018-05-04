#' Adds ISO or CDC datefields to a data set with a "Date" (Date) field
#'
#' @return dataframe (or tibble) with week and year fields added
#'
#' @export
#'

add_datefields <- function(df, type = "ISO"){
  type <- match.arg(type)
  # get morb df (or accepts any df with Date field) and adds year and week numbers
  if (type == "ISO"){
    df <- df %>%
      #add iso wk/yr
      mutate(week_epidemiar = epidemiar::epiweek(Date, system = "who"),
             year_epidemiar = epidemiar::epiyear(Date, system = "who"))
  } else if (type == "CDC"){
    #add CDC epi wk/yr
    mutate(week_epidemiar = epidemiar::epiweek(Date, system = "cdc"),
           year_epidemiar = epidemiar::epiyear(Date, system = "cdc"))
  }
  df
}
