#' Adds ISO or CDC datefields to a data set with a "obs_date" (Date) field
#'
#' @return dataframe (or tibble) with week and year fields added
#'
#'

add_datefields <- function(df, type = "ISO"){
  type <- match.arg(type)
  # get morb df (or accepts any df with obs_date field) and adds year and week numbers
  if (type == "ISO"){
    df <- df %>%
      #add iso wk/yr
      dplyr::mutate(week_epidemiar = epidemiar::epiweek(obs_date, system = "who"),
                    year_epidemiar = epidemiar::epiyear(obs_date, system = "who"))
  } else if (type == "CDC"){
    #add CDC epi wk/yr
    dplyr::mutate(week_epidemiar = epidemiar::epiweek(obs_date, system = "cdc"),
                  year_epidemiar = epidemiar::epiyear(obs_date, system = "cdc"))
  }
  df
}
