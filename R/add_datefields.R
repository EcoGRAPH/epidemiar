#' Adds ISO or CDC datefields to a data set with a "obs_date" (Date) field
#'
#' @return dataframe (or tibble) with week and year fields added
#'
#'

add_datefields <- function(df, type = "ISO"){
  type <- match.arg(type, c("ISO", "CDC"))
  # get morb df (or accepts any df with obs_date field) and adds year and week numbers
  if (type == "ISO"){
    df <- df %>%
      #add iso wk/yr
      dplyr::mutate(week_epidemiar = lubridate::isoweek(obs_date),
                    year_epidemiar = lubridate::isoyear(obs_date))
  } else if (type == "CDC"){
    #add CDC epi wk/yr
    dplyr::mutate(week_epidemiar = lubridate::epiweek(obs_date),
                  year_epidemiar = lubridate::epiyear(obs_date))
  }
  df
}
