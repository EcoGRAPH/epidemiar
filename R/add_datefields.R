add_datefields <- function(df, type = "ISO"){
  match.arg(system, c("ISO", "CDC"))
  # get morb df (or accepts any df with Date field) and adds year and week numbers
  if (type == "ISO"){
    df <- df %>%
      #add iso wk/yr
      mutate(week_epidemiar = epidemiar::epiweek(Date, system = "who"),
             year_epidemiar = epidemiar::epiyear(Date, system = "who"))
  } else if (system == "CDC"){
    #add CDC epi wk/yr
    mutate(week_epidemiar = epidemiar::epiweek(Date, system = "cdc"),
           year_epidemiar = epidemiar::epiyear(Date, system = "cdc"))
  }
  df
}
