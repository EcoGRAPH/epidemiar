#'Create historical reference data from daily environmental datasets.
#'
#'This function takes in a daily environmental data set, formatted the same as
#'env_data in the run_epidemia() function, and turns it into a climatology
#'dataset of historical statistics per week over all years.
#'
#'@param daily_env_data Daily environmental data. Must be in long format - one
#'  row for each date, environmental variable, and geographic grouping. Same
#'  format as the input env_data in run_epidemia().
#'@param groupfield The column name of the field for district or geographic area
#'  unit division names of environmental data (unquoted field name). If there
#'  are no groupings (all one area), user should give a field that contains the
#'  same value throughout.
#'@param obsfield Field name of the environmental data variables (unquoted field
#'  name).
#'@param valuefield Field name of the value of the environmental data variable
#'  observations (unquoted field name).
#'@param week_type The standard (WHO ISO-8601 or CDC epi weeks) that the weeks
#'  of the year that match with the epidemiological data when running
#'  run_epidemia().
#'@param env_info Lookup table for environmental data - reference creation
#'  method (e.g. sum or mean) for daily to weekly values.
#'
#'
#'
#'@return Returns a dataset with mean and other reference statistics of the
#'  environmental variable per geographic grouping per week of year, in the same
#'  format that would be acceptable as env_info in run_epidemia().
#'
#'
#'@examples env_daily_to_ref(daily_env_data = am_env_data,
#'                           groupfield = woreda_name,
#'                           obsfield = environ_var_code,
#'                           valuefield = obs_value,
#'                           week_type = "ISO",
#'                           env_info = am_env_info)
#' See epidemiar-demo: https://github.com/EcoGRAPH/epidemiar-demo
#'
#'@export
#'
#'@importFrom magrittr %>%
#'@importFrom rlang !!
#'

env_daily_to_ref <- function(daily_env_data,
                             groupfield,
                             obsfield,
                             valuefield,
                             week_type = c("ISO", "CDC"),
                             env_info){

  #dplyr programming steps for passing of field names
  quo_groupfield <- rlang::enquo(groupfield)
  quo_obsfield <- rlang::enquo(obsfield)
  quo_valuefield <- rlang::enquo(valuefield)

  #First, turn daily into weekly values, function to summarize as appropriate
  message("Summarizing environmental data to weekly values")
  env_weekly <- daily_env_data %>%
    #get reference/summarizing method from user supplied env_info
    dplyr::inner_join(env_info %>%
                        dplyr::select(!!quo_obsfield, reference_method),
                      by = rlang::set_names(rlang::quo_name(quo_obsfield),
                                            rlang::quo_name(quo_obsfield))) %>%
    #add week, year fields
    epidemiar::add_datefields(week_type) %>%
    #group by grouping, env var, and date week
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield, year_epidemiar, week_epidemiar) %>%
    #calculate with case_when at row level (fx is not vectorized, so can't be used inside summarize)
    dplyr::mutate(val_epidemiar = case_when(
      reference_method == "sum"  ~ sum(!!quo_valuefield, na.rm = TRUE),
      reference_method == "mean" ~ mean(!!quo_valuefield, na.rm = TRUE),
      #default is mean, but since inner_join with info table should not be invoked
      TRUE                       ~ mean(!!quo_valuefield, na.rm = TRUE))) %>%
    #now summarize
    #val_epi is the same for the whole grouped set, so just taking the first value
    dplyr::summarize(val_epidemiar = first(val_epidemiar)) %>%
    #ungroup to end
    dplyr::ungroup()

  #Then, summarize weekly values over the years, multiple reference statistics
  message("Summarizing weekly values over years")
  env_ref <- env_weekly %>%
    #dropping year from grouping list
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield, week_epidemiar) %>%
    #summarize over years
    dplyr::summarize(ref_value = mean(val_epidemiar, na.rm = TRUE),
                     ref_sd = sd(val_epidemiar, na.rm = TRUE),
                     ref_yrcount = n(),
                     ref_max = max(val_epidemiar, na.rm = TRUE),
                     ref_uq = quantile(val_epidemiar, probs = 0.75, na.rm = TRUE),
                     ref_median = median(val_epidemiar, na.rm = TRUE),
                     ref_lq = quantile(val_epidemiar, probs = 0.25, na.rm = TRUE),
                     ref_min = min(val_epidemiar, na.rm = TRUE)) %>%
    #ungroup to end
    dplyr::ungroup()

  env_ref
}
