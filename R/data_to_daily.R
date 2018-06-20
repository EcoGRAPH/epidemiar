## Converts data to daily data
#' With any interval or missing data, converts into daily values
#' with optional linear interpolation
#' @export
#'
data_to_daily <- function(data_notdaily, valuefield, interpolate = TRUE){
  #converts any interval data to daily (with NA missing) & optional linear interpolation of NA
  #Date field as.Date()

  #dplyr programming, NSE
  quo_valuefield <- rlang::enquo(valuefield)

  data_1day <- data_notdaily %>%
    #should handle all grouping/categories, and therefore any ragged data
    dplyr::group_by_at(dplyr::vars(-Date, -!!quo_valuefield)) %>%
    #all explicit missing data - line for every Date for all groupings (above)
    tidyr::complete(Date = tidyr::full_seq(c(min(data_notdaily$Date), max(data_notdaily$Date)), 1)) %>%
    dplyr::ungroup()

  if (interpolate){
    data_1day <- data_1day %>%
      #should handle all grouping/categories, and therefore any ragged data
      dplyr::group_by_at(dplyr::vars(-Date, -!!quo_valuefield)) %>%
      #will not extrapolate beyond last known value, that will happen inside run_epidemia()
      mutate(!!quo_name(quo_valuefield) := epidemiar::na_approx(!!quo_valuefield)) %>%
      #finish by ungrouping
      dplyr::ungroup()
  }

  data_1day
}
