## Converts data to daily data
#' With any interval or missing data, converts into daily values
#' with optional linear interpolation
#' @export
#'
data_to_daily <- function(data_notdaily, valuefield, interpolate = TRUE){
  #converts any interval data to daily (with NA missing) & optional linear interpolation of NA
  #Date field as.Date()

  #dplyr programming, NSE
  quo_valuefield <- enquo(valuefield)

  data_1day <- data_notdaily %>%
    #should handle all grouping/categories, and therefore any ragged data
    group_by_at(vars(-Date, -!!quo_valuefield)) %>%
    #all explicit missing data - line for every Date for all groupings (above)
    complete(Date = full_seq(c(min(data_notdaily$Date), max(data_notdaily$Date)), 1)) %>%
    ungroup()

  #Note: unable to get group_by_at to use -var_names properly inside function (testing line by line runs okay)

  if (interpolate){
    data_1day <- data_1day %>%
      #should handle all grouping/categories, and therefore any ragged data
      group_by_at(vars(-Date, -!!quo_valuefield)) %>%
      #will not extrapolate beyond last known value, that will happen inside run_epidemia()
      mutate(val_epidemiar = epidemiar::na_approx(!!quo_valuefield)) %>%
      #finish by ungrouping
      ungroup()
  }

  data_1day
}
