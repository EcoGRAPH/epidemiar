#Formatting and Report section data calculators

## Environmental Data for report
#' Formats env data for report
#' @export
#'
environ_report_format <- function(env_ext_data, env_ref_data, quo_groupfield,
                                  quo_obsfield, env_used, env_info,
                                  week_type, report_dates){
  #daily env data
  env_data_varused <- env_ext_data %>%
    dplyr::filter(!!quo_obsfield %in% env_used)

  #reference/climatology environmental data
  env_ref_varused <- env_ref_data %>%
    dplyr::filter(!!quo_obsfield %in% env_used)

  ##properly summarize to weekly (from daily)
  env_data_varused_sum <- env_data_varused %>%
    #get reference/summarizing method from user supplied env_info
    dplyr::left_join(env_info %>%
                       dplyr::select(!!quo_obsfield, reference_method),
                     by = rlang::set_names(rlang::quo_name(quo_obsfield),
                                           rlang::quo_name(quo_obsfield))) %>%
    #add week, year fields
    epidemiar::add_datefields(week_type) %>%
    #trim dates to reduce processing (dates are rough, technically just need week prior to start. 8 is not magical)
    dplyr::filter(Date >= report_dates$full$min - 8 & Date <= report_dates$full$max + 8) %>%
    #group by grouping, env var, and date week
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield, year_epidemiar, week_epidemiar) %>%
    #calculate with case_when at row level (fx is not vectorized, so can't be used inside summarize)
    dplyr::mutate(val_epidemiar = case_when(
      reference_method == "sum"  ~ sum(val_epidemiar, na.rm = TRUE),
      reference_method == "mean" ~ mean(val_epidemiar, na.rm = TRUE),
      #default is mean
      TRUE                       ~ mean(val_epidemiar, na.rm = TRUE))) %>%
    #now summarize
    #max Date of that week is how the weekly dates are set up
    dplyr::summarize(Date = max(Date),
                     #val_epi is the same for the whole grouped set, so just taking the first value
                     val_epidemiar = first(val_epidemiar),
                     #will be same throughout week
                     reference_method = first(reference_method),
                     #observed/interpolated/extended -- Mode, whatever source was most often that week.
                     data_source = Mode(data_source, na.rm = TRUE)) %>%
    #ungroup to end
    dplyr::ungroup()

  #filter exact dates
  environ_timeseries <- env_data_varused_sum %>%
    dplyr::filter(Date >= report_dates$full$min & Date <= report_dates$full$max) %>%
    dplyr::arrange(!!quo_groupfield, Date, !!quo_obsfield)

  # add climatology data
  # climatology is based on week number
  #   (hopefully set up with the same type as was selected when ref data was created, will add checks)
  environ_timeseries <- environ_timeseries %>%
    #join
    dplyr::left_join(env_ref_varused %>%
                       dplyr::select(!!quo_obsfield, !!quo_groupfield, week_epidemiar,
                                     ref_value, ref_sd, ref_median, ref_uq, ref_lq),
                     #NSE fun
                     by = rlang::set_names(c(rlang::quo_name(quo_groupfield),
                                             rlang::quo_name(quo_obsfield),
                                             "week_epidemiar"),
                                           c(rlang::quo_name(quo_groupfield),
                                             rlang::quo_name(quo_obsfield),
                                             "week_epidemiar")))
}


## Setting up summary data
#' Creates summary data
#' @export
#'
create_summary_data <- function(ed_res, quo_groupfield, report_dates){

  #levels
  alert_level <- c("Low", "Medium", "High")

  #Early Detection
  ed_summary <- ed_res %>%
    #get the alert series
    dplyr::filter(series == "ed") %>%
    #filter to early detection period
    dplyr::filter(Date %in% report_dates$ed_sum$seq) %>%
    #group (because need to look at period per group level)
    dplyr::group_by(!!quo_groupfield) %>%
    #summarize to 1 obs per grouping
    dplyr::summarize(ed_alert_count = sum(value, na.rm = TRUE)) %>%
    # create 3 levels (0, 1, 2 = >1)
    dplyr::mutate(warning_level = if_else(ed_alert_count > 1, 2, ed_alert_count),
                  #factor to label
                  ed_sum_level = factor(warning_level, levels = 0:2,
                                        labels = alert_level, ordered = TRUE)) %>%
    #ungroup
    dplyr::ungroup() %>%
    #select minimal cols
    dplyr::select(!!quo_groupfield, ed_alert_count, ed_sum_level)


  #Early Warning: ED results on forecast
  ew_summary <- ed_res %>%
    #get the alert series
    dplyr::filter(series == "ew",
                  #get the forecast results ##not needed anymore b/c of new ew series, but just for completeness
                  Date %in% report_dates$forecast$seq) %>%
    #group
    dplyr::group_by(!!quo_groupfield) %>%
    #summarize to 1 obs per grouping
    dplyr::summarize(ew_alert_count = sum(value, na.rm = TRUE)) %>%
    # create 3 levels (0, 1, 2 = >1)
    dplyr::mutate(warning_level = if_else(ew_alert_count > 1, 2, ew_alert_count),
                  #factor to label
                  ew_level = factor(warning_level, levels = 0:2,
                                    labels = alert_level, ordered = TRUE)) %>%
    #ungroup
    dplyr::ungroup() %>%
    #select minimal cols
    dplyr::select(!!quo_groupfield, ew_alert_count, ew_level)

  #join results
  summary_data <- dplyr::inner_join(ed_summary, ew_summary,
                                    by = rlang::set_names(rlang::quo_name(quo_groupfield),
                                                          rlang::quo_name(quo_groupfield)))

  summary_data
}

#' Creates summary of incidence in ED period
#' @export
#'
create_epi_summary <- function(obs_res, quo_groupfield, report_dates){
  #using obs_res - if cases/incidence becomes a user set choice, this might make it easier (value is already what it needs to be)
  #but note that (as of writing this) that obs_res using the original, UNinterpolated values (so that end users are disturbed to see case data where there should not be)

  epi <- obs_res %>%
    #epi data is weekly, get the data for the early detection summary period
    dplyr::filter(Date %in% report_dates$ed_sum$seq) %>%
    #group by groupings
    dplyr::group_by(!!quo_groupfield) %>%
    #get mean incidence
    dplyr::summarize(mean_inc = mean(value, na.rm = TRUE))

}




## Calculate anomalies
#' Calculates anomalies
#' @export
#'
calc_env_anomalies <- function(env_ts, quo_groupfield, quo_obsfield, report_dates){
  # anomalies
  anom_env <- env_ts %>%
    # only mapping those in the early detection period
    dplyr::filter(Date %in% report_dates$ed_sum$seq) %>%
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
    # anomaly value is observed value minus the ref value from env_ref
    dplyr::mutate(anom = val_epidemiar - ref_value) %>%
    # summarized over ED period
    dplyr::summarize(anom_ed_mean = mean(anom, na.rm = TRUE)) %>%
    dplyr::ungroup()
}
