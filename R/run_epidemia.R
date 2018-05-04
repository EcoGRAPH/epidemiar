#' Run EPIDEMIA forecast models and early detection algorithm
#'
#' @param epi_data Epidemiological data with case numbers per week, with date field "Date"
#' @param casefield The field name for the case counts
#' @param populationfield Population field to give population numbers over time. Used to calculated incidence. Also optionally used in Farrington method for populationOffset.
#' @param groupfield The field name for districts or area divisions of epidemiological AND environmental data if no groupings (all one area), user should give a field with the same value throughout
#' @param week_type The standard (WHO ISO-8601 or CDC epi weeks) that the weeks of the year in epidemiological and environmental reference data use (required: dates listed are LAST day of week).
#' @param report_period The number of timesteps (weeks) that the whole report will cover. report_period - forecast_future is the number of weeks of past (known) data that will be included.
#' @param detection_period The number of timesteps (weeks) that the early detection will run over
#' @param ed_method Which method for early detection should be used ("Farrington" is only current) option.
#' @param ed_control All parameters for early detection algorithm, passed through to that subroutine
#' @param env_data Daily environmental data for same groupfields and Date range. Must be in long format, and must start laglen (in fc_control) days before epi_data for forecasting.
#' @param obsfield Field name of the environmental data observation types
#' @param valuefield Field name of the value of the environmental data observations
#' @param forecast_future Number of weeks from the end of the epi_data to extend the forecast out
#' @param fc_control Parameters for forecasting, including model and clusters
# <<>><<>> More documentation here!
#' @param env_ref_data Historical averages by week of year for environmental variables. Used in extended environmental data into the future for long forecast time, to calculate anomalies in early detection period, and to display on timeseries in reports
#' @param env_info Lookup table for environmental data - reference creation method, report label, GA ID, etc.
#'
#'
#' @return Returns a suite of summary and report data.
#'#'
#' @examples
#'
#' @export
#'


## Main Modeling (Early Detection, Forecasting) Function
run_epidemia <- function(epi_data, casefield, populationfield, groupfield, week_type = c("ISO", "CDC"),
                         report_period = 26,
                         ed_summary_period = 4, ed_method = c("Farrington", "EARS"), ed_control = NULL,
                         env_data, obsfield, valuefield, forecast_future = 4,
                         fc_control = NULL, env_ref_data, env_info){
  #temporary argument descriptions, until move into epidemiar and roxygenate appropriately
  #epi_data: epidemiological data with case number of time, currently weekly only, with date field "Date"
  #casefield: the field name for the case counts
  #populationfield: optional field to give population numbers over time
  #                 used in Farrington method for optional populationOffset
  #                 Planned use in report generation function to allow for case or incidence to be displayed
  #groupfield: the field name for districts or area divisions of epidemiological AND environmental data
  #             if no groupings (all one area), user should give a field with the same value throughout
  #week_type: for epidemiological data, whether the user wants to mark weeks of the year via WHO ISO-8601 or CDC epi weeks standards
  #           used internally in Farrington methods (not seen by user)
  #           planned to use in report generation function for display of week values on report
  #report_period: the number of timesteps (weeks) that the whole report will cover
  #               report_period - forecast_period will be the length of time early detection will run over
  #detection_period: the number of timesteps (weeks) that the early detection will run over
  #                   Evaluation range: Last detection_period weeks in epi data
  #ed_method: which method for early detection should be used
  #ed_control: all parameters for early detection algorithm, passed through to that subroutine
  #env_data: daily environmental data for same groupfields and Date range
  #           must be in long format
  #           must start <<>> days before epi_data: check??
  #           <<>>?separate data prep function (to be run first) to be written to take {8|X}-day to daily?
  #obsfield: field name of the environmental data observation types
  #valuefield: field name of the value of the environmental data observations
  #forecast_future: number of timesteps (weeks) from the end of the epi_data to extend the forecast out for
  #fc_control: <<stuff for modeling, tbd, Justin's code>>
  #env_ref_data: historical averages (by week), for anomalies, and report display on timeseries
  #env_info: lookup table for environmental data - reference creation method, future extension method, report label, GA ID, etc.

  #dplyr programming steps for passing of field names
  quo_casefield <- enquo(casefield)
  quo_popfield <- enquo(populationfield)
  quo_groupfield <- enquo(groupfield)
  quo_obsfield <- enquo(obsfield)
  quo_valuefield <- enquo(valuefield)

  #create alphabetical list of unique groups
  groupings <- pull(epi_data, !!quo_groupfield) %>% unique() %>% sort()
  #create alphabetical list of all unique environmental variables
  env_variables <- pull(env_data, !!quo_obsfield) %>% unique() %>% sort()

  ## Create report date information - for passing to interval functions, and report output
  #REM: report_period is full # of weeks of report.  forecast_future is how many of those weeks should be in the future.
  #full report
  report_dates <- list(full = list(min = max(epi_data$Date, na.rm = TRUE) -
                                     as.difftime((report_period - forecast_future - 1),
                                                 unit = "weeks"),
                                   max = max(epi_data$Date, na.rm = TRUE) +
                                     as.difftime(forecast_future,
                                                 units = "weeks")))
  report_dates$full$seq <- report_dates$full %>% {seq.Date(.$min, .$max, "week")}
  #dates with known epidemological data
  report_dates$known <- list(min = report_dates$full$min,
                             max = max(epi_data$Date, na.rm = TRUE))
  report_dates$known$seq <- report_dates$known %>% {seq.Date(.$min, .$max, "week")}
  #forecast period
  report_dates$forecast <- list(min = report_dates$known$max + as.difftime(1, units = "weeks"),
                                #could calculate from forecast_future, but already done so in $full
                                max = report_dates$full$max)
  report_dates$forecast$seq <- report_dates$forecast %>% {seq.Date(.$min, .$max, "week")}
  #early detection summary period (ED runs over full report, this is for summary in defined ED period)
  report_dates$ed_sum <- list(min = report_dates$known$max - as.difftime(ed_summary_period - 1, units = "weeks"),
                              max = report_dates$known$max)
  report_dates$ed_sum$seq <- report_dates$ed_sum %>% {seq.Date(.$min, .$max, "week")}


  ## Data checks and cleaning
  #check for NAs and interpolate as necessary
  #cases_epidemiar field name from data clearning (epi)
  epi_data <- epi_NA_interpolate(epi_data, quo_casefield, quo_groupfield) %>%
    #and sort by alphabetical groupfield
    arrange(!!quo_groupfield, Date)
  #val_epidemiar field name from data cleaning (env)
  env_data <- env_NA_interpolate(env_data, quo_obsfield, quo_valuefield, quo_groupfield) %>%
    #and sort by alphabetical groupfield
    arrange(!!quo_groupfield, !!quo_obsfield, Date)

  ## Set up output report data format
  #create observed data series
  obs_res <- epi_data %>%
    #include only observed data from requested start of report
    filter(Date >= report_dates$full$min) %>%
    mutate(series = "obs",
           #INCIDENCE; also note use of original not interpolated cases
           value = !!quo_casefield / !!quo_popfield * 1000,
           lab = "Observed",
           upper = NA,
           lower = NA) %>%
    select(!!quo_groupfield, Date, series, value, lab, upper, lower)


  ## Forecast
  fc_res_all <- run_forecast(epi_data, quo_popfield, quo_groupfield, groupings,
                             env_data, quo_obsfield, quo_valuefield, env_variables,
                             fc_control, env_ref_data, env_info, report_dates, week_type)


  ## Early detection
  #need to calculate early detection on existing epi data & FUTURE FORECASTED results
  future_fc <- fc_res_all$fc_epi %>%
    #get future forecasted results ONLY
    filter(Date %in% report_dates$forecast$seq)
  #combine existing and future
  obs_fc_epi <- bind_rows(epi_data, future_fc) %>%
    mutate(cases_epidemiar = ifelse(!are_na(cases_epidemiar),
                                    cases_epidemiar,
                                    fc_cases)) %>%
    #will be lost by end, but need for early detection methods using surveillance::sts objects
    epidemiar::add_datefields() %>%
    #arrange (for viewing/checking)
    arrange(!!quo_groupfield, Date)

  #run early detection on combined dataset
  ed_res <- run_early_detection(obs_fc_epi, quo_popfield, quo_groupfield,
                                groupings, ed_method, ed_control, report_dates)


  ## Combine epi datasets
  epi_res <- bind_rows(obs_res, fc_res_all$fc_res, ed_res)
  #add week fields
  modeling_results_data <- epidemiar::add_datefields(epi_res, week_type)


  ## Prep Environmental Data for report
  #using extended environmental data from forecast functions
  environ_timeseries <- environ_report_format(env_ext_data = fc_res_all$env_data_extd, env_ref_data,
                                              quo_groupfield, quo_obsfield,
                                              env_used = fc_res_all$env_variables_used, env_info,
                                              week_type, report_dates)

  ##Environmental Anomaly Data (during ED period)
  environ_anomalies <- calc_env_anomalies(env_ts = environ_timeseries,
                                          quo_groupfield, quo_obsfield, report_dates)

  ## Calculate Early Detection and Early Warning summary data
  summary_data <- create_summary_data(ed_res, quo_groupfield, report_dates)

  ## Calculate epidemiological (incidence) summary
  epi_summary <- create_epi_summary(obs_res, quo_groupfield, report_dates)

  ## Parameters and metadata that might be useful in report generation
  # all of these may not be needed
  fieldnames <- list(casefield = quo_name(quo_casefield),
                     populationfield = quo_name(quo_popfield),
                     groupfield = quo_name(quo_groupfield),
                     obsfield = quo_name(quo_obsfield),
                     valuefield = quo_name(quo_valuefield))
  params_meta <- create_named_list(fieldnames, week_type, ed_method, groupings,
                                   env_variables_used = fc_res_all$env_variables_used,
                                   report_dates, env_info)

  #for now
  all_results <- create_named_list(summary_data,
                                   epi_summary,
                                   modeling_results_data,
                                   environ_timeseries,
                                   environ_anomalies,
                                   params_meta)
  return(all_results)

}
