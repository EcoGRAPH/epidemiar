#'Run EPIDEMIA forecast models and early detection algorithm
#'
#'@param epi_data Epidemiological data with case numbers per week, with date
#'  field "obs_date".
#'@param casefield The column name of the field that contains disease case
#'  counts.
#'@param populationfield Column name of the population field to give population
#'  numbers over time. Used to calculated incidence. Also optionally used in
#'  Farrington method for populationOffset.
#'@param groupfield The column name of the field for district or geographic area
#'  unit division names of epidemiological AND environmental data. If there are
#'  no groupings (all one area), user should give a field that contains the same
#'  value throughout.
#'@param week_type The standard (WHO ISO-8601 or CDC epi weeks) that the weeks
#'  of the year in epidemiological and environmental reference data use
#'  (required: dates listed are LAST day of week).
#'@param report_period The number of weeks that the entire report will cover.
#'  The \code{report_period} minus \code{forecast_future} is the number of weeks
#'  of past (known) data that will be included.
#'@param detection_period The number of weeks that will be considered the "early
#'  detection period". It will count back from the week of last known
#'  epidemiological data.
#'@param ed_method Which method for early detection should be used ("Farrington"
#'  is only current option, or "None").
#'@param ed_control All parameters for early detection algorithm, passed through
#'  to that subroutine.
#'@param env_data Daily environmental data for the same groupfields and date
#'  range as the epidemiological data. It may contain extra data (other
#'  districts or date ranges). The data must be in long format (one row for each
#'  date and environmental variable combination), and must start \code{laglen}
#'  (in \code{fc_control}) days before epi_data for forecasting.
#'@param obsfield Field name of the environmental data variables.
#'@param valuefield Field name of the value of the environmental data variable
#'  observations.
#'@param forecast_future Number of futre weeks from the end of the
#'  \code{epi_data} to produce forecasts.
#'@param fc_control Parameters for forecasting, including which environmental
#'  variable to include and any geographic clusters.
#'@param env_ref_data Historical averages by week of year for environmental
#'  variables. Used in extended environmental data into the future for long
#'  forecast time, to calculate anomalies in early detection period, and to
#'  display on timeseries in reports.
#'@param env_info Lookup table for environmental data - reference creation
#'  method (e.g. sum or mean), report labels, etc.
#'
#'
#'@return Returns a suite of summary and report data.
#'
#'  1. \code{summary_data}: Early detection and early warning alerts levels for
#'  each woreda. Early detection alerts (ed_alert_count) are alerts that are
#'  triggered during the early detection period, which is defined as the 4 most
#'  recent weeks of known epidemiology data. Similarly, early warning alerts
#'  (ew_alert_count) are alerts in the future forecast estimates. “High” level
#'  indicates two or more weeks in this period had incidences greater than the
#'  alert threshold, “Medium” means that one week was in alert status, and “Low”
#'  means no weeks had alerts (ed_sum_level and ew_level, respectively).
#'
#'  2. \code{epi_summary}: Mean disease incidence per geographic group during
#'  the early detection period.
#'
#'  3. \code{modeling_results_data}:These are multiple timeseries values for
#'  observed, forecast, and alert thresholds of disease incidence, over the
#'  report period, for each geographic unit. These data can be used in creating
#'  the individual geographic unit control charts.
#'
#'  4. \code{environ_timeseries}: These are multiple timeseries for the
#'  environmental variables during the report period for each geographic unit.
#'
#'  5. \code{environ_anomalies}: These data are the recent (during the early
#'  detection period) differences (anomalies) of the environmental variable
#'  values from the climatology/reference mean.
#'
#'  6. \code{params_meta}: This lists dates, settings, and parameters that
#'  \code{run_epidemiar()} was called with.
#'
#'  7. \code{regression_object}: This is the regression object from the general
#'  additive model (GAM, parallelized with BAM) regression. This is only for
#'  statistical investigation of the model, and is usually not saved (very large
#'  object).
#'
#'@examples See model_forecast_script.R in epidemiar-demo for full example: https://github.com/EcoGRAPH/epidemiar-demo
#'
#'@export
#'
#'@importFrom magrittr %>%
#'@importFrom rlang !!


## Main Modeling (Early Detection, Forecasting) Function
run_epidemia <- function(epi_data, casefield, populationfield, inc_per = 1000,
                         groupfield, week_type = c("ISO", "CDC"),
                         report_period = 26,
                         ed_summary_period = 4, ed_method = c("Farrington", "None"), ed_control = NULL,
                         env_data, obsfield, valuefield, forecast_future = 4,
                         fc_control = NULL, env_ref_data, env_info){
  #temporary argument descriptions, until move into epidemiar and roxygenate appropriately
  #epi_data: epidemiological data with case number of time, currently weekly only, with date field "obs_date"
  #casefield: the field name for the case counts
  #populationfield: population numbers over time
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
  #           must start 180 days before first day in epi_data
  #obsfield: field name of the environmental data observation types
  #valuefield: field name of the value of the environmental data observations
  #forecast_future: number of timesteps (weeks) from the end of the epi_data to extend the forecast out for
  #fc_control: control options for forecasting
  #env_ref_data: historical averages (by week), for anomalies, and report display on timeseries
  #env_info: lookup table for environmental data - reference creation method, future extension method, report label, GA ID, etc.

  #dplyr programming steps for passing of field names
  quo_casefield <- rlang::enquo(casefield)
  quo_popfield <- rlang::enquo(populationfield)
  quo_groupfield <- rlang::enquo(groupfield)
  quo_obsfield <- rlang::enquo(obsfield)
  quo_valuefield <- rlang::enquo(valuefield)

  #create alphabetical list of unique groups
  #must remain in alpha order for early detection using surveillance package to capture results properly
  groupings <- dplyr::pull(epi_data, !!quo_groupfield) %>% unique() %>% sort()
  #create alphabetical list of all unique environmental variables
  env_variables <- dplyr::pull(env_data, !!quo_obsfield) %>% unique() %>% sort()

  ## Create report date information - for passing to interval functions, and report output
  #REM: report_period is full # of weeks of report.  forecast_future is how many of those weeks should be in the future.
  #full report
  report_dates <- list(full = list(min = max(epi_data$obs_date, na.rm = TRUE) -
                                     lubridate::as.difftime((report_period - forecast_future - 1),
                                                            unit = "weeks"),
                                   max = max(epi_data$obs_date, na.rm = TRUE) +
                                     lubridate::as.difftime(forecast_future,
                                                            units = "weeks")))
  report_dates$full$seq <- report_dates$full %>% {seq.Date(.$min, .$max, "week")}
  #dates with known epidemological data
  report_dates$known <- list(min = report_dates$full$min,
                             max = max(epi_data$obs_date, na.rm = TRUE))
  report_dates$known$seq <- report_dates$known %>% {seq.Date(.$min, .$max, "week")}
  #forecast period
  report_dates$forecast <- list(min = report_dates$known$max + lubridate::as.difftime(1, units = "weeks"),
                                #could calculate from forecast_future, but already done so in $full
                                max = report_dates$full$max)
  report_dates$forecast$seq <- report_dates$forecast %>% {seq.Date(.$min, .$max, "week")}
  #early detection summary period (ED runs over full report, this is for summary in defined ED period)
  report_dates$ed_sum <- list(min = report_dates$known$max - lubridate::as.difftime(ed_summary_period - 1, units = "weeks"),
                              max = report_dates$known$max)
  report_dates$ed_sum$seq <- report_dates$ed_sum %>% {seq.Date(.$min, .$max, "week")}


  ## Data checks and cleaning
  #check for NAs and interpolate as necessary
  #cases_epidemiar field name from data cleaning (epi)
  epi_data <- epi_NA_interpolate(epi_data, quo_casefield, quo_groupfield) %>%
    #and sort by alphabetical groupfield
    dplyr::arrange(!!quo_groupfield, obs_date)
  #val_epidemiar field name from data cleaning (env)
  env_data <- env_NA_interpolate(env_data, quo_obsfield, quo_valuefield, quo_groupfield) %>%
    #and sort by alphabetical groupfield
    dplyr::arrange(!!quo_groupfield, !!quo_obsfield, obs_date)

  ## Set up output report data format
  #create observed data series
  obs_res <- epi_data %>%
    #include only observed data from requested start of report
    dplyr::filter(obs_date >= report_dates$full$min) %>%
    dplyr::mutate(series = "obs",
                  #INCIDENCE; also note use of original not interpolated cases
                  value = !!quo_casefield / !!quo_popfield * inc_per,
                  lab = "Observed",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)


  ## Forecast
  fc_res_all <- run_forecast(epi_data, quo_popfield, inc_per, quo_groupfield, groupings,
                             env_data, quo_obsfield, quo_valuefield, env_variables,
                             fc_control, env_ref_data, env_info, report_dates, week_type)


  ## Early detection
  #need to calculate early detection on existing epi data & FUTURE FORECASTED results
  future_fc <- fc_res_all$fc_epi %>%
    #get future forecasted results ONLY
    dplyr::filter(obs_date %in% report_dates$forecast$seq)
  #combine existing and future
  obs_fc_epi <- dplyr::bind_rows(epi_data, future_fc) %>%
    dplyr::mutate(cases_epidemiar = ifelse(!rlang::are_na(cases_epidemiar),
                                           cases_epidemiar,
                                           fc_cases)) %>%
    #will be lost by end, but need for early detection methods using surveillance::sts objects
    epidemiar::add_datefields() %>%
    #arrange (for viewing/checking)
    dplyr::arrange(!!quo_groupfield, obs_date)

  #run event detection on combined dataset
  ed_res <- run_event_detection(epi_fc_data = obs_fc_epi,
                                quo_popfield, inc_per,
                                quo_groupfield, groupings,
                                ed_method, ed_control, report_dates)


  ## Combine epi datasets
  epi_res <- dplyr::bind_rows(obs_res, fc_res_all$fc_res, ed_res)
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
                                   env_dt_ranges = fc_res_all$env_dt_ranges,
                                   report_dates, env_info)
  #regression object for future other use or troubleshooting
  regression_object <- fc_res_all$reg_obj

  #collect results
  all_results <- create_named_list(summary_data,
                                   epi_summary,
                                   modeling_results_data,
                                   environ_timeseries,
                                   environ_anomalies,
                                   params_meta,
                                   regression_object)
  return(all_results)

}
