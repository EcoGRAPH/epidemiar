#'Run EPIDEMIA forecast models and early detection algorithm.
#'
#'The Epidemic Prognosis Incorporating Disease and Environmental Monitoring for
#'Integrated Assessment (EPIDEMIA) Forecasting System is a set of tools coded in
#'free, open-access software, that integrate surveillance and environmental data
#'to model and create short-term forecasts for environmentally-mediated
#'diseases. This function, `epidemiar::run_epidemia()` is the central function
#'to model and forecast a wide range of environmentally-mediated diseases.
#'
#'For more a longer description of the package, see the overview vignette:
#'\code{vignette("overview-epidemiar", package = "epidemiar")}
#'
#'For more details see the vignette on input data and modeling parameters:
#'\code{vignette("data-modeling", package = "epidemiar")}
#'
#'@param epi_data Epidemiological data with case numbers per week, with date
#'  field "obs_date".
#'@param casefield The column name of the field that contains disease case
#'  counts (unquoted field name).
#'@param populationfield Column name of the population field to give population
#'  numbers over time (unquoted field name). Used to calculated incidence. Also
#'  optionally used in Farrington method for populationOffset.
#'@param inc_per Number for what unit of population the incidence should be
#'    reported in, e.g. incidence rate of 3 per 1000 people.
#'@param groupfield The column name of the field for district or geographic area
#'  unit division names of epidemiological AND environmental data (unquoted
#'  field name). If there are no groupings (all one area), user should give a
#'  field that contains the same value throughout.
#'@param week_type String indicating the standard (WHO ISO-8601 or CDC epi
#'  weeks) that the weeks of the year in epidemiological and environmental
#'  reference data use ["ISO" or "CDC"]. (Required: epidemiological observation
#'  dates listed are LAST day of week).
#'@param report_period The number of weeks that the entire report will cover.
#'  The \code{report_period} minus \code{forecast_future} is the number of weeks
#'  of past (known) data that will be included.
#'@param ed_summary_period The number of weeks that will be considered the "early
#'  detection period". It will count back from the week of last known
#'  epidemiological data.
#'@param ed_method Which method for early detection should be used ("Farrington"
#'  is only current option, or "None").
#'@param ed_control All parameters for early detection algorithm, passed through
#'  to that subroutine.
#'@param env_data Daily environmental data for the same groupfields and date
#'  range as the epidemiological data. It may contain extra data (other
#'  districts or date ranges). The data must be in long format (one row for each
#'  date and environmental variable combination), and must start at absolutel
#'  minimum \code{laglen} (in \code{fc_control}) days before epi_data for
#'  forecasting.
#'@param obsfield Field name of the environmental data variables (unquoted field
#'  name).
#'@param valuefield Field name of the value of the environmental data variable
#'  observations (unquoted field name).
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
#' For more details see the vignette on the output data:
#' \code{vignette("output-report-data", package = "epidemiar")}
#'
#'@examples "See model_forecast_script in epidemiar-demo for full example:
#'https://github.com/EcoGRAPH/epidemiar-demo"
#'
#'@export
#'
#'@importFrom magrittr %>%
#'@importFrom rlang !!


## Main Modeling (Early Detection, Forecasting) Function
run_epidemia <- function(epi_data = NULL,
                         casefield = NULL,
                         populationfield = NULL,
                         inc_per = 1000,
                         groupfield = NULL,
                         week_type = c("ISO", "CDC"),
                         report_period = 26,
                         ed_summary_period = 4,
                         ed_method = c("None", "Farrington"),
                         ed_control = NULL,
                         env_data = NULL,
                         obsfield = NULL,
                         valuefield = NULL,
                         forecast_future = 4,
                         fc_control = NULL,
                         env_ref_data = NULL,
                         env_info = NULL){


  # Non-standard evaluation quosures ----------------------------------------

  # dplyr programming steps for passing of field names
  quo_casefield <- rlang::enquo(casefield)
  quo_popfield <- rlang::enquo(populationfield)
  quo_groupfield <- rlang::enquo(groupfield)
  quo_obsfield <- rlang::enquo(obsfield)
  quo_valuefield <- rlang::enquo(valuefield)

  #Note: if field name does not exist in any dataset, enquo() will throw an error.


  # Preparing: Input checking -----------------------------------------------

  # 1. Test for critical inputs
  # This will not check if they've assigned the right thing to the argument, or got the argument order correct if not explicit argument declarations
  # But, no other checks can really proceed if things are missing
  # NSE is a little tricky.
  #can't test directly on fields-to-be-enquo'd because it'll try to evaluate them, and complain that the object (actually field name) doesn't exist
  #naming the quosures AS the input fields to create more meaningful error messages if the items are missing
  #populationfield eventually to be non necessary, but as of right now, things are reported in incidence, so population is critical
  nec_nse <- list(casefield = quo_casefield, groupfield = quo_groupfield, obsfield = quo_obsfield,
                  valuefield = quo_valuefield, populationfield = quo_popfield)
  necessary <- create_named_list(epi_data, env_data, env_ref_data, env_info, fc_control)
  #ed_control can be NULL if ed_method == None.
  # rest has defaults
  # Note: only checking if control list exists, nothing about what is in the list (later checks)
  #initialize missing info msgs & flag
  missing_msgs <- ""
  missing_flag <- FALSE
  #loop through all necessary fields, checking if argument exists, collecting list of missing
  for (nse in seq_along(nec_nse)){
    #testing if quosure was created on NULL object.
    if(rlang::quo_is_null(nec_nse[[nse]])){
      missing_flag <- TRUE
      missing_msgs <- paste0(missing_msgs, names(nec_nse[nse]), sep = "\n")
    }
  }
  for (arg in seq_along(necessary)){
    if (is.null(necessary[[arg]])){
      missing_flag <- TRUE
      missing_msgs <- paste0(missing_msgs, names(necessary[arg]), sep = "\n")
    }
  }
  #if missing, stop and give error message
  if (missing_flag){
    stop("Missing critical argument(s). Please make sure the following is/are included:\n", missing_msgs)
  }

  # 2. match.arg for arguments with options
  #Note: using message() instead of warning() to get message to appear right away
  ed_method <- tryCatch({
    match.arg(ed_method, c("None", "Farrington"))
  }, error = function(e){
    message("Warning: Given 'ed_method' does not match 'None' or 'Farrington', running as 'None'.")
    "None"
  }, finally = {
    if (length(ed_method) > 1){
      #if ed_method was missing at run_epidemia() call, got assigned c("None", "Farrington")
      message("Warning: 'ed_method' was missing, running as 'None'.")
      #no return, because in match.arg() it will take the first item, which is "None".
    }
  })
  week_type <- tryCatch({
    match.arg(week_type, c("ISO", "CDC"))
  }, error = function(e){
    message("Warning: Given 'week_type' does not match 'ISO' or 'CDC', running as 'ISO'.")
    "ISO"
  }, finally = {
    if (length(week_type) > 1){
      #if week_type was missing at run_epidemia() call, got assigned c("ISO", "CDC")
      message("Warning: 'week_type' was missing, running as 'ISO'.")
      #no return, because in match.arg() it will take the first item, which is "ISO".
    }
  })

  # 3. More input checking
  check_results <- input_check(epi_data,
                               quo_casefield,
                               quo_popfield,
                               inc_per,
                               quo_groupfield,
                               week_type,
                               report_period,
                               ed_summary_period,
                               ed_method,
                               ed_control,
                               env_data,
                               quo_obsfield,
                               quo_valuefield,
                               forecast_future,
                               fc_control,
                               env_ref_data,
                               env_info)
  #if warnings, just give message and continue
  if (check_results$warn_flag){
    message(check_results$warn_msgs)
  }
  #if then if errors, stop and return error messages
  if (check_results$err_flag){
    #prevent possible truncation of all error messages
    options(warning.length = 4000L)
    stop(check_results$err_msgs)
  }



  # Preparing: generating listings and date sets ----------------------------

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



  # Preparing: data checks for NA and interpolation -------------------------

  #check for NAs and interpolate as necessary
  #Note: cases_epidemiar is field name returned (epi)
  epi_data <- epi_NA_interpolate(epi_data, quo_casefield, quo_groupfield) %>%
    #and sort by alphabetical groupfield
    dplyr::arrange(!!quo_groupfield, obs_date)
  #Note: val_epidemiar is field name returned (env)
  env_data <- env_NA_interpolate(env_data, quo_obsfield, quo_valuefield, quo_groupfield) %>%
    #and sort by alphabetical groupfield
    dplyr::arrange(!!quo_groupfield, !!quo_obsfield, obs_date)



  # Set up output report data format ----------------------------------------

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



  # Forecasting -------------------------------------------------------------

  fc_res_all <- run_forecast(epi_data, quo_popfield, inc_per, quo_groupfield, groupings,
                             env_data, quo_obsfield, quo_valuefield, env_variables,
                             fc_control, env_ref_data, env_info, report_dates, week_type)



  # Event detection ---------------------------------------------------------

  #need to calculate event detection on existing epi data & FUTURE FORECASTED results
  future_fc <- fc_res_all$fc_epi %>%
    #get future forecasted results ONLY
    dplyr::filter(obs_date %in% report_dates$forecast$seq)
  #combine existing and future
  obs_fc_epi <- dplyr::bind_rows(epi_data, future_fc) %>%
    dplyr::mutate(cases_epidemiar = ifelse(!rlang::are_na(cases_epidemiar),
                                           cases_epidemiar,
                                           fc_cases)) %>%
    #will be lost by end, but need for event detection methods using surveillance::sts objects
    epidemiar::add_datefields() %>%
    #arrange (for viewing/checking)
    dplyr::arrange(!!quo_groupfield, obs_date)

  #run event detection on combined dataset
  ed_res <- run_event_detection(epi_fc_data = obs_fc_epi,
                                quo_popfield, inc_per,
                                quo_groupfield, groupings,
                                ed_method, ed_control, report_dates)



  # Combining forecast and event detection results --------------------------

  ## Combine epi datasets
  epi_res <- dplyr::bind_rows(obs_res, fc_res_all$fc_res, ed_res)
  #add week fields
  modeling_results_data <- epidemiar::add_datefields(epi_res, week_type)



  # Format other data for report --------------------------------------------

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


  # Collect all results -----------------------------------------------------

  all_results <- create_named_list(summary_data,
                                   epi_summary,
                                   modeling_results_data,
                                   environ_timeseries,
                                   environ_anomalies,
                                   params_meta,
                                   regression_object)
  return(all_results)

}
