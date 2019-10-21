
#'Run EPIDEMIA model validation statistics
#'
#'This function takes a few more arguments than `epidemiar::run_epidemia()` to
#'generate statistics on model validation. The function will evaluate a number
#'of weeks (`total_timesteps`) starting from a specified week (`date_start`) and
#'will look at the n-week ahead forecast (1 to `timesteps_ahead` number of
#'weeks) and compare the values to the observed number of cases. An optional
#'`reporting_lag` argument will censor the last known data back that number of
#'weeks. The validation statistics include Root Mean Squared Error (RMSE) and
#'Mean Absolute Error (MAE), and an R-squared staistic both in total and per
#'geographic grouping (if present).
#'
#'@param date_start Date to start testing for model validation.
#'@param total_timesteps Number of weeks from `week_start` to run validation
#'  tests.
#'@param timesteps_ahead Number of weeks for testing the n-week ahead forecasts.
#'  Results will be generated from 1-week ahead through `weeks_ahead` number of
#'  weeks.
#'@param reporting_lag Number of timesteps to simulate reporting lag. For
#'  instance, if you have weekly data, and a reporting_lag of 1 week, and are
#'  working with a timesteps_ahead of 1 week, then that is functional equivalent
#'  to reporting lag of 0, and timesteps_ahead of 2 weeks. I.e. You are
#'  forecasting next week, but you don't know this week's data yet, you only
#'  know last week's numbers.
#'@param per_timesteps When creating a timeseries of validation results, create
#'  sets with per_timesteps number of steps in each set. Should be a minimum of
#'  10 timesteps. Last grouping may not contain full set.
#'@param skill_test Logical parameter indicating whether or not to run
#'  validations also on two naïve models for a skill test comparison. The naïve
#'  models are "persistence": the last known value (case counts) carried
#'  forward, and "average week" where the predicted value is the average of that
#'  week of the year, as calculated from historical data.
#'@param epi_data See description in `run_epidemia()`.
#'@param casefield See description in `run_epidemia()`.
#'@param populationfield See description in `run_epidemia()`.
#'@param groupfield See description in `run_epidemia()`.
#'@param week_type See description in `run_epidemia()`.
#'@param report_period The number of weeks that the entire report will cover.
#'  The \code{report_period} minus \code{forecast_future} is the number of weeks
#'  of past (known) data that will be included. Overwritten to be `weeks_ahead`
#'  + 1 for validation runs.
#'@param ed_summary_period Overwritten to 1 for validation runs (no-op for no
#'  event detection during validation runs).
#'@param ed_method Overwritten to "none" for validation runs.
#'@param env_data See description in `run_epidemia()`.
#'@param obsfield See description in `run_epidemia()`.
#'@param valuefield See description in `run_epidemia()`.
#'@param forecast_future Number of future weeks from the end of the
#'  \code{epi_data} to produce forecasts, as in `run_epidemia()`, but
#'  overwritten as `weeks_ahead` for validation runs.
#'@param fc_control See description in `run_epidemia()`. Note,
#'  fc_control$value_type is overwritten as "cases" for validation runs.
#'@param env_ref_data See description in `run_epidemia()`.
#'@param env_info See description in `run_epidemia()`.
#'@param model_cached See description in `run_epidemia()`.
#'@param model_choice See description in `run_epidemia()`.
#'@param ... Accepts other arguments that are normally part of `run_epidemia()`,
#'  but ignored for validation runs. For example, `inc_per`, `ed_control`,
#'  `model_run`.
#'
#'
#'@return Returns a list of validation statistics. Statistics are calculated on
#'  the n-week ahead forecast and the actual observed case counts. Statistics
#'  returned are  Mean Absolute Error (MAE), Root Mean Squared Error (RMSE), and
#'  proportion of observed values that were inside the prediction interval
#'  (prop_interval). The first object `validation_overall` is the results
#'  overall, and `validation_grouping` is the results per geographic grouping.
#'
#'@export
#'


run_validation <- function(date_start = NULL,
                           total_timesteps = 12,
                           timesteps_ahead = 2,
                           reporting_lag = 0,
                           per_timesteps = 12,
                           skill_test = FALSE,
                           #for run_epidemia()
                           epi_data = NULL,
                           casefield = NULL,
                           populationfield = NULL,
                           groupfield = NULL,
                           week_type = c("ISO", "CDC"),
                           report_period = 3, #default is timesteps_ahead default + 1
                           ed_summary_period = 1, #0 causes errors, 1 and "none" is no-op equivalent
                           ed_method = "none",
                           env_data = NULL,
                           obsfield = NULL,
                           valuefield = NULL,
                           forecast_future = 2, #default same as weeks_ahead default
                           fc_control = NULL,
                           env_ref_data = NULL,
                           env_info = NULL,
                           model_cached = NULL,
                           model_choice = c("poisson-bam", "negbin"),
                           ...){
  #date_start: week to start reporting of results
  #total_timesteps: number of weeks forward from week_start to gather test results
  #timesteps_ahead: calculate stats on 1 to n week ahead predictions

  #this means that the start of calculations will be date_start minus timesteps_ahead # of weeks
  #then trimmed at the end to start at date_start.

  # Non-standard evaluation quosures ----------------------------------------

  # dplyr programming steps for passing of field names
  quo_casefield <- rlang::enquo(casefield)
  quo_popfield <- rlang::enquo(populationfield)
  quo_groupfield <- rlang::enquo(groupfield)
  quo_obsfield <- rlang::enquo(obsfield)
  quo_valuefield <- rlang::enquo(valuefield)

  #Note: if field name does not exist in any dataset, enquo() will throw an error.


  # Adjust parameters for validation runs -----------------------------------

  #Assumed that run_epidemia() parameters just copied and pasted, so adjust for validation
  #new lengths
  forecast_future <- timesteps_ahead
  report_period <- forecast_future + 1
  #no event detection
  ed_summary_period <- 1
  ed_method <- "none"
  #report out in CASES for validation
  fc_control$value_type <-  "cases"

  #for params accepted by run_epidemia, but are meaningless for validation runs
  # e.g. `inc_per`, `ed_control`, `model_run`
  #captured, but then do nothing with them
  # Also used for hidden raw_data argument for testing/development
  dots <- list(...)


  # All loop prep ------------------------------------------------------

  #Set up for looping
  #preserve full data
  epi_data_orig <- epi_data
  env_data_orig <- env_data

  #Pull obs from original
  # Will have extra dates, but will be trimmed back to user requested dates later
  obs_only <- epi_data_orig %>%
    dplyr::select(!!quo_groupfield, obs_date, !!quo_casefield) %>%
    #rename observation
    dplyr::rename(obs := !!quo_name(quo_casefield))


  #Skill test loop set up
  if (skill_test == TRUE){
    models_to_run = c(model_choice, "naive-persistence", "naive-averageweek")
  } else {
    models_to_run = c(model_choice)
  }

  # Skill test loop ---------------------------------------------------------

  #skill test collection
  all_validations <- vector("list", length = length(models_to_run))
  #add names
  names(all_validations) <- models_to_run

  #model loop
  for (m in seq_along(models_to_run)){

    this_model <- models_to_run[m]

    #If naive-averageweek, timesteps_ahead is meaningless, just use 1 and change to NA later
    if (this_model == "naive-averageweek"){
      this_timesteps_ahead <- 1
      this_forecast_future <- this_timesteps_ahead
      this_report_period <- this_forecast_future + 1
    } else {
      this_timesteps_ahead <- timesteps_ahead
      this_forecast_future <- this_timesteps_ahead
      this_report_period <- this_forecast_future + 1
    }

    # Week loop ---------------------------------------------------------------

    #Create list of dates
    #the start of calculations will be date_start minus timesteps_ahead  # of weeks
    date_list <- date_start + lubridate::weeks(-this_timesteps_ahead:(total_timesteps-1))

    #output will be list of dataframes (forecasts) until we collapse later
    fcs_list <- vector("list", length = length(date_list))

    #loop
    for (i in seq_along(date_list)){
      this_dt <- date_list[i]

      message("Validation run - date: ", this_dt) # for testing for now

      #set up data
      #censoring as appropriate with reporting_lag (used for both epi & env)
      censor_date <- this_dt - lubridate::weeks(reporting_lag)
      epi_data <- epi_data_orig %>%
        dplyr::filter(obs_date <= censor_date)
      env_data <- env_data_orig %>%
        dplyr::filter(obs_date <= censor_date)

      #run_epidemia
      #passing quosures, which will have an escape built into run_epidemia()
      reportdata <- run_epidemia(epi_data = epi_data,
                                 casefield = quo_casefield,
                                 populationfield = quo_popfield,
                                 inc_per = inc_per,
                                 groupfield = quo_groupfield,
                                 week_type = "ISO",
                                 report_period = this_report_period, #this
                                 ed_summary_period = ed_summary_period,
                                 ed_method = ed_method,
                                 ed_control = ed_control,
                                 env_data = env_data,
                                 obsfield = quo_obsfield,
                                 valuefield = quo_valuefield,
                                 forecast_future = this_forecast_future, #this
                                 fc_control = fc_control,
                                 env_ref_data = env_ref_data,
                                 env_info = env_info,
                                 model_cached = model_cached,
                                 model_choice = this_model) ##models_to_run



      #pull needed and reformat
      fcs_list[[i]] <- reportdata$modeling_results_data %>%
        #get forecasts only
        dplyr::filter(series == "fc") %>%
        #get base date of report ('current date' in relation to forecast)
        dplyr::mutate(base_date = this_dt,
                      #due to reporting lag, date of last known data
                      last_known_date = censor_date,
                      #how many weeks ahead is the prediction
                      timestep_ahead = difftime(obs_date, base_date) %>% as.numeric(units = "weeks")) %>%
        #don't need 0 week predictions (same week)
        dplyr::filter(timestep_ahead > 0)

      # #also need to grab observations per run, will merge later
      # fcs_list[[i]] <- dplyr::bind_rows(fcs_list[[i]],
      #                                   reportdata$modeling_results_data %>%
      #                                     dplyr::filter(series == "obs"))


    } #end timestep loop

    #have list of dataframes
    #collapse/bindrows
    fcs_only <- dplyr::bind_rows(fcs_list) %>%
      #nicely arrange
      dplyr::arrange(!!quo_groupfield, obs_date, timestep_ahead)


    # #Split forecasts and observations to join instead
    # fcs_only <- fcs_obs %>%
    #   dplyr::filter(series == "fc")
    # obs_only <- fcs_obs %>%
    #   dplyr::filter(series == "obs") %>%
    #   dplyr::select(!!quo_groupfield, obs_date, value) %>%
    #   dplyr::rename(obs = value)

    #join
    fc_join <- fcs_only %>%
      dplyr::left_join(obs_only,
                       #NSE fun
                       by = rlang::set_names(c(rlang::quo_name(quo_groupfield),
                                               "obs_date"),
                                             c(rlang::quo_name(quo_groupfield),
                                               "obs_date")))

    #Filter to report weeks (trim off edges gathered b/c of weeks_ahead)
    fc_trim <- fc_join %>%
      dplyr::filter(between(obs_date,
                            date_start,
                            date_start + lubridate::weeks(total_timesteps-1))) %>%
    #Add column for showing reporting_lag
      dplyr::mutate(reporting_lag = reporting_lag)

    #timestep_ahead is meaningless for average week.
    # NA may cause unexpected results with grouping, so replace with 0
    if (this_model == "naive-averageweek"){
      fc_trim <- fc_trim %>%
        dplyr::mutate(timestep_ahead = 0)
    }

    ## Calculate statistics
    val_results <- calc_val_stats(fc_trim, quo_groupfield, per_timesteps, dots)

    #add results to list by name
    all_validations[[this_model]] <- val_results

  } #end model loop


  message("Validation run finished.")
  all_validations

} #end run validation



#' Calculate validation statistics from forecast results.
#'
#' Helper function to calculate the validation statistics from each model run.
#' Mean Absolute Error (MAE), Root Mean Square Error (RMSE), Proportion of
#' observations in in prediction interval, and R^2. Calculates it both at a
#' global model level per timestep ahead, and at a geographical grouping level
#' per timestep ahead. Also calculates a timeseries of evaluation metrics at
#' every per_timesteps number of timesteps per grouping (if applicable) and
#' timestep_ahead.
#'
#' @param fc_trim The forecast results of one model type, combined with observed
#'   values, trimmed to user requested date range.
#' @param dots The non-required arguments to run_validation() for developer
#'   testing.
#'
#' @return
#'
calc_val_stats <- function(fc_trim, quo_groupfield, per_timesteps, dots){
  # MAE: mean(|obs - pred|)
  # RMSE: sqrt(mean((obs - pred)^2))
  # Proportion in Interval: 1/T if inside, summed. Over all non-NA entries.
  # R2 (R^2): 1 - SSE/TSS.  SSE = sum((obs-pred)^2). TSS = sum((obs - mean(obs))^2).
  #     B/c involves mean of group of observations, must be calculated after grouping

  #per line stats
  fc_stats <- fc_trim %>%
    dplyr::mutate(diff = obs - value,
                  absdiff = abs(diff),
                  diffsq = diff ^ 2,
                  predinterval = ifelse(obs >= lower & obs <= upper, TRUE, FALSE))


  #overall timestep_ahead
  validation_overall <- fc_stats %>%
    dplyr::group_by(timestep_ahead) %>%
    #Now calc TSS part of R2
    dplyr::mutate(meanobs = mean(obs),
                  total_squares = (obs - meanobs)^2) %>%
    #stat calc
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     RMSE = sqrt(MSE),
                     prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)),
                     SSE = sum(diffsq, na.rm = TRUE),
                     TSS = sum(total_squares, na.rm = TRUE),
                     R2 = 1 - (SSE/TSS))


  #overall timestep_ahead by grouping
  validation_grouping <- fc_stats %>%
    dplyr::group_by(!!quo_groupfield, timestep_ahead) %>%
    #Now calc TSS part of R2
    dplyr::mutate(meanobs = mean(obs),
                  total_squares = (obs - meanobs)^2) %>%
    #stat calc
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     RMSE = sqrt(MSE),
                     prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)),
                     SSE = sum(diffsq, na.rm = TRUE),
                     TSS = sum(total_squares, na.rm = TRUE),
                     R2 = 1 - (SSE/TSS))


  #timeseries calculations
  # minimum of ~10 timesteps per summary
  fc_timeseries <- fc_stats %>%
    dplyr::group_by(!!quo_groupfield, timestep_ahead) %>%
    #create group every per_timesteps
    mutate(ts_group = gl(n(), per_timesteps, n())) %>%
    #group by new timeseries group
    group_by(ts_group, add = TRUE) %>%
    #count number of timesteps in group (last group may not have complete set)
    mutate(n_steps = length(unique(obs_date)))

  validation_timeseries <- fc_timeseries %>%
    dplyr::group_by(!!quo_groupfield, timestep_ahead, ts_group) %>%
    #Now calc TSS part of R2
    dplyr::mutate(meanobs = mean(obs),
                  total_squares = (obs - meanobs)^2) %>%
    #stat calc
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     RMSE = sqrt(MSE),
                     prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)),
                     SSE = sum(diffsq, na.rm = TRUE),
                     TSS = sum(total_squares, na.rm = TRUE),
                     R2 = 1 - (SSE/TSS))




  #return all
  # and raw data with hidden option
  #possibly make "time series" version for clean full data table
  if (!is.null(dots[['raw_data']])){
    if (dots[['raw_data']] == TRUE){
      val_stats <- create_named_list(validation_overall, validation_grouping, validation_timeseries, raw_stats = fc_stats)
    } #end raw data TRUE
  } else {
    #normal return with just results
    val_stats <- create_named_list(validation_overall, validation_grouping, validation_timeseries)
  }
} #end calc_val_stats()




#' View overall model validation statistics
#'
#' Small function to pull out just overall validation statistics.
#'
#' @param validations The set of validation statistics produced by
#'   run_validation()
#'
#' @return A list of tibbles containing only the model overall statistics (and
#'   not including the geographic grouping results, if present).
#'
#' @export
#'
view_overall_validations <- function(validations){
  lapply(validations, `[[`, "validation_overall")
}

#' Save overall model validation statistics
#'
#' Small function to pull out just overall validation statistics and save to
#' csv.
#'
#' @param validations The set of validation statistics produced by
#'   run_validation()
#' @param save_file File name to save results into csv format
#'
#' @return A csv file containing only the model overall statistics (and not
#'   including the geographic grouping results, if present).
#'
#' @export
#'
save_overall_validations <- function(validations, save_file){
  lapply(validations, `[[`, "validation_overall") %>%
    bind_rows(.id = "model") %>%
    write_csv(save_file)
}

#' Save geographic grouping model validation statistics
#'
#' Small function to pull out validation statistics per geographic grouping and
#' save to csv.
#'
#' @param validations The set of validation statistics produced by
#'   run_validation()
#' @param save_file File name to save results into csv format
#'
#' @return A csv file containing the model validation statistics for the
#'   geographic grouping results.
#'
#' @export
#'
save_geog_validations <- function(validations, save_file){
  lapply(validations, `[[`, "validation_grouping") %>%
    bind_rows(.id = "model") %>%
    write_csv(save_file)
}
