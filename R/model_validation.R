
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
#'  a moving window with per_timesteps width number of time points. Should be a
#'  minimum of 10 timesteps.
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
#'@return Returns a nested list of validation statistics. Statistics are
#'  calculated on the n-week ahead forecast and the actual observed case counts.
#'  Statistics returned are  Mean Absolute Error (MAE), Root Mean Squared Error
#'  (RMSE). The first object is `skill_scores`, which contains
#'  `skill_overall` and `skill_grouping`. The second list is `validations`,
#'  which contains lists per model run (the forecast model and then optionally
#'  the naive models). Within each, `validation_overall` is the results overall,
#'  and `validation_grouping` is the results per geographic grouping.
#'
#'@export
#'
run_validation <- function(date_start = NULL,
                           total_timesteps = 26,
                           timesteps_ahead = 2,
                           reporting_lag = 0,
                           per_timesteps = 12,
                           skill_test = TRUE,
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
                           forecast_future = 2, #default same as timesteps_ahead default
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
  forecast_future <- timesteps_ahead + reporting_lag
  report_period <- forecast_future + 1
  #no event detection
  ed_summary_period <- 1
  ed_method <- "none"
  #report out in CASES for validation
  fc_control$value_type <- "cases"

  #for params accepted by run_epidemia, but are meaningless for validation runs
  # e.g. `inc_per`, `ed_control`, `model_run`
  #captured, but then do nothing with them
  # Also used for hidden raw_data argument for testing/development
  dots <- list(...)

  #Create parameter metadata
  metadata <- create_named_list(date_start,
                                total_timesteps,
                                timesteps_ahead,
                                reporting_lag,
                                per_timesteps,
                                skill_test,
                                casefield = quo_name(quo_casefield),
                                date_created = Sys.Date())


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

    #If naive-averageweek, timesteps_ahead is meaningless, just use 1
    if (this_model == "naive-averageweek"){
      this_timesteps_ahead <- 1
      this_forecast_future <- this_timesteps_ahead
      this_report_period <- this_forecast_future + 1
    } else {
      #use modified forecast_future which is timesteps_ahead + reporting_lag
      this_timesteps_ahead <- forecast_future #timesteps_ahead
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
      #censoring as appropriate
      #reporting_lag will be handled with offset timesteps
      epi_data <- epi_data_orig %>%
        dplyr::filter(obs_date <= this_dt)
      env_data <- env_data_orig %>%
        dplyr::filter(obs_date <= this_dt)

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
        dplyr::mutate(preadj_date = this_dt,
                      #how many weeks ahead is the prediction (not adjusting for reporting lag yet)
                      timestep_ahead_orig = difftime(obs_date, preadj_date) %>%
                                          as.numeric(units = "weeks")) %>%
        #don't need 0 week predictions (same week)
        dplyr::filter(timestep_ahead_orig > 0)


    } #end timestep loop

    #have list of dataframes
    #collapse/bindrows
    fcs_only <- dplyr::bind_rows(fcs_list) %>%
      #nicely arrange
      dplyr::arrange(!!quo_groupfield, timestep_ahead_orig, obs_date)


    #join
    fc_join <- fcs_only %>%
      dplyr::left_join(obs_only,
                       #NSE fun
                       by = rlang::set_names(c(rlang::quo_name(quo_groupfield),
                                               "obs_date"),
                                             c(rlang::quo_name(quo_groupfield),
                                               "obs_date")))

    #make all the reporting_lag adjustments
    # basically, we ran extra forecast future steps
    # so we now can simply shift everything backwards except for averageweek
    if (this_model == "naive-averageweek"){
      fc_join <- fc_join %>%
        dplyr::mutate(run_date = preadj_date,
                      #timestep_ahead is meaningless for average week.
                      # NA may cause unexpected results with grouping, so replace with 0
                      timestep_ahead = 0,
                      #Add column for showing reporting_lag
                      reporting_lag = reporting_lag)
    } else {
      fc_join <- fc_join %>%
        dplyr::mutate(run_date = preadj_date - lubridate::weeks(reporting_lag),
                      timestep_ahead = timestep_ahead_orig - reporting_lag,
                      #Add column for showing reporting_lag
                      reporting_lag = reporting_lag) %>%
        #filter out the timesteps that are now less than 1 step
        dplyr::filter(timestep_ahead > 0)
    }


    #Filter to report weeks (trim off edges gathered b/c of weeks_ahead, etc.)
    fc_trim <- fc_join %>%
      dplyr::filter(between(obs_date,
                            date_start,
                            date_start + lubridate::weeks(total_timesteps-1)))


    ## Calculate statistics
    val_results <- calc_val_stats(fc_trim, quo_groupfield, per_timesteps, dots)

    #add results to list by name
    all_validations[[this_model]] <- val_results

  } #end model loop



  #Get skill test list of results
  if (skill_test == TRUE){
    #calc skill comparison statistics
    skill_overall <- calc_skill(get_overall_validations(all_validations))
    skill_grouping <- calc_skill(get_group_validations(all_validations), quo_groupfield)
    skill_scores <- create_named_list(skill_overall, skill_grouping)

    val_return <- create_named_list(skill_scores, validations = all_validations, metadata)
  } else {
    #just the one model validation datasets
    val_return <- create_named_list(all_validations, metadata)
  }




  message("Validation run finished.")
  val_return

} #end run validation



#'Calculate validation statistics from forecast results.
#'
#'Helper function to calculate the validation statistics from each model run.
#'Mean Absolute Error (MAE), Root Mean Square Error (RMSE), Proportion of
#'observations in in prediction interval, and R^2. Calculates it both at a
#'global model level per timestep ahead, and at a geographical grouping level
#'per timestep ahead. Also calculates a timeseries of evaluation metrics at
#'every per_timesteps number of timesteps per grouping (if applicable) and
#'timestep_ahead.
#'
#'@param fc_trim The forecast results of one model type, combined with observed
#'  values, trimmed to user requested date range.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_validation()/run_epidemia().
#'@param per_timesteps When creating a timeseries of validation results, create
#'  a moving window with per_timesteps width number of time points. Should be a
#'  minimum of 10 timesteps.
#'@param dots The non-required arguments to run_validation() for developer
#'  testing.
#'
#'@return A named list of validation statistic results: validation_overall,
#'  validation_grouping, validation_timeseries
#'
calc_val_stats <- function(fc_trim, quo_groupfield, per_timesteps, dots){
  # MAE: mean(|obs - pred|)
  # RMSE: sqrt(mean((obs - pred)^2))
  # R2 (R^2): 1 - SSE/TSS.  SSE = sum((obs-pred)^2). TSS = sum((obs - mean(obs))^2).
  #     B/c involves mean of group of observations, must be calculated after grouping

  #Removed
  # Proportion in Interval: 1/T if inside, summed. Over all non-NA entries.

  #per line stats
  fc_stats <- fc_trim %>%
    dplyr::mutate(diff = obs - value,
                  absdiff = abs(diff),
                  diffsq = diff ^ 2)
                  #,predinterval = ifelse(obs >= lower & obs <= upper, TRUE, FALSE))


  #overall timestep_ahead
  validation_overall <- fc_stats %>%
    dplyr::group_by(timestep_ahead) %>%
    #Now calc TSS part of R2
    dplyr::mutate(meanobs = mean(obs),
                  total_squares = (obs - meanobs)^2) %>%
    #stat calc
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     #prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)),
                     SSE = sum(diffsq, na.rm = TRUE),
                     TSS = sum(total_squares, na.rm = TRUE)) %>%
    #and mutate for final calc
    dplyr::mutate(RMSE = sqrt(MSE),
                  R2 = 1 - (SSE/TSS)) %>%
    #drop unneeded columns
    dplyr::select(-SSE, -TSS, -MSE)



  #overall timestep_ahead by grouping
  validation_grouping <- fc_stats %>%
    dplyr::group_by(!!quo_groupfield, timestep_ahead) %>%
    #Now calc TSS part of R2
    dplyr::mutate(meanobs = mean(obs),
                  total_squares = (obs - meanobs)^2) %>%
    #stat calc
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     #prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)),
                     SSE = sum(diffsq, na.rm = TRUE),
                     TSS = sum(total_squares, na.rm = TRUE)) %>%
    #and mutate for final calc
    dplyr::mutate(RMSE = sqrt(MSE),
                  R2 = 1 - (SSE/TSS)) %>%
    #drop unneeded columns
    dplyr::select(-SSE, -TSS, -MSE)



  #timeseries calculations
  # minimum of ~10 timesteps per summary
  # ROLLING window
  validation_timeseries <- fc_stats %>%
    dplyr::group_by(!!quo_groupfield, timestep_ahead) %>%
    #rollapply for get mean of obs
    dplyr::mutate(meanobs = zoo::rollmeanr(x = obs,
                                           k = per_timesteps,
                                           fill = NA),
                  total_squares = (obs - meanobs)^2,
                  MAE = zoo::rollmeanr(x = absdiff,
                                       k = per_timesteps,
                                       fill = NA),
                  MSE = zoo::rollmeanr(x = diffsq,
                                       k = per_timesteps,
                                       fill = NA),
                  RMSE = sqrt(MSE),
                  #prop_interval = zoo::rollsumr(x = predinterval,
                  #                              k = per_timesteps,
                  #                              fill = NA) /
                  #  zoo::rollsumr(x = !is.na(predinterval),
                  #                k = per_timesteps,
                  #                fill = NA),
                  SSE = zoo::rollsumr(x = diffsq,
                                      k = per_timesteps,
                                      fill = NA),
                  TSS = zoo::rollsumr(x = total_squares,
                                      k = per_timesteps,
                                      fill = NA),
                  R2 = 1 - (SSE/TSS)) %>%
    #rename columns to be clearer
    dplyr::rename(forecast = value,
                  observed = obs) %>%
    # drop unneeded columns
    dplyr::select(-series, -preadj_date, -timestep_ahead_orig, -run_date,
                  -diff, -absdiff, -diffsq, -meanobs, -total_squares, -MSE, -SSE, -TSS) %>%
    # for now, drop R2 until can figure out how to include better
    dplyr::select(-R2)



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


#' Get overall model validation statistics
#'
#' Small function to pull out just overall validation statistics.
#'
#' @param validations The set of validation statistics produced by
#'   run_validation() - only the list of validation data sets, not including the skill metrics.
#'
#' @return A list of tibbles containing only the model overall statistics (and
#'   not including the geographic grouping results, if present).
#'
#' @export
#'
get_overall_validations <- function(validations){
  lapply(validations, `[[`, "validation_overall")
}


#' Get geographic grouping model validation statistics
#'
#' Small function to pull out just the geographic grouping validation statistics.
#'
#' @param validations The set of validation statistics produced by
#'   run_validation() - only the list of validation data sets, not including the skill metrics.
#'
#' @return A list of tibbles containing only the model geographic grouping statistics.
#'
#' @export
#'
get_group_validations <- function(validations){
  lapply(validations, `[[`, "validation_grouping")
}



#' Calculate model skill comparison statistics
#'
#' Helper function to calculate the relative improvement of the forecast over the specified naive model.
#' Skill score = (score_fc - score_naive) / (score_perfect - score_naive)
#' Skill metric has an upper bound of 1. No improvement is 0. Lower bound depends on statistic.
#'
#'@param fc_stat The forecast model statistic value.
#'@param naive_stat The naive model statistic value (same statistic as forecast model).
#'@param perfect_stat The value of a perfect score for that stastistic.
#'
#'@return Skill score: the relative improvement the forecast model has over the naive model.
#'
#'@export
#'
calc_skill_stat <- function(fc_stat, naive_stat, perfect_stat){
  skill_stat <- (fc_stat - naive_stat) / (perfect_stat - naive_stat)
}


#' Calculate the forecast model skill score compared to the naive model predictions.
#'
#'@param val_list A list of 3 datasets of validation results: the first is the forecast model, the following two are the naive model results, as created by binding the results of calc_val_stats() in run_validation().
#'@param grp Optional inclusion of quo_groupfield when calculating skill scores by groupfield.
#'
#'@return Single dataset with skill scores of the main forecast model against each of the naive models, per timestep ahead, and optionally, per geographic grouping
#'
calc_skill <- function(val_list, grp = NULL){

  #separate out, rename columns, and join/crossing
  val_fc <- val_list[[1]] %>%
    dplyr::rename(fc_MAE = MAE,
                  fc_RMSE = RMSE,
                  #fc_prop_interval = prop_interval,
                  fc_R2 = R2) %>%
    dplyr::select(group_cols(), timestep_ahead, starts_with("fc_"))

  val_np <- val_list$`naive-persistence` %>%
    dplyr::rename(np_MAE = MAE,
                  np_RMSE = RMSE,
                  #np_prop_interval = prop_interval,
                  np_R2 = R2) %>%
    dplyr::select(group_cols(), timestep_ahead, starts_with("np_"))

  val_naw <- val_list$`naive-averageweek` %>%
    rename(naw_MAE = MAE,
           naw_RMSE = RMSE,
           #naw_prop_interval = prop_interval,
           naw_R2 = R2) %>%
    #no timestep_ahead for average week, all same
    select(group_cols(), starts_with("naw_"))

  #appropriate joins
  if (is.null(grp)){
    #join together
    val_join <- val_fc %>%
      #join with persistence
      dplyr::left_join(val_np,
                       by = "timestep_ahead") %>%
      #join with average week (1 value to all timesteps ahead)
      tidyr::crossing(val_naw)
  } else {
    #else join with groupfield
    #join together
    val_join <- val_fc %>%
      #join with persistence
      dplyr::left_join(val_np,
                       #NSE fun
                       by = rlang::set_names(c(rlang::quo_name(grp),
                                               "timestep_ahead"),
                                             c(rlang::quo_name(grp),
                                               "timestep_ahead"))) %>%
      #join with average week (1 value to all timesteps ahead)
      dplyr::left_join(val_naw,
                       by = rlang::set_names(rlang::quo_name(grp),
                                             rlang::quo_name(grp)))
  } #end joinings

  #perfect skill metrics
  perfect_MAE <- 0
  perfect_RMSE <- 0
  #perfect_prop_interval <- 1
  perfect_R2 <- 1

  #calc skill metrics of fc model to each of naive models
  val_skill <- val_join %>%
    mutate(skill_MAE_persistence = calc_skill_stat(fc_MAE, np_MAE, perfect_MAE),
           skill_RMSE_persistence = calc_skill_stat(fc_RMSE, np_RMSE, perfect_RMSE),
           #skill_interval_persistence = calc_skill_stat(fc_prop_interval, np_prop_interval,
           #                                             perfect_prop_interval),
           skill_R2_persistence =  calc_skill_stat(fc_R2, np_R2, perfect_R2),
           skill_MAE_averageweek = calc_skill_stat(fc_MAE, naw_MAE, perfect_MAE),
           skill_RMSE_averageweek = calc_skill_stat(fc_RMSE, naw_RMSE, perfect_RMSE),
           #skill_interval_averageweek = calc_skill_stat(fc_prop_interval, naw_prop_interval,
           #                                             perfect_prop_interval),
           skill_R2_averageweek =  calc_skill_stat(fc_R2, naw_R2, perfect_R2)) %>%
    #select final stats only
    select(group_cols(), timestep_ahead, starts_with("skill_"))

  val_skill
}



#' Save overall model validation statistics
#'
#' Small function to pull out just overall validation statistics and save to
#' csv.
#'
#' @param validations The set of validation statistics produced by
#'   run_validation() - only the list of validation data sets, not including the skill metrics.
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
#'   run_validation() - only the list of validation data sets, not including the skill metrics.
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
