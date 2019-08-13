
#'Run EPIDEMIA model validation statistics
#'
#'This function takes a few more arguments than `epidemiar::run_epidemia()` to
#'generate statistics on model validation. The function will evaluate a number
#'of weeks (`total_weeks`) starting from a specified week (`week_start`) and
#'will look at the n-week ahead forecast (1 to `weeks_ahead` number of weeks)
#'and compare the values to the observed number of cases. The validation
#'statistics include Mean Squared Error (MSE) and Mean Absolute Error (MAE),
#'both in total and per geographic grouping (if present).
#'
#'@param week_start Date to start testing for model validation.
#'@param total_weeks Number of weeks from `week_start` to run validation tests.
#'@param weeks_ahead Number of weeks for testing the n-week ahead forecasts.
#'  Results will be generated from 1-week ahead through `weeks_ahead` number of
#'  weeks.
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
#'  returned are Mean Squared Error (MSE), Mean Absolute Error (MAE), and
#'  proportion of observed values that were inside the prediction interval
#'  (prop_interval). The first object `validation_overall` is the results
#'  overall, and `validation_grouping` is the results per geographic grouping.
#'
#'@export
#'


run_validation <- function(week_start = NULL,
                           total_weeks = 12,
                           weeks_ahead = 2,
                           skill_test = FALSE,
                           #for run_epidemia()
                           epi_data = NULL,
                           casefield = NULL,
                           populationfield = NULL,
                           groupfield = NULL,
                           week_type = c("ISO", "CDC"),
                           report_period = 3, #default is weeks_ahead default + 1
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
  #week_start: week to start reporting of results
  #total_weeks: number of weeks forward from week_start to gather test results
  #weeks_ahead: calculate stats on 1 to n week ahead predictions

  #this means that the start of calculations will be week_start minus weeks_ahead # of weeks
  #then trimmed at the end to start at week_start.

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
  forecast_future <- weeks_ahead
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

    # Week loop ---------------------------------------------------------------

    #create list of dates
    #the start of calculations will be week_start minus weeks_ahead # of weeks
    week_list <- week_start + lubridate::weeks(-weeks_ahead:(total_weeks-1))

    #output will be list of dataframes until we collapse later
    fcs_list <- vector("list", length = length(week_list))

    #loop
    for (i in seq_along(week_list)){
      this_wk <- week_list[i]

      message("Validation run - date: ", this_wk) # for testing for now

      #set up data
      epi_data <- epi_data_orig %>%
        dplyr::filter(obs_date <= this_wk)
      env_data <- env_data_orig %>%
        dplyr::filter(obs_date <= this_wk)

      #run_epidemia
      #passing quosures, which will have an escape built into run_epidemia()
      reportdata <- run_epidemia(epi_data = epi_data,
                                 casefield = quo_casefield,
                                 populationfield = quo_popfield,
                                 inc_per = inc_per,
                                 groupfield = quo_groupfield,
                                 week_type = "ISO",
                                 report_period = report_period,
                                 ed_summary_period = ed_summary_period,
                                 ed_method = ed_method,
                                 ed_control = ed_control,
                                 env_data = env_data,
                                 obsfield = quo_obsfield,
                                 valuefield = quo_valuefield,
                                 forecast_future = forecast_future,
                                 fc_control = fc_control,
                                 env_ref_data = env_ref_data,
                                 env_info = env_info,
                                 model_cached = model_cached,
                                 model_choice = this_model) ##models_to_run



      #pull needed and reformat
      fcs_list[[i]] <- reportdata$modeling_results_data %>%
        #get forecasts only
        dplyr::filter(series == "fc") %>%
        #get date of last known epi data
        dplyr::mutate(known_week = this_wk,
                      #how many weeks ahead is the prediction
                      week_ahead = difftime(obs_date, known_week) %>% as.numeric(units = "weeks")) %>%
        #don't need 0 week predictions (same week)
        dplyr::filter(week_ahead > 0)

      #also need to grab observations per run, will merge later
      fcs_list[[i]] <- dplyr::bind_rows(fcs_list[[i]],
                                        reportdata$modeling_results_data %>%
                                          dplyr::filter(series == "obs"))


    } #end week loop

    #have list of dataframes
    #collapse/bindrows
    fcs_obs <- dplyr::bind_rows(fcs_list)

    #Split forecasts and observations to join instead
    fcs_only <- fcs_obs %>%
      dplyr::filter(series == "fc")

    obs_only <- fcs_obs %>%
      dplyr::filter(series == "obs") %>%
      dplyr::select(!!quo_groupfield, obs_date, value) %>%
      dplyr::rename(obs = value)

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
                            week_start,
                            week_start + lubridate::weeks(total_weeks-1)))


    ## Calculate statistics
    #per line stats
    fc_stats <- fc_trim %>%
      dplyr::mutate(diff = value - obs,
                    absdiff = abs(diff),
                    diffsq = diff ^ 2,
                    predinterval = ifelse(obs >= lower & obs <= upper, TRUE, FALSE))


    #overall week_ahead
    validation_overall <- fc_stats %>%
      dplyr::group_by(week_ahead) %>%
      dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                       MSE = mean(diffsq, na.rm = TRUE),
                       prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)))


    #overall week_ahead by grouping
    validation_grouping <- fc_stats %>%
      dplyr::group_by(week_ahead, !!quo_groupfield) %>%
      dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                       MSE = mean(diffsq, na.rm = TRUE),
                       prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval))) %>%
      dplyr::arrange(!!quo_groupfield, week_ahead)


    #return both
    # and raw data with hidden option
    if (!is.null(dots[['raw_data']])){
      if (dots[['raw_data']] == TRUE){
        val_results <- create_named_list(validation_overall, validation_grouping, raw_stats = fc_stats)
      } #end raw data TRUE
    } else {
      #normal return with just results
      val_results <- create_named_list(validation_overall, validation_grouping)
    }

    #add results to list by name
    all_validations[[this_model]] <- val_results

  }#end model loop


  message("Validation run finished.")
  all_validations

} #end run validation


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
