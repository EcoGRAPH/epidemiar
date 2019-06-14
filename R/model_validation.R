#validation function

run_validation <- function(week_start = NULL,
                           total_weeks = 12,
                           weeks_ahead = 2,
                           #for run_epidemia()
                           epi_data = NULL,
                           casefield = NULL,
                           populationfield = NULL,
                           inc_per = 1000,
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
                           model_choice = c("poisson-bam", "negbin")){
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


  # Week loop and prep ------------------------------------------------------

  #Set up for looping
  #preserve full data
  epi_data_orig <- epi_data
  env_data_orig <- env_data
  #create list of dates
  #the start of calculations will be week_start minus weeks_ahead # of weeks
  week_list <- week_start + lubridate::weeks(-weeks_ahead:total_weeks)

  #output will be list of dataframes until we collapse later
  fcs_list <- vector("list", length = length(week_list))

  #loop
  for (i in seq_along(week_list)){
    this_wk <- week_list[i]

    print(paste0("Validation date running: ", this_wk)) # for testing for now

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
                               model_choice = model_choice)



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
  fc_join <- fc_join %>%
    dplyr::mutate(diff = value - obs,
                  absdiff = abs(diff),
                  diffsq = diff ^ 2,
                  predinterval = ifelse(obs >= lower & obs <= upper, TRUE, FALSE))


  #overall week_ahead
  validation_overall <- fc_join %>%
    dplyr::group_by(week_ahead) %>%
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval)))


  #overall week_ahead by grouping
  validation_grouping <- fc_join %>%
    dplyr::group_by(week_ahead, !!quo_groupfield) %>%
    dplyr::summarize(MAE = mean(absdiff, na.rm = TRUE),
                     MSE = mean(diffsq, na.rm = TRUE),
                     prop_interval = sum(predinterval, na.rm = TRUE) / sum(!is.na(predinterval))) %>%
    dplyr::arrange(!!quo_groupfield, week_ahead)


  #return both
  val_results <- create_named_list(validation_overall, validation_grouping)


} #end run validation
