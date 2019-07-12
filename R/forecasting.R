# All run_epidemiar() subfunctions related to forecasting
## Forecasting

#' Runs the forecast modeling
#'
#'@param epi_data Epidemiological data with case numbers per week, with date
#'  field "obs_date".
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param inc_per Number for what unit of population the incidence should be
#'  reported in, e.g. incidence rate of 3 per 1000 people.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param env_data Daily environmental data for the same groupfields and date
#'  range as the epidemiological data. It may contain extra data (other
#'  districts or date ranges). The data must be in long format (one row for each
#'  date and environmental variable combination), and must start at absolutel
#'  minimum \code{laglen} (in \code{fc_control}) days before epi_data for
#'  forecasting.
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables
#'@param quo_valuefield Quosure of user given field name of the value of the
#'  environmental data variable observations.
#'@param env_variables alphabetical list of all unique environmental variables
#'  present in the original env_data dataset.
#'@param fc_control Parameters for forecasting, including which environmental
#'  variable to include and any geographic clusters.
#'@param env_ref_data Historical averages by week of year for environmental
#'  variables. Used in extended environmental data into the future for long
#'  forecast time, to calculate anomalies in early detection period, and to
#'  display on timeseries in reports.
#'@param env_info Lookup table for environmental data - reference creation
#'  method (e.g. sum or mean), report labels, etc.
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param week_type String indicating the standard (WHO ISO-8601 or CDC epi
#'  weeks) that the weeks of the year in epidemiological and environmental
#'  reference data use ["ISO" or "CDC"].
#'@param model_run TRUE/FALSE flag for whether to only generate the model
#'  regression object plus metadata. This model can be cached and used later on
#'  its own, skipping a large portion of the slow calculations for future runs.
#'@param model_cached The output of a previous model_run = TRUE run of
#'  run_epidemia() that produces a model (regression object) and metadata. The
#'  metadata will be used for input checking and validation. Using a prebuilt
#'  model saves on processing time, but will need to be updated periodically.
#'@param model_choice Critical argument to choose the type of model to generate.
#'  The options are versions that the EPIDEMIA team has used for forecasting.
#'  The first supported options is "poisson-gam" ("p") which is the original
#'  epidemiar model: a Poisson regression using bam (for large data GAMs), with
#'  a smoothed cyclical for seasonality. The default for fc_control$anom_env is
#'  TRUE for using the anomalies of environmental variables rather than their
#'  raw values. The second option is "negbin" ("n") which is a negative binomial
#'  regression using glm, with no external seasonality terms - letting the
#'  natural cyclical behavior of the environmental variables fill that role. The
#'  default for fc_control$anom_env is FALSE and uses the actual observation
#'  values in the modeling. The fc_control$anom_env can be overruled by the user
#'  providing a value, but this is not recommended unless you are doing
#'  comparisons.
#'
#'
#'@return Named list containing:
#'fc_epi: Full forecasted resulting dataset.
#'fc_res: The forecasted series in report format.
#'env_data_extd: Data set of the environmental data variables extended into the
#'  unknown/future.
#'env_variables_used: list of environmental variables that were used in the
#'  modeling (had to be both listed in model variables input file and present the
#'  env_data dataset)
#'env_dt_ranges: Date ranges of the input environmental data.
#'reg_obj: The regression object from modeling.
#'Unless model_run is TRUE, in which case only the regression object is returned.
#'
#'
run_forecast <- function(epi_data,
                         quo_popfield,
                         inc_per,
                         quo_groupfield,
                         groupings,
                         env_data,
                         quo_obsfield,
                         quo_valuefield,
                         env_variables,
                         fc_control,
                         env_ref_data,
                         env_info,
                         report_dates,
                         week_type,
                         model_run,
                         model_cached = NULL,
                         model_choice){

  message("Preparing for forecasting...")

  #set up default parallel processing number of cores to use number
  #if user-supplied, use that, otherwise create a default number
  #used in anomalize_env() and forecast_regression()
  if (!is.null(fc_control$ncores)) {
    ncores <- fc_control$ncores
  } else ncores <- max(parallel::detectCores(logical=FALSE) - 1, 1)

  # create the modeling variable
  # epi_data <- mutate(epi_data, logcase = log(cases_epidemiar + 1))
  epi_data <- dplyr::mutate(epi_data, modeledvar = floor(cases_epidemiar))

  # trim to the needed env variables as dictated by the model
  env_data <- pull_model_envvars(env_data, quo_obsfield, fc_control)
  #create alphabetical list of ONLY USED unique environmental variables
  env_variables_used <- dplyr::pull(env_data, !!quo_obsfield) %>% unique() %>% sort()

  # extract start & end dates for each variable for log file
  env_dt_ranges <- dplyr::group_by(env_data, !!quo_obsfield) %>%
    dplyr::summarize(start_dt = min(obs_date), end_dt = max(obs_date))

  # extend data into future, for future forecast portion
  env_data_extd <- extend_env_future(env_data,
                                     quo_groupfield,
                                     groupings,
                                     quo_obsfield,
                                     quo_valuefield,
                                     env_ref_data,
                                     env_info,
                                     env_variables_used,
                                     report_dates,
                                     week_type)
  epi_data_extd <- extend_epi_future(epi_data,
                                     quo_popfield,
                                     quo_groupfield,
                                     groupings,
                                     report_dates)

  # format the data for forecasting algorithm
  env_fc <- env_format_fc(env_data_extd,
                          quo_groupfield,
                          quo_obsfield)
  epi_fc <- epi_format_fc(epi_data_extd,
                          quo_groupfield,
                          fc_control)

  # anomalizing the environ data, if requested.
  #For backwards compatibility, poisson-gam default is TRUE.
  #Default for negbin is FALSE
  if (is.null(fc_control[["anom_env"]])){
    #fc_control$anom_env <- TRUE
    fc_control$anom_env <- dplyr::case_when(
      model_choice == "poisson-gam" ~ TRUE,
      model_choice == "negbin" ~ FALSE,
      model_choice == "null-persistence" ~ FALSE,
      model_choice == "null-weekaverage" ~ FALSE,
      #should never occur, but if it does, default to FALSE
      TRUE ~ FALSE)
  }
  if (fc_control$anom_env){
    message("Anomalizing the environmental variables...")
    env_fc <- anomalize_env(env_fc,
                            quo_groupfield,
                            env_variables_used,
                            ncores)
  }

  # create the lags
  epi_lag <- lag_environ_to_epi(epi_fc,
                                quo_groupfield,
                                groupings,
                                env_fc,
                                env_variables_used,
                                laglen = fc_control$lag_length)

  # add week of year, needed for null-weekaverage model
  # here with week_type, else need to pass that in to further down functions
  epi_lag <- add_datefields(epi_lag, week_type)


  # If only model_run, then return to run_epidemia() here
  if (model_run){
    model_run_result <- forecast_regression(epi_lag,
                                            quo_groupfield,
                                            groupings,
                                            env_variables_used,
                                            report_dates,
                                            req_date = report_dates$full$max,
                                            ncores,
                                            fit_freq = "once",
                                            model_run,
                                            model_cached,
                                            model_choice,
                                            theta = fc_control$theta)

    model_run_only <- create_named_list(env_variables_used,
                                        env_dt_ranges,
                                        reg_obj = model_run_result)
    return(model_run_only)
  }


  #Split regression call depending on {once|week} model fit frequency
  # default "once"
  if (!is.null(fc_control$fit_freq)) {
    fit_freq <- fc_control$fit_freq
  } else fit_freq <- "once"

  if (fit_freq == "once"){
    message("Generating forecasts...")
    #for single fit, call with last week (and subfunction has switch to return all)
    forereg_return <- forecast_regression(epi_lag,
                                          quo_groupfield,
                                          groupings,
                                          env_variables_used,
                                          report_dates,
                                          req_date = report_dates$full$max,
                                          ncores,
                                          fit_freq,
                                          model_run,
                                          model_cached,
                                          model_choice,
                                          theta = fc_control$theta)
    preds_catch <- forereg_return$date_preds
    reg_obj <- forereg_return$regress

  } else if (fit_freq == "week") {
    # for each week of report, run forecast
    # initialize: prediction returns 4 columns
    preds_catch <- data.frame()
    #loop by week
    for (w in seq_along(report_dates$full$seq)){
      message("Forecasting week ", w, " starting at ", Sys.time())
      dt <- report_dates$full$seq[w]
      forereg_return <- forecast_regression(epi_lag,
                                            quo_groupfield,
                                            groupings,
                                            env_variables_used,
                                            report_dates,
                                            req_date = dt,
                                            ncores,
                                            fit_freq,
                                            model_run,
                                            model_cached,
                                            model_choice,
                                            theta = fc_control$theta)

      dt_preds <- forereg_return$date_preds
      preds_catch <- rbind(preds_catch, as.data.frame(dt_preds))

      #taking advantage that only result will be of the last loop through
      reg_obj <- forereg_return$regress
    }

  } else stop("Model fit frequency unknown") #shouldn't happen with default "once"


  # Interval calculation
  preds_catch <- preds_catch %>%
    dplyr::mutate(fc_cases = fit,
                  fc_cases_upr = fit+1.96*sqrt(fit),
                  fc_cases_lwr = fit-1.96*sqrt(fit))

  # extract fc series into report format
  fc_res <- preds_catch %>%
    dplyr::mutate(series = "fc",
                  value = dplyr::case_when(
                    #if reporting in case counts
                    fc_control$value_type == "cases" ~ fc_cases,
                    #if incidence
                    fc_control$value_type == "incidence" ~ fc_cases / !!quo_popfield * inc_per,
                    #otherwise
                    TRUE ~ NA_real_),
                  lab = "Forecast Trend",
                  upper = dplyr::case_when(
                    #if reporting in case counts
                    fc_control$value_type == "cases" ~ fc_cases_upr,
                    #if incidence
                    fc_control$value_type == "incidence" ~ fc_cases_upr / !!quo_popfield * inc_per,
                    #otherwise
                    TRUE ~ NA_real_),
                  lower = dplyr::case_when(
                    #if reporting in case counts
                    fc_control$value_type == "cases" ~ fc_cases_lwr,
                    #if incidence
                    fc_control$value_type == "incidence" ~ fc_cases_lwr / !!quo_popfield * inc_per,
                    #otherwise
                    TRUE ~ NA_real_)
                  #value = fc_cases / !!quo_popfield * inc_per,
                  #upper = fc_cases_upr / !!quo_popfield * inc_per,
                  #lower = fc_cases_lwr / !!quo_popfield * inc_per
    ) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  # return list with res and other needed items
  fc_res_full <- create_named_list(fc_epi = preds_catch,
                                   fc_res,
                                   env_data_extd,
                                   env_variables_used,
                                   env_dt_ranges,
                                   reg_obj)
}

#forecasting helper functions
# this creates a modified b-spline basis (which is a piecewise polynomial)

#' Truncates poly. Creates a modified b-spline basis.
#'
#' The modified basis splines are used to capture any long term trends per
#' geographic group.
#'
#'@param x Vector of weekly observation dates.
#'@param degree Degree passed to splines::bs().
#'@param maxobs Date of the last known value.
#'@param minobs Date of the first known value.
#'
#'@returna A modified b-spline basis with the last basis splines reversed and
#'  the second to last basis spline function removed (to reduce the edge effects
#'  of using splines).
#'
truncpoly <- function(x = NULL, degree = 6, maxobs = NULL, minobs = NULL){

  # Some of the later functions will convert date to type spline, so
  # it's best to go ahead and convert now. Left_join doesn't convert
  # dates to numeric on the fly.
  x <- as.numeric(x)

  # create frame to hold modified b-spline basis
  xdf  <- data.frame(x=x)

  # figure out where we apparently have data
  apparentminobs <- min(xdf$x, na.rm=TRUE)
  apparentmaxobs <- max(xdf$x, na.rm=TRUE)

  # figure out where the bspline basis will have support
  if (!is.null(minobs)) {

    actualminobs <- max(apparentminobs, minobs)

  } else { actualminobs <- apparentminobs }
  if (!is.null(maxobs)) {

    actualmaxobs <- min(apparentmaxobs, maxobs)

  } else { actualmaxobs <- apparentmaxobs }

  # create a full frame to hold this basis before truncation
  xdf2     <- data.frame(x = actualminobs:actualmaxobs)
  xdf2$bas <- splines::bs(x=xdf2$x, degree=degree)

  # reverse the last spline basis function
  xdf2$bas[,degree] <- rev(xdf2$bas[,degree])

  # delete the next to last spline basis function
  xdf2$bas <- xdf2$bas[,-(degree-1)]

  # merge with original frame
  xdf <- dplyr::left_join(xdf, xdf2, by="x")

  # make values 0 where we extend beyond actualmax/minobs
  for (colnum in 1:ncol(xdf$bas)) {

    xdf$bas[is.na(xdf$bas[,colnum]),colnum] <- 0

  }

  tempdf <- data.frame(xdf$bas)
  names(tempdf) <- paste("modbs_reserved_", names(tempdf), sep="")

  return(tempdf)

}

#' Pull only model env variables.
#'
#'@param env_data Daily environmental data for the same groupfields and date
#'  range as the epidemiological data. It may contain extra data (other
#'  districts or date ranges). The data must be in long format (one row for each
#'  date and environmental variable combination), and must start at absolutel
#'  minimum \code{laglen} (in \code{fc_control}) days before epi_data for
#'  forecasting.
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables
#'@param fc_control Parameters for forecasting, including which environmental
#'  variable to include and any geographic clusters.
#'
#'@return List of environmental variables that were used in the
#'  modeling (had to be both listed in model variables input file and present the
#'  env_data dataset).
#'
pull_model_envvars <- function(env_data, quo_obsfield, fc_control){

  #pull variables from model info input
  model_vars <- fc_control$env_vars %>% dplyr::pull(!!quo_obsfield)

  #filter env_data for those model_vars
  env_data <- env_data %>%
    dplyr::filter(!!quo_obsfield %in% model_vars)
}

#' Extend environmental data into the future.
#'
#'@param env_data Daily environmental data for the same groupfields and date
#'  range as the epidemiological data. It may contain extra data (other
#'  districts or date ranges). The data must be in long format (one row for each
#'  date and environmental variable combination), and must start at absolutel
#'  minimum \code{laglen} (in \code{fc_control}) days before epi_data for
#'  forecasting.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables
#'@param quo_valuefield Quosure of user given field name of the value of the
#'  environmental data variable observations.
#'@param env_ref_data Historical averages by week of year for environmental
#'  variables. Used in extended environmental data into the future for long
#'  forecast time, to calculate anomalies in early detection period, and to
#'  display on timeseries in reports.
#'@param env_info Lookup table for environmental data - reference creation
#'  method (e.g. sum or mean), report labels, etc.
#'@param env_variables_used List of environmental variables that were used in
#'  the modeling
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param week_type String indicating the standard (WHO ISO-8601 or CDC epi
#'  weeks) that the weeks of the year in epidemiological and environmental
#'  reference data use ["ISO" or "CDC"].
#'
#'@return Environmental dataset, with data extended into the future forecast
#'  period. Unknown environmental data with runs of < 2 weeks is
#'  filled in with last known data (i.e. "persistence" method, using the mean of
#'  the previous week of known data). For missing data runs more than 2 weeks, the
#'  values are filled in using a progressive blend of the the mean of the last
#'  known week and the historical means.
#'
extend_env_future <- function(env_data,
                              quo_groupfield,
                              groupings,
                              quo_obsfield,
                              quo_valuefield,
                              env_ref_data,
                              env_info,
                              env_variables_used,
                              report_dates,
                              week_type){

  # Extend environmental data into the future forecast time period, while
  # also dealing with any other missing environmental data.

  # Progressive blend (over # of missing weeks) of
  # mean of last known week & historical/climatological means (weekly value)
  # but only if missing run is larger than 2 weeks * 7 = 14 days
  # if less than, just use persistence/carry forward/last known value
  # E.g. if 20 missing in a run:
  # 1 was filled in with previous week mean (recent value)
  # 2: 19/20 recent + 1/20 historical, 3: 18/20 recent + 2/20 historical, ... 20: 1/20 recent + 19/20 historical.
  # Will ALWAYS include part of recent known data (relevant if recent patterns are departure from climate averages)


  #Do not need data past end of forecast period
  env_trim <- env_data %>%
    dplyr::filter(obs_date <= report_dates$forecast$max)

  #Calculate the earliest of the latest known data dates
  # per env var, per geographic grouping
  earliest_end_known <- env_trim %>%
    #per geographic grouping, per environmental variable
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
    #the last known date for each
    dplyr::summarize(max_dates = max(obs_date, na.rm = TRUE)) %>%
    #the earliest of the last known
    dplyr::pull(max_dates) %>% min()


  #If earliest_end_known is end of forecast period, then no missing data
  if (earliest_end_known >= report_dates$forecast$max){

    env_extended_final <- env_trim

  } else {
    #Some amount of missing data exists

    #Calculate full/complete data table
    #combination of all groups, env vars, and dates (DAILY)
    #from earliest_end_known through the end of the forecast period
    env_future_complete <- tidyr::crossing(obs_date = seq.Date(earliest_end_known + 1,
                                                               report_dates$forecast$max, 1),
                                           group_temp = groupings,
                                           obs_temp = env_variables_used)
    #and fix names with NSE
    env_future_complete <- env_future_complete %>%
      dplyr::rename(!!rlang::quo_name(quo_groupfield) := group_temp,
                    !!rlang::quo_name(quo_obsfield) := obs_temp)

    #could have ragged env data per variable per grouping
    #so, antijoin with env_known_fill first to get the actually missing rows
    env_future_missing <- env_future_complete %>%
      dplyr::anti_join(env_trim, by = rlang::set_names(c(rlang::quo_name(quo_groupfield),
                                                         rlang::quo_name(quo_obsfield),
                                                         "obs_date"),
                                                       c(rlang::quo_name(quo_groupfield),
                                                         rlang::quo_name(quo_obsfield),
                                                         "obs_date")))



    #bind with existing data (NAs for everything else)
    env_future <- dplyr::bind_rows(env_trim, env_future_missing) %>%
      #mark which are about to be filled in
      dplyr::mutate(data_source = ifelse(is.na(val_epidemiar), "Extended", data_source))


    #function for rle needed to honor group_bys
    # computes an rle based on if value is NA or not
    # returns the number of rows in the run
    get_rle_na_info <- function(x){
      x_na_rle <- rle(is.na(x))
      run_id = rep(seq_along(x_na_rle$lengths), times = x_na_rle$lengths)
      run_tot <- rep(x_na_rle$lengths, times = x_na_rle$lengths)
      as_tibble(create_named_list(run_id, run_tot))
    }


    env_na_rle <- env_future %>%
      dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
      #make doubly sure in sorted date order
      arrange(obs_date) %>%
      #since adding multiple columns, use do instead of mutate
      do(cbind(., get_rle_na_info(.$val_epidemiar))) %>%
      #mutate(run_total = rle_na_count(val_epidemiar)) %>%
      #add a groupby with the new run ID
      group_by(!!quo_groupfield, !!quo_obsfield, run_id) %>%
      #creates an index of where that row is in the run
      mutate(id_in_run = seq_along(val_epidemiar))

    #find 1st NA, then take mean of previous week, input for that day
    #first NA now can be found with is.na(val_epidemiar) & id_in_run == 1
    #use zoo::rollapply for mean

    #Fill in first NA of a run with the mean of previous week
    env_na1fill <- env_na_rle %>%
      #confirm proper grouping
      dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
      #create a 1 day lag variable since need previous 7 days not including current
      mutate(val_lag1 = dplyr::lag(val_epidemiar, n = 1),
             #if_else to find the first NA
             val_epidemiar = ifelse(is.na(val_epidemiar) & id_in_run == 1,
                                    #zoo:rollapply to calculate mean of last 7 days (week) on lagged var
                                    zoo::rollapply(data = val_lag1,
                                                   width = 7,
                                                   FUN = mean,
                                                   align = "right",
                                                   na.rm = TRUE),
                                    #if not first NA, then contine with original val_epidemiar value
                                    val_epidemiar)) %>%
      #drop unneeded lag column
      select(-val_lag1)

    ##Prep for blending previous week mean & historical averages for other missing

    #Prep ref data - get only used vars
    env_ref_varused <- env_ref_data %>%
      dplyr::filter(!!quo_obsfield %in% env_variables_used)


    #joins for ref summary type, and summary for week
    env_join_ref <- env_na1fill %>%
      #add week, year fields
      epidemiar::add_datefields(week_type) %>%
      #get reference/summarizing method from user supplied env_info
      dplyr::left_join(env_info %>%
                         dplyr::select(!!quo_obsfield, reference_method),
                       by = rlang::set_names(rlang::quo_name(quo_obsfield),
                                             rlang::quo_name(quo_obsfield))) %>%
      #get weekly ref value
      dplyr::left_join(env_ref_varused %>%
                         dplyr::select(!!quo_obsfield, !!quo_groupfield, week_epidemiar, ref_value),
                       #NSE fun
                       by = rlang::set_names(c(rlang::quo_name(quo_groupfield),
                                               rlang::quo_name(quo_obsfield),
                                               "week_epidemiar"),
                                             c(rlang::quo_name(quo_groupfield),
                                               rlang::quo_name(quo_obsfield),
                                               "week_epidemiar")))

    #calculate NA missing values using carry|blend
    env_filled <- env_join_ref %>%
      #order very important for filling next step
      dplyr::arrange(!!quo_groupfield, !!quo_obsfield, obs_date) %>%
      dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
      #propagate last known value down rows
      dplyr::mutate(last_known = val_epidemiar) %>%
      #fill down, so missing weeks has "last known value" IN row for calculations
      tidyr::fill(last_known, .direction = "down") %>%
      #calculate parts (for all, will only use when needed)
      # with progressive blending based on id_in_run and run_tot
      mutate(recent_modifier = (run_tot - id_in_run - 1) / run_tot,
             recent_part = recent_modifier * last_known,
             historical_modifier = (id_in_run - 1) / run_tot,
             #historical is by week, so get pseudo-daily value depending on reference method,
             # i.e. how to summarize a week of data
             historical_value = dplyr::case_when(
               reference_method == "mean" ~ ref_value,
               reference_method == "sum"  ~ ref_value / 7,
               #default as if mean
               TRUE             ~ ref_value),
             historical_part = historical_modifier * historical_value,
             #testing
             val_orig = val_epidemiar,
             #only fill NA values
             val_epidemiar = ifelse(is.na(val_epidemiar),
                                    #persist if <15 days, blend if greater
                                    ifelse(run_tot < 15,
                                           last_known,
                                           recent_part + historical_part),
                                    #if notNA, then use existing val_epidemiar value
                                    val_epidemiar))

    #clean up
    env_extended_final <- env_filled %>%
      #remove all added columns to match original format
      select(-c(run_id, run_tot, id_in_run,
                week_epidemiar, year_epidemiar,
                last_known,
                reference_method, ref_value,
                recent_modifier, recent_part,
                historical_modifier, historical_value, historical_part,
                val_orig)) %>%
      #fill everything except original value field
      #for any other column that got vanished during crossing, etc.
      tidyr::fill(dplyr::everything(), -!!quo_valuefield, .direction = "down") %>%
      #ungroup to end
      ungroup()

  } #end else, meaning some missing data

  return(env_extended_final)

} # end extend_env_future





#' Extend epidemiology dataframe into future.
#'
#'@param epi_data Epidemiological data with case numbers per week, with date
#'  field "obs_date".
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'
#'@return Epidemiological dataset extended past the known epi data time range
#'  and into the future/forecast period. Case numbers are filled in the NA (to
#'  be forecasted), and the population is estimated in a persistence method.
#'
extend_epi_future <- function(epi_data, quo_popfield, quo_groupfield, groupings, report_dates){
  #extended epi data into future dates
  #for use in modeling later (results will be put elsewhere), this is for env and lags and modeling dataset
  epi_future <- tidyr::crossing(obs_date = report_dates$forecast$seq,
                                group_temp = groupings)
  #and fix names with NSE
  epi_future <- epi_future %>%
    dplyr::rename(!!quo_name(quo_groupfield) := group_temp)

  #bind with exisiting data (NAs for everything else in epi_future)
  extended_epi <- dplyr::bind_rows(epi_data, epi_future) %>%
    dplyr::arrange(!!quo_groupfield, obs_date)

  #fill population down
  extended_epi <- tidyr::fill(extended_epi, !!quo_popfield, .direction = "down")

  extended_epi
}

#' Format env data for modeling
#'
#'@param env_data_extd An environmental dataset extended into the
#'  future/forecast period with estimated values for the environmental
#'  variables.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables.
#'
#'@return An environmental dataset formatted to pass over to BAM/GAM modeling.
#'
env_format_fc <- function(env_data_extd, quo_groupfield, quo_obsfield){
  #turns long format into wide format - one entry per day per group
  #1: groupfield, 2: Date, 3: numericdate, 4+: env var (column name is env name)
  env_spread <- env_data_extd %>%
    dplyr::mutate(numericdate = as.numeric(obs_date)) %>%
    dplyr::select(!!quo_groupfield, !!quo_obsfield, obs_date, numericdate, val_epidemiar) %>%
    tidyr::spread(key = !!quo_obsfield, value = val_epidemiar)

  env_spread
}

#' Format epi data for modeling
#'
#'@param epi_data_extd An epidemiological dataset extended into the
#'  future/forecast period with NA values for to-be-forecasted case numbers.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param fc_control Parameters for forecasting, including which environmental
#'  variable to include and any geographic clusters.
#'
#'@return An epidemiological dataset formatted to pass over to BAM/GAM modeling.
#'
epi_format_fc <- function(epi_data_extd, quo_groupfield, fc_control){

  #Get cluster information from model
  epi_format <- epi_data_extd %>%
    #join with cluster info
    dplyr::left_join(fc_control$clusters,
                     #NSE
                     by = rlang::set_names(rlang::quo_name(quo_groupfield),
                                           rlang::quo_name(quo_groupfield))) %>%
    #set cluster id as factor, must be for regression later
    dplyr::mutate(cluster_id = as.factor(cluster_id),
                  #need numeric date for regression
                  numericdate = as.numeric(obs_date))

  epi_format
}

#' Convert environmental data into anomalies.
#'
#' Raw environmental values are not used in modeling, but rather their
#' anomalies, departures for the historical "normal". We are looking at the
#' influence of deviation from normal in the environmental factors to help
#' explain deviations from normal in the human cases.
#'
#'@param env_fc Environmental data formatted for forecasting by env_format_fc().
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param env_variables_used List of environmental variables that were used in
#'  the modeling.
#'@param ncores The number of physical cores to use in parallel processing, set
#'  in fc_control$ncores, else the max of the number of physical core available
#'  minus 1, or 1 core.
#'
#'@return Environmental dataset in same format as env_fc but with the residuals
#'  from a GAM with geographic unit and cyclical cubic regression spline on day
#'  of year per geographic group in place of the original raw values.
#'
anomalize_env <- function(env_fc, quo_groupfield, env_variables_used, ncores) {

  # Loop through each environmental variable replacing non-NA observations
  # with residuals from a gam with only geographic area (group) and day of year

  #Originally written with dataframes. Tibble/NSE/dplyr conversion not yet fully done.
  # #new tibble
  # env_fc_anom <- env_fc %>%
  #   mutate(group_factor = factor(!!quo_groupfield),
  #          doy = format(obs_date, "%j")) %>% as.numeric()

  # needed data for gam
  group_factor <- env_fc %>% pull(!!quo_groupfield) %>% factor()
  doy <- env_fc %>% pull(obs_date) %>% format("%j") %>% as.numeric()

  env_fc <- as.data.frame(env_fc)

  # loop through environmental columns
  # note: brittle on column position until rewrite
  for (curcol in 4:ncol(env_fc)) {

    #if more than one geographic area
    if (nlevels(group_factor) > 1){
      tempbam <- mgcv::bam(env_fc[,curcol] ~ group_factor + s(doy, bs="cc", by=group_factor),
                           data=env_fc,
                           discrete = TRUE,
                           nthreads = ncores)
    } else {
      #if only 1 geographic area, then run without group_factor
      tempbam <- mgcv::bam(env_fc[,curcol] ~ s(doy, bs="cc"),
                           data=env_fc,
                           discrete = TRUE,
                           nthreads = ncores)
    }

    # could perhaps more cleverly be figured out by understanding the na.options of bam,
    # but for the moment just replace non-NA observations with their residuals
    env_fc[!is.na(env_fc[,curcol]),curcol] <- tempbam$residuals

  }

  env_fc <- tibble::as_tibble(env_fc)
  return(env_fc)

}

#' Lag the environmental data.
#'
#'@param epi_fc An epidemiological dataset extended into the
#'  future/forecast period with NA values for to-be-forecasted case numbers,
#'  as formatted for forecasting by epi_format_fc().
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param env_fc Environmental data formatted for forecasting by env_format_fc().
#'@param env_variables_used List of environmental variables that were used in
#'  the modeling.
#'@param laglen The maximum number of days in the past to consider interactions
#'  between the environmental variable anomalies and the disease case counts.
#'
#'@return Wide dataset based on epidemiological data dates with five
#'  bandsummaries per environmental variable, from the basis spline summaries of
#'  the lagged environmental variable.
#'
lag_environ_to_epi <- function(epi_fc, quo_groupfield, groupings,
                               env_fc, env_variables_used, laglen){

  #create lag frame
  datalagger <- tidyr::crossing(group_temp = groupings,
                                obs_date = unique(epi_fc$obs_date),
                                lag = seq(from = 0, to = laglen - 1, by = 1)) %>%
    # #same order from originally written expand.grid
    # arrange(lag, Date, group_temp) %>%
    #add lagging date
    dplyr::mutate(laggeddate = obs_date - as.difftime(lag, units = "days"))

  #and fix names with NSE
  datalagger <- datalagger %>%
    dplyr::rename(!!quo_name(quo_groupfield) := group_temp)

  #add env data
  datalagger <- dplyr::left_join(datalagger, env_fc,
                                 #because dplyr NSE, notice flip order
                                 by = rlang::set_names(c(rlang::quo_name(quo_groupfield), "obs_date"),
                                                       c(rlang::quo_name(quo_groupfield), "laggeddate")))

  # pivot lagged environmental data to epi data
  epi_lagged <- epi_fc #to more easily debug and rerun
  for (curcol in which(colnames(env_fc) %in% env_variables_used)){
    valuevar <- colnames(env_fc)[curcol]
    #wide data for all lags of that env var
    meandat <- datalagger %>%
      dplyr::select(!!quo_groupfield, obs_date, lag, valuevar) %>%
      tidyr::spread(key = lag, value = valuevar)
    #rename lag columns (but not groupfield or Date)
    names(meandat)[-(1:2)] <- paste0(valuevar, "_", names(meandat)[-(1:2)])

    #join cur var wide data to epi data
    epi_lagged <- dplyr::left_join(epi_lagged, meandat,
                                   #dplyr NSE
                                   by = rlang::set_names(c(rlang::quo_name(quo_groupfield), "obs_date"),
                                                         c(rlang::quo_name(quo_groupfield), "obs_date")))
  } #end pivot loop

  # # set up distributed lag basis functions (creates 5 basis functions)
  # lagframe <- data.frame(x = seq(from = 1, to = laglen, by = 1))
  # alpha <- 1/4
  # distlagfunc <- splines::ns(lagframe$x, intercept = TRUE,
  #                            knots = quantile(lagframe$x,
  #                                             probs=seq(from = alpha, to = 1 - alpha,
  #                                                       by = alpha),
  #                                             na.rm = TRUE))
  # dlagdeg <- pracma::size(distlagfunc)[2]

  # set up distributed lag basis functions (creates 7 basis functions)
  alpha <- 1/4
  distlagfunc <- splines::bs(x=seq(from=1, to=laglen, by=1), intercept=TRUE,
                             knots=quantile(seq(from=1, to=laglen, by=1),
                                            probs=seq(from=alpha, to=1-alpha, by=alpha),
                                            na.rm=TRUE))
  dlagdeg <- ncol(distlagfunc)


  # create actual distributed lag summaries
  for (curvar in env_variables_used){
    bandsum <- matrix(data = rep(0, nrow(epi_lagged) * dlagdeg),
                      nrow = nrow(epi_lagged), ncol = dlagdeg)
    #first column of that variable (0 lag)
    mindex <- which(colnames(epi_lagged) == paste0(curvar, "_0"))
    #temp working matrix
    bandtemp <- as.matrix(epi_lagged[, (mindex:(mindex+laglen-1))])
    #distributed lag summaries
    for (j in 1:dlagdeg){
      bandsum[, j] <- bandtemp %*% distlagfunc[,j]
    }
    bandsum <- data.frame(bandsum)
    names(bandsum) <- paste0("bandsum_", curvar, "_", 1:dlagdeg)

    # we used to do a submatrix here so that the regression formulae would
    # be more easily written, but this was incompatible with dplyr
    epi_lagged <- dplyr::bind_cols(epi_lagged, bandsum)

    #created summary value for each basis function (5) per env variable per group per week (based on epidemiological data time unit)

  } #end distr lag summary loop

  #only keep bandsummaries (daily lags can be removed to free up a lot of space)
  #  note: ^ matches beginning of string, otherwise we'd get the bandsummaries too, which we want to keep
  for (cvar in env_variables_used){
    epi_lagged[, which(grepl(paste0("^", cvar, "_"), colnames(epi_lagged)))] <- NULL
  }

  epi_lagged
}

#' Run forecast regression
#'
#'@param epi_lag Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), as output by lag_environ_to_epi().
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param env_variables_used List of environmental variables that were used in
#'  the modeling.
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param req_date The end date of requested forecast regression. When fit_freq
#'  == "once", this is the last date of the full report, the end date of the
#'  forecast period.
#'@param ncores The number of physical cores to use in parallel processing, set
#'  in fc_control$ncores, else the max of the number of physical core available
#'  minus 1, or 1 core.
#'@param fit_freq String indicating "once" or "weekly" on how often to fit the
#'  model - once for the whole report, or every week of the report. Unless
#'  otherwise needed, the value should be "once", as weekly drastically
#'  increases processing time.
#'@param model_run TRUE/FALSE flag for whether to only generate the model
#'  regression object plus metadata. This model can be cached and used later on
#'  its own, skipping a large portion of the slow calculations for future runs.
#'@param model_cached The output of a previous model_run = TRUE run of
#'  run_epidemia() that produces a model (regression object) and metadata. The
#'  metadata will be used for input checking and validation. Using a prebuilt
#'  model saves on processing time, but will need to be updated periodically.
#'@param model_choice Critical argument to choose the type of model to generate.
#'  The options are versions that the EPIDEMIA team has used for forecasting.
#'  The first supported options is "poisson-bam" ("p") which is the original
#'  epidemiar model: a Poisson regression using bam (for large data GAMs), with
#'  a smoothed cyclical for seasonality. The default for fc_control$anom_env is
#'  TRUE for using the anomalies of environmental variables rather than their
#'  raw values. The second option is "negbin" ("n") which is a negative binomial
#'  regression using glm, with no external seasonality terms - letting the
#'  natural cyclical behavior of the environmental variables fill that role. The
#'  default for fc_control$anom_env is FALSE and uses the actual observation
#'  values in the modeling. The fc_control$anom_env can be overruled by the user
#'  providing a value, but this is not recommended unless you are doing
#'  comparisons.
#'@param theta From fc_control$theta, the value of theta for a "negbin" model.
#'  If present, will use glm(..., family = MASS::negative.binomial(theta)).  If
#'  missing, will use MASS::glm.nb().
#'
#'@return Named list containing:
#'date_preds: Full forecasted resulting dataset.
#'reg_obj: The regression object from modeling.
#'Unless model_run is TRUE, in which case only the regression object is returned.
#'
#'
forecast_regression <- function(epi_lag,
                                quo_groupfield,
                                groupings,
                                env_variables_used,
                                report_dates,
                                req_date,
                                ncores,
                                fit_freq,
                                model_run,
                                model_cached = NULL,
                                model_choice,
                                theta){

  if (fit_freq == "once"){
    #single fits use all the data available
    last_known_date <-  report_dates$known$max
  } else if (fit_freq == "week"){
    # for "week" model fits, forecasts are done knowing up to just before that date
    last_known_date <- req_date - lubridate::as.difftime(1, units = "days")
  }

  ## Set up data

  #mark known or not
  epi_lag <- epi_lag %>%
    dplyr::mutate(known = ifelse(obs_date <= last_known_date, 1, 0))

  # ensure that quo_name(quo_groupfield) is a factor - gam/bam will fail if given a character,
  # which is unusual among regression functions, which typically just coerce into factors.
  epi_lag <- epi_lag %>% dplyr::mutate(!!rlang::quo_name(quo_groupfield) := factor(!!quo_groupfield))
  #number of geographic area groupings
  n_groupings <- epi_lag %>% pull(!!quo_groupfield) %>% nlevels()

  #number of clusters
  n_clusters <- nlevels(epi_lag$cluster_id)

  # create a doy field so that we can use a cyclical spline
  epi_lag <- dplyr::mutate(epi_lag, doy = as.numeric(format(obs_date, "%j")))

  # create modified bspline basis in epi_lag file to model longterm trends
  epi_lag <- cbind(epi_lag, truncpoly(x=epi_lag$obs_date,
                                      degree=6,
                                      maxobs=max(epi_lag$obs_date[epi_lag$known==1], na.rm=TRUE)))



  ## If model_cached is NOT given, then create model / run regression
  if (is.null(model_cached)){

    #create variable bandsummaries equation piece
    #  e.g. 'bandsummaries_{var1} * cluster_id' for however many env var bandsummaries there are
    bandsums_list <- grep("bandsum_*", colnames(epi_lag), value = TRUE)
    bandsums_cl_list <- paste(bandsums_list, ": cluster_id")
    #need variant without known multiplication if <= 1 clusters
    if (n_clusters > 1) {
      bandsums_eq <- glue::glue_collapse(bandsums_cl_list, sep =" + ")
    } else {
      bandsums_eq <- glue::glue_collapse(bandsums_list, sep = " + ")
    }

    # get list of modbspline reserved variables and format for inclusion into model
    modb_list <- grep("modbs_reserved_*", colnames(epi_lag), value = TRUE)
    # variant depending on >1 geographic area groupings
    if (n_groupings > 1){
      modb_list_grp <- paste(modb_list, ":", rlang::quo_name(quo_groupfield))
      modb_eq <- glue::glue_collapse(modb_list_grp, sep = " + ")
    } else {
      modb_eq <- glue::glue_collapse(modb_list, sep = " + ")
    }

    #filter to known
    epi_known <- epi_lag %>% dplyr::filter(known == 1)


    # Model building switching point

    regress <- build_model(model_choice,
                           n_groupings,
                           quo_groupfield,
                           modb_eq,
                           bandsums_eq,
                           epi_known,
                           ncores,
                           theta)

  } else {
    #if model_cached given, then use that as regress instead of building a new one (above)

    #message with model input
    message("Using given cached ", model_cached$model_info$model_choice, " model, created ",
            model_cached$model_info$date_created, ", with epidemiological data up through ",
            model_cached$model_info$known_epi_range$max, ".")

    regress <- model_cached$model_obj
  }

  ## If model run, return regression object to run_forecast() at this point
  if (model_run){
    return(regress)
  }

  ## Creating predictions switching point on model choice
  preds <- create_predictions(model_choice,
                              regress,
                              epi_lag,
                              req_date,
                              ncores)


  ## Clean up
  #remove distributed lag summaries and bspline basis, which are no longer helpful
  band_names <- grep("bandsum_*", colnames(epi_lag), value = TRUE)
  bspl_names <- grep("modbs_reserved_*", colnames(epi_lag), value = TRUE)
  #remove
  epi_lag_trim <- dplyr::select(epi_lag, -dplyr::one_of(band_names))
  epi_lag_trim <- dplyr::select(epi_lag_trim, -dplyr::one_of(bspl_names))


  #now cbind to get ready to return
  epi_preds <- cbind(epi_lag_trim %>%
                       filter(obs_date <= req_date),
                     as.data.frame(preds)) %>%
    #and convert factor back to character for the groupings (originally converted b/c of bam/gam requirements)
    dplyr::mutate(!!rlang::quo_name(quo_groupfield) := as.character(!!quo_groupfield))

  if (fit_freq == "once"){
    #for single model fit, this has all the data we need, just trim to report dates
    date_preds <- epi_preds %>%
      filter(obs_date >= report_dates$full$min)
  } else if (fit_freq == "week"){
    #prediction of interest are last ones (equiv to req_date) per groupfield
    date_preds <- epi_preds %>%
      dplyr::group_by(!!quo_groupfield) %>%
      dplyr::filter(obs_date == req_date)
  }

  forecast_reg_results <- create_named_list(date_preds,
                                            regress)
}


#'Build the appropriate model
#'
#'@param model_choice Critical argument to choose the type of model to generate.
#'  The options are versions that the EPIDEMIA team has used for forecasting.
#'  The first supported options is "poisson-bam" ("p") which is the original
#'  epidemiar model: a Poisson regression using bam (for large data GAMs), with
#'  a smoothed cyclical for seasonality. The default for fc_control$anom_env is
#'  TRUE for using the anomalies of environmental variables rather than their
#'  raw values. The second option is "negbin" ("n") which is a negative binomial
#'  regression using glm, with no external seasonality terms - letting the
#'  natural cyclical behavior of the environmental variables fill that role. The
#'  default for fc_control$anom_env is FALSE and uses the actual observation
#'  values in the modeling. The fc_control$anom_env can be overruled by the user
#'  providing a value, but this is not recommended unless you are doing
#'  comparisons.
#'@param n_groupings Count of the number of geographic groupings in the model.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param modb_eq Pieces of the regression formula that include the modified
#'  basis functions to account for long term trend (with or without groupings,
#'  as appropriate).
#'@param bandsums_eq Pieces of the regression formula that include the b-spline
#'  bandsummaries of the environmental factors.
#'@param epi_known Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), with column marking if "known"
#'  data and groupings converted to factors.
#'@param ncores The number of physical cores to use in parallel processing, set
#'  in fc_control$ncores, else the max of the number of physical core available
#'  minus 1, or 1 core.
#'@param theta From fc_control$theta, the value of theta for a "negbin" model.
#'  If present, will use glm(..., family = MASS::negative.binomial(theta)).  If
#'  missing, will use MASS::glm.nb().

#'
#'@return Regression object
#'
#'
build_model <- function(model_choice,
                        n_groupings,
                        quo_groupfield,
                        modb_eq,
                        bandsums_eq,
                        epi_known,
                        ncores,
                        theta){

  #POISSON-BAM (set as default in first round input checking)
  if (model_choice == "poisson-bam"){

    message("Building Poisson model using bam() and forced cyclical...")

    #due to dplyr NSE and bandsum eq and modb_eq pieces, easier to create
    #expression to give to modeling function
    #different versions if multiple geographic area groupings or not
    if (n_groupings > 1){
      reg_eq <- stats::as.formula(paste("modeledvar ~ ",
                                        rlang::quo_name(quo_groupfield),
                                        " + s(doy, bs=\"cc\", by=",
                                        rlang::quo_name(quo_groupfield),
                                        ") + ",
                                        modb_eq, " + ",
                                        bandsums_eq))
    } else {
      reg_eq <- stats::as.formula(paste("modeledvar ~ ",
                                        "s(doy, bs=\"cc\") + ",
                                        modb_eq, " + ",
                                        bandsums_eq))
    }

    # run bam
    # Using discrete = TRUE was much faster than using parallel with bam.
    regress <- mgcv::bam(reg_eq, data = epi_known,
                         family=poisson(),
                         control=mgcv::gam.control(trace=FALSE),
                         discrete = TRUE,
                         nthreads = ncores)


  } else if (model_choice == "negbin"){
    #NEGATIVE BINOMIAL using GLM

    message("Building negative binomial model...")

    #due to dplyr NSE and bandsum eq and modb_eq pieces, easier to create
    #expression to give to modeling function
    #different versions if multiple geographic area groupings or not
    #No cycical (as opposed to bam with s())
    if (n_groupings > 1){
      reg_eq <- stats::as.formula(paste("modeledvar ~ ",
                                        rlang::quo_name(quo_groupfield), " + ",
                                        modb_eq, " + ",
                                        bandsums_eq))
    } else {
      reg_eq <- stats::as.formula(paste("modeledvar ~ ",
                                        modb_eq, " + ",
                                        bandsums_eq))
    }

    # run glm
    # Which negative binomial function depends on if fc_control$theta exists
    if(!is.null(theta)){
      message("Theta value provided. Running with glm(..., family = MASS::negative.binomial(theta = ", theta, "))...")
      regress <- stats::glm(reg_eq,
                            data = epi_known,
                            #theta value REQUIRED
                            #family = MASS::negative.binomial(theta=2.31),
                            family = MASS::negative.binomial(theta = theta))
    } else {
      message("Theta estimate (fc_control$theta) not provided, running with MASS::glm.nb()...")
      regress <- MASS::glm.nb(reg_eq,
                              data = epi_known)
    }


  } else if (model_choice == "null-persistence"){

    #null model
    #persistence (carry forward)
    #no regression object

    #create "model" using known data.
    #Will fill down in create_predictions
    regress <- epi_known %>%
      #grouping by geographical unit
      dplyr::group_by(!!quo_groupfield) %>%
      #prediction is 1 lag (previous week)
      #fit is name of value from regression models
      dplyr::mutate(fit = dplyr::lag(modeledvar, n = 1)) %>%
      #cleaning up as not needed, and for bug hunting
      dplyr::select(-dplyr::starts_with("band")) %>%
      dplyr::select(-dplyr::starts_with("modbs"))



  } else if (model_choice == "null-averageweek"){

    #null model
    #average of week of year (from historical data)
    #not a regression object

    #create "model" (averages) using known data.
    regress <- epi_known %>%
      #calculate averages per geographic group per week of year
      dplyr::group_by(!!quo_groupfield, week_epidemiar) %>%
      dplyr::summarize(fit = mean(modeledvar, na.rm = TRUE))


  } else {
    #Shouldn't happen, just in case.
    stop("Error in selecting model choice.")
  }
} # end build_model()



#'Create the appropriate predictions/forecasts.
#'
#'@param model_choice Critical argument to choose the type of model to generate.
#'  The options are versions that the EPIDEMIA team has used for forecasting.
#'  The first supported options is "poisson-gam" ("p") which is the original
#'  epidemiar model: a Poisson regression using bam (for large data GAMs), with
#'  a smoothed cyclical for seasonality. The default for fc_control$anom_env is
#'  TRUE for using the anomalies of environmental variables rather than their
#'  raw values. The second option is "negbin" ("n") which is a negative binomial
#'  regression using glm, with no external seasonality terms - letting the
#'  natural cyclical behavior of the environmental variables fill that role. The
#'  default for fc_control$anom_env is FALSE and uses the actual observation
#'  values in the modeling. The fc_control$anom_env can be overruled by the user
#'  providing a value, but this is not recommended unless you are doing
#'  comparisons.
#'@param regress The regression object, either the user-supplied one, or
#'  the one just generated.
#'@param epi_lag Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), with groupings as a factor.
#'@param req_date The end date of requested forecast regression. When fit_freq
#'  == "once", this is the last date of the full report, the end date of the
#'  forecast period.
#'@param ncores The number of physical cores to use in parallel processing, set
#'  in fc_control$ncores, else the max of the number of physical core available
#'  minus 1, or 1 core.
#'
#'@return A dataset from predict() using the regression object generated in
#'  build_model or a newly created one. The dataset includes the
#'  predicted/forecast values through the end of the report requested.
#'
#'
create_predictions <- function(model_choice,
                               regress,
                               epi_lag,
                               req_date,
                               ncores){

  #POISSON-BAM (set as default in first round input checking)
  if (model_choice == "poisson-bam"){

    message("Creating Poisson predictions...")


    ## Create predictions from either newly generated model, or given one

    #output prediction (through req_date)
    preds <- mgcv::predict.bam(regress,
                               newdata = epi_lag %>% dplyr::filter(obs_date <= req_date),
                               se.fit = TRUE,       # included for backwards compatibility
                               type="response",
                               discrete = TRUE,
                               n.threads = ncores)



  } else if (model_choice == "negbin"){
    #NEGATIVE BINOMIAL using GLM

    message("Creating negative binomial predictions...")


    ## Create predictions from either newly generated model, or given one

    #output prediction (through req_date)
    preds <- stats::predict.glm(regress,
                                newdata = epi_lag %>% dplyr::filter(obs_date <= req_date),
                                se.fit = TRUE,       # included for backwards compatibility
                                type="response")


  } else if (model_choice == "null-persistence"){

    message("Creating predictions using persistence null model...")

    #persistence model just carries forward the last known value
    #the important part is the forecast / trailing end part
    #manipulating to be in quasi-same format as the other models return

    #cleaning up as not needed, and for bug hunting
    epi_lag <- epi_lag %>%
      dplyr::select(-dplyr::starts_with("band")) %>%
      dplyr::select(-dplyr::starts_with("modbs"))

    #regress is a tibble not regression object here
    # has a variable fit with lag of 1 on known data
    #epi_lag has the newer rows
    preds <- epi_lag %>%
      #filter to requested date
      dplyr::filter(obs_date <= req_date) %>%
      #join to get "fit" values from "model"
      #join on all shared columns (i.e. everything in regress not "fit") to prevent renaming
      dplyr::left_join(regress, by = names(regress)[!names(regress) %in% c("fit")]) %>%
      #important at end/fc section, when we fill down
      tidyr::fill(fit, .direction = "down") %>%
      #format into nominal regression predict output
      dplyr::select(fit) %>%
      as.data.frame()

  } else if (model_choice == "null-averageweek"){

    message("Creating predictions using average week of year null model...")

    #average week null model calculates the average cases of that
    # week of year from historical data
    #manipulating to be in quasi-same format as the other models return

    #regress is the averages per week of year from known data

    epi_lag <- epi_lag %>%
      #filter to requested date
      dplyr::filter(obs_date <= req_date)

    #join back
    preds <- epi_lag %>%
      #join to get average values
      #join on all shared columns (i.e. everything in regress not "fit") to prevent renaming
      # and so don't need column names not passed into this function
      dplyr::left_join(regress, by = names(regress)[!names(regress) %in% c("fit")]) %>%
      #format into nominal regression output
      dplyr::select(fit) %>%
      as.data.frame()


  } else {
    #Shouldn't happen, just in case.
    stop("Error in selecting model choice.")
  }

} #end create_predictions()
