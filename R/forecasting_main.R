# Main run_epidemiar() subfunctions related to forecasting

#' Runs the forecast modeling
#'
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables.
#'@param quo_valuefield Quosure of user given field name of the value of the
#'  environmental data variable observations.
#'@param env_variables List of environmental variables that exist in env_data.
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param valid_run Internal TRUE/FALSE for whether this is part of a validation run.
#'
#'@inheritParams run_epidemia
#'
#'@return Named list containing:
#'fc_epi: Full forecasted resulting dataset.
#'fc_res: The forecasted series in report format.
#'env_data_extd: Data set of the environmental data variables extended into the
#'  unknown/future.
#'env_variables_used: list of environmental variables that were used in the
#'  modeling (had to be listed in model variables input file and present the
#'  env_data and env_info datasets)
#'env_dt_ranges: Date ranges of the input environmental data.
#'reg_obj: The regression object from modeling.
#'Unless model_run is TRUE, in which case only the regression object is returned.
#'
#'
run_forecast <- function(epi_data,
                          quo_popfield,
                          quo_groupfield,
                          env_data,
                          quo_obsfield,
                          quo_valuefield,
                          env_ref_data,
                          env_info,
                          fc_model_family,
                          report_settings,
                          #internal/calculated
                          valid_run,
                          groupings,
                          env_variables,
                          report_dates){

  message("Preparing for forecasting...")


  # trim to the needed env variables as dictated by the model
  env_data <- pull_model_envvars(env_data = env_data,
                                 quo_obsfield = quo_obsfield,
                                 env_var = report_settings[["env_var"]])
  #create alphabetical list of ONLY USED unique environmental variables
  env_variables_used <- dplyr::pull(env_data, !!quo_obsfield) %>% unique() %>% sort()

  # extract start & end dates for each variable for log file
  env_dt_ranges <- dplyr::group_by(env_data, !!quo_obsfield) %>%
    dplyr::summarize(start_dt = min(.data$obs_date), end_dt = max(.data$obs_date))

  # extend data into future, for future forecast portion
  env_data_extd <- extend_env_future(env_data,
                                     quo_groupfield,
                                     quo_obsfield,
                                     quo_valuefield,
                                     env_ref_data,
                                     env_info,
                                     fc_model_family, #reduced processing for naive models
                                     #pull from report_settings
                                     epi_date_type = report_settings[["epi_date_type"]],
                                     #calculated/internal
                                     valid_run,
                                     groupings,
                                     env_variables_used,
                                     report_dates)

  #extend into future and/or gaps in requested report dates & known data
  epi_data_extd <- extend_epi_future(epi_data,
                                     quo_popfield,
                                     quo_groupfield,
                                     #calculated/internal
                                     groupings,
                                     report_dates)

  # format the data for forecasting algorithm
  env_fc <- env_format_fc(env_data_extd,
                          quo_groupfield,
                          quo_obsfield)
  epi_fc <- epi_format_fc(epi_data_extd,
                          quo_groupfield,
                          fc_clusters = report_settings[["fc_clusters"]])

  # anomalizing the environ data, if requested.
  # note: brittle on format from env_format_fc(), edit with caution
  # AND not a naive model run

  if (!fc_model_family == "naive-persistence" & !fc_model_family == "naive-averageweek"){
    if (report_settings[["env_anomalies"]]){
      message("Anomalizing the environmental variables...")
      env_fc <- anomalize_env(env_fc,
                              quo_groupfield,
                              nthreads = report_settings[["fc_nthreads"]],
                              #calculated/internal
                              env_variables_used)
    }
  } #end not if naive model run


  # create the lags
  epi_lag <- lag_environ_to_epi(epi_fc,
                                env_fc,
                                quo_groupfield,
                                report_settings,
                                #calculated/internal
                                groupings,
                                env_variables_used)

  # add week of year, needed for null-weekaverage model
  # switch epi_date_type to week_type needed for add_datefields()
  week_type <- dplyr::case_when(
    report_settings[["epi_date_type"]] == "weekISO" ~ "ISO",
    report_settings[["epi_date_type"]] == "weekCDC"  ~ "CDC",
    #default as if mean
    TRUE             ~ NA_character_)
  epi_lag <- add_datefields(epi_lag, week_type)


  # If only model_run, then return to run_epidemia() here
  if (report_settings[["model_run"]]){
  model_run_result <- forecast_regression(epi_lag,
                                          quo_groupfield,
                                          fc_model_family,
                                          report_settings,
                                          #internal calculated
                                          groupings,
                                          env_variables_used,
                                          report_dates,
                                          req_date = report_dates$full$max)

    model_run_only <- create_named_list(env_variables_used,
                                        env_dt_ranges,
                                        reg_obj = model_run_result)
    return(model_run_only)
  }


  #Split regression call depending on {once|week} model fit frequency

  if (report_settings[["dev_fc_fit_freq"]] == "once"){
    #for single fit, call with last week (and subfunction has switch to return all)
    forereg_return <- forecast_regression(epi_lag,
                                          quo_groupfield,
                                          fc_model_family,
                                          report_settings,
                                          #internal calculated
                                          groupings,
                                          env_variables_used,
                                          report_dates,
                                          req_date = report_dates$full$max)
    preds_catch <- forereg_return$date_preds
    reg_obj <- forereg_return$regress

  } else if (report_settings[["dev_fc_fit_freq"]] == "week") {
    message("DEVELOPER: Fitting model per week...")
    # for each week of report, run forecast
    # initialize: prediction returns 4 columns
    preds_catch <- data.frame()
    #loop by week
    for (w in seq_along(report_dates$full$seq)){
      message("Forecasting week ", w, " starting at ", Sys.time())
      dt <- report_dates$full$seq[w]
      forereg_return <- forecast_regression(epi_lag,
                                            quo_groupfield,
                                            fc_model_family,
                                            report_settings,
                                            #internal calculated
                                            groupings,
                                            env_variables_used,
                                            report_dates,
                                            req_date = dt)

      dt_preds <- forereg_return$date_preds
      preds_catch <- rbind(preds_catch, as.data.frame(dt_preds))

      #taking advantage that only result will be of the last loop through
      reg_obj <- forereg_return$regress
    }

  } else stop("Developer setting model fit frequency unknown, please review/remove report_settings$dev_fc_fit_freq parameter.")


  # Interval calculation
  preds_catch <- preds_catch %>%
    dplyr::mutate(fc_cases = .data$preds,
                  fc_cases_upr = .data$preds+1.96*sqrt(.data$preds),
                  fc_cases_lwr = .data$preds-1.96*sqrt(.data$preds))

  #Trim fc report results ONLY (not full epi dataset) to report period
  preds_catch_trim <- preds_catch %>%
    dplyr::filter(.data$obs_date >= report_dates$full$min)

  # extract fc series into report format
  # if else off of report_value_type of reporting in terms of cases or incidence
  # using full if else blocks to do all 3 at once, rather than if_elses in each variable
  if (report_settings[["report_value_type"]] == "cases"){
    fc_res <- preds_catch_trim %>%
      dplyr::mutate(series = "fc",
                    lab = "Forecast Trend",
                    value = .data$fc_cases,
                    upper = .data$fc_cases_upr,
                    lower = .data$fc_cases_lwr)
  } else if (report_settings[["report_value_type"]] == "incidence"){
    fc_res <- preds_catch_trim %>%
      dplyr::mutate(series = "fc",
                    lab = "Forecast Trend",
                    value = .data$fc_cases / !!quo_popfield * report_settings[["report_inc_per"]],
                    upper = .data$fc_cases_upr / !!quo_popfield * report_settings[["report_inc_per"]],
                    lower = .data$fc_cases_lwr / !!quo_popfield * report_settings[["report_inc_per"]])

  } else { #shouldn't happen
    fc_res <- preds_catch_trim %>%
      dplyr::mutate(series = "fc",
                    lab = "Forecast Trend",
                    value = NA_real_,
                    upper = NA_real_,
                    lower = NA_real_)

  }

  #select only the needed columns
  fc_res <- fc_res %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)

  # case_when evaluates ALL RHS. Not appropriate to use here as populationfield may not exist.
  # value = dplyr::case_when(
  #   #if reporting in case counts
  #   report_settings[["report_value_type"]] == "cases" ~ fc_cases,
  #   #if incidence
  #   report_settings[["report_value_type"]] == "incidence" ~ fc_cases / !!quo_popfield * report_settings[["report_inc_per"]],
  #   #otherwise
  #   TRUE ~ NA_real_),
  #
  # upper = dplyr::case_when(
  #   #if reporting in case counts
  #   report_settings[["report_value_type"]] == "cases" ~ fc_cases_upr,
  #   #if incidence
  #   report_settings[["report_value_type"]] == "incidence" ~ fc_cases_upr / !!quo_popfield * report_settings[["report_inc_per"]],
  #   #otherwise
  #   TRUE ~ NA_real_),
  #
  # lower = dplyr::case_when(
  #   #if reporting in case counts
  #   report_settings[["report_value_type"]] == "cases" ~ fc_cases_lwr,
  #   #if incidence
  #   report_settings[["report_value_type"]] == "incidence" ~ fc_cases_lwr / !!quo_popfield * report_settings[["report_inc_per"]],
  #   #otherwise
  #   TRUE ~ NA_real_)


  # return list with res and other needed items
  fc_res_full <- create_named_list(fc_epi = preds_catch,
                                   fc_res,
                                   env_data_extd,
                                   env_variables_used,
                                   env_dt_ranges,
                                   reg_obj)
}


#' Run forecast regression
#'
#'@param epi_lag Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), as output by lag_environ_to_epi().
#'@param env_variables_used List of environmental variables that were used in
#'  the modeling.
#'@param req_date The end date of requested forecast regression. When fit_freq
#'  == "once", this is the last date of the full report, the end date of the
#'  forecast period.
#'
#'@inheritParams run_epidemia
#'@inheritParams run_forecast
#'
#'@return Named list containing:
#'date_preds: Full forecasted resulting dataset.
#'reg_obj: The regression object from modeling.
#'Unless model_run is TRUE, in which case only the regression object is returned.
#'
#'
forecast_regression <- function(epi_lag,
                                quo_groupfield,
                                fc_model_family,
                                report_settings,
                                #internal calculated
                                groupings,
                                env_variables_used,
                                report_dates,
                                req_date){


  if (report_settings[["dev_fc_fit_freq"]] == "once"){
    #single fits use all the data available
    last_known_date <-  report_dates$known$max
  } else if (report_settings[["dev_fc_fit_freq"]] == "week"){
    # for "week" model fits, forecasts are done knowing up to just before that date
    last_known_date <- req_date - lubridate::as.difftime(1, units = "days")
  }

  ## Set up data

  #mark within prior known range or not
  epi_lag <- epi_lag %>%
    dplyr::mutate(input = ifelse(.data$obs_date <= last_known_date, 1, 0))

  # ensure that quo_name(quo_groupfield) is a factor - gam/bam will fail if given a character,
  # which is unusual among regression functions, which typically just coerce into factors.
  epi_lag <- epi_lag %>% dplyr::mutate(!!rlang::quo_name(quo_groupfield) := factor(!!quo_groupfield))


  # if (report_settings[["fc_cyclicals"]] == TRUE){
  #   # create a doy field so that we can use a cyclical spline
  #   epi_lag <- dplyr::mutate(epi_lag, doy = as.numeric(format(.data$obs_date, "%j")))
  # }

  if (report_settings[["fc_splines"]] == "modbs"){
    # create modified bspline basis in epi_lag file to model longterm trends
    epi_lag <- cbind(epi_lag, truncpoly(x=epi_lag$obs_date,
                                        degree=6,
                                        maxobs=max(epi_lag$obs_date[epi_lag$input==1], na.rm=TRUE)))
  }

  #filter to input data for model building
  epi_input <- epi_lag %>% dplyr::filter(.data$input == 1)


  ## If model_cached is NOT given, then create model / run regression
  if (is.null(report_settings[["model_cached"]])){

    # Model building switching point

    regress <- build_model(fc_model_family,
                           quo_groupfield,
                           epi_lag,
                           report_settings,
                           #calc/internal
                           env_variables_used)

  } else {
    #if model_cached given, then use that as regress instead of building a new one (above)

    model_cached <- report_settings[["model_cached"]]

    #message with model input
    message("Using given cached ", model_cached$model_info$fc_model_family, " model, created ",
            model_cached$model_info$date_created, ", with epidemiological data up through ",
            model_cached$model_info$known_epi_range$max, ".")

    regress <- model_cached$model_obj
  }

  ## If model run, return regression object to run_forecast() at this point
  if (report_settings[["model_run"]]){
    return(regress)
  }

  # ## Error check all model results if using batch_bam/tp
  # if (report_settings[["fc_splines"]] == "tp"){
  #   check_bb_models(regress)
  # }


  ## Creating predictions switching point on model choice
  preds <- create_predictions(fc_model_family,
                              report_settings,
                              regress,
                              epi_lag,
                              env_variables_used,
                              req_date)


  #now cbind to get ready to return
  epi_preds <- cbind(epi_lag %>%
                       dplyr::filter(.data$obs_date <= req_date),
                     #column will be named preds
                     as.data.frame(preds)) %>%
    #and convert factor back to character for the groupings (originally converted b/c of bam/gam requirements)
    dplyr::mutate(!!rlang::quo_name(quo_groupfield) := as.character(!!quo_groupfield)) %>%
    #remake into tibble
    tibble::as_tibble()


  if (report_settings[["dev_fc_fit_freq"]] == "once"){
    #for single model fit, this has all the data we need,
    # trimming to report dates will happen later AFTER event detection
    date_preds <- epi_preds
      #%>% dplyr::filter(.data$obs_date >= report_dates$full$min)
  } else if (report_settings[["dev_fc_fit_freq"]] == "week"){
    #prediction of interest are last ones (equiv to req_date) per groupfield
    date_preds <- epi_preds %>%
      dplyr::group_by(!!quo_groupfield) %>%
      dplyr::filter(.data$obs_date == req_date)
    #note 'week' fits are limit to report period only, and workaround for farrington spin up will not work
  }

  forecast_reg_results <- create_named_list(date_preds,
                                            regress)
}


#'Build the appropriate model
#'
#'@param epi_input Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), with column marking if "known"
#'  data and groupings converted to factors.
#'@param env_variables_used a list of environmental variables that will be used in the
#'  modeling (had to be listed in model variables input file and present the
#'  env_data and env_info datasets)
#'
#'@inheritParams run_epidemia
#'@inheritParams run_forecast
#'
#'@return Regression object
#'
#'
build_model <- function(fc_model_family,
                        quo_groupfield,
                        epi_input,
                        report_settings,
                        #calc/internal
                        env_variables_used){

  #1. check and handle naive models
  # else is the user supplied model family
  #2. call build_equation to handle all the different equation pieces
  #3. call mgcv::bam or batchapply::bam_batch() as relevant

  #number of geographic area groupings
  n_groupings <- epi_input %>% dplyr::pull(!!quo_groupfield) %>% nlevels()
  #number of clusters
  n_clusters <- nlevels(epi_input$cluster_id)


  if (fc_model_family == "naive-persistence"){

    #naive model
    #persistence (carry forward)
    #no regression object

    #create "model" using known data.
    #Will fill down in create_predictions
    regress <- epi_input %>%
      #grouping by geographical unit
      dplyr::group_by(!!quo_groupfield) %>%
      #prediction is 1 lag (previous week)
      #preds is name of value from regression models
      dplyr::mutate(preds = dplyr::lag(.data$cases_epidemiar, n = 1))

  } else if (fc_model_family == "naive-averageweek"){

    #naive model
    #average of week of year (from historical data)
    #not a regression object

    #create "model" (averages) using known data.
    regress <- epi_input %>%
      #calculate averages per geographic group per week of year
      dplyr::group_by(!!quo_groupfield, .data$week_epidemiar) %>%
      dplyr::summarize(preds = mean(.data$cases_epidemiar, na.rm = TRUE))


  } else {
    #user supplied model family

    #Formula override: developer mode
    if (!is.null(report_settings[["dev_fc_formula"]])){
      message("DEVELOPER: Using user-supplied formula: ", report_settings[["dev_fc_formula"]])
      reg_eq <- report_settings[["dev_fc_formula"]]
      #note, if using formula override AND cyclicals,
      # dev users should put fc_cyclicals = TRUE, else message about discrete ignored.
      # dev users also need to set fc_splines appropriately

    } else {
      #build equation
      reg_eq <- build_equation(quo_groupfield,
                               epi_input,
                               report_settings,
                               n_groupings,
                               n_clusters,
                               env_variables_used)
    }

  #run the regression
    message("Creating regression model...")

  if (report_settings[["fc_splines"]] == "modbs"){
    if (report_settings[["fc_cyclicals"]]) {
      #yes cyclicals
      regress <- mgcv::bam(reg_eq,
                           data = epi_input,
                           family = fc_model_family,
                           control = mgcv::gam.control(trace=FALSE),
                           discrete = TRUE,
                           nthreads = report_settings[["fc_nthreads"]])
    } else {
      #no cyclicals
      regress <- mgcv::bam(reg_eq,
                           data = epi_input,
                           family = fc_model_family,
                           control = mgcv::gam.control(trace=FALSE))
    }
    #end modbs
  } else if (report_settings[["fc_splines"]] == "tp"){

    #tibble to dataframe, and turn all env wide data into each own sub matrix
    epi_input_tp <- format_lag_ca(epi_input,
                                  env_variables_used,
                                  report_settings)

    # create a cluster for clusterapply to use
    bb_cluster <- parallel::makeCluster(max(1, (report_settings[["ncores"]]-1), na.rm = TRUE))

    regress <- clusterapply::batch_bam(data = epi_input_tp,
                                       bamargs = list("formula" = reg_eq,
                                                      "family" = fc_model_family,
                                                      "discrete" = TRUE),
                                       over = "cluster_id",
                                       cluster = bb_cluster)

    #stop the cluster (if model run, won't use again,
    #  so starts and ends for modeling building or predictions)
    parallel::stopCluster(bb_cluster)

  } #end thin plate

  } #end else user supplied family

return(regress)

} # end build_model()


#'Create the appropriate regression equation.
#'
#'@param epi_input Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), with groupings as a factor,
#'  trimmed to data being used to create the model
#'@param n_groupings Count of the number of geographic groups (groupfield) in
#'  total.
#'@param n_clusters Count of the number of clusters in total
#'@param env_variables_used a list of environmental variables that will be used in the
#'  modeling (had to be listed in model variables input file and present the
#'  env_data and env_info datasets)
#'
#'@inheritParams run_epidemia
#'@inheritParams run_forecast
#'
#'@return A formula to be used in the regression call, built based on settings
#'  for cyclicals, spline type, and the number of geographic groupings present.
#'
#'
build_equation <- function(quo_groupfield,
                           epi_input,
                           report_settings,
                           n_groupings,
                           n_clusters,
                           env_variables_used){

  #switch split between modbs and tp spline options
  # equation depends on spline choice, cyclical choice, # (>1 or not) geogroups, # (>1 or not) clusters

  if (report_settings[["fc_splines"]] == "modbs"){

    #message("Creating equation for modified b-splines....")

    #create variable bandsummaries equation piece
    #  e.g. 'bandsummaries_{var1} * cluster_id' for however many env var bandsummaries there are
    bandsums_list <- grep("bandsum_*", colnames(epi_input), value = TRUE)
    bandsums_cl_list <- paste(bandsums_list, ": cluster_id")
    #need variant without known multiplication if <= 1 clusters
    if (n_clusters > 1) {
      bandsums_eq <- glue::glue_collapse(bandsums_cl_list, sep =" + ")
    } else {
      bandsums_eq <- glue::glue_collapse(bandsums_list, sep = " + ")
    }

    # get list of modbspline reserved variables and format for inclusion into model
    modb_list <- grep("modbs_reserved_*", colnames(epi_input), value = TRUE)
    # variant depending on >1 geographic area groupings
    if (n_groupings > 1){
      modb_list_grp <- paste(modb_list, ":", rlang::quo_name(quo_groupfield))
      modb_eq <- glue::glue_collapse(modb_list_grp, sep = " + ")
    } else {
      modb_eq <- glue::glue_collapse(modb_list, sep = " + ")
    }

    #cyclical or not
    if (report_settings[["fc_cyclicals"]]) {
      #TRUE, include cyclicals

      message("Including seasonal cyclicals into model...")

      #build equation

      #need different formulas if 1+ or only 1 geographic grouping
      if (n_groupings > 1){
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          rlang::quo_name(quo_groupfield),
                                          " + s(doy, bs=\"cc\", by=",
                                          rlang::quo_name(quo_groupfield),
                                          ") + ",
                                          modb_eq, " + ",
                                          bandsums_eq))
      } else {
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          "s(doy, bs=\"cc\") + ",
                                          modb_eq, " + ",
                                          bandsums_eq))
      }
    } else {
      # FALSE, no cyclicals

      #build equation

      #need different formulas if 1+ or only 1 geographic grouping
      if (n_groupings > 1){
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          rlang::quo_name(quo_groupfield), " + ",
                                          modb_eq, " + ",
                                          bandsums_eq))
      } else {
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          modb_eq, " + ",
                                          bandsums_eq))
      }
    } #end cyclicals if else

    #end if modbs
  } else if (report_settings[["fc_splines"]] == "tp"){

    message("Creating equation using thin plate splines...")

    #create s({}, by = {}, bs = 'tp', id = {unique})
    # numericdate column for long-term trend
    # lag column-matrix for environmental variables
    # ids need to be unique for each, but do not have to be sequential
      # id = 1 reserved for cyclicals which may or may not be present
      # id = 2 reserved for long-term trend
      # id = 3+ for each of the environmental variables

    #for geogroup
    #need different formulas if 1+ or only 1 geographic grouping
    tp_geo_eq <- if (n_groupings > 1){
      paste0("s(numericdate, by = ", rlang::quo_name(quo_groupfield),
             ", bs = \'tp\', id = 2)")
    } else {
      paste0("s(numericdate, ", "bs = \'tp\', id = 2)")
    }

    #for each env var
    idn_var <- seq(from = 3, to = (3-1+length(env_variables_used)))
    tp_env_eq_list <- paste0("s(lag, by = ", env_variables_used,
                             ", bs = \'tp\', id = ", idn_var, ")")
    tp_env_eq <- glue::glue_collapse(tp_env_eq_list, sep = " + ")


    if (report_settings[["fc_cyclicals"]]) {
      #TRUE, include cyclicals

      message("Including seasonal cyclicals into model...")

      #build formula

      #need different formulas if 1+ or only 1 geographic grouping
      if (n_groupings > 1){
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          rlang::quo_name(quo_groupfield),
                                          #cyclical
                                          " + s(doy, bs=\"cc\", by=",
                                          rlang::quo_name(quo_groupfield),
                                          ", id = 1) + ",
                                          #tp
                                          tp_geo_eq, " + ",
                                          tp_env_eq))
      } else {
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          "s(doy, bs=\"cc\", id = 1) + ",
                                          tp_geo_eq, " + ",
                                          tp_env_eq))
      }
    } else {
      # FALSE, no cyclicals

      #build formula

      #need different formulas if 1+ or only 1 geographic grouping
      if (n_groupings > 1){
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          rlang::quo_name(quo_groupfield), " + ",
                                          tp_geo_eq, " + ",
                                          tp_env_eq))
      } else {
        reg_eq <- stats::as.formula(paste("cases_epidemiar ~ ",
                                          tp_geo_eq, " + ",
                                          tp_env_eq))
      }
    } #end else cyclicals


  } #end splines tp

  #return
  reg_eq

} #end build_equation


#'Create the appropriate predictions/forecasts.
#'
#'@param regress The regression object, either the user-supplied one from
#'  `report_settings$model_cached`, or the one just generated.
#'@param epi_lag Epidemiological dataset with basis spline summaries of the
#'  lagged environmental data (or anomalies), with groupings as a factor.
#'@param env_variables_used a list of environmental variables that will be used in the
#'  modeling (had to be listed in model variables input file and present the
#'  env_data and env_info datasets)
#'@param req_date The end date of requested forecast regression. When fit_freq
#'  == "once", this is the last date of the full report, the end date of the
#'  forecast period.
#'
#'@inheritParams run_epidemia
#'
#'@return A dataset from predict() using the regression object generated in
#'  build_model or a newly created one. The dataset includes the
#'  predicted/forecast values through the end of the report requested.
#'
#'
create_predictions <- function(fc_model_family,
                               report_settings,
                               regress,
                               epi_lag,
                               env_variables_used,
                               req_date){


  #handle naive models
  if (fc_model_family == "naive-persistence"){

    message("Creating predictions using persistence naive model...")

    #persistence model just carries forward the last known value
    #the important part is the forecast / trailing end part
    #manipulating to be in quasi-same format as the other models return

    #regress is a tibble not regression object here
    # has a variable preds with lag of 1 on known data
    #epi_lag has the newer rows
    preds <- epi_lag %>%
      #filter to requested date
      dplyr::filter(.data$obs_date <= req_date) %>%
      #join to get "preds" values from "model"
      #join on all shared columns (i.e. everything in regress not "preds") to prevent renaming
      dplyr::left_join(regress, by = names(regress)[!names(regress) %in% c("preds")]) %>%
      #important at end/fc section, when we fill down
      tidyr::fill(.data$preds, .direction = "down") %>%
      #format into nominal regression predict output
      dplyr::select(.data$preds) %>%
      as.data.frame()

  } else if (fc_model_family == "naive-averageweek"){

    message("Creating predictions using average week of year naive model...")

    #average week null model calculates the average cases of that
    # week of year from historical data
    #manipulating to be in quasi-same format as the other models return

    #regress is the averages per week of year from known data

    epi_lag <- epi_lag %>%
      #filter to requested date
      dplyr::filter(.data$obs_date <= req_date)

    #join back
    preds <- epi_lag %>%
      #join to get average values
      #join on all shared columns (i.e. everything in regress not "preds") to prevent renaming
      # and so don't need column names not passed into this function
      dplyr::left_join(regress, by = names(regress)[!names(regress) %in% c("preds")]) %>%
      #format into nominal regression output
      dplyr::select(.data$preds) %>%
      as.data.frame()


  } else {
    #user supplied family, use predict.bam on regression object (regress)

    message("Creating predictions...")

    if (report_settings[["fc_splines"]] == "modbs"){
      #output prediction (through req_date)
      preds <- mgcv::predict.bam(regress,
                                 newdata = epi_lag %>% dplyr::filter(.data$obs_date <= req_date),
                                 type = "response",
                                 discrete = TRUE,
                                 n.threads = report_settings[["fc_nthreads"]],
                                 #default, and environmental predictors should not be NA
                                 #but setting explicit since it is assumed to return
                                 # the same number of rows as in newdata
                                 na.action = stats::na.pass)

    } else if (report_settings[["fc_splines"]] == "tp"){

      pred_input <- epi_lag %>%
        dplyr::filter(.data$obs_date <= req_date)

      #tibble to dataframe, and turn all env wide data into each own sub matrix
      pred_input_tp <- format_lag_ca(pred_input,
                                    env_variables_used,
                                    report_settings)


      # create a cluster for clusterapply to use
      pred_cluster <- parallel::makeCluster(max(1, (report_settings[["ncores"]]-1), na.rm = TRUE))

      preds <- clusterapply::predict.batch_bam(models = regress,
                                               predictargs = list("type"="response"),
                                               over = "cluster_id",
                                               newdata = pred_input_tp,
                                               cluster = pred_cluster)

      #stop the cluster
      parallel::stopCluster(pred_cluster)


    } #end else if fc_splines

  } #end else user supplied fc_family

  return(preds)

} #end create_predictions()
