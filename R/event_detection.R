# All run_epidemiar() subfunctions related to early detection

#'Main subfunction for running event detection algorithm.
#'
#'@param epi_fc_data Internal pass of epidemiological data complete with future
#'  forecast values.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param ed_method An extract of report_settings$ed_method after defaults have
#'  been applied - which method for early detection should be used ("farrington"
#'  or default "none", currently).
#'@param ed_control An extract of report_settings$ed_control - all parameters
#'  for early detection algorithm, passed through to that subroutine.
#'@param val_type An extract of report_settings$report_value_type after defaults
#'  applies - whether to return epidemiological report values in "incidence" or
#'  "cases" (default).
#'@param inc_per An extract of report_settings$report_inc_per after defaults
#'  applies - number for what unit of population the incidence should be
#'  reported in, e.g. incidence rate per 1000 people. Ignored when
#'  report_settings$report_value_type is 'cases'.
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param valid_run Internal TRUE/FALSE for whether this is part of a validation
#'  run.
#'
#'@return Returns a list of three generated series: "ed" : early detection
#'  alerts (ed period of most recent epi data) "ew" : early warning alerts
#'  (forecast/future portion) "thresh" : threshold values per week
#'
run_event_detection <- function(epi_fc_data,
                                quo_groupfield,
                                quo_popfield,
                                #rpt settings items
                                ed_method,
                                ed_control,
                                val_type,
                                inc_per,
                                #internal/calc
                                groupings,
                                report_dates,
                                valid_run){
  #message("Running early detection...")

  #only supporting Farrington Improved method from Surveillance right now,
  #leaving option open for expanding later
  #Note: using exact matches, which can because of match.arg() at beginning on run_epidemia()

  if (ed_method == "farrington") {

    message("Running early detection: Farrington...")
    ed_far_res <- run_farrington(epi_fc_data,
                                 quo_groupfield,
                                 quo_popfield,
                                 ed_control,
                                 val_type,
                                 inc_per,
                                 groupings,
                                 report_dates)
    return(ed_far_res)

  } else if (ed_method == "none") {

    if(!valid_run){
      message("Skipping early detection...")
    }
    ed_far_res <- run_no_detection(epi_fc_data,
                                   quo_groupfield,
                                   report_dates)

  }

}

#' Run the Farrington early detection algorithm
#'
#'@inheritParams run_event_detection
#'
#'@return Returns a list of three generated series from the Farrington algorithm:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
run_farrington <- function(epi_fc_data,
                           quo_groupfield,
                           quo_popfield,
                           ed_control,
                           val_type,
                           inc_per,
                           groupings,
                           report_dates){
  ## Make sts objects
  #check about population offset
  # did the user set population offset
  if (!is.null(ed_control[["populationOffset"]])){
    #is it set to true
    if (ed_control[["populationOffset"]] == TRUE){
      #if so, did they give the population field
      if (!is.null(quo_popfield)){
        epi_stss <- make_stss(epi_fc_data,
                              quo_groupfield,
                              quo_popfield,
                              groupings)
      } else stop("Population offset is TRUE, but population field not given")
      #<<>> add to earlier input checks so fails early rather than later?
    } else epi_stss <- make_stss(epi_fc_data,
                                 quo_groupfield,
                                 quo_popfield = NULL,
                                 groupings) #popoffset is FALSE, so no pop to sts
  } else epi_stss <- make_stss(epi_fc_data,
                               quo_groupfield,
                               quo_popfield = NULL,
                               groupings) #if null, default is false, so pop = NULL

  ## Set up new control list for Farrington (using their names)
  far_control <- list()

  #get evaluation period (range of row numbers)
  far_control[["range"]] <- seq(nrow(epi_stss[[1]]) - length(report_dates$full$seq) + 1,
                                nrow(epi_stss[[1]]))

  #test for all other parameters that can be passed onto Farrington flexible method
  # if not null, use user parameter, otherwise leave as null to use its defaults
  if (!is.null(ed_control[["w"]])){
    far_control[["w"]] <- ed_control[["w"]]
  }
  if (!is.null(ed_control[["reweight"]])){
    far_control[["reweight"]] <- ed_control[["reweight"]]
  }
  if (!is.null(ed_control[["weightsThreshold"]])){
    far_control[["weightsThreshold"]] <- ed_control[["weightsThreshold"]]
  }
  if (!is.null(ed_control[["alpha"]])){
    far_control[["alpha"]] <- ed_control[["alpha"]]
  }
  if (!is.null(ed_control[["trend"]])){
    far_control[["trend"]] <- ed_control[["trend"]]
  }
  if (!is.null(ed_control[["pThresholdTrend"]])){
    far_control[["pThresholdTrend"]] <- ed_control[["pThresholdTrend"]]
  }
  if (!is.null(ed_control[["limit54"]])){
    far_control[["limit54"]] <- ed_control[["limit54"]]
  }
  if (!is.null(ed_control[["powertrans"]])){
    far_control[["powertrans"]] <- ed_control[["powertrans"]]
  }
  if (!is.null(ed_control[["fitFun"]])){
    far_control[["fitFun"]] <- ed_control[["fitFun"]]
  }
  if (!is.null(ed_control[["populationOffset"]])){
    far_control[["populationOffset"]] <- ed_control[["populationOffset"]]
  }
  if (!is.null(ed_control[["noPeriods"]])){
    far_control[["noPeriods"]] <- ed_control[["noPeriods"]]
  }
  if (!is.null(ed_control[["pastWeeksNotIncluded"]])){
    far_control[["pastWeeksNotIncluded"]] <- ed_control[["pastWeeksNotIncluded"]]
  }
  if (!is.null(ed_control[["thresholdMethod"]])){
    far_control[["thresholdMethod"]] <- ed_control[["thresholdMethod"]]
  }

  #set number of years to go back in time
  # allow user set b, else calculate maximum number of years previous data available
  # includes allowance for window value, w # of weeks
  if (is.null(ed_control[["b"]])){

    #subtract window to earliest report date
    #    this will appropriately increase the amount of time needed without altering the actual available data date
    #    allow default w=3
    if (!is.null(ed_control[["w"]])){
      adjdt <- report_dates$full$min - lubridate::weeks(ed_control[["w"]])
    } else adjdt <- report_dates$full$min - lubridate::weeks(3)

    #calculate number of years difference between earliest available data date and adjusted report date "start"
    #    using interval(), because this allows time_lenth() to deal with leap years, etc.
    yrdiff <- lubridate::interval(min(epi_fc_data$obs_date), adjdt) %>%
      lubridate::time_length(unit = "years")

    #get the minimum integer year value to feed to Farrington control
    #(cannot round up, must only request data that exists)
    far_control[["b"]] <- floor(yrdiff)

  } else far_control[["b"]] <- ed_control[["b"]]


  ## input check / overrides:
    # #This test is useless now that implicit missing weeks are handled.
    # #do all groups have the same number of weeks? Farrington will error otherwise.
    # wks_diff_grps <- epi_fc_data %>%
    #   dplyr::group_by(!!quo_groupfield) %>%
    #   dplyr::count() %>%
    #   dplyr::pull(.data$n) %>%
    #   range() %>% diff()

  #only run if b > 0. If 0 full years available (or b=0 requested), then "none"
  if (far_control[["b"]] < 1) {

    message("Warning: Less than 1 full year of epidemiological data available or requested (ed_control$b). Cannot run Farrington, skipping event detection.")
    far_res <- run_no_detection(epi_fc_data,
                                quo_groupfield,
                                report_dates)


  }
  # else if (wks_diff_grps > .Machine$double.eps ^ 0.5){
  #   #do all groups have the same number of weeks? Farrington will error otherwise.
  #   #using small tolerance rather than == 0.
  #   message("Warning: Groups do not have the same number of weeks of epidemiological data. Cannot run Farrington, skipping event detection.")
  #   far_res <- run_no_detection(epi_fc_data,
  #                               quo_groupfield,
  #                               report_dates)
  # }
  else {
    #if all okay, then run Farrington

    #run Farringtons
    far_res_list <- vector('list', length(epi_stss))
    for (i in 1:length(epi_stss)){
      #far_res_list[[i]] <- surveillance::farringtonFlexible(epi_stss[[i]], control = far_control)
      far_res_list[[i]] <- tryCatch({
        #successful run will have Farrington results
        surveillance::farringtonFlexible(epi_stss[[i]], control = far_control)},
        error = function(e){
          #failed run
          #use pre-farrington sts for the appropriate evaluation weeks (range)
          #will have the correct format and matches the rest of the results
          #but will not have thresholds or alert values
          message(paste0("Farrington model failure on ", groupings[i],
                         "and will not have thresholds values or alerts for this group.",
                         "Continuing with remaining groups.",
                         "Error from Farrington:", e))
            #print(i); print(groupings[i]); print(e)
          epi_stss[[i]][far_control$range,]})
    }


    #results into output report data form
    far_res <- stss_res_to_output_data(stss_res_list = far_res_list,
                                       epi_fc_data,
                                       quo_groupfield,
                                       quo_popfield,
                                       val_type,
                                       inc_per,
                                       groupings,
                                       report_dates)

  }

  far_res
}

#' Make the list of sts objects
#'
#'@inheritParams run_event_detection
#'
#'@return A list of surveillance time series (sts) objects,
#'one for each geographic grouping.
#'
make_stss <- function(epi_fc_data,
                      quo_groupfield,
                      quo_popfield,
                      groupings){
  #create a list of surveillance::sts objects, one for each group
  stss <- vector('list', length(groupings))
  for (i in 1:length(groupings)){
    g <- groupings[i]
    g_data <- dplyr::filter(epi_fc_data, !!quo_groupfield == g) %>%
      #confirming sorting by date
      dplyr::arrange(.data$obs_date)
    #Surveillance::sts() expecting a dataframe
    g_df <- as.data.frame(g_data)
    #get NA interpolated case field
    g_cases <- dplyr::select(g_df, .data$cases_epidemiar) %>%
      #sts() likes matrices
      as.matrix()
    #if population field given, get population
    #only is passed in when popoffset = TRUE & population field is given #<<pop>>
    if (!is.null(quo_popfield)){
      g_pop <- dplyr::select(g_df, !!quo_popfield) %>%
        #sts() likes matrices
        as.matrix()
    } else g_pop <- NULL

    #start year and week values
    g_start <- c(g_df[ which(g_df$obs_date == min(g_df$obs_date)), "year_epidemiar"],
                 g_df[ which(g_df$obs_date == min(g_df$obs_date)), "week_epidemiar"])

    #make that group's sts object
    # due to R 3.6+ now verbose on S3 methods overwriting
    # wrapping first surveillance package call in suppressMesssages()
    # because somewhere it calls spatstats with print method for boxx class that overwrites cli's
    stss[[i]] <- suppressMessages(surveillance::sts(observed = g_cases,
                                                    start = g_start,
                                                    frequency = 52,  #weekly
                                                    population = g_pop,
                                                    epochAsDate = TRUE,
                                                    epoch = as.numeric(g_data$obs_date)))
  } # end for loop
  stss
}

#' Formats output data from sts result objects
#'
#'@param stss_res_list List of sts output object from Farrington algorithm.
#'@inheritParams run_event_detection
#'
#'@return Returns a list of three series from the Farrington sts result output:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
stss_res_to_output_data <- function(stss_res_list,
                                    epi_fc_data,
                                    quo_groupfield,
                                    quo_popfield,
                                    val_type,
                                    inc_per,
                                    groupings,
                                    report_dates){
  #take results of a surveillance event detection and reshape to output data format
  #stss to dfs
  stss_res_dfs <- lapply(stss_res_list, surveillance::as.data.frame)

  #add groupings back in (lost when creating sts object)
  #note importance of alphabetical order for groupings & initial sort at beginning of main function
  stss_res_grp <- mapply(cbind, stss_res_dfs, group_temp = groupings, SIMPLIFY = FALSE)
  #flatten out of list (now that we have the grouping labels)
  stss_res_flat <- do.call(rbind, stss_res_grp) %>%
    #fix group name field with dplyr programming
    dplyr::rename(!!rlang::as_name(quo_groupfield) := .data$group_temp) %>%
    #and convert to character for joining
    dplyr::mutate(!!rlang::as_name(quo_groupfield) := as.character(!!quo_groupfield))

  #recover population (for incidence calculations), not present if popoffset was FALSE
  #only if optional population field was given
  if (!rlang::quo_is_null(quo_popfield)) {
    stss_res_flat <- stss_res_flat %>%
      dplyr::left_join(epi_fc_data %>%
                         dplyr::select(!!quo_groupfield, !!quo_popfield, .data$obs_date),
                       by = rlang::set_names(c(rlang::as_name(quo_groupfield),
                                               "obs_date"),
                                             c(rlang::as_name(quo_groupfield),
                                               "epoch")))
  }


  #recover alert censor flag (when observed was NA in 'prev' report period)
  stss_res_flat <- stss_res_flat %>%
    dplyr::left_join(epi_fc_data %>%
                       dplyr::select(!!quo_groupfield, .data$obs_date, .data$censor_flag),
                     by = rlang::set_names(c(rlang::as_name(quo_groupfield),
                                             "obs_date"),
                                           c(rlang::as_name(quo_groupfield),
                                             "epoch")))


  #gather early detection (pre-forecast) event detection alert series
  #early detection alerts show for all time previous and including early detection period
  #"historical" alerts were wanted
  ed_alert_res <- stss_res_flat %>%
    dplyr::filter(.data$epoch %in% report_dates$prev$seq) %>%
    dplyr::mutate(series = "ed",
                  obs_date = .data$epoch,
                  value = .data$alarm,
                  lab = "Early Detection Alert",
                  upper = NA,
                  lower = NA) %>%
    #censor alarms to NA for when observed value was actually NA in 'prev' period
    dplyr::mutate(value = ifelse(.data$censor_flag == TRUE,
                                 NA_integer_,
                                 .data$value)) %>%
    # #surveillance returns an alarm value (0) for when observed is NA, we want NA in this case
    #   #this should no longer happen with blending of modelled values to force threshold generation
    # dplyr::mutate(value = ifelse(is.na(.data$observed),
    #                              NA_integer_,
    #                              .data$value)) %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)

  #gather early WARNING event detection alert series
  ew_alert_res <- stss_res_flat %>%
    dplyr::filter(.data$epoch %in% report_dates$forecast$seq) %>%
    dplyr::mutate(series = "ew",
                  obs_date = .data$epoch,
                  value = .data$alarm,
                  lab = "Early Warning Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)

  #gather event detection threshold series
  ed_thresh_res <- stss_res_flat %>%
    dplyr::filter(.data$epoch %in% report_dates$full$seq) %>%
    dplyr::mutate(series = "thresh",
                  obs_date = .data$epoch,
                  #value calculations change depending on report_value_type
                  #case_when is not viable because it evaluates ALL RHS
                  value = if(val_type == "cases"){
                      .data$upperbound
                    } else if (val_type == "incidence"){
                      .data$upperbound / !!quo_popfield * inc_per
                    } else {NA_real_},
                  # value = dplyr::case_when(
                  #   #if reporting in case counts
                  #   val_type == "cases" ~ upperbound,
                  #   #if incidence
                  #   val_type == "incidence" ~ upperbound / !!quo_popfield * inc_per,
                  #   #otherwise
                  #   TRUE ~ NA_real_),
                  lab = "Alert Threshold",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)


  #combine ed results
  ed <- rbind(ed_alert_res, ew_alert_res, ed_thresh_res)

  ed
}

#' Run No outbreak detection algorithm
#'
#'@inheritParams run_event_detection
#'
#'@return Returns a list of three generated series with all NAs:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
run_no_detection <- function(epi_fc_data,
                             quo_groupfield,
                             report_dates){


  #early detection (pre-forecast, obstensibly though not nec. known) event detection alert series
  ed_alert_res <- epi_fc_data %>%
    dplyr::filter(.data$obs_date %in% report_dates$prev$seq) %>%
    dplyr::mutate(series = "ed",
                  value = NA_integer_,
                  lab = "Early Detection Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)

  #gather early WARNING event detection alert series
  ew_alert_res <- epi_fc_data %>%
    dplyr::filter(.data$obs_date %in% report_dates$forecast$seq) %>%
    dplyr::mutate(series = "ew",
                  value = NA_integer_,
                  lab = "Early Warning Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)

  #gather event detection threshold series
  ed_thresh_res <- epi_fc_data %>%
    dplyr::filter(.data$obs_date %in% report_dates$full$seq) %>%
    dplyr::mutate(series = "thresh",
                  value = NA_real_,
                  lab = "Alert Threshold",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, .data$obs_date, .data$series, .data$value, .data$lab, .data$upper, .data$lower)

  #combine ed results
  ed <- rbind(ed_alert_res, ew_alert_res, ed_thresh_res)

  ed


}
