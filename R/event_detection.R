# All run_epidemiar() subfunctions related to early detection

## Early Detection
#'Main subfunction for running event detection algorithm.
#'
#'@param epi_fc_data Internal pass of epidemiological data complete with future
#'  forecast values.
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param inc_per Number for what unit of population the incidence should be
#'  reported in, e.g. incidence rate of 3 per 1000 people.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param ed_method Which method for early detection should be used ("Farrington"
#'  is only current option, or "None").
#'@param ed_control All parameters for early detection algorithm, passed through
#'  to that subroutine.
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param vt From match.arg evaluation of fc_control$value_type, whether to return
#'  epidemiological report values in "incidence" (default) or "cases".
#'@param mc From match.arg evaluation of model_choice. Reserved for future overrides on value_type depending on
#'  model choice selection.

#'
#'@return Returns a list of three generated series:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
run_event_detection <- function(epi_fc_data,
                                quo_popfield,
                                inc_per,
                                quo_groupfield,
                                groupings,
                                ed_method,
                                ed_control,
                                report_dates,
                                vt,
                                mc){
  message("Running early detection")

  #only supporting Farrington Improved method from Surveillance right now,
  #leaving option open for expanding later
  #Note: using exact matches, which can because of match.arg() at beginning on run_epidemia()

  if (ed_method == "farrington") {

    ed_far_res <- run_farrington(epi_fc_data,
                                 quo_popfield,
                                 inc_per,
                                 quo_groupfield,
                                 groupings,
                                 ed_control,
                                 report_dates,
                                 vt,
                                 mc)
    return(ed_far_res)

  } else if (ed_method == "none") {

    ed_far_res <- run_no_detection(epi_fc_data,
                                   quo_groupfield,
                                   report_dates)

  }

}

#' Run the Farrington early detection algorithm
#'
#'@param epi_fc_data Internal pass of epidemiological data complete with future
#'  forecast values.
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param inc_per Number for what unit of population the incidence should be
#'  reported in, e.g. incidence rate of 3 per 1000 people.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param ed_control All parameters for early detection algorithm, passed through
#'  to that subroutine.
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param vt From match.arg evaluation of fc_control$value_type, whether to return
#'  epidemiological report values in "incidence" (default) or "cases".
#'@param mc From match.arg evaluation of model_choice. Reserved for future overrides on value_type depending on
#'  model choice selection.

#'
#'@return Returns a list of three generated series from the Farrington algorithm:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
run_farrington <- function(epi_fc_data,
                           quo_popfield,
                           inc_per,
                           quo_groupfield,
                           groupings,
                           ed_control,
                           report_dates,
                           vt,
                           mc){
  ## Make sts objects
  #check about population offset
  # did the user set population offset
  if (!is.null(ed_control[["populationOffset"]])){
    #is it set to true
    if (ed_control[["populationOffset"]] == TRUE){
      #if so, did they give the population field
      if (!is.null(quo_popfield)){
        epi_stss <- make_stss(epi_fc_data,
                              quo_popfield,
                              quo_groupfield,
                              groupings)
      } else stop("Population offset is TRUE, but population field not given")
      #<<>> add to earlier input checks so fails early rather than later
    } else epi_stss <- make_stss(epi_fc_data,
                                 quo_popfield = NULL,
                                 quo_groupfield,
                                 groupings) #popoffset is FALSE, so no pop to sts
  } else epi_stss <- make_stss(epi_fc_data,
                               quo_popfield = NULL,
                               quo_groupfield,
                               groupings) #if null, default is false, so pop = NULL

  ## Set up new control list for Farrington (using their names)
  far_control <- list()

  #get evaluation period (range of row numbers)
  far_control[["range"]] <- seq(nrow(epi_stss[[1]]) - length(report_dates$full$seq) + 1, nrow(epi_stss[[1]]))

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
    ed_control[["trend"]] <- ed_control[["trend"]]
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

    #get the minimum integer year value to feed to Farrington control (cannot round up, must only request data that exists)
    far_control[["b"]] <- floor(yrdiff)

  } else far_control[["b"]] <- ed_control[["b"]]


  #run Farringtons
  far_res_list <- vector('list', length(epi_stss))
  for (i in 1:length(epi_stss)){
    far_res_list[[i]] <- surveillance::farringtonFlexible(epi_stss[[i]], control = far_control)
  }

  #results into output report data form
  far_res <- stss_res_to_output_data(stss_res_list = far_res_list,
                                     epi_fc_data,
                                     quo_popfield,
                                     inc_per,
                                     quo_groupfield,
                                     groupings,
                                     report_dates,
                                     vt,
                                     mc)

  far_res
}

#' Make the list of sts objects
#'
#'@param epi_fc_data Internal pass of epidemiological data complete with future
#'  forecast values.
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'
#'@return A list of surveillance time series (sts) objects,
#'one for each geographic grouping.
#'
make_stss <- function(epi_fc_data, quo_popfield, quo_groupfield, groupings){
  #create a list of surveillance::sts objects, one for each group
  stss <- vector('list', length(groupings))
  for (i in 1:length(groupings)){
    g <- groupings[i]
    g_data <- dplyr::filter(epi_fc_data, !!quo_groupfield == g) %>%
      #confirming sorting by date
      dplyr::arrange(obs_date)
    #Surveillance::sts() expecting a dataframe
    g_df <- as.data.frame(g_data)
    #get NA interpolated case field
    g_cases <- dplyr::select(g_df, cases_epidemiar) %>%
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
    stss[[i]] <- surveillance::sts(observed = g_cases,
                                   start = g_start,
                                   frequency = 52,  #weekly <<>>
                                   population = g_pop,
                                   epochAsDate = TRUE,
                                   epoch = as.numeric(g_data$obs_date))
  } # end for loop
  stss
}

#' Formats output data from sts result objects
#'
#'@param stss_res_list List of sts output object from Farrington algorithm.
#'@param epi_fc_data Internal pass of epidemiological data complete with future
#'  forecast values.
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param inc_per Number for what unit of population the incidence should be
#'  reported in, e.g. incidence rate of 3 per 1000 people.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param groupings A unique list of the geographic groupings (from groupfield).
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'@param vt From match.arg evaluation of fc_control$value_type, whether to return
#'  epidemiological report values in "incidence" (default) or "cases".
#'@param mc From match.arg evaluation of model_choice. Reserved for future overrides on value_type depending on
#'  model choice selection.

#'
#'@return Returns a list of three series from the Farrington sts result output:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
stss_res_to_output_data <- function(stss_res_list,
                                    epi_fc_data,
                                    quo_popfield,
                                    inc_per,
                                    quo_groupfield,
                                    groupings,
                                    report_dates,
                                    vt,
                                    mc){
  #take results of a surveillance event detection and reshape to output data format
  #stss to dfs
  stss_res_dfs <- lapply(stss_res_list, surveillance::as.data.frame)

  #add groupings back in (lost when creating sts object)
  #note importance of alphabetical order for groupings & initial sort at beginning of main function
  stss_res_grp <- mapply(cbind, stss_res_dfs, group_temp = groupings, SIMPLIFY = FALSE)
  #flatten out of list (now that we have the grouping labels)
  stss_res_flat <- do.call(rbind, stss_res_grp) %>%
    #fix group name field with dplyr programming
    dplyr::rename(!!quo_name(quo_groupfield) := group_temp) %>%
    #and convert to character for joining
    dplyr::mutate(!!rlang::quo_name(quo_groupfield) := as.character(!!quo_groupfield))

  #recover population (for incidence calculations), not present if popoffset was FALSE #<<pop>>
  stss_res_flat <- stss_res_flat %>%
    dplyr::left_join(epi_fc_data %>%
                       dplyr::select(!!quo_groupfield, !!quo_popfield, obs_date),
                     by = rlang::set_names(c(rlang::quo_name(quo_groupfield),
                                      "obs_date"),
                                    c(rlang::quo_name(quo_groupfield),
                                      "epoch")))

  #gather early detection (KNOWN - pre-forecast) event detection alert series
  ed_alert_res <- stss_res_flat %>%
    dplyr::filter(epoch %in% report_dates$known$seq) %>%
    dplyr::mutate(series = "ed",
                  obs_date = epoch,
                  value = alarm,
                  lab = "Early Detection Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #gather early WARNING event detection alert series
  ew_alert_res <- stss_res_flat %>%
    dplyr::filter(epoch %in% report_dates$forecast$seq) %>%
    dplyr::mutate(series = "ew",
                  obs_date = epoch,
                  value = alarm,
                  lab = "Early Warning Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #gather event detection threshold series
  ed_thresh_res <- stss_res_flat %>%
    dplyr::mutate(series = "thresh",
                  obs_date = epoch,
                  value = calc_return_value(cases = upperbound,
                                            c_quo_tf = FALSE,
                                            q_pop = quo_popfield,
                                            inc_per,
                                            vt,
                                            mc),
                  #value = upperbound / !!quo_popfield * inc_per, #Incidence, from stss & epi_fc_data
                  lab = "Alert Threshold",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #combine ed results
  ed <- rbind(ed_alert_res, ew_alert_res, ed_thresh_res)

  ed
}

#' Run No outbreak detection algorithm
#'
#'@param epi_fc_data Internal pass of epidemiological data complete with future
#'  forecast values.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param report_dates Internally generated set of report date information: min,
#'  max, list of dates for full report, known epidemiological data period,
#'  forecast period, and early detection period.
#'
#'@return Returns a list of three generated series with all NAs:
#' "ed" : early detection alerts (ed period of most recent epi data)
#' "ew" : early warning alerts (forecast/future portion)
#' "thresh" : threshold values per week
#'
run_no_detection <- function(epi_fc_data, quo_groupfield, report_dates){


  #early detection (KNOWN - pre-forecast) event detection alert series
  ed_alert_res <- epi_fc_data %>%
    dplyr::filter(obs_date %in% report_dates$known$seq) %>%
    dplyr::mutate(series = "ed",
                  value = NA_integer_,
                  lab = "Early Detection Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #gather early WARNING event detection alert series
  ew_alert_res <- epi_fc_data %>%
    dplyr::filter(obs_date %in% report_dates$forecast$seq) %>%
    dplyr::mutate(series = "ew",
                  value = NA_integer_,
                  lab = "Early Warning Alert",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #gather event detection threshold series
  ed_thresh_res <- epi_fc_data %>%
    dplyr::filter(obs_date %in% report_dates$full$seq) %>%
    dplyr::mutate(series = "thresh",
                  value = NA_real_,
                  lab = "Alert Threshold",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #combine ed results
  ed <- rbind(ed_alert_res, ew_alert_res, ed_thresh_res)

  ed


}
