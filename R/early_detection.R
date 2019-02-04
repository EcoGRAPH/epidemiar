# All run_epidemiar() subfunctions related to early detection

## Early Detection
#' Run early detection algorithm
#'
#'
run_early_detection <- function(epi_fc_data, quo_popfield, inc_per,
                                quo_groupfield, groupings,
                                ed_method, ed_control, report_dates){
  message("Running early detection")

  #only supporting Farrington Improved method from Surveillance right now,
  #leaving option open for expanding later
  if (ed_method == "Farrington") {
    ed_far_res <- run_farrington(epi_fc_data, quo_popfield, inc_per,
                                 quo_groupfield, groupings,
                                 ed_control, report_dates)
    return(ed_far_res)
  } else stop("Early Detection method not supported")
}

#' Run the Farrington early detection algorithm
#'
run_farrington <- function(epi_fc_data, quo_popfield, inc_per,
                           quo_groupfield, groupings,
                           ed_control, report_dates){
  ## Make sts objects
  #check about population offset
  # did the user set population offset
  if (!is.null(ed_control[["populationOffset"]])){
    #is it set to true
    if (ed_control[["populationOffset"]] == TRUE){
      #if so, did they give the population field
      if (!is.null(quo_popfield)){
        epi_stss <- make_stss(epi_fc_data, quo_popfield, quo_groupfield, groupings)
      } else stop("Population offset is TRUE, but population field not given")
    } else epi_stss <- make_stss(epi_fc_data, quo_popfield = NULL, quo_groupfield, groupings) #popoffset is FALSE, so no pop to sts
  } else epi_stss <- make_stss(epi_fc_data, quo_popfield = NULL, quo_groupfield, groupings) #if null, default is false, so pop = NULL
  #though note that pop is still a required field atm, so this path will fail later in early detection

  ## Set up new control list for Farrington (using their names)
  far_control <- list()

  #get evaluation period (range of row numbers)
  far_control[["range"]] <- seq(nrow(epi_stss[[1]]) - length(report_dates$full$seq) + 1, nrow(epi_stss[[1]]))

  #set number of years to go back in time
  # allow user set b, else calculate maximum number of years previous data available
  if (is.null(ed_control[["b"]])){
    #probably more properly done with isoyears and isoweeks, honestly.  <<>>
    daydiff <- difftime(report_dates$full$min, min(epi_fc_data$obs_date), "days") %>% as.numeric()
    far_control[["b"]] <- floor(daydiff / 365.242)

  } else far_control[["b"]] <- ed_control[["b"]]

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

  #run Farringtons
  far_res_list <- vector('list', length(epi_stss))
  for (i in 1:length(epi_stss)){
    far_res_list[[i]] <- surveillance::farringtonFlexible(epi_stss[[i]], control = far_control)
  }

  #results into output report data form
  far_res <- stss_res_to_output_data(stss_res_list = far_res_list, epi_fc_data,
                                     quo_popfield, inc_per,
                                     quo_groupfield, groupings, report_dates)

  far_res
}

#' Make the list of sts objects
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
    #only is passed in when popoffset = TRUE & population field is given
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
stss_res_to_output_data <- function(stss_res_list, epi_fc_data,
                                    quo_popfield, inc_per,
                                    quo_groupfield, groupings, report_dates){
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

  #recover population (for incidence calculations), not present if popoffset was FALSE
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
                  value = upperbound / !!quo_popfield * inc_per, #Incidence, from stss & epi_fc_data
                  lab = "Alert Threshold",
                  upper = NA,
                  lower = NA) %>%
    dplyr::select(!!quo_groupfield, obs_date, series, value, lab, upper, lower)

  #combine ed results
  ed <- rbind(ed_alert_res, ew_alert_res, ed_thresh_res)

  ed
}
