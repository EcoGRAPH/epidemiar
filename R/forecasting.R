# All run_epidemiar() subfunctions related to forecasting


## Forecasting
#' Runs the forecast modeling
#' @export
#'
run_forecast <- function(epi_data, quo_popfield, quo_groupfield, groupings,
                         env_data, quo_obsfield, quo_valuefield, env_variables,
                         fc_control, env_ref_data, env_info, report_dates, week_type){
  message("Running forecasts")

  #set up default parallel processing number of cores to use number
  #if user-supplied, use that, otherwise create a default number
  #used in anomalize_env() and forecast_regression()
  if (!is.null(fc_control[["ncores"]])){
    ncores <- fc_control[["ncores"]]
  } else ncores <- max(detectCores(logical=FALSE) - 1, 1)

  # create the modeling variable
  # epi_data <- mutate(epi_data, logcase = log(cases_epidemiar + 1))
  epi_data <- dplyr::mutate(epi_data, modeledvar = floor(cases_epidemiar))

  # trim to the needed env variables as dictated by the model
  env_data <- pull_model_envvars(env_data, quo_obsfield, fc_control)
  #create alphabetical list of ONLY USED unique environmental variables
  env_variables_used <- dplyr::pull(env_data, !!quo_obsfield) %>% unique() %>% sort()

  # extend data into future, for future forecast portion
  env_data_extd <- extend_env_future(env_data, quo_groupfield, groupings, quo_obsfield, quo_valuefield,
                                     env_ref_data, env_info, env_variables_used, report_dates, week_type)
  epi_data_extd <- extend_epi_future(epi_data, quo_popfield, quo_groupfield,
                                     groupings, report_dates)

  # format the data for forecasting algorithm
  env_fc <- env_format_fc(env_data_extd, quo_groupfield, quo_obsfield)
  epi_fc <- epi_format_fc(epi_data_extd, quo_groupfield, fc_control)

  # anomalizing the environ data
  env_fc <- anomalize_env(env_fc, quo_groupfield, quo_obsfield, ncores)

  # create the lags
  epi_lag <- lag_environ_to_epi(epi_fc, quo_groupfield, groupings,
                                env_fc, env_variables_used, laglen = fc_control$lag_length)

  # for each week of report, run forecast
  # initialize: prediction returns 4 columns
  preds_catch <- data.frame()
  #loop by week
  for (w in seq_along(report_dates$full$seq)){
    message("Forecasting week ", w, " starting at ", Sys.time())
    dt <- report_dates$full$seq[w]
    dt_preds <- forecast_regression(epi_lag, quo_groupfield, groupings,
                                    env_variables_used,
                                    req_date = dt,
                                    ncores)
    preds_catch <- rbind(preds_catch, as.data.frame(dt_preds))
  }

  #create cases value from log
  # preds_catch <- preds_catch %>%
  #   mutate(fc_cases = exp(fit.fit) + 1,
  #          fc_cases_upr = exp(fit.upr) + 1,
  #          fc_cases_lwr = exp(fit.lwr) + 1)

  # Since we're not doing prediction intervals and since we're modeling untransformed data, this is
  # just an identity transformation, but we retain the variables for compatibility and perhaps further

  # debugging
  #print(head(preds_catch))

  # expansion. This is just a guess at how it might work.
  preds_catch <- preds_catch %>%
    dplyr::mutate(fc_cases = fit,
                  fc_cases_upr = fit+1.96*sqrt(fit),
                  fc_cases_lwr = fit-1.96*sqrt(fit))

  # extract fc series into report format
  fc_res <- preds_catch %>%
    dplyr::mutate(series = "fc",
                  value = fc_cases / !!quo_popfield * 1000, #Incidence per 1000
                  lab = "Forecast Trend",
                  upper = fc_cases_upr / !!quo_popfield * 1000,
                  lower = fc_cases_lwr / !!quo_popfield * 1000) %>%
    dplyr::select(!!quo_groupfield, Date, series, value, lab, upper, lower)

  # return list with res and other needed items
  fc_res_full <- create_named_list(fc_epi = preds_catch, fc_res,
                                   env_data_extd, env_variables_used)
}

#forecasting helper functions
# this creates a modified b-spline basis, which is still (piecewise) polynomial, so
# we will keep this name
#' Truncates poly
#' @export
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

#' Pull only model env variables
#' @export
#'
pull_model_envvars <- function(env_data, quo_obsfield, fc_control){
  #extract values (environment ids) of vars.1 to vars.x from GA output
  model_var_ids <- fc_control$model[, which(grepl("vars.", names(fc_control$model)))] %>%
    t() %>% as.data.frame()

  #lookup ids against envcodes and get environmental names as in obsfield (will need to add checks) <<>>
  model_vars <- model_var_ids %>%
    dplyr::left_join(fc_control$envcodes, by = c("V1" = "environ_var_id")) %>%
    dplyr::pull(!!quo_obsfield)

  #filter env_data for those model_vars
  env_data <- env_data %>%
    dplyr::filter(!!quo_obsfield %in% model_vars)
}

#' Extend environmental data into the future
#' @export
#'
extend_env_future <- function(env_data, quo_groupfield, groupings, quo_obsfield, quo_valuefield,
                              env_ref_data, env_info, env_variables_used, report_dates, week_type){
  #extending time into future forecast
  #if <=4 weeks, use persistence of mean of previous week before NA
  #if >4 weeks, use progressive blend of mean of previous week & historical means
  #Note: env data is DAILY, env ref is WEEKLY

  ## 1. Complete given env_data, as could have ragged env data -
  # these were marked interpolation, though technically should be considered extrapolation instead
  env_data <- env_data %>%
    #mark which ones still have val_epidemiar as NA as extended series
    dplyr::mutate(data_source = ifelse(is.na(val_epidemiar), "Extended", data_source))

  #find 1st NA, then take mean of previous week, input for that day
  env_known <- env_last_week_mean(env_df = env_data, env_variables_used, quo_groupfield, quo_obsfield, groupings)
  #fill in any following days with mean of previous week value
  env_known_fill <- env_fill_down(env_df = env_known, quo_groupfield, quo_obsfield, quo_valuefield)

  #trim env var to end of forecast period, don't need data past forecast period if running historical.
  env_known_fill <- env_known_fill %>%
    dplyr::filter(Date <= report_dates$forecast$max)

  ## 2. if need future data, now add value-empty future dataframe (to fill in)
  #the minimum date of the maximum dates for each env var data in each groupings
  min_known_env <- env_known_fill %>%
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
    dplyr::summarize(max_dates = max(Date, na.rm = TRUE)) %>%
    dplyr::pull(max_dates) %>% min()
  #if missing at least some of the env data for the forecast period
  if (report_dates$forecast$max > min_known_env) {
    #create combination of all groups, env vars, and future dates (DAILY)
    env_future <- tidyr::crossing(Date = seq.Date(min_known_env + 1, report_dates$forecast$max, 1),
                                  group_temp = groupings,
                                  obs_temp = env_variables_used)
    #and fix names with NSE
    env_future <- env_future %>%
      dplyr::rename(!!quo_name(quo_groupfield) := group_temp,
                    !!quo_name(quo_obsfield) := obs_temp)
    #could have ragged env data per variable per grouping
    #so, antijoin with env_known_fill first to get actually missing future entries
    env_future_missing <- env_future %>%
      dplyr::anti_join(env_known_fill, by = set_names(c(quo_name(quo_groupfield),
                                                        quo_name(quo_obsfield),
                                                        "Date"),
                                                      c(quo_name(quo_groupfield),
                                                        quo_name(quo_obsfield),
                                                        "Date")))

    #bind with existing data (NAs for everything else in env_future)
    extended_env <- dplyr::bind_rows(env_known_fill, env_future_missing) %>%
      #mark which are about to be filled in
      dplyr::mutate(data_source = ifelse(is.na(val_epidemiar), "Extended", data_source))

    ## 2.1. if forecast <= 4 weeks, then just use persistance and carry down
    if (length(report_dates$forecast$seq) <= 4){

      #find 1st NA (1st day in future here), then take mean of previous week, input for that day
      extended_env_mean <- env_last_week_mean(env_df = extended_env, env_variables_used, quo_groupfield, quo_obsfield, groupings)
      #fill down all NAs with last known non-NA value
      extended_env_fill <- env_fill_down(env_df = extended_env_mean, quo_groupfield, quo_obsfield, quo_valuefield)

    } #end 2.1 if fc wks <= 4

    ## 2.2. if forecast > 4 weeks, calculate blend of persistance and historical week to daily env
    #last week mean
    if (length(report_dates$forecast$seq) > 4){

      ## For forceast of daily values per week
      #are a blend of that 'last week mean' in report_dates$known$max+1 & historical ref during that week
      #rem: env data is DAILY, env ref is WEEKLY
      #rem: In the env_ref, some var may be mean, some may be sum, etc. for the week
      #Plan: get value for that week (week_epidemiar last date of week, +1 to start of next week), then use fill down for rest of daily values in week

      #set up first week kick off value (mean last week of known data)
      #find 1st NA (1st day in future here), then take mean of previous week, input for that day
      extended_env_mean <- env_last_week_mean(env_df = extended_env, env_variables_used, quo_groupfield, quo_obsfield, groupings)

      #prep ref data - get only used vars
      env_ref_varused <- env_ref_data %>%
        dplyr::filter(!!quo_obsfield %in% env_variables_used)

      #joins for ref summary type, and summary for week
      extended_env_join <- extended_env_mean %>%
        #add week, year fields
        epidemiar::add_datefields(week_type) %>%
        #get reference/summarizing method from user supplied env_info
        dplyr::left_join(env_info %>%
                           dplyr::select(!!quo_obsfield, reference_method),
                         by = set_names(rlang::quo_name(quo_obsfield),
                                        rlang::quo_name(quo_obsfield))) %>%
        #get weekly ref value
        dplyr::left_join(env_ref_varused %>%
                           dplyr::select(!!quo_obsfield, !!quo_groupfield, week_epidemiar, ref_value),
                         #NSE fun
                         by = set_names(c(rlang::quo_name(quo_groupfield),
                                          rlang::quo_name(quo_obsfield),
                                          "week_epidemiar"),
                                        c(rlang::quo_name(quo_groupfield),
                                          rlang::quo_name(quo_obsfield),
                                          "week_epidemiar")))

      #Monster loop
      #first day of week (so can use fill down, later)
      d1_wks <- report_dates$forecast$seq - 6
      #number of weeks in forecast
      nwks <- length(report_dates$forecast$seq)
      #result dataframe initialize
      extended_env_calc <- data.frame()
      #loop b/c too complicated for just dplyr group_by
      for (g in groupings){
        g_data <- extended_env_join %>% dplyr::filter(!!quo_groupfield == g)
        for (e in env_variables_used){
          g_e_data <- g_data %>% dplyr::filter(!!quo_obsfield == e)
          last_mean <- g_e_data[g_e_data$Date == report_dates$known$max + 1, "val_epidemiar"] %>% dplyr::pull()
          for (w in seq_along(d1_wks)){
            #date value
            d <- d1_wks[w]
            ref <- g_e_data[g_e_data$Date == d, "ref_value"] %>% dplyr::pull()
            persist_part <- ((nwks - w) / nwks) * last_mean
            #historical part also depends on env_ref summary method - if mean, then use ref as in. if sum, then ref/7.
            r_meth <- g_e_data[g_e_data$Date == d, "reference_method"]
            hx_part <- dplyr::case_when(
              r_meth == "mean" ~ (w / nwks) * ref,
              r_meth == "sum"  ~ (w / nwks) * ref/7,
              TRUE             ~ (w / nwks) * ref) #default as if mean
            #only fill NAs, i.e. don't fill any existing data
            if (is.na(g_e_data[g_e_data$Date == d, "val_epidemiar"])){
              g_e_data[g_e_data$Date == d, "val_epidemiar"] <- persist_part + hx_part
            } #end if NA
          } #end x loop

          #create updated tbl (either adding original, or modified depending on if above)
          extended_env_calc <- rbind(extended_env_calc, g_e_data)
        } #end e loop
      } #g loop

      #remove extra columns to keep consistency from other paths
      extended_env_calc <- extended_env_calc %>%
        dplyr::select(-reference_method, -ref_value, -week_epidemiar, -year_epidemiar)

      #then fill down to get the rest of the columns filled appropriately
      extended_env_fill <- env_fill_down(env_df = extended_env_calc, quo_groupfield, quo_obsfield, quo_valuefield)
    } #end if fc wks >4

  } else {
    #end if report_dates$forecast$max > min_known_env
    #if already have all needed 'future' env data, filled just in case of ragged env data
    extended_env_fill <- env_known_fill
  }
  extended_env_fill
}

#' Calculate mean of last week env values
#' @export
#'
env_last_week_mean <- function(env_df, env_variables_used, quo_groupfield, quo_obsfield, groupings){
  #gets mean of previous week of daily env data, puts in first NA slot
  #per grouping, per env variable

  #find 1st NA, then take mean of previous week, input for that week
  env_mean <- data.frame()
  #loop b/c too complicated for just dplyr group_by
  for (g in groupings){
    g_data <- env_df %>% dplyr::filter(!!quo_groupfield == g)
    for (e in env_variables_used){
      g_e_data <- g_data %>% dplyr::filter(!!quo_obsfield == e)
      #check if NA value or not in set
      if (any(is.na(g_e_data$val_epidemiar))){
        #index of first NA value
        #in this step, this is often a known calculable date depending on when called,
        #but since this works, no need to change the algorithm
        indx <- min(which(is.na(g_e_data$val_epidemiar)), na.rm = TRUE)
        #calc mean of previous 7 row (days) -- will be towards the end, so don't have to worry about not having 7 days prior
        mean_prev_wk <- g_e_data[(indx-1):(indx-7), "val_epidemiar"] %>%
          dplyr::pull() %>% mean(na.rm = TRUE)
        #add value into appropriate place
        g_e_data[indx, "val_epidemiar"] <- mean_prev_wk
      }
      #create updated tbl (either adding original, or modified depending on if above)
      env_mean <- rbind(env_mean, g_e_data)
    } #e loop
  } #g loop
  env_mean
}

#' Fill env data down
#' @export
#'
env_fill_down <- function(env_df, quo_groupfield, quo_obsfield, quo_valuefield){
  #to fill down values (except for original value field) for remaining length of dataset given
  env_filled <- env_df %>%
    #order very important for filling next step
    dplyr::arrange(!!quo_groupfield, !!quo_obsfield, Date) %>%
    dplyr::group_by(!!quo_groupfield, !!quo_obsfield) %>%
    #fill everything except original value field
    tidyr::fill(dplyr::everything(), -!!quo_valuefield, .direction = "down") %>%
    ungroup()
  env_filled
}

#' Extend epidemiology dataframe into future
#' @export
#'
extend_epi_future <- function(epi_data, quo_popfield, quo_groupfield, groupings, report_dates){
  #extended epi data into future dates
  #for use in modeling later (results will be put elsewhere), this is for env and lags and modeling dataset
  epi_future <- tidyr::crossing(Date = report_dates$forecast$seq,
                                group_temp = groupings)
  #and fix names with NSE
  epi_future <- epi_future %>%
    dplyr::rename(!!quo_name(quo_groupfield) := group_temp)

  #bind with exisiting data (NAs for everything else in epi_future)
  extended_epi <- dplyr::bind_rows(epi_data, epi_future) %>%
    dplyr::arrange(!!quo_groupfield, Date)

  #fill population down
  extended_epi <- tidyr::fill(extended_epi, !!quo_popfield, .direction = "down")

  extended_epi
}

#' Format env data for modeling
#' @export
#'
env_format_fc <- function(env_data_extd, quo_groupfield, quo_obsfield){
  env_spread <- env_data_extd %>%
    dplyr::mutate(numericdate = as.numeric(Date)) %>%
    dplyr::select(!!quo_groupfield, !!quo_obsfield, Date, numericdate, val_epidemiar) %>%
    tidyr::spread(key = !!quo_obsfield, value = val_epidemiar)

  #no longer needed with updated extend_env_future()
  #env_spread <- fill(env_spread, everything(), .direction = "down")

  env_spread
}

#' Format epi data for modeling
#' @export
#'
epi_format_fc <- function(epi_data_extd, quo_groupfield, fc_control){
  #cluster information from model
  cluster_groups <- fc_control$model[, which(grepl("clusters.", names(fc_control$model)))] %>%
    t() %>% as.data.frame()
  #finangle to get group id values
  cluster_groups <- cbind(rownames(cluster_groups), data.frame(cluster_groups, row.names = NULL))
  colnames(cluster_groups) <- c("cl_string", "cluster_id")
  cluster_groups <- cluster_groups %>%
    dplyr::mutate(cl_string = as.character(cl_string)) %>%
    tidyr::separate(cl_string, c("cl_text", "group_id"), sep = "[.]") %>%
    dplyr::select(-cl_text) %>%
    dplyr::mutate(group_id = as.numeric(group_id),
                  #must be factor for regression later
                  cluster_id = as.factor(cluster_id))

  epi_format <- epi_data_extd %>%
    dplyr::mutate(numericdate = as.numeric(Date)) %>%
    #get group_id
    dplyr::left_join(fc_control$groupcodes,
                     #NSE
                     by = set_names(rlang::quo_name(quo_groupfield),
                                    rlang::quo_name(quo_groupfield))) %>%
    #get cluster_id
    dplyr::left_join(cluster_groups, by = "group_id")

  epi_format
}

#' Convert environmental data into anomalies
#' @export
#'
anomalize_env <- function(env_fc, quo_groupfield, quo_obsfield, ncores) {

  # the following code only works on dataframes (e.g. env_fc[,curcol] instead of env_fc[[curcol]])
  # but I believe there's no use fixing it into tibble format since we're probably going to do it all with NSE,
  # and this conversion won't be necessary anyway
  env_fc <- as.data.frame(env_fc)

  # loop through environmental columns - sorry, can't figure out how to do this with NSE today
  regionfactor <- factor(env_fc[,1])
  doy          <- as.numeric(format(env_fc$Date, "%j"))
  for (curcol in 4:ncol(env_fc)) {

    cl <- parallel::makeCluster(ncores)

    tempbam <- mgcv::bam(env_fc[,curcol] ~ regionfactor + mgcv::s(doy, bs="cc", by=regionfactor),
                         data=env_fc,
                         chunk.size=1000,
                         cluster=cl)

    # could perhaps more cleverly be figured out by understanding the na.options of bam,
    # but for the moment just replace non-NA observations with their residuals
    env_fc[!is.na(env_fc[,curcol]),curcol] <- tempbam$residuals

    parallel::stopCluster(cl)

  }

  env_fc <- as.tibble(env_fc)
  return(env_fc)

}

#' Lag the env data
#' @export
#'
lag_environ_to_epi <- function(epi_fc, quo_groupfield, groupings,
                               env_fc, env_variables_used, laglen){

  #create lag frame
  datalagger <- tidyr::crossing(group_temp = groupings,
                                Date = unique(epi_fc$Date),
                                lag = seq(from = 0, to = laglen - 1, by = 1)) %>%
    # #same order from originally written expand.grid
    # arrange(lag, Date, group_temp) %>%
    #add lagging date
    dplyr::mutate(laggeddate = Date - as.difftime(lag, units = "days"))

  #and fix names with NSE
  datalagger <- datalagger %>%
    dplyr::rename(!!quo_name(quo_groupfield) := group_temp)

  #add env data
  datalagger <- dplyr::left_join(datalagger, env_fc,
                                 #because dplyr NSE, notice flip order
                                 by = set_names(c(rlang::quo_name(quo_groupfield), "Date"),
                                                c(rlang::quo_name(quo_groupfield), "laggeddate")))

  # pivot lagged environmental data to epi data
  epi_lagged <- epi_fc #to more easily debug and rerun
  for (curcol in which(colnames(env_fc) %in% env_variables_used)){
    valuevar <- colnames(env_fc)[curcol]
    #wide data for all lags of that env var
    meandat <- datalagger %>%
      dplyr::select(!!quo_groupfield, Date, lag, valuevar) %>%
      tidyr::spread(key = lag, value = valuevar)
    #rename lag columns (but not groupfield or Date)
    names(meandat)[-(1:2)] <- paste0(valuevar, "_", names(meandat)[-(1:2)])

    #join cur var wide data to epi data
    epi_lagged <- dplyr::left_join(epi_lagged, meandat,
                                   #dplyr NSE
                                   by = set_names(c(rlang::quo_name(quo_groupfield), "Date"),
                                                  c(rlang::quo_name(quo_groupfield), "Date")))
  } #end pivot loop

  # set up distributed lag basis functions
  lagframe <- data.frame(x = seq(from = 1, to = laglen, by = 1))
  alpha <- 1/4
  distlagfunc <- splines::ns(lagframe$x, intercept = TRUE,
                             knots = quantile(lagframe$x,
                                              probs=seq(from = alpha, to = 1 - alpha,
                                                        by = alpha),
                                              na.rm = TRUE))
  dlagdeg <- pracma::size(distlagfunc)[2]

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

  } #end distr lag summary loop

  #only keep bandsummaries (daily lags can be removed to free up a lot of space)
  #  note: ^ matches beginning of string, otherwise we'd get the bandsummaries too, which we want to keep
  for (cvar in env_variables_used){
    epi_lagged[, which(grepl(paste0("^", cvar, "_"), colnames(epi_lagged)))] <- NULL
  }

  epi_lagged
}

#' Run forecast regression
#' @export
#'
forecast_regression <- function(epi_lag, quo_groupfield, groupings,
                                env_variables_used, req_date, ncores){

  #forecasts are always done knowing up to just before that date
  # <<>> need to test what happens with daily vs 8day to daily data
  last_known_date <- req_date - lubridate::as.difftime(1, units = "days")

  #mark known or not
  epi_lag <- epi_lag %>%
    dplyr::mutate(known = ifelse(Date <= last_known_date, 1, 0))

  #number of clusters (different function for lm() needed if <= 1)
  n_clusters <- nlevels(epi_lag$cluster_id)

  #create variable bandsummaries equation piece
  #  e.g. 'bandsummaries_{var1} * cluster_id' for however many env var bandsummaries there are
  bandsums_list <- grep("bandsum_*", colnames(epi_lag), value = TRUE)
  bandsums_cl_list <- paste0(bandsums_list, "*cluster_id")
  #note, glue:: to distinguish b/t very different dplyr::collapse
  #need variant without known multiplication if <= 1 clusters
  if (n_clusters > 1) {
    bandsums_eq <- glue::collapse(bandsums_cl_list, sep =" + ")
  } else {
    bandsums_eq <- glue::collapse(bandsums_list, sep = " + ")
  }

  # create a doy field so that we can use a cyclical spline
  epi_lag <- dplyr::mutate(epi_lag, doy = as.numeric(format(Date, "%j")))

  # create modified bspline basis in epi_lag file to model longterm trends
  epi_lag <- cbind(epi_lag, truncpoly(x=epi_lag$Date,
                                      degree=6,
                                      maxobs=max(epi_lag$Date[epi_lag$known==1], na.rm=TRUE)))

  # get list of modbspline reserved variables and format for inclusion into model
  modb_list <- grep("modbs_reserved_*", colnames(epi_lag), value = TRUE)
  modb_list <- paste(modb_list, "*", rlang::quo_name(quo_groupfield))
  modb_eq <- glue::collapse(modb_list, sep = " + ")

  # ensure that quo_name(quo_groupfield) is a factor - gam/bam will fail if given a character,
  # which is unusual among regression functions, which typically just coerce into factors.
  epi_lag <- epi_lag %>% dplyr::mutate(!!rlang::quo_name(quo_groupfield) := factor(!!quo_groupfield))

  #filter to known
  epi_known <- epi_lag %>%
    dplyr::filter(known == 1)

  #due to dplyr NSE and bandsum eq and modb_eq pieces, easier to create expression to give to lm()
  reg_eq <- as.formula(paste("modeledvar ~ ", modb_eq,
                             "+s(doy, bs=\"cc\", by=",
                             rlang::quo_name(quo_groupfield),
                             ") + ",
                             rlang::quo_name(quo_groupfield), "+",
                             bandsums_eq))

  # set up clusters for parallel gam
  cl <- parallel::makeCluster(ncores)

  # run bam
  #message("Beginning bam on historical epi data")
  cluster_regress <- mgcv::bam(reg_eq, data = epi_known,
                               family=poisson(),
                               chunk.size=1000,
                               cluster=cl,
                               control=gam.control(trace=FALSE))

  # shut down cluster
  parallel::stopCluster(cl)

  #output prediction (through req_date now)
  # cluster_preds <- predict(cluster_regress,
  #                          newdata = epi_lag %>% filter(Date <= req_date),
  #                          se.fit = TRUE,
  #                          interval = "prediction",
  #                          level = 0)
  cluster_preds <- mgcv::predict(cluster_regress,
                                 newdata = epi_lag %>% dplyr::filter(Date <= req_date),
                                 se.fit = TRUE,                # included for backwards compatibility
                                 # interval = "prediction",    # cannot include for backwards compatibility, will probably break something
                                 type="response")

  #remove distributed lag summaries and bspline basis, which are no longer helpful
  band_names <- grep("bandsum_*", colnames(epi_lag), value = TRUE)
  bspl_names <- grep("modbs_reserved_*", colnames(epi_lag), value = TRUE)
  #remove
  epi_lag_trim <- dplyr::select(epi_lag, -dplyr::one_of(band_names))
  epi_lag_trim <- dplyr::select(epi_lag_trim, -dplyr::one_of(bspl_names))

  #now cbind to get ready to return
  epi_preds <- cbind(epi_lag_trim %>% filter(Date <= req_date),
                     as.data.frame(cluster_preds))

  #prediction of interest are last ones (equiv to req_date) per groupfield
  date_preds <- epi_preds %>%
    dplyr::group_by(!!quo_groupfield) %>%
    dplyr::filter(Date == req_date)
}

