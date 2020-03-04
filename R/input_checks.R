
#'Set defaults of any missing report_settings parameters
#'
#'Function sets defaults to report_settings parameters.
#'
#'@param raw_settings The report_settings object as given by the user.
#'@param env_variables List of all unique environmental variables in env_data.
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables.
#'@param groupings List of all unique geographical groupings in epi_data.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'
#'@inheritParams run_epidemia
#'
#'@return Returns a full report_settings object, using user supplied values or
#'  defaults is option was missing.
#'

set_report_defaults <- function(raw_settings,
                                env_info,
                                env_ref_data,
                                env_variables,
                                quo_obsfield,
                                groupings,
                                quo_groupfield){

  #set up list in case no report_settings were given
  if (is.null(raw_settings)){
    new_settings <- list()
  } else {
    #copy over to begin before editing/updating below
    new_settings <- raw_settings
  }

  if (is.null(raw_settings[["report_period"]])){
    new_settings[["report_period"]] <- 26
  }

  if (is.null(raw_settings[["report_inc_per"]])){
    new_settings[["report_inc_per"]] <- 1000
    #okay if not used, if report_value_type is cases instead of incidence
  }

  if (is.null(raw_settings[["epi_interpolate"]])){
    new_settings[["epi_interpolate"]] <- FALSE
  }

  if (is.null(raw_settings[["ed_summary_period"]])){
    new_settings[["ed_summary_period"]] <- 4
  }

  if (is.null(raw_settings[["model_run"]])){
    new_settings[["model_run"]] <- FALSE
  }

  if (is.null(raw_settings[["model_cached"]])){
    new_settings[["model_cached"]] <- NULL
  }

  if (is.null(raw_settings[["env_lag_length"]])){
    #maybe make default based on data length, but for now
    new_settings[["env_lag_length"]] <- 180
  }

  if (is.null(raw_settings[["fc_cyclicals"]])){
    new_settings[["fc_cyclicals"]] <- FALSE
  }

  if (is.null(raw_settings[["fc_future_period"]])){
    new_settings[["fc_future_period"]] <- 8
  }

  #default false, with explicit false for naive models (probably ok w/out, just being careful)
  if (is.null(raw_settings[["env_anomalies"]])){
    new_settings[["env_anomalies"]] <- dplyr::case_when(
      fc_model_family == "naive-persistence" ~ FALSE,
      fc_model_family == "naive-weekaverage" ~ FALSE,
      #default to FALSE
      TRUE ~ FALSE)
  }


  # For things that are being string matched:
  # tolower to capture upper and lower case user-input variations since match.arg is case sensitive
  # but must only try function if ed_method is not null (i.e. was given)

  #report_value_type
  # if provided, prepare for matching
  if (!is.null(raw_settings[["report_value_type"]])){
    new_settings[["report_value_type"]] <- tolower(raw_settings[["report_value_type"]])
  } else {
    #if not provided/missing/null
    message("Note: 'report_value_type' was not provided, returning results in case counts ('cases').")
    new_settings[["report_value_type"]] <- "cases"
  }
  #try match
  new_settings[["report_value_type"]] <- tryCatch({
    match.arg(new_settings[["report_value_type"]], c("cases", "incidence"))
  }, error = function(e){
    message("Warning: Given 'report_value_type' does not match 'cases' or 'incidence', running as 'cases'.")
    "cases"
  }, finally = {
    #failsafe default
    "cases"
  })

  # epi_date_type
  # if provided, prepare for matching
  if (!is.null(raw_settings[["epi_date_type"]])){
    #want to keep ISO and CDC capitalized, but drop 'Week' to 'week' if had been entered that way
    first_char <- substr(raw_settings[["epi_date_type"]], 1, 1) %>%
      tolower()
    #remainder of user entry
    rest_char <- substr(raw_settings[["epi_date_type"]], 2, nchar(raw_settings[["epi_date_type"]]))
    #paste back together
    new_settings[["epi_date_type"]] <- paste0(first_char, rest_char)
  } else {
    #if not provided/missing/null
    message("Note: 'epi_date_type' was not provided, running as weekly, ISO/WHO standard ('weekISO').")
    new_settings[["epi_date_type"]] <- "weekISO"
  }
  #try match
  new_settings[["epi_date_type"]] <- tryCatch({
    match.arg(new_settings[["epi_date_type"]], c("weekISO", "weekCDC")) #"monthly" reserved for future
  }, error = function(e){
    message("Warning: Given 'epi_date_type' does not match 'weekISO' or 'weekCDC', running as 'weekISO' (weekly, ISO/WHO standard).")
    "weekISO"
  }, finally = {
    #failsafe default
    "weekISO"
  })


  # ed_method
  # if provided, prepare for matching
  if (!is.null(raw_settings[["ed_method"]])){
    new_settings[["ed_method"]] <- tolower(raw_settings[["ed_method"]])
  } else {
    #if not provided/missing/null
    message("Note: 'ed_method' was not provided, running as 'none'.")
    new_settings[["ed_method"]] <- "none"
  }
  #try match
  new_settings[["ed_method"]] <- tryCatch({
    match.arg(new_settings[["ed_method"]], c("none", "farrington"))
  }, error = function(e){
    message("Warning: Given 'ed_method' does not match 'none' or 'farrington', running as 'none'.")
    "none"
  }, finally = {
    #failsafe default to no event detection
    "none"
  })


  # For more complicated defaults

  #env_var -- what is listed in env_data, env_ref_data, & env_info
  if (is.null(raw_settings[["env_var"]])){

    #create list of all environmental variables in env_info
    env_info_variables <- dplyr::pull(env_info, !!quo_obsfield)

    #create list of all environmental variables in env_ref_data
    env_ref_variables <- dplyr::pull(env_ref_data, !!quo_obsfield)

    #env_variables already gen list of env_data

    #Two sets of intersection to create list that are present in all three
    env_data_info <- dplyr::intersect(env_variables, env_info_variables)
    default_env_var <- dplyr::intersect(env_data_info, env_ref_variables)
    new_settings[["env_var"]] <- dplyr::tibble(obs_temp = default_env_var) %>%
                                 #rename NSE fun
                                 dplyr::rename(!!rlang::quo_name(quo_obsfield) := .data$obs_temp)

    #message result
    message("No user supplied list of environmetal variables to use. Using: ", paste(default_env_var, ""),
            " based on presence in env_data, env_ref_data, and env_info.\n")
  }

  #nthreads
  #default value is 1 for 1 core machines, 2 for multi-core (testing shows no additional value past 2)
  #if user-supplied, use that cap at 2, otherwise create a default number
  #used to decide if run anomalize_env() prior to forecasting
  if (!is.null(raw_settings[["fc_nthreads"]])) {
    # nthreads above 2 is not actually helpful
    new_settings[["fc_nthreads"]] <- ifelse(raw_settings[["fc_nthreads"]] > 1, 2, 1)
  } else {
    #no value fed in, so test and determine
    new_settings[["fc_nthreads"]] <- ifelse(parallel::detectCores(logical=FALSE) > 1, 2, 1)
  } #end else for ncores not given


  #fc_clusters
  #default is one cluster, probably not what you actually want for any type of large system
  if (is.null(raw_settings[["fc_clusters"]])){
    #create tbl of only one cluster
    #groupings already exist as list of geographic groups
    cluster_tbl <- tibble::tibble(group_temp = groupings, cluster_id = 1) %>%
      #and fix names with NSE
      dplyr::rename(!!rlang::quo_name(quo_groupfield) := .data$group_temp)
    #assign
    new_settings[["fc_clusters"]] <- cluster_tbl
  }


  # Developer options
  if (is.null(raw_settings[["dev_fc_fit_freq"]])){
    new_settings[["dev_fc_fit_freq"]] <- "once"
  }
  if (is.null(raw_settings[["dev_fc_modbsplines"]])){
    new_settings[["dev_fc_modbsplines"]] <- FALSE
  }
  if (is.null(raw_settings[["dev_fc_formula"]])){
    new_settings[["dev_fc_formula"]] <- NULL
  }


  new_settings

}



#'Functions to check input to epidemiar
#'
#'Function does basic existance checks and variety of logic checks on input data
#'to run_epidemia().
#'
#'@param quo_casefield Quosure of user given field containing the disease case
#'  counts.
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param quo_obsfield Quosure of user given field name of the environmental data
#'  variables
#'@param quo_valuefield Quosure of user given field name of the value of the
#'  environmental data variable observations.
#'
#'@inheritParams run_epidemia
#'
#'@return Returns a flag if there were any errors, plus accompanying error
#'  messages. Also returns a flag and messages for warnings, as well.
#'
#'
#'

input_check <- function(epi_data,
                        env_data,
                        env_ref_data,
                        env_info,
                        quo_casefield,
                        quo_popfield,
                        quo_groupfield,
                        quo_obsfield,
                        quo_valuefield,
                        fc_model_family,
                        report_settings){

  # Want ALL data checks to happen, whether or not error happen before the end of the tests.
  # Want to collect all errors, and return all of them to console
  # Exception: if some early/critical errors happen, it will need to be returned early or skip certain later checks that would fail to check properly.
  # Note: does not test for integer value (versus simply numeric), which is non-trivial, for various items.

  # Create err_flag (binary if any error) and err_msgs (all error messages) variables
  # these will be passed to each sub check function and returned
  err_flag <- FALSE
  err_msgs <- ""

  # Also collect any warning messages to display
  warn_flag <- FALSE
  warn_msgs <- "Warning messages:\n"


  # Existing & Types --------------------------------------------------------

  # Quick test for some simple settings
  if (!is.numeric(report_settings[["report_inc_per"]]) || report_settings[["report_inc_per"]] <= 0){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'inc_per' must be numeric and a positive number.\n")
  }


  # epi_data tests
  # has obs_date as Date
  if (!"obs_date" %in% colnames(epi_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column 'obs_date' in the epidemiological dataset, 'epi_data'.\n")
  } else if (!class(epi_data$obs_date) == "Date"){
    #has obs_date, now check type
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'obs_date' in the epidemiological dataset, 'epi_data', must be type Date.\n")
  }
  # has casefield
  if (!rlang::quo_name(quo_casefield) %in% colnames(epi_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_casefield), ", in the epidemiological dataset, 'epi_data'.\n")
  }
  # has groupfield
  if(!rlang::quo_name(quo_groupfield) %in% colnames(epi_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_groupfield), ", in the epidemiological dataset, 'epi_data'.\n")
  }
  # has populationfield, but only if given as it is optional
  #testing if quosure was created on NULL object.
  if(!rlang::quo_is_null(quo_popfield)){
    if(!rlang::quo_name(quo_popfield) %in% colnames(epi_data)){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "The specified column ", rlang::quo_name(quo_popfield), ", for population must be in the in the epidemiological dataset, 'epi_data'.\n")
    }
  }


  # env_data tests
  # has obs_date as Date
  if(!"obs_date" %in% colnames(env_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column 'obs_date' in the environmental dataset, 'env_data'.\n")
  } else if(!class(env_data$obs_date) == "Date"){
    #has obs_date, now check type
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'obs_date' in the environmental dataset, 'env_data', must be type Date.\n")
  }
  # has groupfield
  if(!rlang::quo_name(quo_groupfield) %in% colnames(env_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_groupfield), ", in the environmental dataset, 'env_data'.\n")
  }
  # has obsfield
  if(!rlang::quo_name(quo_obsfield) %in% colnames(env_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_obsfield), ", in the environmental dataset, 'env_data'.\n")
  }
  # has valuefield
  if(!rlang::quo_name(quo_valuefield) %in% colnames(env_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_valuefield), ", in the environmental dataset, 'env_data'.\n")
  }


  # env_ref tests
  # has groupfield
  if(!rlang::quo_name(quo_groupfield) %in% colnames(env_ref_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_groupfield), ", in the environmental reference dataset, 'env_ref_data'.\n")
  }
  # has obsfield
  if(!rlang::quo_name(quo_obsfield) %in% colnames(env_ref_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_obsfield), ", in the environmental reference dataset, 'env_ref_data'.\n")
  }
  #has week_epidemiar
  if(!"week_epidemiar" %in% colnames(env_ref_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column 'week_epidemiar' for the week of the year in the environmental reference dataset, 'env_ref_data'.\n")
  } else if(!(is.numeric(env_ref_data$week_epidemiar) | is.integer(env_ref_data$week_epidemiar))){
    #week_epidemiar exists, now check class/type
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "The column 'week_epidemiar' in 'env_ref_data' must be numeric or integer type (integer values only).\n")
  }
  #has ref_value
  if(!"ref_value" %in% colnames(env_ref_data)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column 'ref_value' for the historical reference value in the dataset, 'env_ref_data'.\n")
  }

  # env_info
  # has obsfield
  if(!rlang::quo_name(quo_obsfield) %in% colnames(env_info)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_obsfield), ", in the environmental metadata file, 'env_info'.\n")
  }
  # has reference_method
  if(!"reference_method" %in% colnames(env_info)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column 'reference_method' in 'env_info' for how to summarize values from daily to weekly ('sum' or 'mean').\n")
  }

  # Lengths of Report Sections --------------------------------------

  #special flag for all dates
  rpt_len_flag <- FALSE

  # check data types
  if (!(is.numeric(report_settings[["report_period"]]) | is.integer(report_settings[["report_period"]]))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'report_settings$report_period' must be numeric or integer type - integer number of weeks only.\n")
    rpt_len_flag <- TRUE
  }
  if (!(is.numeric(report_settings[["fc_future_period"]]) | is.integer(report_settings[["fc_future_period"]]))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'report_settings$forecast_future' must be numeric or integer type - integer number of weeks only.\n")
    rpt_len_flag <- TRUE
  } else if (report_settings[["fc_future_period"]] > 13){
    # warn on long forecasts
    warn_flag <- TRUE
    warn_msgs <- paste(warn_msgs, "Warning: It is not recommended to forecast more than 12 weeks into the future. You are forecasting for ", report_settings[["fc_future_period"]], " weeks.\n")
  }
  if (!(is.numeric(report_settings[["ed_summary_period"]]) | is.integer(report_settings[["ed_summary_period"]]))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'report_settings$ed_summary_period' must be numeric or integer type - integer number of weeks only.\n")
    rpt_len_flag <- TRUE
  }

  # report length must be equal to or larger than forecast and ED period together
  if (!rpt_len_flag){
    if (report_settings[["report_period"]] < report_settings[["ed_summary_period"]] + report_settings[["fc_future_period"]]){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "The report length ", report_settings[["report_period"]], " must be longer than the early detection period ", report_settings[["ed_summary_period"]], " plus the forecast ", report_settings[["fc_future_period"]], ".\n")
    }
  }


  # Models & Caching --------------------------------------------------------

  #check if fc_model_family is cached that a cached model was given, else fail with error
  if (fc_model_family == "cached"){
    if (is.null(report_settings[["model_cached"]])){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "If 'fc_model_family' == 'cached', a cached model must be supplied in 'report_settings$model_cached'.\n")
    }
  }


  #if given a full model
  if (!is.null(report_settings[["model_cached"]])){

    #check that $model_info and $model_obj exists in model_cached
    if (all(c("model_obj", "model_info") %in% names(report_settings[["model_cached"]]))){

      #if model looks okay so far, then check further

      #make sure given model (if given) is a regression object (using basic "lm" as test)
      #model_cached$model_obj
      classes <- class(report_settings[["model_cached"]][["model_obj"]])
      if (!"lm" %in% classes){
        err_flag <- TRUE
        err_msgs <- paste(err_msgs, "The object in 'report_settings$model_cached$model_obj' is not a regression object, found classes are: ", classes, ".\n")
      } #end lm check

      #if using a cached model, the model family from the cached model will be used
      #warn about overriding any user input family
      if (fc_model_family != "cached"){
        warn_flag <- TRUE
        warn_msgs <- paste(warn_msgs, "Warning: the cached model family ", report_settings$model_cached$model_info$fc_model_family, " will override any user input.",
                           "Found 'fc_model_family' set to ",  fc_model_family, "instead of 'cached'.\n")
      }

      #end if names
    } else {
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "The given cached model is missing $model_obj and/or $model_info.\n")
    } #end else on if names

  } #end if !is.null model_cached

  # things that must exist in model_cached$model_info
  # model_cached$model_info$fc_model_family
  # model_cached$model_info$date_created
  # model_cached$model_info$known_epi_range$max
  #but will probably give decent error messages on their own if missing.




  # Control lists -----------------------------------------------------------

  # Forecasting

  # model_env
  # has obsfield
  if(!rlang::quo_name(quo_obsfield) %in% colnames(report_settings[["env_var"]])){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_obsfield), ", to indicate the list of model environmental variables in 'report_settings$env_vars'.\n")
  } else {
    #does have obsfield,
    #check that model variables exist in env data and env ref data

    #but only if no other problems so far, since that could cause errors in the checks below
    if (!err_flag){

      #pull variables from model info input
      model_vars <- report_settings[["env_var"]] %>% dplyr::pull(!!quo_obsfield)
      #pull variables in env data
      env_in_data <- env_data %>% dplyr::pull(!!quo_obsfield) %>% unique()
      #pull variables in env ref data
      env_in_ref <- env_ref_data %>% dplyr::pull(!!quo_obsfield) %>% unique()

      if (!all(model_vars %in% env_in_data)){
        err_flag <- TRUE
        err_msgs <- paste(err_msgs, "Model variable(s) given in 'report_settings$env_var' is/are missing from 'env_data':\n",
                          model_vars[which(!model_vars %in% env_in_data)], "\n")
      }
      if (!all(model_vars %in% env_in_ref)){
        err_flag <- TRUE
        err_msgs <- paste(err_msgs, "Model variable(s) given in 'report_settings$env_var' is/are missing from 'env_ref_data':\n",
                          model_vars[which(!model_vars %in% env_in_ref)], "\n")
      }
    } #end err_flag
  } #end else obsfield

  #clusters
  # has groupfield
  if(!rlang::quo_name(quo_groupfield) %in% colnames(report_settings[["fc_clusters"]])){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_groupfield), ", in 'report_settings$clusters'.\n")
  }
  # has cluster_id
  if(!"cluster_id" %in% colnames(report_settings[["fc_clusters"]])){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "There must be a column 'cluster_id' in 'report_settings$clusters'.\n")
  }
  #now check that all geographic groupings from epi data have a cluster assigned
  #as long as no previous errors
  if (!err_flag){
    #groupings in cluster info
    model_cl <- report_settings[["fc_clusters"]] %>% dplyr::pull(!!quo_groupfield)
    #groupings in epidemiological data
    groups_epi <- dplyr::pull(epi_data, !!quo_groupfield) %>% unique()
    #check all in cluster list
    if (!all(groups_epi %in% model_cl)){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "Geographic groupings present in the epidemiological data are missing in 'report_settings$clusters':\n",
                        groups_epi[which(!groups_epi %in% model_cl)])
    }
    #Don't need to check environmental data. Extra env data for other groupings not in epidemiological data are just ignored.
  }

  #lag_length
  #already checked existance and numeric/integer type
  #check that enough environmental data exists for lag length selected

  #but only if no other problems so far, since that could cause errors in the checks below
  if (!err_flag){
    #subset to env variables as dictated by the model
    env_model_data <- pull_model_envvars(env_data, quo_obsfield, env_var = report_settings$env_var)
    #get earliest dates available
    env_start_dts <- env_model_data %>% dplyr::group_by(!!quo_obsfield) %>% dplyr::summarize(start_dt = min(.data$obs_date))
    #date needed by laglength and first epidemiological data date
    need_dt <- min(epi_data$obs_date) - as.difftime(report_settings[["env_lag_length"]], units = "days")
    #all env dates equal or before needed date?
    if (!all(env_start_dts$start_dt <= need_dt)){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "Not enough environmental data for a lag length of ", report_settings[["env_lag_length"]],
                        "days.\n Epidemiological start is", min(epi_data$obs_date),
                        "therefore environmental data is needed starting", need_dt, "for variables:\n",
                        env_start_dts[which(!env_start_dts$start_dt <= need_dt),1])

    }
  } #end err_flag


  # ed_method & ed_control

  if (report_settings[["ed_method"]] == "farrington"){

    #controls for Farrington all have defaults in farringtonFlexible() and can be missing, just warn
    if (is.null(report_settings[["ed_control"]])){
      #warning if missing though
      warn_flag <- TRUE
      warn_msgs <- paste(warn_msgs, "Warning: Early Detection controls not found, running with surveillance package defaults.\n")
    }
  }



  ## Return
  create_named_list(err_flag, err_msgs, warn_flag, warn_msgs)

} #end input_check()







