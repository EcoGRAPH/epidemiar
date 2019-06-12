#'Functions to check input to epidemiar
#'
#'Function does basic existance checks and variety of logic checks on input data
#'to run_epidemia().
#'
#'@param epi_data Epidemiological data with case numbers per week, with date
#'  field "obs_date".
#'@param quo_casefield Quosure of user given field containing the disease case
#'  counts.
#'@param quo_popfield Quosure of user-given field containing population values.
#'@param inc_per Number for what unit of population the incidence should be
#'  reported in, e.g. incidence rate of 3 per 1000 people.
#'@param quo_groupfield Quosure of the user given geographic grouping field to
#'  run_epidemia().
#'@param week_type String indicating the standard (WHO ISO-8601 or CDC epi
#'  weeks) that the weeks of the year in epidemiological and environmental
#'  reference data use ["ISO" or "CDC"]. (Required: epidemiological observation
#'  dates listed are LAST day of week).
#'@param report_period The number of weeks that the entire report will cover.
#'  The \code{report_period} minus \code{forecast_future} is the number of weeks
#'  of past (known) data that will be included.
#'@param ed_summary_period The number of weeks that will be considered the
#'  "early detection period". It will count back from the week of last known
#'  epidemiological data.
#'@param ed_method Which method for early detection should be used ("Farrington"
#'  is only current option, or "None").
#'@param ed_control All parameters for early detection algorithm, passed through
#'  to that subroutine.
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
#'@param forecast_future Number of futre weeks from the end of the
#'  \code{epi_data} to produce forecasts.
#'@param fc_control Parameters for forecasting, including which environmental
#'  variable to include and any geographic clusters.
#'@param env_ref_data Historical averages by week of year for environmental
#'  variables. Used in extended environmental data into the future for long
#'  forecast time, to calculate anomalies in early detection period, and to
#'  display on timeseries in reports.
#'@param env_info Lookup table for environmental data - reference creation
#'  method (e.g. sum or mean), report labels, etc.
#'@param model_obj Deprecated, use model_cached.
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
#'@return Returns a flag if there were any errors, plus accompanying error
#'  messages. Also returns a flag and messages for warnings, as well.
#'
#'
#'

input_check <- function(epi_data,
                        quo_casefield,
                        quo_popfield,
                        inc_per,
                        quo_groupfield,
                        week_type,
                        report_period,
                        ed_summary_period,
                        ed_method,
                        ed_control = NULL,
                        env_data,
                        quo_obsfield,
                        quo_valuefield,
                        forecast_future,
                        fc_control,
                        env_ref_data,
                        env_info,
                        model_obj = NULL,
                        model_cached = NULL,
                        model_choice = NULL){

  #NULL defaults same as run_epidemia(), but excluding the necessary items already checked

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
  # inc_per (0 could be "cases" so allow when that is built)
  if (!is.numeric(inc_per) || inc_per <= 0){
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
  if (!(is.numeric(report_period) | is.integer(report_period))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'report_period' must be numeric or integer type - integer number of weeks only.\n")
    rpt_len_flag <- TRUE
  }
  if (!(is.numeric(forecast_future) | is.integer(forecast_future))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'forecast_future' must be numeric or integer type - integer number of weeks only.\n")
    rpt_len_flag <- TRUE
  } else if (forecast_future > 12){
    # warn on long forecasts
    warn_flag <- TRUE
    warn_msgs <- paste(warn_msgs, "Warning: It is not recommended to forecast more than 12 weeks into the future. You are forecasting for ", forecast_future, " weeks.\n")
  }
  if (!(is.numeric(ed_summary_period) | is.integer(ed_summary_period))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'ed_summary_period' must be numeric or integer type - integer number of weeks only.\n")
    rpt_len_flag <- TRUE
  }

  # report length must be equal to or larger than forecast and ED period together
  if (!rpt_len_flag){
    if (report_period < ed_summary_period + forecast_future){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "The report length ", report_period, " must be longer than the early detection period ", ed_summary_period, " plus the forecast ", forecast_future, ".\n")
    }
  }


# Models & Caching --------------------------------------------------------

  #use model_cached not old model_obj
  if (!is.null(model_obj) & is.null(model_cached)){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "Please use the new 'model_cached' argument, and not deprecated 'model_obj'.\n")
  }

  #make sure model_choice matches between cached model and settings.
  if (!model_cached$model_info$model_choice == model_choice){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "The model choice of the given cached model, ", model_cached$model_info$model_choice, " does not match the current setting of 'model_choice', ", model_choice, ".\n")
  }

  #make sure given model (if given) is a regression object (using basic "lm" as test)
  #model_cached$model_obj
  if (!is.na(model_cached)){
    classes <- class(model_cached$model_obj)
    if(!"lm" %in% classes){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "The object in 'model_cached$model_obj' is not a regression object, found classes are: ", classes, ".\n")
    }
  }

  # things that must exist in model_cached$model_info
  # model_cached$model_info$model_choice
  # model_cached$model_info$date_created
  # model_cached$model_info$known_epi_range$max
  #but will probably give decent error messages on their own if missing.

  # Control lists -----------------------------------------------------------

  # Forecasting
  #special flag for initial forecasting checks
  fc_flag <- FALSE

  #fc_control
  # env_vars, clusters, lag_length
  if (is.null(fc_control[["env_vars"]])){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "List of environmental variables to use for modeling is missing. Check 'fc_control$env_vars'.\n")
    fc_flag <- TRUE
  }
  if (is.null(fc_control[["clusters"]])){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "Cluster information to use for modeling is missing. Check 'fc_control$clusters'.\n")
    fc_flag <- TRUE
  }
  if (is.null(fc_control[["lag_length"]])){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "Length of maximum lag in days, 'lag_length', is missing. Check 'fc_control$lag_length'.\n")
    fc_flag <- TRUE
  } else if (!(is.numeric(fc_control[["lag_length"]]) | is.integer(fc_control[["lag_length"]]))){
    err_flag <- TRUE
    err_msgs <- paste(err_msgs, "'lag_length' must be an integer number of days. Check 'fc_control$lag_length'.\n")
    fc_flag <- TRUE
  }


  #if no initial fc errors, then continue
  if (!fc_flag){

    # model_env
    # has obsfield
    if(!rlang::quo_name(quo_obsfield) %in% colnames(fc_control$env_vars)){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_obsfield), ", to indicate the list of model environmental variables in 'fc_control$env_vars'.\n")
    } else {
      #does have obsfield,
      #check that model variables exist in env data and env ref data

      #but only if no other problems so far, since that could cause errors in the checks below
      if (!err_flag){

        #pull variables from model info input
        model_vars <- fc_control$env_vars %>% dplyr::pull(!!quo_obsfield)
        #pull variables in env data
        env_in_data <- env_data %>% dplyr::pull(!!quo_obsfield) %>% unique()
        #pull variables in env ref data
        env_in_ref <- env_ref_data %>% dplyr::pull(!!quo_obsfield) %>% unique()

        if (!all(model_vars %in% env_in_data)){
          err_flag <- TRUE
          err_msgs <- paste(err_msgs, "Model variable(s) given in 'fc_control$env_vars' is/are missing from 'env_data':\n",
                            model_vars[which(!model_vars %in% env_in_data)], "\n")
        }
        if (!all(model_vars %in% env_in_ref)){
          err_flag <- TRUE
          err_msgs <- paste(err_msgs, "Model variable(s) given in 'fc_control$env_vars' is/are missing from 'env_ref_data':\n",
                            model_vars[which(!model_vars %in% env_in_ref)], "\n")
        }
      } #end err_flag
    } #end else obsfield

    #clusters
    # has groupfield
    if(!rlang::quo_name(quo_groupfield) %in% colnames(fc_control$clusters)){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "There must be a column ", rlang::quo_name(quo_groupfield), ", in 'fc_control$clusters'.\n")
    }
    # has cluster_id
    if(!"cluster_id" %in% colnames(fc_control$clusters)){
      err_flag <- TRUE
      err_msgs <- paste(err_msgs, "There must be a column 'cluster_id' in 'fc_control$clusters'.\n")
    }
    #now check that all geographic groupings from epi data have a cluster assigned
    #as long as no previous errors
    if (!err_flag){
      #groupings in cluster info
      model_cl <- fc_control$clusters %>% dplyr::pull(!!quo_groupfield)
      #groupings in epidemiological data
      groups_epi <- dplyr::pull(epi_data, !!quo_groupfield) %>% unique()
      #check all in cluster list
      if (!all(groups_epi %in% model_cl)){
        err_flag <- TRUE
        err_msgs <- paste(err_msgs, "Geographic groupings present in the epidemiological data are missing in 'fc_control$clusters':\n",
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
      env_model_data <- pull_model_envvars(env_data, quo_obsfield, fc_control)
      #get earliest dates available
      env_start_dts <- env_model_data %>% dplyr::group_by(!!quo_obsfield) %>% summarize(start_dt = min(obs_date))
      #date needed by laglength and first epidemiological data date
      need_dt <- min(epi_data$obs_date) - as.difftime(fc_control$lag_length, units = "days")
      #all env dates equal or before needed date?
      if (!all(env_start_dts$start_dt <= need_dt)){
        err_flag <- TRUE
        err_msgs <- paste(err_msgs, "Not enough environmental data for a lag length of ", fc_control$lag_length,
                          "days.\n Epidemiological start is", min(epi_data$obs_date),
                          "therefore environmental data is needed starting", need_dt, "for variables:\n",
                          env_start_dts[which(!env_start_dts$start_dt <= need_dt),1])

      }
    } #end err_flag

  } #end fc flag

  # ed_method & ed_control

  if (ed_method == "Farrington"){

    # if Farrington, then check for controls for Farrington
    #w = 4, reweight = TRUE, weightsThreshold = 2.58, trend = TRUE, pThresholdTrend = 0, populationOffset = TRUE, noPeriods = 10, pastWeeksNotIncluded = 4, thresholdMethod = "nbPlugin"
    #allow b

    #actually, all have defaults in farringtonFlexible() and can be missing
    if (is.null(ed_control)){
      #warning if missing though
      warn_flag <- TRUE
      warn_msgs <- paste(warn_msgs, "Warning: Early Detection controls not found, running with surveillance package defaults.\n")
    }
  }



  ## Return
  create_named_list(err_flag, err_msgs, warn_flag, warn_msgs)

} #end input_check()







