# Functions to check input to epidemiar

input_check_exists <- function(){

  # First, before any other input checks or processing, we need to make sure all the necessary fields, data, and arguments were given
  # This will not check if they've assigned the right thing to the argument, or got the argument order correct if not explicit argument declarations
  # But, no other checks can really proceed if things are missing

  nec_flds <- c(casefield, groupfield, obsfield, valuefield) #populationfield

  nec_data <- c(epi_data, env_data, env_ref_data, env_info)

  nec_cntls <- c(week_type, ed_control, fc_control) #rest has defaults
  # Note: only checking if control list exists, nothing about what is in the list (later check)

  #combine all
  necessary <- c(nec_flds, nec_data, nec_cntls)

  #initialize missing info msgs & flag
  missing_msgs <- "Missing critical arguments. Please make sure the following are included:"
  missing_flag <- FALSE

  #loop through all necessary fields, checking if argument exists, collecting list of missing
  for (arg_name in necessary){
    if (!hasArg(arg_name)){
      missing_flag <- TRUE
      missing_list <- paste(missing_msgs, arg_name, sep = "\n")
    }
  }


}


input_check <- function(epi_data,
                        quo_casefield,
                        quo_populationfield,
                        inc_per,
                        quo_groupfield,
                        week_type,
                        report_period,
                        ed_summary_period,
                        ed_method,
                        ed_control,
                        env_data,
                        quo_obsfield,
                        quo_valuefield,
                        forecast_future,
                        fc_control,
                        env_ref_data,
                        env_info){

  # Want ALL data checks to happen, whether or not error happen before the end of the tests.
  # Want to collect all errors, and return all of them to console
  # Exception: if some early/critical errors happen, it will need to be returned early, otherwise later checks won't work at all

  # Create err_flag (binary if any error) and err_msgs (all error messages) variables
  # these will be passed to each sub check function and returned
  err_flag <- FALSE
  err_msgs <- ""

  # Check if all important field names where given
  # Note: these are the ENQUO fields, so will be "" if argument not given
  check_givenfields()



}

# Existing & Types --------------------------------------------------------

# fields given
# casefield, populationfield, groupfield, obsfield, valuefield


# simple settings
# inc_per (0 could be "cases" so allow)
# week_type


# epi_data
# has obs_date as Date
# has casefield, populationfield, groupfield


# clusters
# has groupfield
# rem to allow for 1


# env_data
# has obs_date as Date
# has groupfield, obsfield, valuefield


# env_ref


# env_info
# has obsfield
# has reference_method
#? has report_label


# model_env
# has obsfield



# Control lIsts -----------------------------------------------------------

## ED
# if Farrington, then controls for Farrington
#w = 4, reweight = TRUE, weightsThreshold = 2.58,
#trend = TRUE, pThresholdTrend = 0,
#populationOffset = TRUE,
#noPeriods = 10, pastWeeksNotIncluded = 4,
#thresholdMethod = "nbPlugin
#allow b
#check which ones have default is missing, allowed to be missing


# forecasting controls list
#env_vars = model_env
#clusters
#lag_length
#fit_freq
#ncores


# Dates & Lengths of Report Sections --------------------------------------

#report_period
#ed_summary_period
#forecast_future

# Dependencies & Calculations ---------------------------------------------


# all of geog groupings in epi_data in clusters

# all of env var in model_env in env_info

# all of geog groupings * env var in model_env existing in env_ref_data


## Forecasting
# laglength of env data avail
