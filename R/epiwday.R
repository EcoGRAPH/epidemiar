
#' Get epidemiological weekday number of a date-time
#'
#' Returns the epidemiological weekday number using either the WHO or CDC
#' system.
#'
#' The WHO system uses the ISO 8601 standard in which weeks start on Monday,
#' while in the CDC system weeks start on Sunday.
#'
#' @inheritParams lubridate::isoweek
#' @inheritParams epiweek
#'
#' @return The weekday number (1--7) as an integer vector.
#'
#' @inherit lubridate::isoweek references
#'
#' @examples
#' epiwday(as.Date("2005-01-01")) # 6
#' epiwday(as.Date("2005-01-01"), system = "cdc") # 7
#'
#' @export
#'
epiwday <- function(x, system = "who") {
  match.arg(system, c("who", "cdc"))
  if (system == "who") {
    as.integer(lubridate::wday(x - 1))
  } else if (system == "cdc") {
    as.integer(lubridate::wday(x))
  }
}
