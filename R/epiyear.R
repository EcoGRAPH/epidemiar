
#' Get epidemiological year of a date-time
#'
#' Returns the epidemiological year using either the WHO or CDC system, both of
#' which use a reoccuring leap year resulting in week numbers from 1--53.
#'
#' The WHO system uses the ISO 8601 standard in which weeks start on Monday,
#' while in the CDC system weeks start on Sunday.
#'
#' In both systems, week 1 of a given year is identified as being the first week
#' of the year with at least four days in that calendar year. As a result, the
#' epidemiological year may differ from the Gregorian calendar year for some
#' dates near January 1.
#'
#' Internally, this function uses [lubridate::isoyear()] for WHO years and a
#' modified version of that functions code for CDC years
#'
#' @inheritParams lubridate::isoweek
#' @inheritParams epiweek
#'
#' @inherit lubridate::isoweek references
#' @inherit lubridate::isoyear return
#'
#' @examples
#' epiyear(as.Date("2005-01-01")) # 2004
#' epiyear(as.Date("2005-01-01"), system = "cdc") # 2005
#'
#' @export
#'

epiyear <- function(x, system = "who") {
  match.arg(system, c("who", "cdc"))
  if (system == "who") {
    as.integer(lubridate::isoyear(x))
  } else if (system == "cdc") {
    # code modified from lubridate::isoweek()
    # replaced dn calculation with wday(x) to treat Sundays as day 1
    xday <- lubridate::make_datetime(lubridate::year(x), lubridate::month(x),
                                     lubridate::day(x))
    dn <- lubridate::wday(x)
    nth <- xday + lubridate::ddays(4 - dn)
    as.integer(lubridate::year(nth))
  }
}