
#' Get epidemiological week number of a date-time
#'
#' Returns the epidemiological week number using either the WHO or CDC system,
#' both of which use a reoccuring leap year resulting in week numbers from
#' 1--53.
#'
#' The WHO system uses the ISO 8601 standard in which weeks start on Monday,
#' while in the CDC system weeks start on Sunday.
#'
#' In both systems, week 1 of a given year is identified as being the first week
#' of the year with at least four days in that calendar year. As a result, the
#' epidemiological year may differ from the Gregorian calendar year for some
#' dates near January 1.
#'
#' Internally, this function uses [lubridate::isoweek()] for WHO weeks and a
#' modified version of that functions code for CDC weeks.
#'
#' @inheritParams lubridate::isoweek
#' @param system Either "who" or "cdc". WHO epidemiological weeks start on
#'   Monday. CDC epidemiological weeks (MMWR weeks) start on Sunday. The default
#'   is "who".
#'
#' @inherit lubridate::isoweek return references
#'
#' @examples
#' epiweek(as.Date("2005-01-01")) # 52
#' epiweek(as.Date("2005-01-01"), system = "cdc") # 1
#'
#' @export
#'
epiweek <- function(x, system = "who") {
  match.arg(system, c("who", "cdc"))
  if (system == "who") {
    lubridate::isoweek(x)
  } else if (system == "cdc") {
    # code modified from lubridate::isoweek()
    # replaced dn calculation with wday(x) to treat Sundays as day 1
    xday <- lubridate::make_datetime(lubridate::year(x), lubridate::month(x),
                                     lubridate::day(x))
    dn <- lubridate::wday(x)
    nth <- xday + lubridate::ddays(4 - dn)
    jan1 <- lubridate::make_datetime(lubridate::year(nth), 1, 1)
    1L + as.integer(difftime(nth, jan1, units = "days"))%/%7L
  }
}
