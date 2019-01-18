
#' Create dates from epidemiological weeks and years
#'
#' Calculate Gregorian calendar dates based on epidemiological years, week
#' numbers, and weekday numbers.
#'
#' Arguments `year`, `week`, and `weekday` are recycled as necessary.
#'
#' @param year epidemiological year
#' @param week eidemiological week number (1--53).
#' @param weekday epidemiological weekday number (1--7). Day 1 is a Monday in
#'   the WHO system  and a Sunday in the CDC system.
#' @inheritParams epiweek
#'
#' @return A vector of class `Date`.
#' @export
#'
#' @examples
#' make_date_yw(2017, 1)
#' make_date_yw(2017, 1, 2)
#' make_date_yw(2017, 1, system = "cdc")
#' make_date_yw(2017, 1, 2, system = "cdc")
#'
#' # arguments are recycled
#' make_date_yw(2017, 1:10)
#' make_date_yw(2017, 1, 1:7)
#' make_date_yw(2010:2017, 1)
#'
make_date_yw <- function(year = 1970L, week = 1L, weekday = 1L, system = "who") {

  match.arg(system, c("who", "cdc"))

  lengths <- vapply(list(year, week, weekday), length, 1, USE.NAMES = FALSE)
  if (min(lengths) == 0L) as.Date(integer(), lubridate::origin)

  # recycle arguments
  N <- max(lengths)
  y <- rep_len(as.integer(year), N)
  w <- rep_len(as.integer(week), N)
  d <- rep_len(as.integer(weekday), N)

  out <-
    ifelse(
      is.na(y) | is.na(w) | is.na(d), NA,
      {
        jan1 <- lubridate::make_date(y, 1, 1)
        wday <- epiwday(jan1, system)
        to_add <- ifelse(wday <= 4, 1, 8) - wday
        wk1 <- jan1 + to_add
        day1 <- wk1 + (w - 1) * 7
        day1 + d - 1
      }
    )
  as.Date(out, origin)

}

#' @export
#' @rdname make_date_yw
#' @usage NULL
parse_epiweek <- function(...) {
  .Deprecated("make_date_yw")
  make_date_yw(...)
}
