adjust_hsv <- function(col, h = NULL, s = NULL, v = NULL, mode = "relative") {
  if (!is.null(h)) assertthat::assert_that(is.numeric(h), h >= -1, h <= 1)
  if (!is.null(s)) assertthat::assert_that(is.numeric(s), s >= -1, s <= 1)
  if (!is.null(v)) assertthat::assert_that(is.numeric(v), v >= -1, v <= 1)
  assertthat::assert_that(tolower(mode) %in% c("relative", "absolute"))
  
  x <- col2hsv(col)
  
  if (tolower(mode) == "relative") {
    if (!is.null(h)) x["h",] <- clamp(x["h",] + h, 0, 1)
    #if (!is.null(h)) x["h",] <- (x["h",] + 0.5) %% 1
    if (!is.null(s)) x["s",] <- clamp(x["s",] + s, 0, 1)
    if (!is.null(v)) x["v",] <- clamp(x["v",] + v, 0, 1)
  }
  else if (tolower(mode) == "absolute") {
    if (!is.null(h)) x["h",] <- clamp(h, 0, 1)
    #if (!is.null(h)) x["h",] <- (x["h",] + 0.5) %% 1
    if (!is.null(s)) x["s",] <- clamp(s, 0, 1)
    if (!is.null(v)) x["v",] <- clamp(v, 0, 1)
  }
  
  hsv2hex(x)
}
