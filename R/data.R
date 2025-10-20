#' @title Realized variance of Amazon and SPY
#'
#' @description Intraday realized covariances of the returns between the Amazon stock (\code{rvarAMZN}) and the SPDR S&P 500 ETF (\code{rvarSPY}) using five minutes data, for the period of September 13th, 2023 to September 12, 2024.
#' @format A 2 by 2 by 250 array
#' @source Anne MacKay
#' @examples
#' data(realvar, package = "ksm")
#' bopt <- bandwidth_optim(
#'  x = realvar,
#'  criterion = "lscv",
#'  kernel = "Wishart",
#'  h = 4L
#'  )
"realvar"
