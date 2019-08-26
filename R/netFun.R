#' netFun
#'
#' This is the network reconstruction function that will be used to build aggregate networks.
#' @param x Numeric matrix with samples in columns.
#' @keywords netFun
#' @export
#' @examples
#' netFun()


netFun <- function(x) {
    stats::cor(t(x), method="pearson") 
}
