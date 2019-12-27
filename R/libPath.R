#' libPath
#' This function provides the path of the compiled library
#' @return path of the compiled library
#'
#' @importFrom utils setTxtProgressBar
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_polygon
#' @importFrom utils installed.packages
#'
#' @details This function provides the path of the compiled library so that the c++ code can be directly called from c++ routines in other other packages. You need to specify
#' something like PKG_LIBS += '${R_HOME}/bin/Rscript -e "cat(cpgsR::libPath())"'
#' @examples
#' libPath()
#'
#' @export
#'
#'
libPath <- function() {
  cat(sprintf(
    "%s/cpgsR/libs/cpgsR%s",
    installed.packages()["cpgsR", "LibPath"][1],
    .Platform$dynlib.ext
  ))
}
