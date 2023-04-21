#' Get path to example csvs (heavily inspired by readr::readr_example)
#'
#' archnetan includes csvs used during the analysis of corresponding publications included in
#' inst/examples which this function helps make easier to access
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @export
#' @examples
#' archnetan_example()
#' archnetan_example("originalDataset.csv")
archnetan_example <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("examples", package = "archnetan"))
  } else {
    system.file("examples", file, package = "archnetan", mustWork = TRUE)
  }
}
