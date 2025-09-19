.onAttach <- function(libname, pkgname) {
  local_ver <- utils::packageVersion(pkgname)
  cran_ver  <- "9.2.0"  # latest CRAN release

  if (local_ver < cran_ver) {
    packageStartupMessage(
      sprintf("%s %s: \n", pkgname, local_ver),
      "Please note this version differs slightly from that published in ",
      "Behavior Research Methods. Please install the new version from CRAN. ",
      "See the README for how to update: https://github.com/JimGrange/flankr"
    )
  }
}
