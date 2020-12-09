.onAttach <- function(libname, pkgname) {
  msg <- paste0("Functions related to Stewart's potential are deprecated.\n",
               "Please use the `potential` package instead.\n",
               "https://riatelab.github.io/potential/")
  packageStartupMessage(msg)
}