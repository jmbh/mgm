# I copied this piece of code from the bootnet package where it was copied from Lavaan mainly:

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  pkgname <- paste0(pkgname, "-parallel")
  packageStartupMessage("This is ", pkgname)
  packageStartupMessage("A parallel version of the mgm package")
}
