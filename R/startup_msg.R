# I copied this piece of code from the bootnet package where it was copied from Lavaan mainly:

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage("Please report issues on Github: https://github.com/jmbh/mgm/issues")
}
