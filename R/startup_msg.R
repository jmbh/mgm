# I copied this piece of code from the bootnet package where it was copied from Lavaan mainly:

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage("Note that the syntax changed considerably from version 1.1-8 to 1.2-0.")
  packageStartupMessage("Please report any bugs: https://github.com/jmbh/mgm/issues")
}
