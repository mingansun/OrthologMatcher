
# show the message when the package is attached
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to use OrthologMatcher!")
}

