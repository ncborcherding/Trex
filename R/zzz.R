.onLoad <- function(libname, pkgname) {
  .split_dependencies(pkgname)
}

TrexEnv <- NULL

.onLoad <- function(libname, pkgname) {
  TrexEnv <<- basilisk::BasiliskEnvironment(
    envname = "TrexEnv",
    pkgname = "Trex",
    packages = c(
        "pandas==2.2.2",
        "tensorflow==2.16.2",
        "keras==3.4.1"
    )
  )
}
