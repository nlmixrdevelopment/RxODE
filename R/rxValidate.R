##' Validate RxODE
##'
##' This allows easy validation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param type Type of test or fitler of test type
##' @param check Use devtools::check instead
##' @author Matthew L. Fidler
##' @export
rxValidate <- function(type = NULL, check = FALSE) {
  ## rxVersion(" Validation", TRUE);
  .tests <- c(
    "cran", "norm", "demo", "lvl2", "parsing",
    "focei", "indLin", "lincmt",
    "plot", "print",
    "parseLincmt"
  )
  if (is.character(type)) {
    if (type == "covr") {
      Sys.setenv("NOT_CRAN" = "true", "covr" = "true")
      on.exit({
        setwd(old.wd)
        Sys.unsetenv("rxCran")
      })
      covr::report()
    } else {
      if (any(type == .tests)) {
        if (check) {
          devtools::check(env_vars = c(
            NOT_CRAN = "true",
            rxCran = type
          ))
        } else {
          old.wd <- getwd()
          on.exit({
            setwd(old.wd)
            Sys.unsetenv("rxCran")
          })
          Sys.setenv("rxCran" = type)
          path <- file.path(system.file("tests", package = "RxODE"), "testthat")
          setwd(path)
          pt <- proc.time()
          testthat::test_dir(path)
          message("================================================================================")
          print(proc.time() - pt)
          message("================================================================================")
        }
      } else {
        old.wd <- getwd()
        on.exit({
          setwd(old.wd)
          Sys.unsetenv("rxCran")
        })
        Sys.setenv("rxCran" = "true")
        path <- file.path(system.file("tests", package = "RxODE"), "testthat")
        setwd(path)
        pt <- proc.time()
        testthat::test_dir(path, filter = type)
        message("================================================================================")
        print(proc.time() - pt)
        message("================================================================================")
      }
    }
  } else {
    old.wd <- getwd()
    on.exit({
      setwd(old.wd)
      Sys.unsetenv("rxCran")
    })
    path <- file.path(system.file("tests", package = "RxODE"), "testthat")
    setwd(path)
    for (t in .tests) {
      Sys.setenv("rxCran" = t)
      message(sprintf("%s only tests", t))
      message("================================================================================")
      pt <- proc.time()
      testthat::test_dir(path)
      message("================================================================================")
      print(proc.time() - pt)
      message("================================================================================")
    }
  }
}

##' @rdname rxValidate
##' @export
rxTest <- rxValidate
