.pkgModelCurrent <- TRUE

.setPkgModels <- function(value) { ## For testing
  assignInMyNamespace(".pkgModelCurrent", value)
}

.norm2 <- function(obj) {
  if (inherits(obj, "RxODE")) {
    if (exists(".linCmtM", obj)) {
      return(get(".linCmtM", obj))
    }
  }
  setNames(RxODE::rxModelVars(obj)$model["normModel"], NULL)
}

.isWritable <- function(...) {
  .ret <- try(assertthat::is.writeable(...), silent = TRUE)
  if (inherits(.ret, "try-error")) {
    .ret <- FALSE
  }
  .ret
}

.rxPkgInst <- function(obj) {
  .wd <- getwd()
  if (regexpr(obj$package, .wd) != -1) {
    .inst <- gsub(paste0("(", obj$package, ").*"), "\\1", .wd)
  } else {
    .inst <- system.file(package = obj$package)
  }
  if (.isWritable(.inst)) {
    if (regexpr("inst$", .inst) != -1) {
      return(.inst)
    }
    .inst2 <- file.path(.inst, "inst")
    if (file.exists(.inst2)) {
      return(.inst2)
    }
    .html <- file.path(.inst, "html")
    if (file.exists(.html)) {
      return(.inst)
    }
    return(.inst2)
  } else {
    .inst <- "~/.rxCache/"
    if (.isWritable(.inst)) {
      return(.inst)
    }
    return(rxTempDir())
  }
}

.rxPkgDir <- function(obj) {
  return(file.path(.rxPkgInst(obj), "rx"))
}

.rxPkgDll <- function(obj) {
  obj$mdir <- .rxPkgDir(obj)
  .pkgInfo <- getLoadedDLLs()[[obj$package]]
  if (!all(is.null(.pkgInfo))) {
    if (obj$isValid()) {
      .tmp <- .pkgInfo
      class(.tmp) <- "list"
      return(.tmp$path)
    } else {
      return(file.path(obj$mdir, basename(obj$rxDll$dll)))
    }
  } else {
    return(file.path(obj$mdir, basename(obj$rxDll$dll)))
  }
}

.rxNewMvStr <- function(obj) {
  gsub("[.].*", paste0("_new_", .Platform$r_arch, "_model_vars"), basename(obj$rxDll$dll))
}

.rxPkgLoaded <- function(pkg) {
  .si <- sessionInfo()
  return(length(intersect(pkg, c(
    names(.si$otherPkgs) ## ,names(.si$loadedOnly)
  ))) != 0)
}

.rxUseI <- new.env(parent = emptyenv())
# 1
.rxUseI$i <- 1L
.rxUseCdir <- ""

#' Use model object in your package
#' @param obj model to save.
#' @param internal If this is run internally.  By default this is FALSE
#' @inheritParams usethis::use_data
#' @return Nothing; This is used for its side effects and shouldn't be called by a user
#' @export
rxUse <- function(obj, overwrite = TRUE, compress = "bzip2",
                  internal = FALSE) {
  rxReq("usethis")
  rxReq("devtools")
  internal <- internal
  if (missing(obj)) {
    .env <- new.env()
    assign("internal", internal, .env)
    assign("overwrite", overwrite, .env)
    assign("compress", compress, .env)
    sapply(list.files(devtools::package_file("inst/rx"), full.names = TRUE),
      unlink,
      force = TRUE, recursive = TRUE
    )
    .models <- NULL
    for (.f in list.files(
      path = devtools::package_file("data"),
      pattern = "\\.rda$", full.names = TRUE
    )) {
      load(.f, envir = .env)
      .f2 <- basename(.f)
      .f2 <- substr(.f2, 0, nchar(.f2) - 4)
      if (is(.env[[.f2]], "RxODE")) {
        .env[[.f2]]$package <- NULL
        .minfo(sprintf("recompile '%s'", .f2))
        .models <- c(.models, .f2)
        eval(parse(text = sprintf("rxUse(%s, internal=internal, overwrite=overwrite, compress=compress)", .f2)),
          envir = .env
        )
        .docFile <- file.path(devtools::package_file("R"), paste0(.f2, "-doc.R"))
        if (!file.exists(.docFile)) {
          (sprintf("creating documentation '%s'", .docFile))
          sink(.docFile)
          .tmp <- .env[[.f2]]
          .mv <- rxModelVars(.tmp)
          cat(sprintf("#' %s RxODE model\n", .f2))
          cat("#'\n")
          cat(sprintf(
            "#' @format An \\emph{RxODE} model with %s parameters, %s ODE states, and %s calc vars.\n",
            length(.tmp$params), length(.tmp$state) + .mv$extraCmt, length(.tmp$lhs)
          ))
          cat("#'\n")
          cat(sprintf("#'\\emph{Parameters (%s$params)}\n", .f2))
          cat("#'\n")
          cat("#' \\describe{\n")
          .def <- rxInits(.tmp)
          .defs <- paste0(" (default=", .def, ")")
          .defs[is.na(.def)] <- ""
          cat(paste(paste0("#'   \\item{", .tmp$params, "}{", .defs, "}\n"), collapse = ""))
          cat("#'}\n")
          .state <- .tmp$state
          ##
          if (.mv$extraCmt == 2) {
            .state <- c(.state, "depot", "central")
          } else if (.mv$extraCmt == 1) {
            .state <- c(.state, "central")
          }
          if (length(.state) > 0) {
            cat("#'\n")
            cat(sprintf("#' \\emph{State %s$state}\n", .f2))
            cat("#'\n")
            cat("#' \\describe{\n")
            cat(paste(paste0("#'   \\item{", .state, "}{ (=", seq_along(.state), ")}\n"), collapse = ""))
            cat("#' }\n")
          }
          .lhs <- .tmp$lhs
          if (length(.lhs) > 0) {
            cat("#'\n")
            cat(sprintf("#' \\emph{Calculated Variables %s$lhs}\n", .f2))
            cat("#'\n")
            cat("#' \\describe{\n")
            cat(paste(paste0("#'   \\item{", .lhs, "}{}\n"), collapse = ""))
            cat("#' }\n")
          }
          cat("#'\n")
          cat("#' \\emph{Model Code}\n") # sprintf(,.f2)
          cat("#'\n")
          .code <- deparse(body(eval(parse(text = paste("function(){", .norm2(.tmp), "}")))))
          .code[1] <- "RxODE({"
          .code[length(.code)] <- "})"
          cat(paste(paste0("#' ", .code, "\n"), collapse = ""))
          cat("#'\n")
          cat(paste(paste0("#' @seealso \\code{\\link[RxODE]{eventTable}}, \\code{\\link[RxODE]{et}}, \\code{\\link[RxODE]{rxSolve}}, \\code{\\link[RxODE]{RxODE}}\n")))
          cat("#' \n")
          cat("#' @examples\n")
          cat("#' ## Showing the model code\n")
          cat(sprintf("#' summary(%s)\n", .f2))
          cat("#'\n")
          cat(sprintf('"%s"\n', .f2))
          sink()
        }
      }
    }
    if (!dir.exists(devtools::package_file("src"))) {
      dir.create(devtools::package_file("src"), recursive = TRUE)
    }
    .pkg <- basename(usethis::proj_get())
    .rx <- loadNamespace("RxODE")
    sapply(
      list.files(.rxUseCdir, pattern = "[.]c", full.names = TRUE),
      function(x) {
        .minfo(sprintf("copy '%s'", basename(x)))
        .env <- .rx$.rxUseI
        .rxUseI <- .env$i
        .f0 <- gsub(
          "^#define (.*) _rx(.*)$",
          paste0("#define \\1 _rxp", .rxUseI, "\\2"), readLines(x)
        )
        assign("i", .rxUseI + 1, envir = .env)
        .f0 <- c("#include <RxODE.h>\n#include <RxODE_model_shared.h>", .f0)
        .w <- which(.f0 == "#include \"extraC.h\"")
        if (length(.w) > 0) .f0 <- .f0[-.w[1]]
        writeLines(text = .f0, con = file.path(devtools::package_file("src"), basename(x)))
      }
    )
    .inits <- paste0("R_init0_", .pkg, "_", .models)
    .tmp <- paste0("{\"", .pkg, "_", .models, "_model_vars\", (DL_FUNC) &", .pkg, "_", .models, "_model_vars, 0},\\")
    .tmp[length(.tmp)] <- substr(.tmp[length(.tmp)], 0, nchar(.tmp[length(.tmp)]) - 1)
    .extraC <- c(
      "#define compiledModelCall \\",
      .tmp,
      paste0("SEXP ", .pkg, "_", .models, "_model_vars();"),
      paste0("void ", .inits, "();"),
      paste0("void R_init0_", .pkg, "_RxODE_models(){"),
      paste0("  ", .inits, "();"),
      "}"
    )
    sink(file.path(devtools::package_file("src"), paste0(.pkg, "_compiled.h")))
    if (.pkg == "RxODE") {
      cat("#include <R.h>\n#include <Rinternals.h>\n#include <stdlib.h> // for NULL\n#include <R_ext/Rdynload.h>\n#include \"../inst/include/RxODE.h\"\n#include \"../inst/include/RxODE_model_shared.h\"\n")
    } else {
      cat("#include <R.h>\n#include <Rinternals.h>\n#include <stdlib.h> // for NULL\n#include <R_ext/Rdynload.h>\n#include <RxODE.h>\n#include <RxODE_model_shared.h>\n")
    }
    cat(paste(.extraC, collapse = "\n"))
    cat("\n")
    sink()
    .files <- list.files(devtools::package_file("src"))
    .files <- .files[regexpr("RxODE_model_shared", .files) == -1]
    if (all(regexpr(paste0("^", .pkg), .files) != -1)) {
      .minfo(sprintf("only compiled models in this package, creating '%s_init.c'", .pkg))
      sink(file.path(devtools::package_file("src"), paste0(.pkg, "_init.c")))
      cat("#include <R.h>\n#include <Rinternals.h>\n#include <stdlib.h> // for NULL\n#include <R_ext/Rdynload.h>\n")
      cat("#include <RxODE.h>\n")
      cat("#include <RxODE_model_shared.h>\n")
      cat(paste0('#include "', .pkg, '_compiled.h"\n'))
      cat(sprintf("void R_init_%s(DllInfo *info){\n", .pkg))
      cat(sprintf("  R_init0_%s_RxODE_models();\n", .pkg))
      cat("  static const R_CallMethodDef callMethods[]  = {\n  compiledModelCall\n  {NULL, NULL, 0}\n  };\n")
      cat("  R_registerRoutines(info, NULL, callMethods, NULL, NULL);\n")
      cat("  R_useDynamicSymbols(info,FALSE);\n")
      cat("}\n")
      cat(paste(paste0("void R_unload_", .pkg, "_", .models, "(DllInfo *info);\n"), collapse = ""))
      cat(sprintf("void R_unload_%s(DllInfo *info){\n", .pkg))
      cat(paste(paste0("  R_unload_", .pkg, "_", .models, "(info);\n"), collapse = ""))
      cat("}\n")
      sink()
    }
    if (!file.exists(devtools::package_file("R/rxUpdated.R")) && .pkg != "RxODE") {
      sink(devtools::package_file("R/rxUpdated.R"))
      cat(".rxUpdated <- new.env(parent=emptyenv())\n")
      sink()
    }
    unlink(devtools::package_file("inst/rx"), recursive = TRUE, force = TRUE)
    if (length(list.files(devtools::package_file("inst"))) == 0) {
      unlink(devtools::package_file("inst"), recursive = TRUE, force = TRUE)
    }
    return(invisible(TRUE))
  } else {
    .modName <- as.character(substitute(obj))
    .pkg <- basename(usethis::proj_get())
    .env <- new.env(parent = baseenv())
    assign(.modName, RxODE(.norm2(obj), package = .pkg, modName = .modName), .env)
    assignInMyNamespace(".rxUseCdir", dirname(rxC(.env[[.modName]])))
    assign("internal", internal, .env)
    assign("overwrite", overwrite, .env)
    assign("compress", compress, .env)
    eval(parse(text = sprintf("usethis::use_data(%s, internal=internal, overwrite=overwrite, compress=compress)", .modName)), envir = .env)
  }
}

#' Creates a package from compiled RxODE models
#'
#' @param ... Models to build a package from
#' @param package String of the package name to create
#' @param action Type of action to take after package is created
#' @param name Full name of author
#' @param license is the type of license for the package.
#' @inheritParams usethis::create_package
#' @inheritParams RxODE
#' @author Matthew Fidler
#' @return this function returns nothing and is used for its side effects
#' @export
rxPkg <- function(..., package,
                  wd = getwd(),
                  action = c("install", "build", "binary", "create"),
                  license = c("gpl3", "lgpl", "mit", "agpl3"),
                  name = "Firstname Lastname",
                  fields = list()) {
  if (missing(package)) {
    stop("'package' needs to be specified")
  }
  action <- match.arg(action)
  license <- match.arg(license)
  .owd <- getwd()
  .op <- options()
  on.exit({
    setwd(.owd)
    options(.op)
  })
  .dir <- wd
  if (!dir.exists(.dir)) {
    dir.create(.dir)
  }
  setwd(.dir)
  options(
    usethis.description = list(`Title` = "This is generated from RxODE"),
    usethis.full_name = ifelse(missing(name), getOption("usethis.full_name", "Firstname Lastname"), name)
  )
  .dir2 <- file.path(.dir, package)
  usethis::create_package(.dir2,
    fields = fields,
    rstudio = FALSE,
    roxygen = TRUE,
    check_name = TRUE,
    open = FALSE
  )
  setwd(.dir2)
  usethis::use_package("RxODE", "LinkingTo")
  usethis::use_package("RxODE", "Depends")
  if (license == "gpl3") {
    usethis::use_gpl3_license()
  } else if (license == "lgpl") {
    usethis::use_lgpl_license()
  } else if (license == "agpl3") {
    usethis::use_agpl3_license()
  } else if (license == "mit") {
    usethis::use_mit_license()
  }
  .p <- devtools::package_file("DESCRIPTION")
  writeLines(c(
    readLines(.p),
    "NeedsCompilation: yes",
    "Biarch: true"
  ), .p)
  ## Now use rxUse for each item
  .env <- new.env()
  .lst <- as.list(match.call()[-1])
  .w <- which(names(.lst) == "")
  .lst <- .lst[.w]

  for (.i in seq_along(.lst)) {
    .v <- as.character(deparse(.lst[[.i]]))
    assign(.v, eval(.lst[[.i]], envir = parent.frame(1)), .env)
    print(.env[[.v]])
    eval(parse(text = sprintf("rxUse(%s)", .v)), envir = .env)
  }

  ## Final rxUse to generate all code
  rxUse()

  .p <- file.path(devtools::package_file("R"), "rxUpdated.R")
  .f <- readLines(.p)
  .w <- which(regexpr("@useDynLib", .f) != -1)

  if (length(.w) == 0) {
    .f <- c(
      paste0("#' @useDynLib ", package, ", .registration=TRUE"),
      "#' @import RxODE",
      .f
    )
    writeLines(.f, .p)
  }
  devtools::document()
  if (!file.exists("configure.win")) {
    writeLines(c(
      "#!/bin/sh",
      "echo \"unlink('src', recursive=TRUE);RxODE::rxUse()\" > build.R",
      "${R_HOME}/bin/Rscript build.R",
      "rm build.R"
    ), "configure.win")
  }
  if (!file.exists("configure")) {
    writeLines(c(
      "#!/bin/sh",
      "echo \"unlink('src', recursive=TRUE);RxODE::rxUse()\" > build.R",
      "${R_HOME}/bin/Rscript build.R",
      "rm build.R"
    ), "configure")
    if (!file.exists("configure.ac")) {
      writeLines(
        "## dummy autoconf script",
        "configure.ac"
      )
    }
  }
  if (action == "install") {
    devtools::install()
  } else if (action == "build") {
    devtools::build()
  } else if (action == "binary") {
    devtools::build(binary = TRUE)
  }
  invisible()
}
