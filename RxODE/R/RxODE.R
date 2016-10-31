"RxODE" <-
function(model, modName = basename(wd), wd = getwd(), 
   filename = NULL, do.compile = TRUE, ...)
{
   if(!missing(model) && !missing(filename))
      stop("must specify exactly one of 'model' or 'filename'")

   if(!is.null(filename) && missing(model))
      model <- paste0(readLines(filename, ...), collapse="\n")

   # RxODE compilation manager (location of parsed code, generated C, 
   # shared libs, etc.)

   cmpMgr <- rx.initCmpMgr(model, modName, wd) 
   # NB: the set of model variables (modelVars) is only available 
   # after parsing, thus it needs to be dynamically computed in cmpMgr
   get.modelVars <- cmpMgr$get.modelVars

   .version <- "0.5"          # object version
   .last.solve.args <- NULL   # to be populated by solve()

   solve <- 
   function(params, events, inits = NULL, stiff = TRUE, transit_abs = FALSE, 
        atol = 1.0e-8, rtol = 1.0e-6, ...)
   {
      event.table <- events$get.EventTable()
      modelVars <- get.modelVars()

      # preserve input arguments. 
      .last.solve.args <<-
         list(params = params, events = events$copy(),
              inits = inits, stiff = stiff, 
              transit_abs = transit_abs, atol = atol, rtol = rtol, ...)

      # check that starting values for all needed parameters are 
      # specified in the input "params" 
      if (length(setdiff(modelVars$params, names(params)))) {
         msg <- 
            paste("var(s) not found in input pars.\n", 
               paste(setdiff(modelVars$params, names(params)), collapse=" "))
         stop(msg)
      }
      params <- params[modelVars$params]
      s <- as.list(match.call(expand.dots = TRUE)) 
      wh <- grep(pattern="S\\d+$", names(s))[1]
      # HACK: fishing scaling variables "S1 S2 S3 ..." from params call
      # to solve(). Maybe define a "scale=c(central=7.6, ...)" argument
      # similar to "params="?
      scaler.ix <- 0
      if (!is.na(wh)) {
         if (s[[wh]] %in% names(params)) {
            scaler <- params[s[[wh]]]
            scaler.ix <- as.numeric(substring(names(s)[wh], 2))
         } else {
            warning(paste("scaler variable not found:", s[[wh]]))
         }
      }

      state_vars <- modelVars$state
      neq  <- length(state_vars)
      lhs_vars <- modelVars$lhs
      nlhs <- length(lhs_vars)

      ntime <- dim(event.table)[1]
      ret   <- rep(0.0, ntime*neq)
      lhs   <- rep(0.0, ntime*nlhs)
      rc <- as.integer(0)  # return code 0 (success) or IDID in call_dvode.c
      if (is.null(inits)) inits <- rep(0.0, neq)

      # may need to reload (e.g., when we re-start R and
      # re-instantiate the RxODE object from a save.image.
      cmpMgr$dynLoad()
      xx <- .C(cmpMgr$ode_solver,
            as.integer(neq),
            as.double(params),
            as.double(event.table$time),
            as.integer(event.table$evid),
            length(event.table$time),
            as.double(inits),
            as.double(event.table$amt[event.table$evid>0]),
            as.double(ret),
            as.double(atol),
            as.double(rtol),
            as.integer(stiff),
            as.integer(transit_abs),
            as.integer(nlhs),
            as.double(lhs),
            rc
            )

      rc <- xx[[length(xx)]]
      if(rc!=0)
         stop(sprintf("could not solve ODE, IDID=%d (see further messages)", rc))
      x <- cbind(
            matrix(xx[[8]], ncol=neq, byrow=T),
            if(nlhs) matrix(xx[[14]], ncol=nlhs, byrow=T) else NULL
         )
      colnames(x) <- c(state_vars, lhs_vars)

      if (scaler.ix) {
         x[, scaler.ix] <- x[, scaler.ix]/scaler
      }
      
      
      cbind(time=event.table$time, x)[events$get.obs.rec(),]
   }

   if (do.compile) {
      cmpMgr$parse()
      cmpMgr$compile()
      cmpMgr$dynLoad()
   }

   out <- 
      list(modName = modName, 
         model = model,           # actual model code
         get.modelVars = get.modelVars,  # extract model variables (pars, lhs, etc)
         solve = solve, 
         cmpMgr = cmpMgr, 
         dynLoad = cmpMgr$dynLoad, 
         dynUnload = cmpMgr$dynUnload,
         isValid = cmpMgr$isValid, 
         version = .version,
         delete = cmpMgr$delete, 
         # the next is for backward compatibility and will be deprecated
         parse = cmpMgr$parse, compile = cmpMgr$compile, 
         get.index = cmpMgr$get.index,
         run = solve, 
         getObj = function(obj) get(obj, envir = environment(solve)))
   class(out) <- "RxODE"
   out
}

"print.RxODE" <-
function(x, ...)
{
   valid <- x$cmpMgr$isValid()
   if(!valid){
      .msg <- "invalid object, needs to be re-created"
   } else {
      .ready <- x$cmpMgr$getObj(".compiled")
      .msg <- if(.ready) "ready to run" else "needs compilation"
   }
   cat(sprintf('RxODE model named "%s" (%s)\n', x$modName, .msg))
   invisible(x)
}

"rx.initCmpMgr" <-
function(model, modName, wd)
{
   # Initialize the RxODE compilation manager (all parsing,
   # compilation, and loading of dynamic shared objects is 
   # done through this object. It also keeps tracks of 
   # filenames for the C file, dll, and the translator 
   # model file, parameter file, state variables, etc.

   .modName <- modName  
   .wd <- wd  

   .parsed <- FALSE
   .compiled <- FALSE
   .modelVars <- list()            # params, state, LHS in the model

   is.win <- .Platform$OS.type=="windows"
   win.path <- function(p) gsub("\\\\", "/", utils::shortPathName(p))
   R_ARCH <- .Platform$r_arch

   # locations for the various libs and header files, plus 
   # commands needed to parse, compile, solve ODE-defined models
   .bin  <- system.file(file.path("bin",  R_ARCH), package = "RxODE")
   .libs <- system.file(file.path("libs", R_ARCH), package = "RxODE")
   .incl <- system.file("include", package = "RxODE")

   # model-specific directory under .md (default current dir)
   .mdir <- file.path(.wd, sprintf("%s.d", .modName))

   if(is.win){  
      # on Windows ensure there are no spaces in path names and 
      # no backslashes to avoid tripping system("R CMD SHLIB")
      .wd <- win.path(.wd)
      .bin <- win.path(.bin)
      .libs <- win.path(.libs)
      .incl <- win.path(.incl)
      .mdir <- win.path(.mdir)
   }

   # filenames and shell command required for parsing (these are unique
   # within a given model directory
   .modfile <- file.path(.mdir, "model.txt")  # copy of user-specified model
   .errfile <- file.path(.mdir, "error.txt")

   # files needed for compiling the C output of the parsed ODE
   .cfile <- file.path(.mdir, sprintf("%s.c", .modName))
   .ofile <- file.path(.mdir, sprintf("%s.o", .modName))
   .dllfile.0 <- sprintf("%s%s", .modName, .Platform$dynlib.ext) 
   .dllfile <- file.path(.mdir, .dllfile.0)
   .parfile <- file.path(.mdir, "ODE_PARS.txt")
   .stvfile <- file.path(.mdir, "STATE_VARS.txt")
   .lhsfile <- file.path(.mdir, "LHS_VARS.txt")

   # shell calls to parsing and compiling (via R CMD SHLIB)
   .gflibs <- if (is.win) "-lRblas -lgfortran" else "-lgfortran" # TODO: check do we need this
   #.sh <- if(is.win) "shell" else "system"
   .sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths
   .dvode <- file.path(.mdir, "call_dvode.o")
   .prefix <- paste(.mdir, "/", sep="")
   .parse.cmd <- 
      sprintf("%s/tran.exe %s %s %s 2>%s", 
         .bin, .modfile, .cfile, .prefix, .errfile)
   .shlib <- 
      sprintf("%s/bin/R CMD SHLIB %s %s", 
         Sys.getenv("R_HOME"), .cfile, .dvode)
   if(!file.exists(.mdir))
      dir.create(.mdir, recursive = TRUE)

   cat(model, file = .modfile, "\n")   

   # Hack: copy "call_dvode.c" to .mdir to avoid dyn.load() errors
   common <- system.file("common", package = "RxODE")
   if(is.win) 
      common <- win.path(common)
   # replace "dydt" and "calc_lhs" by Rx_ODE_mod_<modName>_* inside the code 
   # of call_dvode.c to avoid symbol conflicts
   safe_name <- gsub("\\W", "_", .modName)  # replace non-alphanumeric by "_"
   .dydt <- paste0("RxODE_mod_", safe_name, "_dydt")
   .calc_lhs <- paste0("RxODE_mod_", safe_name, "_calc_lhs")
   .ode_solver <- paste0("RxODE_mod_", safe_name, "_ode_solver")
   src <- readLines(file.path(common, "call_dvode.c"))
   src <- gsub("dydt", .dydt, src)
   src <- gsub("calc_lhs", .calc_lhs, src)
   src <- gsub("ode_solver", .ode_solver, src)
   writeLines(src, file.path(.mdir, "call_dvode.c"))
   .objName <- .dydt

   parse <- function(force = FALSE){
      do.it <- force || !.parsed 
      if(!do.it)
         return(invisible(.parsed))

      rc <- do.call(.sh, list(.parse.cmd))  # parse command (shell)
      if(file.exists(.errfile))
         err <- readLines(.errfile)
      else
         err <- NULL
      if(rc!=0 && length(err))
         stop(sprintf("error parsing ODE system in %s: %s", .modfile, err))

      # HACK: replace the common dydt by RxODE_mod_<mod>_dydt
      src <- gsub("dydt", .dydt, readLines(.cfile))
      src <- gsub("calc_lhs", .calc_lhs, src)
      writeLines(src, .cfile)    # overwrite with replaced symbol
      .parsed <<- TRUE
      .compiled <<- FALSE
      # we now know the params and set of variables defined in the model
      # TODO: What about the builtin variables?
      # buitin <- c("t", "tlast", "podo")
      .modelVars <<- 
         list(
            params = scan(.parfile, what  = "", quiet = TRUE),
            state = scan(.stvfile, what = "", quiet = TRUE),
            lhs = scan(.lhsfile, what = "", quiet = TRUE)
         )

      invisible(.parsed)
   }

   compile <- function(force = FALSE){
      do.it <- force || !.parsed || !.compiled
      if(!do.it)
         return(invisible(.compiled))
      
      if(!.parsed)
         parse(force = force)

      # may need to unload previous model object code
      if (is.loaded(.objName)) try(dyn.unload(.dllfile), silent = TRUE)

      #on.exit(unlink("Makevars"))
      cat(sprintf("PKG_CPPFLAGS=-I%s\n",.incl), file="Makevars")
      cat(
         sprintf("PKG_LIBS=-L%s -lodeaux $(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)", .libs),
         file="Makevars", append=TRUE
      )

      # create SHLIB
      rc <- try(do.call(.sh, list(.shlib)), silent = FALSE)
      if(inherits(rc, "try-error"))
         stop(sprintf("error compiling %s", .cfile))

      # dyn load it
      rc <- try(dyn.load(.dllfile), silent = TRUE)
      if(inherits(rc, "try-error"))
         stop(sprintf("error loading dll file %s", .dllfile))
      
      .compiled <<- TRUE

      invisible(.compiled)
   }

   dynLoad <- function(force = FALSE){
      # NB: we may need to reload (e.g., when we re-start R and
      # re-instantiate the RxODE object from a save.image.
      if(is.loaded(.objName)){
         if(!force){
            return()
         } else {
            dyn.unload(.dllfile)
         }
      } 
      rc <- try(dyn.load(.dllfile), silent = TRUE)
      if(inherits(rc, "try-error"))
         stop(sprintf("error loading dll file %s", .dllfile))
   }

   dynUnload <- function(){
      if(!is.loaded(.objName))
         return()
      rc <- try(dyn.unload(.dllfile), silent = TRUE)
      if(inherits(rc, "try-error"))
         stop(sprintf("error unloading dll file %s", .dllfile))
   }

   delete <- function(){
      if (is.loaded(.objName))
         try(dyn.unload(.dllfile), silent = TRUE)
      .parsed <<- FALSE
      .compiled <<- FALSE
      .modelVars <<- NULL
      # TODO: should we remove all objects in the closure?
      # as the object is no longer valid. Need a valid.object()
      unlink(.mdir, recursive = TRUE) # leave dir
   }

   isValid <- function(){
      # need a better valid.object()
      file.exists(.mdir) 
   }

   get.index = function(s) {
      # return the (one) state varible index
      if(length(s)!=1) 
         warning("only one state variable should be input", immediate = TRUE)
      ix <- match(s, scan(.stvfile, what = "", quiet = T), nomatch=0)
      if(!ix) stop(paste("state var not found:", s))
      ix
   }

   out <- 
      list(parse = parse, compile = compile, 
         dynLoad = dynLoad, dynUnload = dynUnload,
         ode_solver = .ode_solver,   # name of C function
         modelDir = .mdir,               # model directory
         dllfile = .dllfile,
         get.modelVars = function() .modelVars,   
         isValid = isValid, delete = delete,
         get.index = get.index, 
         getObj = function(obj) get(obj, envir = environment(parse))
      )
   class(out) <- "RxCompilationManager"
   out
}

"print.RxCompilationManager" <-
function(x, ...)
{
   modName <- x$getObj(".modName")
   cat(sprintf("RxCompilationManager for RxODE model '%s'\n", modName))
   invisible(x)
}

