"RxODE" <-
    function(model, modName = basename(wd), wd = getwd(),
             filename = NULL, do.compile = NULL, flat = FALSE,...)
{
   if(!missing(model) && !missing(filename))
      stop("must specify exactly one of 'model' or 'filename'")

   if (missing(modName) && missing(wd) & missing(flat)){
       flat <- TRUE;
       modName <- "RxODE";
   }
   
   if(!is.null(filename) && missing(model))
      model <- paste0(readLines(filename, ...), collapse="\n")
   
   # RxODE compilation manager (location of parsed code, generated C, 
   # shared libs, etc.)

   cmpMgr <- rx.initCmpMgr(model, modName, wd, flat) 
   # NB: the set of model variables (modelVars) is only available 
   # after parsing, thus it needs to be dynamically computed in cmpMgr
   get.modelVars <- cmpMgr$get.modelVars

   .version <- "0.6"          # object version
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
      inits <- RxODE.inits(inits,modelVars$state);
      params <- params[modelVars$params];
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

   if (is.null(do.compile) & !cmpMgr$isValid()){
       do.compile <- TRUE
   } else if (is.null(do.compile)) {
       do.compile <- FALSE
   }
   if (do.compile) {
      cmpMgr$parse()
      cmpMgr$compile()
   }
   if (cmpMgr$isValid()){
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

RxODE.inits <- function(vec,names,default = 0){
    ret <- vec;
    nv <- names(vec)
    if (!is.null(nv)){
        ret <- ret[nv %in% names];
        missing <- names[!(names %in%  names(ret))];
        if (is.na(default) & length(missing) > 1){
            stop(sprintf("Missing the following parameter(s): %s.",paste(missing,collapse=", ")))
        }
        if (length(missing) > 1){
            ret[missing] <- default;
            warning(sprintf("Assiged %s to %s.",paste(missing,collapse=", "),default))
        }
        ret <- ret[names];
    }
    return(ret);
}

"solve.RxODE" <- function(x,...){
    x$solve(...);
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

igraph <- function(obj,...){
    ## In case anyone else wants to use the method...
    UseMethod("igraph");
}

RxODE.nodeInfo <- function(x, # RxODE object
                           ...){
    ## RxODE.nodeInfo returns a list containing edgeList, biList and nodes
    mod0 <- strsplit(strsplit(gsub("\n","",gsub(" *=","=",gsub("\n *","\n",gsub("#.*","",x$cmpMgr$model)))),";")[[1]],"=");
    mod <- eval(parse(text=sprintf("c(%s)",paste(unlist(lapply(mod0,function(x){
        if (length(x) == 2){
            return(sprintf("\"%s\"=\"%s\"",x[1],gsub(" *$","",gsub("^ *","",x[2])))) 
        } else {
            return(NULL);
        }
    })),collapse=","))))
    modVars <- x$get.modelVars();
    de <- c()
    for (v in modVars$state){
        de <- c(de,mod[sprintf("d/dt(%s)",v)])
    }
    names(de) <- modVars$state;
    for (v in modVars$lhs){
        de <- gsub(sprintf("\\b%s\\b",v),mod[v],de,perl=TRUE);
    }
    fullDe <- de;
    de <- de[which(regexpr("[(]",de) == -1)];
    parDe <- de;
    negDe <- gsub(" *$","",gsub("^ *","",gsub("^ *[^-][^-+]*","",gsub("[+] *[^-+]*","",de))))
    de <- gsub(" *$","",gsub("^ *","",gsub("- *[^-+]*","",de)));
    nde <- names(de)
    de <- de[de != ""];
    de <- strsplit(de," *[+] *");
    negDe <- negDe[negDe != ""];
    negDe <- lapply(strsplit(negDe," *[-] *"),function(x){return(x[x != ""])});
    for (v in names(negDe)){
        negDe[[v]] <- gsub(sprintf("[*] *%s",v),"",negDe[[v]])
    }
    nodes <- c();
    edges <- NULL;
    edgeList <- list();
    biList <- list();
    for (n in names(de)){
        for (n2 in nde){
            w <- which(regexpr(sprintf("[*] *%s",n2),de[[n]]) != -1);
            if (length(w) == 1){
                nodes <- c(nodes,n2,n);
                var <- gsub(sprintf("[*] *%s",n2),"",de[[n]][w]);
                for (e in edgeList){
                    if (e[1] == n & e[2] == n2){
                        biList[[length(biList) + 1]] <- c(n,n2);
                    }
                }
                edgeList[[var]] <- c(n2,n);
                negDe[[n2]] <- negDe[[n2]][negDe[[n2]] != var];
            }
        }
    }
    nodes <- unique(nodes);
    for (i in 1:length(negDe)){
        if (length(negDe[[i]]) == 1){
            edgeList[[negDe[[i]]]] <- c(names(negDe)[i],".out");
            nodes <- c(nodes,".out");
        }
    }
    w <- which(sapply(fullDe,function(x){
        lx <- tolower(x)
        return(regexpr("kin",lx)  != -1 && regexpr("kout",lx) != -1)
    }))
    if (length(w) == 1){
        idr <- fullDe[w];
        idrl <- tolower(idr);
        isKout <- regexpr("^ *kin *[-]",idrl) != -1;
        isInb <- regexpr("[(] *1 *[-]",idrl) != -1;
        idr.from <- gsub(sprintf("^.*(%s).*",paste(names(fullDe)[-w],collapse="|")),"\\1",idr)
        ef <- names(fullDe)[w]
        if (isKout){
            if (isInb){
                Kout <- "Kout\n\u25BC IC50"
                Kin <- "Kin\n";
                idr <- "Indirect\nEffect (II)"
            } else {
                Kout <- "Kout\n\u25B3 EC50";
                Kin <- "Kin\n";
                idr <- "Indirect\nEffect (IV)"
            }
        } else {
            if (isInb){
                Kout <- "Kout"
                Kin <- "Kin\n\u25BC IC50";
                idr <- "Indirect\nEffect (I)"
            } else {
                Kout <- "Kout\n";
                Kin <- "Kin\n\u25B3 EC50";
                idr <- "Indirect\nEffect (III)"
            }
        }
        edgeList[[idr]] <- c(idr.from,ef);
        edgeList[[Kin]] <- c(".Kin",ef);
        edgeList[[Kout]] <- c(ef,".Kout");
        nodes <- c(nodes,ef,".Kin",".Kout");
    }
    return(list(nodes    = nodes,
                edgeList = edgeList,
                biList   = biList));
} # end function RxODE.nodeInfo

igraph.RxODE <- function(x,                                   # RxODE object
                         shape      = c("square","circle","csquare","rectangle","crectangle","vrectangle","sphere","none"),
                         size       = 30,                     # Size of square
                         colors     = c("accent1"="#0460A9"), # Colors
                         fillColor  = "accent1",              # Color to fill 
                         family     = "sans",                 # font Family
                         font       = 2,                      # Font (1: plain, 2: bold, 3: italic, 4: bold/italic, 5:symbol)
                         labelColor = "white",                # Label Color
                         lineColor  = "accent1",              # Line Color
                         shapeEnd   = c("none","sphere","circle","square","csquare","rectangle","crectangle","vrectangle"),
                         sizeEnd    = 10,                     # Size of End
                         tk         = FALSE,                  # For tkplot; transparancy not supported.
                         ...){
    ## igraph.RxODE returns igraph object from RxODE
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package igraph needed for this function to work. Please install it.",
             call. = FALSE)
    }
    with(RxODE.nodeInfo(x),{
        ret <- eval(parse(text=sprintf("igraph::graph_from_literal(%s);",paste(unlist(lapply(edgeList,function(x){sprintf("\"%s\" -+ \"%s\"",x[1],x[2])})),collapse=","))));
        if (length(shape) > 1){
            shape <- shape[1];
        }
        if (length(shapeEnd) > 1){
            shapeEnd <- shapeEnd[1];
        }
        getColor <- function(x){
            if (any(names(colors) == x)){
                return(colors[x]);
            } else {
                return(x);
            }
        }
        ret <- igraph::set_vertex_attr(ret,"shape",value=shape);
        ret <- igraph::set_vertex_attr(ret,"size",value=size);
        ret <- igraph::set_vertex_attr(ret,"color",value=getColor(fillColor));
        ret <- igraph::set_vertex_attr(ret,"label.family",value=family)
        ret <- igraph::set_vertex_attr(ret,"label.font",value=font);
        ret <- igraph::set_vertex_attr(ret,"label.color",value=getColor(labelColor));
        for (n in nodes){
            if (substring(n,0,1) == "."){
                if (!tk){
                    ret <- igraph::set_vertex_attr(ret,"label.color",igraph::V(ret)[n],value="transparent");
                    ret <- igraph::set_vertex_attr(ret,"shape",igraph::V(ret)[n],value=shapeEnd);
                    ret <- igraph::set_vertex_attr(ret,"size",igraph::V(ret)[n],value=sizeEnd);
                } else{
                    ret <- igraph::set_vertex_attr(ret,"size",igraph::V(ret)[n],value=sizeEnd);
                }
            }
        }
        ret <- igraph::set_edge_attr(ret,"color",value=getColor(lineColor))
        ret <- igraph::set_edge_attr(ret,"curved",value=FALSE)
        ret <- igraph::set_edge_attr(ret,"label.font",value=font)
        ## Set bidirectional to curved.
        eval(parse(text=paste(lapply(biList,function(x){sprintf("ret <- igraph::set_edge_attr(ret,\"curved\",igraph::E(ret)[\"%s\" %%--%% \"%s\"],value=TRUE)",x[1],x[2])}),collapse=";")));
        ## add the labels
        for (l in names(edgeList)){
            v <- edgeList[[l]];
            eg <- igraph::E(ret)[v[1] %->% v[2]];
            ret <- igraph::set_edge_attr(ret,"label",eg,value=l);
        }
        return(ret);
    });
} # end function igraph.RxODE

plot.RxODE <- function(x,
                       family="sans",
                       interactive = FALSE,
                       ...)
{
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package igraph needed for this function to work. Please install it.",
             call. = FALSE)
    }
    ig <- igraph(x,family=family,tk=interactive);
    layout <- -igraph::layout.grid(ig);
    
    if (!interactive){
        plot(ig,edge.label.family=family, layout = layout,...)
    } else {
        igraph::tkplot(ig,edge.label.family=family, layout = layout, ...)
    }
    
    ## plot.RxODE returns nothing, but plots the ode graph diagram
    ## Sort of works with the Rgraphviz package...
    ## But Rgraphviz was removed from CRAN, so I wont use it...
    
    ## oldGraphPar <- graph.par()
    ## if (length(engine) > 1){
    ##     engine <- engine[1];
    ## }
    ## if (length(nodeShape) > 1){
    ##     nodeShape <- nodeShape[1];
    ## } 
    ## ret <- new("graphNEL",nodes=nodes,edgemode="directed");
    ## for (i in 1:length(edgeList)){
    ##     n2 <- edgeList[[i]][1];
    ##     n <- edgeList[[i]][2];
    ##     ret <- addEdge(n2,n,ret,1);
    ##     labs[sprintf("%s~%s",n2,n)] <- var;
    ## }
    ## if (!labels){
    ##     nl <- names(labs)
    ##     labs <- rep("",length(nl));
    ##     names(labs) <- nl;
    ## }
    ## attrs <- getDefaultAttrs();
    ## attrs$node$fillcolor <- colors["accent1"];
    ## attrs$node$color <- colors["accent1"];
    ## attrs$node$fontsize <- nodeFontSize;
    ## attrs$node$fixedsize <- FALSE;
    ## attrs$node$fontcolor <- "white"
    ## attrs$node$shape <- nodeShape;
    ## attrs$node$width <- nodeWidth;
    ## attrs$node$height <- nodeHeight;
    ## attrs$edge$color <- colors["accent1"];
    ## attrs$edge$lwd <- edgeLwd;
    ## attrs$edge$fontsize <- edgeFontSize;
    ## attrs$edge$fontcolor <- colors["accent1"];
    ## ret <- agopen(ret,"RxODE",layoutType=engine,recipEdges="distinct",
    ##               edgeAttrs=list(label=labs),attrs=attrs);
    ## for (i in 1:length(AgEdge(ret))){
    ##     AgEdge(ret)[[i]]@txtLabel@labelJust <- labelJust;
    ##     str(AgEdge(ret)[[i]])
    ## }
    ## plot(ret);
    ## graph.par(oldGraphPar);
    invisible(NULL)
} # end function plot.RxODE

"rx.initCmpMgr" <-
function(model, modName, wd, flat)
{
   # Initialize the RxODE compilation manager (all parsing,
   # compilation, and loading of dynamic shared objects is 
   # done through this object. It also keeps tracks of 
   # filenames for the C file, dll, and the translator 
   # model file, parameter file, state variables, etc.
   .mod.parsed <- gsub("[*][*]","^",model); #Support ** operator since R does
   .mod.parsed <- gsub("[<][-]","=",.mod.parsed); #Support <- operator since R does
   .mod.parsed <- gsub("#[^\n]*\n","\n",.mod.parsed); # Strip comments
   .mod.parsed <- gsub("=([^\n;]*);* *\n","=\\1;\n",.mod.parsed); # Don't require semicolons.
   .mod.parsed <- gsub("[.]","_",.mod.parsed); # Allow [.] notation because R does
   .digest <- digest::digest(.mod.parsed);
   .modName <- modName
   .flat <- flat
   .flatOut <- NULL
   if (.flat){
       .wd <- tempfile("RxODE-");
        .flatOut <- wd;
   } else {
       .wd <- wd
   }
   .parsed <- FALSE
   .compiled <- FALSE

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
      if (!is.null(.flatOut)){
           .flatOut <- win.path(.flatOut)
      }
   }

   # filenames and shell command required for parsing (these are unique
   # within a given model directory
   .modfile <- file.path(.mdir, "model.txt")  # copy of user-specified model
   .errfile <- file.path(.mdir, "error.txt")

   # files needed for compiling the C output of the parsed ODE
   .cfile <- file.path(.mdir, sprintf("%s.c", .modName))
   .ofile <- file.path(.mdir, sprintf("%s.o", .modName))
   .dllfile.0 <- file.path(.mdir,sprintf("%s%s", .modName, .Platform$dynlib.ext))
   if (.flat) {
        .dllfile <- file.path(.flatOut,sprintf("%s-%s-%s%s", .modName, .digest, R_ARCH,.Platform$dynlib.ext)) 
   } else {
       .dllfile <- file.path(.mdir,sprintf("%s-%s%s", .modName, R_ARCH,.Platform$dynlib.ext)) 
   }
   .parfile <- file.path(.mdir, "ODE_PARS.txt")
   .stvfile <- file.path(.mdir, "STATE_VARS.txt")
   .lhsfile <- file.path(.mdir, "LHS_VARS.txt")

    if (.flat){
        .modelVarsFile <- file.path(.flatOut,sprintf("%s-%s.Rdata", .modName, .digest))
   } else {
       .modelVarsFile <- file.path(.mdir,sprintf("%s.Rdata", .modName))
   }
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
      sprintf("%s/bin/R CMD SHLIB %s %s -L%s -lodeaux %s", 
         Sys.getenv("R_HOME"), .cfile, .dvode, .libs, .gflibs)

    .dydt <- .calc_lhs <- .ode_solver <- NULL;
    .md5file <- file.path(.mdir,"model_md5");
    if (file.exists(.modelVarsFile) && (flat || (file.exists(.md5file) && readLines(.md5file) == .digest))){
        .objName <- NULL;
        load(.modelVarsFile);
    } else {
        if (file.exists(.mdir) & !flat){
            ## Unload the .dllfile first...
            try(dyn.unload(.dllfile), silent = TRUE)
            unlink(.mdir, recursive = TRUE)
        }
        .modelVars <- list()            # params, state, LHS in the model
   }

   parse <- function(force = FALSE){
      do.it <- force || !.parsed 
      if(!do.it)
          return(invisible(.parsed))
      if(!file.exists(.mdir))
          dir.create(.mdir, recursive = TRUE)
      cat(.mod.parsed, file = .modfile, "\n")   
      
      ## Hack: copy "call_dvode.c" to .mdir to avoid dyn.load() errors
      common <- system.file("common", package = "RxODE")
      if(is.win) 
        common <- win.path(common)
      # replace "dydt" and "calc_lhs" by Rx_ODE_mod_<modName>_* inside the code 
      # of call_dvode.c to avoid symbol conflicts
      safe_name <- paste0(gsub("\\W", "_", .modName),.digest)  # replace non-alphanumeric by "_"
      .dydt <<- paste0("RxODE_mod_", safe_name, "_dydt")
      .calc_lhs <<- paste0("RxODE_mod_", safe_name, "_calc_lhs")
      .ode_solver <<- paste0("RxODE_mod_", safe_name, "_ode_solver")
      src <- readLines(file.path(common, "call_dvode.c"))
      src <- gsub("dydt", .dydt, src)
      src <- gsub("calc_lhs", .calc_lhs, src)
      src <- gsub("ode_solver", .ode_solver, src)
      writeLines(src, file.path(.mdir, "call_dvode.c"))
      .objName <<- .dydt
      writeLines(.digest,.md5file)


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
      save(.modelVars,.parsed,.dydt,.objName,.calc_lhs,.ode_solver,file=.modelVarsFile)
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
         sprintf("PKG_LIBS=-L%s -lodeaux $(BLAS_LIBS) $(FLIBS)", .libs),
         file="Makevars", append=TRUE
      )

      # create SHLIB
      rc <- try(do.call(.sh, list(.shlib)), silent = FALSE)
      if(inherits(rc, "try-error"))
          stop(sprintf("error compiling %s", .cfile))
      if (file.exists(.dllfile.0)){
          file.copy(.dllfile.0,.dllfile)
      }
      if (.flat){
          unlink(.wd, recursive = TRUE)
      }
      
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
        valid <- file.exists(.dllfile) & file.exists(.modelVarsFile);
        if (valid){
            if (!.flat){
                valid <- file.exists(.md5file);
                if (valid){
                    valid <- readLines(.md5file) == .digest;
                }
            }
        }
        return(valid);
    }

   get.index = function(s) {
      # return the (one) state varible index
      if(length(s)!=1) 
          warning("only one state variable should be input", immediate = TRUE)
       ix <- match(s, scan(.stvfile, what = "", quiet = T), nomatch=0)
      if(!ix) stop(paste("state var not found:", s))
      ix
   }

    if (isValid()){
        .compiled <- TRUE
    }

   out <- 
       list(parse = parse, compile = compile,
            model = .mod.parsed,
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

