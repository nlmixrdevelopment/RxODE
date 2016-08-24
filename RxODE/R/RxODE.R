"RxODE" <-
    function(model, modName = basename(wd), wd = getwd(),
             filename = NULL, do.compile = NULL, flat = FALSE,
             extra.c=NULL,debug = FALSE,...)
{
    if(!missing(model) && !missing(filename))
        stop("must specify exactly one of 'model' or 'filename'")
    
   if (missing(modName) && missing(wd) & missing(flat)){
       flat <- TRUE;
       if (!is.null(filename) && missing(model)){
           modName <- basename(filename);
       } else if (missing(filename)&& file.exists(model)){
           modName <- basename(model)
       } else {
           modName <- "RxODE";
       }
   }
    
   if(!is.null(filename) && missing(model))
       model <- paste0(readLines(filename, ...), collapse="\n")

   if (missing(filename) && file.exists(model)){
       model <- paste0(readLines(model, ...), collapse="\n")
   }
   
   # RxODE compilation manager (location of parsed code, generated C, 
                                        # shared libs, etc.)

    cmpMgr <- rx.initCmpMgr(model, modName, wd, flat, extra.c,debug);
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
      
      
      ret <- cbind(time=event.table$time, x)[events$get.obs.rec(),];
      return(ret)
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

RxODE.quoteVar <- function(var,word=TRUE){
    ret <- c();
    for (i in 1:nchar(var)){
        ret <- c(ret,sprintf("one_of(\"%s\")",substr(var,i,i)));
    }
    ret <- parse(text=sprintf(ifelse(word,"rex::rex(boundary,%s,boundary)",
                                     "rex::rex(%s)"),
                              paste(ret,collapse=",")));
    ret <- eval(ret)
    return(ret);
}

solve.RxODE <- function(obj,...){
    obj$solve(...);
}

predict.RxODE <- solve.RxODE;

RxODE.inits <- function(vec,names,default = 0){
    ret <- vec;
    nv <- names(vec)
    if (!is.null(nv)){
        ret <- ret[nv %in% names];
        missing <- names[!(names %in%  names(ret))];
        if (is.na(default) & length(missing) > 0){
            stop(sprintf("Missing the following parameter(s): %s.",paste(missing,collapse=", ")))
        }
        if (length(missing) > 0){
            ret[missing] <- default;
            warning(sprintf("Assiged %s to %s.",paste(missing,collapse=", "),default))
        }
        ret <- ret[names];
    }
    return(ret);
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

"summary.RxODE" <-
function(x, ...)
{
    print.RxODE(x);
    modelvars <- x$get.modelVars();
    cat(sprintf("dll: %s\n",x$cmpMgr$dllfile))
    cat(sprintf("Jacobian: %s\n",ifelse(modelvars$jac == "fulluser","Full User Specified","Full Internally Caluclated")))
    if (length(modelvars$params) > 0){
        cat("\nUser Supplied Parameters:\n");
        print(modelvars$params);
    }
    cat("\nCompartents:\n");
    tmp <- modelvars$state;
    names(tmp) <- paste0("cmt=",1:length(tmp));
    print(tmp);
    if (length(modelvars$lhs) > 0){
        cat("\nCalculated Variables:\n");
        print(modelvars$lhs);
    }
    cat("\nModel:\n")
    cat(x$model)
    cat("\n");

    invisible(x)
}

igraph <- function(obj,...){
    ## In case anyone else wants to use the method...
    UseMethod("igraph");
}

nodeInfo <- function(x,       # RxODE normalized model
                     modVars, # model Vars
                     ...){
    ## nodeInfo returns a list containing edgeList, biList and nodes
    
    ##  The nodes object is the compartments in the model that will be
    ##  drawn.  Compartments that should be hidden are prefixed with a
    ##  "."

    ## The edgeList has a list of nodes and the directions they
    ## connect like c("A","B") is an arrow from compartent A to B.

    ## The names of the edgeList is the label applied to the arrow.
    ##
    
    ## Move any of the equals logical operators to ~ operators so they
    ## wont be split...
    mod0 <- gsub(rex::rex(capture(one_of(">!<")),"="),
                 "\\1~",x,perl=TRUE);
    mod0 <- gsub(rex::rex("=="),"~~",mod0,perl=TRUE);
    ## Make sure there are no spaces around remaining = sign
    mod0 <- gsub(rex::rex(any_spaces,"=",any_spaces),"=",mod0,perl=TRUE)
    ## remove the newlines
    mod0 <- gsub(rex::rex(newlines),"",mod0,perl=TRUE)
    ## Remove any if statements
    mod0 <- gsub(rex::rex(any_spaces,"if",any_spaces,
                          one_of("("),except_any_of(")"),one_of(")"),
                          any_spaces,one_of("{"),
                          except_any_of("}"),one_of("}")),
                 "",mod0,perl=TRUE);
    ## Remove any else statements
    mod0 <- gsub(rex::rex(any_spaces,"else",any_spaces,one_of("{"),except_any_of("}"),one_of("}")),
                 "",mod0,perl=TRUE);
    ## Split by semicolon, and =.
    mod0 <- strsplit(strsplit(mod0,rex::rex(some_of(";")))[[1]],
                     rex::rex(one_of("=")));
    mod <- eval(parse(text=sprintf("c(%s)",paste(unlist(lapply(mod0,function(x){
        if (length(x) == 2){
            x1 <- gsub(rex::rex(start,any_spaces,capture(anything),any_spaces,end),
                       "\\1",x[1]);
            x2 <- gsub(rex::rex(start,any_spaces,capture(anything),any_spaces,end),
                       "\\1",x[2]);
            return(sprintf("\"%s\"=\"%s\"",x1,x2));
        } else {
            return(NULL);
        }
    })),collapse=","))))
    
    ## Expand  - abc*() and - ()*abc expressions
    expandFactor <- function(var = "F*KA*depot-C2*(CL+Q)+Q*C3" ){
        ## Expand any simple factored expressions.
        pr <- rex::rex(anything,
                       capture(one_of("-+"),any_spaces,
                               or(
                                   ## - abc * (a+b+c+d+e...)
                                   group(except_some_of("*"),
                                         any_spaces,
                                         one_of("*"),
                                         any_spaces,
                                         one_of("("),
                                         except_any_of("()"),
                                         one_of(")")),
                                   ## - (a+b+c+d+e) * abc
                                   group(one_of("("),
                                         except_any_of("()"),
                                         one_of(")"),
                                         any_spaces,
                                         one_of("*"),
                                         or(except_some_of("*-+;"),
                                            newline)))),
                       anything);
        neg <- rex::rex(start,any_spaces,one_of("-+"))
        ## Add initial + if needed
        if (regexpr(neg,var,perl=TRUE) == -1){
            var <- sprintf("+%s",var);
        }
        while(regexpr(pr,var,perl=TRUE) != -1){
            ## Get expression to expand
            totalExpr <- gsub(pr,"\\1",var,perl=TRUE);
            v <- gsub(rex::rex(any_of("-+"),any_spaces),"",
                      gsub(rex::rex(any_of("*"),
                                    any_spaces,
                                    one_of("("),except_any_of("()"),one_of(")"),
                                    any_spaces,
                                    any_of("*")),
                           "",totalExpr,perl=TRUE),perl=TRUE);
            ## Is this an negative expression?
            negVar <- FALSE
            if (regexpr(neg,totalExpr,perl=TRUE) != -1){
                negVar <- TRUE;
            }
            ## Capture the parenthetical expression
            p <- gsub(rex::rex(anything,
                               one_of("("),capture(except_any_of("()")),one_of(")"),
                               anything),"\\1",totalExpr,perl=TRUE)
            ## Add a positive to the beginning of the expression, if needed.
            if (regexpr(rex::rex(start,any_spaces,one_of("-")),
                       p,perl=TRUE) == -1){
                p <- paste0("+",p);
            }
            ## Protect negative expressions for later use
            p <- gsub(rex::rex(one_of("-")),"~~~~",p,perl=TRUE);
            ## Expand the positive expression
            p <- gsub(rex::rex(one_of("+")),sprintf("%s%s*",ifelse(negVar,"-","+"),v),p,perl=TRUE)
            ## Expand the negative expression
            p <- gsub(rex::rex(n_times(one_of("~"),4)),sprintf("%s%s*",ifelse(negVar,"+","-"),v),p,perl=TRUE)
            ## replace in original
            var <- gsub(RxODE.quoteVar(totalExpr,FALSE),p,var,perl=TRUE);
            print(var);
        }
        ## Replace initial +
        var <- gsub(rex::rex(start,any_spaces,one_of("+"),any_spaces),"",var);
        return(var)
    }
    ## Expand factored expressions
    for (i in 1:length(mod)){
        mod[i] <- expandFactor(mod[i]);
    }
    ## now replace defined variables in overall expressions.
    for (i in 1:length(mod)){
        if (i > 1){
            for (j in seq(1,i-1)){
                mod[i] <- gsub(RxODE.quoteVar(names(mod)[j]),
                               mod[j],mod[i],perl=TRUE);
            }
        }
    }
    ## Get the variables from the parser
    
    ## Subset to the ODE equations.
    de <- c()
    for (v in modVars$state){
        de <- c(de,mod[sprintf("d/dt(%s)",v)])
    }
    w <- which(regexpr(rex::rex(start,any_spaces,one_of("-")),de,perl=TRUE) == -1);
    de[w] <- paste0("+",de[w]);
    names(de) <- modVars$state;

    fullDe <- de;
    
    ## Currently only parses simple expressions.
    ## FIXME: Add michelis-menton / Hill kinetics
    de <- de[which(regexpr("[(]",de) == -1)];
    parDe <- de;

    sortDe <- function(x =c("centr/V2*CL","centr/V2*Q")){
        ## sortDe returns
        ## First protect the / operators.
        x <- unlist(lapply(strsplit(gsub(rex::rex(any_spaces,one_of("/"),any_spaces),
                                  "/zzzzzzz",x[x != ""],perl=TRUE),
                             rex::rex(any_spaces,one_of("/*"),any_spaces),
                             perl=TRUE),
                    function(x){
                        x <- paste(gsub(rex::rex(one_of("*"),n_times("z",7)),"/",paste0("*",sort(x))),collapse="");
                        return(x)
                    }));
        return(x);
    } # end function sortDe
    
    ## Find the negative exressions in the ODEs
    negDe <- gsub(rex::rex(one_of("+"),any_spaces,except_any_of("+-")),"",de,perl=TRUE);
    negDe <- negDe[negDe != ""];
    negDe <- lapply(strsplit(negDe,rex::rex(any_spaces,one_of("-"),any_spaces)),sortDe);
    
    
    ## Find the positive expressions in the ODEs
    de <- gsub(rex::rex(one_of("-"),any_spaces,except_any_of("+-")),"",de,perl=TRUE) 
    nde <- names(de)
    de <- de[de != ""];
    de <- lapply(strsplit(de,rex::rex(any_spaces,one_of("+"),any_spaces)),sortDe);
    for (v in names(negDe)){
        negDe[[v]] <- gsub(rex::rex(any_spaces,one_of("*"),any_spaces,RxODE.quoteVar(v)),"",negDe[[v]])
    }
    ## Look for connections between compartments from the positive direction
    nodes <- c();
    edges <- NULL;
    edgeList <- list();
    biList <- list();
    for (n in names(de)){
        for (n2 in nde){
            reg <- rex::rex(any_spaces,one_of("*"),any_spaces,RxODE.quoteVar(n2))
            w <- which(regexpr(reg,de[[n]],perl=TRUE) != -1);
            if (length(w) == 1){
                nodes <- c(nodes,n2,n);
                var <- gsub(reg,"",de[[n]][w],perl=TRUE);
                for (e in edgeList){
                    if (e[1] == n & e[2] == n2){
                        biList[[length(biList) + 1]] <- c(n,n2);
                    }
                }
                edgeList[[var]] <- c(n2,n);
                ## If the corresponding negative direction exists,
                ## remove it from negDe.  The negDe should have overall body clearances.
                noFvar <- gsub(rex::rex(any_spaces,one_of("*"),any_spaces,one_of("Ff"),any_alnums),"",var,perl=TRUE)
                if (any(names(negDe) == n2)){
                    negDe[[n2]] <- negDe[[n2]][negDe[[n2]] != var];
                    negDe[[n2]] <- negDe[[n2]][negDe[[n2]] != noFvar];

                }
            }
        }
    }
    ## Assign overall body clerances
    nodes <- unique(nodes);
    outn <- 1;
    if (length(negDe) > 0){
        for (i in 1:length(negDe)){
            if (length(negDe[[i]]) == 1){
                outNode <- sprintf(".out%s",outn)
                edgeList[[negDe[[i]]]] <- c(names(negDe)[i],outNode);
                nodes <- c(nodes,outNode);
                outn <- outn +1;
            }
        }
    }
    reg <- rex::rex(capture(anything),any_spaces,one_of("*"),any_spaces,capture(one_of("Ff"),any_alnums),capture(anything));
    w <- which(regexpr(reg,names(edgeList),perl=TRUE) != -1);
    if (length(w) > 0){
        names(edgeList)[w] <- gsub(reg,"\\1\\3\n\\2",names(edgeList)[w],perl=TRUE);
    }
    names(edgeList) <- gsub(rex::rex(any_spaces,one_of("*"),any_spaces),"",names(edgeList),perl=TRUE)
    ## Now see if we can recognize kin/kout from indirect response models
    
    ## FIXME: more robust parsing.  Currently requires kin/kout combination
    ## Also assigns suppression by the presence of "(1-" in either the kin or kout term.
    ## w <- which(sapply(fullDe,function(x){
    ##     lx <- tolower(x)
    ##     return(regexpr("kin",lx)  != -1 && regexpr("kout",lx) != -1)
    ## }))
    ## if (length(w) == 1){
    ##     idr <- fullDe[w];
    ##     idrl <- tolower(idr);
    ##     isKout <- regexpr("^ *kin *[-]",idrl) != -1;
    ##     isInb <- regexpr("[(] *1 *[-]",idrl) != -1;
    ##     idr.from <- gsub(sprintf("^.*(%s).*",paste(names(fullDe)[-w],collapse="|")),"\\1",idr)
    ##     ef <- names(fullDe)[w]
    ##     ## FIXME: Get the actual IC50

    ##     ## FIXME: Indicate if this is a Hill equation, or any other
    ##     ## sort of dose-response equations.
    ##     if (isKout){
    ##         if (isInb){
    ##             Kout <- "Kout\n\u25BC IC50"
    ##             Kin <- "Kin\n";
    ##             idr <- "Indirect\nEffect (II)"
    ##         } else {
    ##             Kout <- "Kout\n\u25B3 EC50";
    ##             Kin <- "Kin\n";
    ##             idr <- "Indirect\nEffect (IV)"
    ##         }
    ##     } else {
    ##         if (isInb){
    ##             Kout <- "Kout"
    ##             Kin <- "Kin\n\u25BC IC50";
    ##             idr <- "Indirect\nEffect (I)"
    ##         } else {
    ##             Kout <- "Kout\n";
    ##             Kin <- "Kin\n\u25B3 EC50";
    ##             idr <- "Indirect\nEffect (III)"
    ##         }
    ##     }
    ##     edgeList[[idr]] <- c(idr.from,ef);
    ##     edgeList[[Kin]] <- c(".Kin",ef);
    ##     edgeList[[Kout]] <- c(ef,".Kout");
    ##     nodes <- c(nodes,ef,".Kin",".Kout");
    ## }
    return(list(nodes    = nodes,
                edgeList = edgeList,
                biList   = biList));
} # end function nodeInfo

igraph.rxDll <- function(x,                                   #  object
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
    ## igraph.rxDll returns igraph object from rxDll
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package igraph needed for this function to work. Please install it.",
             call. = FALSE)
    }
    with(nodeInfo(rxModelVars(m)$model["normModel"],rxModelVars(m)),{
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
} # end function igraph.rxDll

plot.rxDll <- function(x,
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
    invisible(NULL)
} # end function plot.rxDll

"rx.initCmpMgr" <-
    function(model, modName, wd, flat, extra.c,debug)
{
   ## Initialize the RxODE compilation manager (all parsing,
   ## compilation, and loading of dynamic shared objects is 
    ## done through this object. It also keeps tracks of 
   ## filenames for the C file, dll, and the translator 
    ## model file, parameter file, state variables, etc.
    if (is.null(extra.c)){
        .extra.c <- "";
    } else {
        .extra.c <- extra.c;
    }
    .debug <- "";
    if (debug){
        .debug <- " -D__DEBUG__";
    }
    .mod.parsed <- model;
    .mod.md5 <- gsub(rex::rex(any_spaces,one_of("#"),anything,newline),"\n",model,perl=TRUE); # Strip comments
    .mod.md5 <- gsub(rex::rex(any_spaces)," ",.mod.md5,perl=TRUE); # strip spurrios spaecs
    .digest <- digest::digest(sprintf("%s%s",.mod.md5,.debug));
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
   .libs <- 
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
       .dllfile <- file.path(.flatOut,sprintf("%s-%s%s", .modName, .digest, R_ARCH,.Platform$dynlib.ext))
   } else {
       .dllfile <- file.path(.mdir,sprintf("%s-%s%s", .modName, R_ARCH,.Platform$dynlib.ext)) 
   }
   # shell calls to parsing and compiling (via R CMD SHLIB)
    .gflibs <- if (is.win) "-lRblas -lgfortran" else "-lgfortran" # TODO: check do we need this
   #.sh <- if(is.win) "shell" else "system"
   .sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths
   .dvode <- file.path(.mdir, "call_dvode.o");

   ## .parse.dll <- sprintf("%s/tran.dll", .bin);
   .shlib <- 
       sprintf("%s/bin/R CMD SHLIB %s %s", 
               Sys.getenv("R_HOME"), .cfile, .dvode)

    safe_name <- sprintf("RxODE_mod_%s_",gsub("\\W", "_", .modName))  # replace non-alphanumeric by "_"
    .trans <-c();
    
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
      file.copy(file.path(common, "call_dvode.c"),
                file.path(.mdir, "call_dvode.c"));
      ## dyn.load(.parse.dll);
      ## on.exit(dyn.unload(.parse.dll));
      .trans <<- .Call(trans,.modfile, .cfile, .extra.c, safe_name);
      print(.trans);
      .parsed <<- TRUE
      .compiled <<- FALSE
      ## buitin <- c("t", "tlast", "podo")
      invisible(.parsed)
   }

   compile <- function(force = FALSE){
      do.it <- force || !.parsed || !.compiled
      if(!do.it)
         return(invisible(.compiled))
      
      if(!.parsed)
         parse(force = force)

      # may need to unload previous model object code
      if (is.loaded(.trans["ode_solver"])) try(dyn.unload(.dllfile), silent = TRUE)
      
      if (.trans["jac"] == "fulluser"){
          .jac <- " -D__JT__=1 -D__MF__=21";
      } else if (.trans["jac"] == "fullint"){
          .jac <- " -D__JT__=2 -D__MF__=22";
      }
      .extra <- .trans[-(1:2)];
      .extra <- paste(sprintf("-D__%s__=%s",toupper(names(.extra)),.extra),collapse=" ")
      ##on.exit(unlink("Makevars"))
      cat(sprintf("PKG_CPPFLAGS=-I%s%s%s %s\n",.incl,.debug,.jac,.extra), file="Makevars")
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
      if(is.loaded(.trans["ode_solver"])){
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
      if(!is.loaded(.trans["ode_solver"]))
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
        return(file.exists(.dllfile));
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
            trans = .trans,   # name of C function
            modelDir = .mdir,               # model directory
            dllfile = .dllfile,
            get.modelVars = function() .modelVars,   
            isValid = isValid, delete = delete,
            get.index = get.index, 
            getObj = function(obj) get(obj, envir = environment(parse)),
            extra.c = .extra.c
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

rxShortPath <- function(p, # Path
                        ...){
    ## rxShortPath returns a short path name under windows
    return(gsub("\\\\", "/", utils::shortPathName(p)))
} # end function rxShortPath

rxDvode <- file.path(system.file("common", package = "RxODE"), "call_dvode.c");
rxInclude <- system.file("include", package = "RxODE");
rxLibs <- system.file(file.path("libs", .Platform$r_arch), package = "RxODE");
if (.Platform$OS.type=="windows"){
    rxDvode <- rxShortPath(rxDvode);
    rxInclude <- rxShortPath(rxInclude);
    rxLibs <- rxShortPath(rxLibs);
}

rxPrefix <- function(model,          # Model or file name of model
                     modName = NULL, # Model name, overrides calculated model name.
                     ...){
    ## rxPrefix returns a prefix for a model.
    if (!is.null(modName)){
        modelPrefix <- sprintf("%s_",gsub("\\W", "_", modName));
    } else if (file.exists(model)){
        modelPrefix <- sprintf("%s_",gsub("\\W", "_", gsub("[.].*$","",base::basename(model))));
    } else {
        parseModel <- tempfile();
        cFile <- tempfile();
        on.exit({unlink(parseModel); unlink(cFile)});
        sink(parseModel);
        cat(model);
        cat("\n");
        sink();
        trans <- rxTrans(parseModel,cFile)
        modelPrefix <- sprintf("rx_%s_",trans["parsed_md5"]);
    }
    modelPrefix <- sprintf("%s%s_",modelPrefix,.Platform$r_arch);
    return(modelPrefix);
} # end function rxPrefix

rxMd5 <- function(model,         # Model File
                  extraC = NULL, # Extra C
                  ...){
    ## rxMd5 returns MD5 of model file.
    ## digest(file=TRUE) includes file times, so it doesn't work for this needs.
    if (class(model) == "character"){
        if (file.exists(model)){
            ret <- readLines(model);
            mod <- paste(ret,collapse="\n");
        } else {
            stop("Requires model to be a file.");
        }
        if (class(extraC) == "character"){
            if (file.exists(extraC)){
                ret <- c(ret,gsub(rex::rex(or(any_spaces,any_newlines)),"",readLines(extraC),perl=TRUE));
            }
        }
        return(list(text=mod,
                    digest=digest::digest(ret)));
    } else {
        rxModelVars(model)$md5;
    }
} # end function rxMd5

rxTrans <- function(model,
                    cFile       = sprintf("%s.c",gsub("[.][^.]*$","",model)), # C file output
                    extraC      = NULL,                                       # Extra C file(s)
                    modelPrefix = "",                                         # Model Prefix
                    md5         = "",                                         # Md5 of model
                    modName     = NULL,                                       # Model name for dll
                    ...){
    ## rxTrans returns a list of compiled properties
    if (class(model) == "character"){
        if (missing(modelPrefix)){
            modelPrefix <- rxPrefix(model,modName);
        }
        if (file.exists(model)){
            if (missing(md5)){
                md5 <- rxMd5(model,extraC)$digest;
            }
        } else {
            stop("This only translates a file (currently; Try rxCompile).");
        }
        parseModel <- tempfile();
        ret <- .Call(trans, model, cFile, extraC, modelPrefix,md5,parseModel);
        ret <- c(ret,parsed_md5 = rxMd5(parseModel,extraC)$digest);
        on.exit(unlink(parseModel));
        return(ret);
    } else {
        return(rxModelVars(model)$trans);
    }
} # end function rxTrans

rxTransMakevars <- function(rxProps,                                                                              # rxTrans translation properties
                            compileFlags =c("parsed_md5","ode_solver","model_vars","calc_lhs","calc_jac","dydt"), # List of compile flags
                            debug        = FALSE,                                                                 # Debug compile?
                            ...){
    ## rxTransCompileFlags returns a string for the compiler options
    neededProps <- c("jac",compileFlags);
    if (all(neededProps %in% names(rxProps))){
        if (rxProps["jac"] == "fulluser"){
            ret <- " -D__JT__=1 -D__MF__=21";
        } else if (rxProps["jac"] == "fullint"){
            ret <- " -D__JT__=2 -D__MF__=22";
        }
        tmp <- rxProps[compileFlags];
        tmp["parsed_md5_str"] <- sprintf("\"\\\"%s\\\"\"",tmp["parsed_md5"]);
        ret <- paste(c(ret,sprintf("-D__%s__=%s",toupper(names(tmp)),tmp)),collapse=" ");
        if (debug){
            ret <- sprintf("%s -D__DEBUG__",ret);
        }
        ret <- sprintf("PKG_CPPFLAGS=-I%s %s\nPKG_LIBS=-L%s -lodeaux $(BLAS_LIBS) $(FLIBS)",rxInclude,ret,rxLibs);
        cat(ret);
        return(ret);
    } else {
        stop("Cannot figure out what needs to be specified in the compiler.")
    }
} # end function rxTransCompileFlags


rxCompile <-  function(model,           # Model
                       dir,             # Directory
                       prefix,          # Prefix
                       extraC  = NULL,  # Extra C File.
                       force   = FALSE, # Force compile
                       modName = NULL,  # Model Name 
                       ...){
    ## rxCompile returns the dll name that was created.
    if (class(model) == "character"){
        dllCopy <- FALSE;
        if (missing(dir)){
            dir <- tempfile();
            dllCopy <-  TRUE;
            on.exit(unlink(dir, recursive = TRUE))
        }
        if (missing(prefix)){
            prefix <- rxPrefix(model, modName);
        }
        if(!file.exists(dir))
            dir.create(dir, recursive = TRUE)
        cFile <- file.path(dir,sprintf("%s.c",substr(prefix,0,nchar(prefix)-1)));
        cDllFile <- file.path(dir,sprintf("%s%s",substr(prefix,0,nchar(prefix)-1), .Platform$dynlib.ext));
        if (dllCopy){
            finalDll <- file.path(getwd(),basename(cDllFile));
        } else {
            finalDll <-  cDllFile;
        }
        if (!file.exists(model)){
            mFile <- sprintf("%s.rx",substr(cFile,0,nchar(cFile)-2));
            sink(mFile);
            cat(model);
            cat("\n");
            sink();
        } else {
            mFile <- model;
        }
        md5 <- rxMd5(mFile,extraC);
        allModVars <- NULL;
        needCompile <- TRUE
        if (file.exists(finalDll)){
            try(dyn.load(finalDll, local = FALSE), silent=TRUE);
            modVars <- sprintf("%smodel_vars",prefix);
            if (is.loaded(modVars)){
                allModVars <- eval(parse(text=sprintf(".Call(\"%s\")",modVars)),envir=.GlobalEnv)
                modVars <- allModVars$md5;
            }
            if (modVars["file_md5"] == md5$digest){
                needCompile <- FALSE;
            }
        }
        if (force || needCompile){
            dvODE <- file.path(dir,"call_dvode.c");
            dvODEo <- file.path(dir,"call_dvode.o");
            file.copy(rxDvode, dvODE);
            Makevars <- file.path(dir,"Makevars");
            trans <- rxTrans(mFile,cFile = cFile,md5=md5$digest,extraC=extraC,...,modelPrefix=prefix);
            if (file.exists(finalDll)){
                if (modVars["parsed_md5"] == trans["parsed_md5"]){
                    cat("Don't need to recompile, minimal change to model detected.\n");
                    needCompile <- FALSE;
                }
            }
            if (force || needCompile){
                ## Setup Makevars
                sink(Makevars);
                cat(rxTransMakevars(trans));
                sink();
                sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths    
                ## Change working directory
                owd <- getwd();
                on.exit(setwd(owd));
                setwd(dir);
                try(dyn.unload(finalDll), silent=TRUE);
                try(unlink(finalDll));
                cmd <- sprintf("%s/bin/R CMD SHLIB %s %s", 
                               Sys.getenv("R_HOME"), base::basename(cFile),
                               base::basename(dvODEo));
                cat(sprintf("%s\n",cmd));
                rc <- try(do.call(sh, list(cmd)), silent = FALSE)
                if(inherits(rc, "try-error"))
                    stop(sprintf("error compiling %s", cFile));
                if (dllCopy){
                    file.copy(cDllFile,finalDll);
                }
                try(dyn.load(finalDll, local = FALSE), silent=TRUE);
                modVars <- sprintf("%smodel_vars",prefix);
                if (is.loaded(modVars)){
                    allModVars <- eval(parse(text=sprintf(".Call(\"%s\")",modVars)),envir=.GlobalEnv)
                }
            }
        }
        .c <- function(...){return(.C(...))};
        .call <- function(...){return(.Call(...))};
        args <- list(model = model, dir = dir, prefix = prefix,
                     extraC = extraC, force = force, modName = modName,
                     ...);
        ret <- list(dll     = finalDll,
                    model   = md5$text,
                    extra   = extraC,
                    modVars = allModVars,
                    .c      = .c,
                    .call   = .call,
                    args    = args);
        class(ret) <- "rxDll";
        return(ret);
    } else if (class(model) == "rxDll") {
        args <- model$args;
        if (!missing(dir)){
            args$dir <- dir;
        }
        if (!missing(prefix)){
            args$prefix <- prefix;
        }
        if (!missing(extraC)){
            args$extraC <- extraC;
        }
        if (!missing(force)){
            args$force <- force;
        }
        if (!missing(modName)){
            args$modName <- modName;
        }
        ret <- rxCompile(model  = args$model,  dir     = args$dir,
                         prefix = args$prefix, extraC  = args$extra,
                         force  = args$force,  modName = args$modName,
                         ...);
        return(ret);
    }
} # end function rxCompile

rxDllLoaded <- function(x,retry = TRUE){
    m <- rxTrans(x);
    if (any(names(m) == "ode_solver")){
        return(is.loaded(m["ode_solver"]));
    } else if (retry) {
        m <- rxCompile(x,force=FALSE);
        return(rxDllLoaded(m,retry = FALSE))
    } else {
        stop("Can't figure out if the object...");
    }
}

rxDll <- function(obj,...){
    UseMethod("rxDll");
}

rxDll.character <- function(model,...){
    return(rxCompile(model,...))
}

rxDll.rxDll <- function(obj){
    return(obj$dll)
}


rxLoad <- function(obj){
    if (!(rxDllLoaded(obj))){
        dll <- rxDll(obj);
        rc <- try(dyn.load(dll), silent = TRUE)
        if(inherits(rc, "try-error"))
            stop(sprintf("error loading dll file %s", dll));
    }
    return(invisible());
}

rxUnLoad <- function(obj){
    if ((rxDllLoaded(obj))){
        dll <- rxDll(obj);
        rc <- try(dyn.unload(dll), silent = TRUE)
        if(inherits(rc, "try-error"))
            stop(sprintf("error unloading dll file %s", dll));
    }
    return(invisible());
}

rxDelete <- function(obj){
    dll <- rxDll(obj);
    rxUnLoad(obj)
    unlink(dll);
}


rxParams <- function(obj,...){
    return(rxModelVars(obj)$params);
}

rxState <- function(obj,state,...){
    if (missing(state)){
        return(rxModelVars(obj)$state);
    } else {
        objState <- rxState(obj);
        if (length(objState) == 1)
            warning("only one state variable should be input", immediate = TRUE);
        w <- which(objState == state)
        if (length(w) != 1){
            stop(sprintf("Cannot locate state \"%s\"",state));
        }
        return(w);
    }
}

rxLhs <- function(obj,...){
    return(rxModelVars(obj)$lhs);
}

rxModelVars <- function(obj,...){
    UseMethod("rxModelVars");
}

rxModelVars.rxDll <- function(obj,...){
    return(obj$modVars)
}

print.rxDll <- function(x,...){
    if (file.exists(x$dll)){
        cat(sprintf("RxODE dll named \"%s\"",basename(x$dll)));
        if (rxDllLoaded(x)){
            cat(" is loaded and ready to use.\n");
        } else {
            cat(" is not loaded now.\n");
        }
    } else {
        cat(sprintf("RxODE dll named \"%s\" has been deleted.\n",basename(x$dll)));
    }
    invisible(x);
}

summary.rxDll <- function(x,...){
    print(x);
    cat(sprintf("dll: %s\n",rxDll(x)));
    cat(sprintf("Jacobian: %s\n",ifelse(rxModelVars(x)$jac == "fulluser","Full User Specified","Full Internally Caluclated")));
    if (length(rxParams(x)) > 0){
        cat("\nUser Supplied Parameters:\n");
        print(rxParams(x));
    }
    cat("\nCompartents:\n");
    tmp <- rxState(x);
    names(tmp) <- paste0("cmt=",1:length(tmp));
    print(tmp);
    if (length(rxLhs(x)) > 0){
        cat("\nCalculated Variables:\n");
        print(rxLhs(x));
    }
    cat("\nModel:\n")
    cat(rxModelVars(x)$model["model"]);
    cat("\n");

    return(invisible(x))
}

solve.rxDll <- function(rxDllObj,               # rxDll object
                        params,              # Parmaeter 
                        events,              # Events
                        inits       = NULL,  # Initial Events
                        stiff       = TRUE,  # Is the system stiff
                        transit_abs = FALSE, # Transit compartment absorption?
                        atol        = 1.0e-8, # Absoltue Tolerance for LSODA solver
                        rtol        = 1.0e-6, # Relative Tolerance for LSODA solver
                        ...) {
    ## solve.rxDll returns a solved object
    event.table <- events$get.EventTable()
    
    ## preserve input arguments. 
    last.solve.args <-
        list(params = params, events = events$copy(),
             inits = inits, stiff = stiff, 
             transit_abs = transit_abs, atol = atol, rtol = rtol, ...);

    ## check that starting values for all needed parameters are
    ## specified in the input "params"
    if (length(setdiff(rxParams(rxDllObj), names(params)))) {
        msg <- paste("var(s) not found in input pars.\n", 
                     paste(setdiff(rxParams(rxDllObj), names(params)), collapse=" "))
        stop(msg)
    }
    inits <- RxODE.inits(inits,rxState(rxDllObj));
    params <- params[rxParams(rxDllObj)];
    s <- as.list(match.call(expand.dots = TRUE)) 
    wh <- grep(pattern="S\\d+$", names(s))[1]
    ## HACK: fishing scaling variables "S1 S2 S3 ..." from params call
    ## to solve(). Maybe define a "scale=c(central=7.6, ...)" argument
    ## similar to "params="?
    scaler.ix <- 0
    if (!is.na(wh)) {
        if (s[[wh]] %in% names(params)) {
            scaler <- params[s[[wh]]]
            scaler.ix <- as.numeric(substring(names(s)[wh], 2))
        } else {
            warning(paste("scaler variable not found:", s[[wh]]))
        }
    }
    state_vars <- rxState(rxDllObj);
    neq  <- length(state_vars)
    lhs_vars <- rxLhs(rxDllObj);
    nlhs <- length(lhs_vars)

    ntime <- dim(event.table)[1]
    ret   <- rep(0.0, ntime*neq)
    lhs   <- rep(0.0, ntime*nlhs)
    rc <- as.integer(0)  # return code 0 (success) or IDID in call_dvode.c
    if (is.null(inits))
        inits <- rep(0.0, neq)

    ## may need to reload (e.g., when we re-start R and
    ## re-instantiate the RxODE object from a save.image.
    ## cmpMgr$dynLoad()
    rxLoad(rxDllObj);
    xx <- m$.c(rxTrans(rxDllObj)["ode_solver"],
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
               );

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
    
    ret <- cbind(time=event.table$time, x)[events$get.obs.rec(),];
    return(ret)
} # end function solve.rxDll
