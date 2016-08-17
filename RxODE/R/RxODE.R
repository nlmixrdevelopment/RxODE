"RxODE" <-
    function(model, modName = basename(wd), wd = getwd(),
             filename = NULL, do.compile = NULL, flat = FALSE,
             strict = FALSE,## reduce.rounding = TRUE,
             clean.flat.dlls=NULL,
             extra.c=NULL,debug = FALSE,...)
{
    if(!missing(model) && !missing(filename))
        stop("must specify exactly one of 'model' or 'filename'")

    if (!missing(model) & class(extra.c) ==  "character"){
        if (!file.exists(extra.c)){
            extra.c <- FALSE;
        }
    }

   if (missing(modName) && missing(wd) & missing(flat)){
       flat <- TRUE;
       if (!is.null(filename) && missing(model)){
           modName <- basename(filename);
           if (missing(clean.flat.dlls)){
               clean.flat.dlls <- TRUE
           }
       } else if (missing(filename)&& file.exists(model)){
           modName <- basename(model)
           if (missing(clean.flat.dlls)){
               clean.flat.dlls <- TRUE
           }
       } else {
           modName <- "RxODE";
           if (missing(clean.flat.dlls)){
               clean.flat.dlls <- FALSE
           }
       }
   }
    if  (!missing(flat) & flat & !missing(modName)& missing(clean.flat.dlls)){
        clean.flat.dlls <- TRUE;
    }
   
   if(!is.null(filename) && missing(model))
       model <- paste0(readLines(filename, ...), collapse="\n")

   if (missing(filename) && file.exists(model)){
       model <- paste0(readLines(model, ...), collapse="\n")
   }
   
   # RxODE compilation manager (location of parsed code, generated C, 
                                        # shared libs, etc.)

    cmpMgr <- rx.initCmpMgr(model, modName, wd, flat, strict, ## reduce.rounding,
                            clean.flat.dlls,extra.c,debug);
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

RxODE.nodeInfo <- function(x, # RxODE object
                           ...){
    ## RxODE.nodeInfo returns a list containing edgeList, biList and nodes
    
    ##  The nodes object is the compartments in the model that will be
    ##  drawn.  Compartments that should be hidden are prefixed with a
    ##  "."

    ## The edgeList has a list of nodes and the directions they
    ## connect like c("A","B") is an arrow from compartent A to B.

    ## The names of the edgeList is the label applied to the arrow.
    ##

    ## First change the totally parsed  model specifciation to a named R vector.
    ## Take out any comments.
    mod0 <- gsub(rex::rex(any_spaces,"#",anything,newline),
                 "",x$cmpMgr$model,perl=TRUE);
    ## make sure that spaces before and after newlines are chomped out
    ## of model.
    mod0 <- gsub(rex::rex(any_spaces,newline,any_spaces),"\n",mod0,
                 perl=TRUE);
    ## Move any of the equals logical operators to ~ operators so they
    ## wont be split...
    mod0 <- gsub(rex::rex(capture(one_of(">!<")),"="),
                 "\\1~",mod0,perl=TRUE);
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
    modVars <- x$get.modelVars();

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

    sortDe <- function(x=c("centr/V2*CL","centr/V2*Q")){
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
    }
    
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
    function(model, modName, wd, flat, strict, ## reduce.rounding,
             clean.flat.dlls,extra.c,debug)
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
       if (clean.flat.dlls){
           for( f in list.files(.flatOut,sprintf("%s-.*",.modName))){
               if (regexpr(.digest,f) == -1){
                   if (regexpr(.Platform$dynlib.ext,f) != -1){
                       try(dyn.unload(f), silent = TRUE)
                   }
                   unlink(f)
               }
           }
       }
   } else {
       .dllfile <- file.path(.mdir,sprintf("%s-%s%s", .modName, R_ARCH,.Platform$dynlib.ext)) 
   }
   .parfile <- file.path(.mdir, "ODE_PARS.txt")
   .stvfile <- file.path(.mdir, "STATE_VARS.txt")
   .lhsfile <- file.path(.mdir, "LHS_VARS.txt")
   .jacfile <- file.path(.mdir, "JAC_TYPE.txt")
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
       sprintf("%s/bin/R CMD SHLIB %s %s", 
         Sys.getenv("R_HOME"), .cfile, .dvode)

    .dydt <- .calc_lhs <- .jac <- .ode_solver <- NULL;
    safe_name <- paste0(gsub("\\W", "_", .modName),.digest)  # replace non-alphanumeric by "_"
    .dydt <- paste0("RxODE_mod_", safe_name, "_dydt")
    .calc_lhs <- paste0("RxODE_mod_", safe_name, "_calc_lhs")
    .jac <- paste0("RxODE_mod_", safe_name, "_calc_jac")
    .ode_solver <- paste0("RxODE_mod_", safe_name, "_ode_solver")
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
      
      src <- readLines(file.path(common, "call_dvode.c"))
      src <- gsub("dydt", .dydt, src)
      src <- gsub("calc_lhs", .calc_lhs, src)
      src <- gsub("calc_jac", .jac, src);
      src <- gsub("ode_solver", .ode_solver, src);
      
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
      src <- gsub("calc_jac", .jac, src)
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
            lhs = scan(.lhsfile, what = "", quiet = TRUE),
            jac = scan(.jacfile, what = "", quiet = TRUE)
         )
      save(.modelVars,.parsed,.dydt,.objName,file=.modelVarsFile)
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

      if (.modelVars$jac == "fulluser"){
          .jac <- " -D__JT__=1 -D__MF__=21";
      } else if (.modelVars$jac == "fullint"){
          .jac <- " -D__JT__=2 -D__MF__=22";
      }
      #on.exit(unlink("Makevars"))
      cat(sprintf("PKG_CPPFLAGS=-I%s%s%s\n",.incl,.debug,.jac), file="Makevars")
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

