rex::register_shortcuts("RxODE");
regIni <- rex::rex(or(group(one_of("_."), "0"), "0", "(0)", "[0]", "{0}"), end);
RxODE.quoteVar <- function(var, word = TRUE){
    ret <- c();
    for (i in 1:nchar(var)){
        ret <- c(ret, sprintf("one_of(\"%s\")", substr(var, i, i)));
    }
    ret <- parse(text = sprintf(ifelse(word, "rex::rex(boundary, %s, boundary)",
                                       "rex::rex(%s)"),
                                paste(ret, collapse = ", ")));
    ret <- eval(ret)
    return(ret);
}
igraph <- function(obj, ...){
    ## In case anyone else wants to use the method...
    UseMethod("igraph");
}

rxSub <- function(regexp, # regular expression with named groups
                  what,   # replacement with \g{named_group} replacement
                  text,   # Text to replace
                  ...){
    ## rxSub returns a replaced character string, like gsub with named groups
    reg <- regexpr(regexp, text, perl = TRUE);
    if (reg != -1){
        n <- apply(rbind(attr(reg, "capture.start"), attr(reg, "capture.length")), 2, function(x){
            return(substr(text, x[1], x[1]+x[2]-1))
        })
        ## Support \g{} backreference in replacement
        names(n) <- paste0("\\g{", attr(reg, "capture.names"), "}");
        m <- what;
        for (i in 1:length(n)){
            m <- gsub(names(n)[i], n[i], m, fixed = TRUE);
        }
        text <- paste0(substr(text, 0, reg-1), m,
                       substr(text, reg+attr(reg, "match.length"),
                              nchar(text)));
        return(text);
    } else{
        return(text);
    }
} # end function rxSub

##' IDR information parsing routine
##' @title IDR information
##' @param cmt Compartment that idr information is being assesed on.
##' @param text Text that is parsed for indirect response information.
##' @return IDR information (if any)
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
idrInfo <- function(cmt = "eff", text = "+Kin-Kout*(1-centr/V2/(EC50+centr/V2))*eff"){
    text <- rxSigmoidInfo(text);
    id <- rex::rex(one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9"));
    resp <- rex::rex("Hill(", capture(except_some_of(", ", newline), name = "emax"), ## Emax
                     ", ", capture(except_some_of(",",newline), name = "cp"),## Cp
                     ",",capture(except_some_of(",",newline), name = "e50"),## e50
                     ",",capture(except_some_of(",",newline), name = "gamma"),## gamma
                     ")")
    kout <- rex::rex(capture(id, name = "kout"), any_spaces);
    kin <- rex::rex(capture(id, name = "kin"), any_spaces);
    R <- rex::rex(capture(cmt, name = "R"), any_spaces);
    inb <- rex::rex("(", any_spaces, or("1", "1.0"), any_spaces, "-", any_spaces, resp, any_spaces, ")", any_spaces);
    stim <- rex::rex("(", any_spaces, or("1", "1.0"), any_spaces, "+", any_spaces, resp, any_spaces, ")", any_spaces);
    stim2 <- rex::rex("(", any_spaces, resp, any_spaces, "+", any_spaces, or("1", "1.0"), any_spaces, ")", any_spaces);
    regs <- list(
        ## IDR1 kin*(1-sig)-kout*R
        c(rex::rex(kin, "*", any_spaces, inb, "-", any_spaces, R, "*", any_spaces, kout),
          "IDR1"),
        c(rex::rex(kin, "*", any_spaces, inb, "-", any_spaces, kout, "*", any_spaces, R),
          "IDR1"),
        c(rex::rex(inb, "*", any_spaces, kin, "-", any_spaces, R, "*", any_spaces, kout),
          "IDR1"),
        c(rex::rex(inb, "*", any_spaces, kin, "-", any_spaces, kout, "*", any_spaces, R),
          "IDR1"),
        ## IDR2
        c(rex::rex(kin, "-", any_spaces, kout, "*", any_spaces, inb, "*", any_spaces, R),
          "IDR2"),
        c(rex::rex(kin, "-", any_spaces, kout, "*", any_spaces, R, "*", any_spaces, inb),
          "IDR2"),
        c(rex::rex(kin, "-", any_spaces, inb, "*", any_spaces, kout, "*", any_spaces, R),
          "IDR2"),
        c(rex::rex(kin, "-", any_spaces, inb, "*", any_spaces, R, "*", any_spaces, kout),
          "IDR2"),
        c(rex::rex(kin, "-", any_spaces, R, "*", any_spaces, kout, "*", any_spaces, inb),
          "IDR2"),
        c(rex::rex(kin, "-", any_spaces, R, "*", any_spaces, inb, "*", any_spaces, kout),
          "IDR2"),
        ## IDR3
        c(rex::rex(kin, "*", any_spaces, stim, "-", any_spaces, kout, "*", any_spaces, R),
          "IDR3"),
        c(rex::rex(stim, "*", any_spaces, kin, "-", any_spaces, kout, "*", any_spaces, R),
          "IDR3"),
        c(rex::rex(kin, "*", any_spaces, stim, "-", any_spaces, R, "*", any_spaces, kout),
          "IDR3"),
        c(rex::rex(stim, "*", any_spaces, kin, "-", any_spaces, R, "*", any_spaces, kout),
          "IDR3"),
        c(rex::rex(kin, "*", any_spaces, stim2, "-", any_spaces, kout, "*", any_spaces, R),
          "IDR3"),
        c(rex::rex(stim2, "*", any_spaces, kin, "-", any_spaces, kout, "*", any_spaces, R),
          "IDR3"),
        c(rex::rex(kin, "*", any_spaces, stim2, "-", any_spaces, R, "*", any_spaces, kout),
          "IDR3"),
        c(rex::rex(stim2, "*", any_spaces, kin, "-", any_spaces, R, "*", any_spaces, kout),
          "IDR3"),
        ## IDR4
        c(rex::rex(kin, "-", any_spaces, kout, "*", any_spaces, stim, "*", any_spaces, R),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, kout, "*", any_spaces, R, "*", any_spaces, stim),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, stim, "*", any_spaces, kout, "*", any_spaces, R),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, stim, "*", any_spaces, R, "*", any_spaces, kout),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, R, "*", any_spaces, stim, "*", any_spaces, kout),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, R, "*", any_spaces, kout, "*", any_spaces, stim),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, kout, "*", any_spaces, stim2, "*", any_spaces, R),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, kout, "*", any_spaces, R, "*", any_spaces, stim2),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, stim2, "*", any_spaces, kout, "*", any_spaces, R),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, stim2, "*", any_spaces, R, "*", any_spaces, kout),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, R, "*", any_spaces, stim2, "*", any_spaces, kout),
          "IDR4"),
        c(rex::rex(kin, "-", any_spaces, R, "*", any_spaces, kout, "*", any_spaces, stim2),
          "IDR4")
    );
    ret <- text
    for (i in 1:length(regs)){
        ret <- rxSub(regs[[i]][1], sprintf("%s(\"\\g{kin}\", \"\\g{kout}\", \\g{e50}, \"\\g{R}\", \\g{cp}, edgeList, nodes)", regs[[i]][2]), ret);
    }
    ret <- gsub(rex::rex(start, "+"), "", ret);
    return(ret);
}

##' Internal function to detect sigmoid functions
##'
##' Used in plotting to figure out if there is a sigmoid function
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
rxSigmoidInfo <- function(text = "+Kin-Kout*(1-centr/V2/(EC50+centr/V2))*eff"){
    ## Emax models
    id <- rex::rex(one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9"));

    cp <- rex::rex(or(id, ## Cp
                      group(id, any_spaces, "/", any_spaces, id), ## AMT/V
                      group("(", any_spaces, id, any_spaces, "/", any_spaces, id, any_spaces, ")") ## (AMT/V)
                      ))
    cpg <- function(first = FALSE, paren = FALSE){
        id <- rex::rex(one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9"));
        cp <- rex::rex(or(id, ## Cp
                          group("(", any_spaces, id, any_spaces, ")"),
                          group(id, any_spaces, "/", any_spaces, id), ## AMT/V
                          group("(", any_spaces, id, any_spaces, "/", any_spaces, id, any_spaces, ")") ## (AMT/V)
                          ));
        if (first){
            g <- rex::rex(capture(id, name = "gamma"));
            cp <- rex::rex(capture(cp, name = "cp"));
        } else {
            g <- rex::rex(capture_group("gamma"));
            cp <- rex::rex(capture_group("cp"));
        }
        if (paren){
            return(rex::rex(cp, any_spaces, or("^", "**"), any_spaces, "(", any_spaces, g, any_spaces, ")"))
        } else {
            return(rex::rex(cp, any_spaces, or("^", "**"), any_spaces, g))
        }
    }

    regs <- list(
        ## Emax*C/(C+E50)
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## Emax*C/(E50+C)
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## C*Emax/(C+E50)
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*Emax/(E50+C)
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C/(C+E50)*Emax
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C/(E50+C)*Emax
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## (C/(C+E50))*Emax
        c(rex::rex("(", any_spaces, capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## (C/(E50+C))*Emax
        c(rex::rex("(", any_spaces, capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C/(C+E50)
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")"), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C/(E50+C)
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")"), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*(E50+c)^-1 variants
        ## Emax*C*(C+E50)^-1
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## Emax*C*(E50+C)^-1
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## C*Emax*(C+E50)^-1
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*Emax*(E50+C)^-1
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*(C+E50)^-1*Emax
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")",
                   any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"),
                   any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*(E50+C)^-1*Emax
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")",
                   any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"), any_spaces, "*",
                   any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## (C*(C+E50))*Emax
        c(rex::rex("(", any_spaces, capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, ")",
                   any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"),
                   any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## (C*(E50+C))*Emax
        c(rex::rex("(", any_spaces, capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"),
                   any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*(C+E50)
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture_group("cp"), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## C*(E50+C)
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        c(rex::rex(capture(cp, name = "cp"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "+", any_spaces, capture_group("cp"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ##  Emax/(1+EC50/Cp) variants
        c(rex::rex(capture(id, name = "emax"), any_spaces, "/", any_spaces,
                   "(", any_spaces, or("1", "1.0"), any_spaces, "+", any_spaces, capture(id, name = "e50"), any_spaces, "/", any_spaces,
                   capture(cp, name = "cp"), any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        c(rex::rex(capture(id, name = "emax"), any_spaces, "/", any_spaces,
                   "(", any_spaces, or("1", "1.0"), any_spaces, "+", any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "/", any_spaces,
                   capture(cp, name = "cp"), any_spaces, ")", any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        c(rex::rex(capture(id, name = "emax"), any_spaces, "/", any_spaces,
                   "(", any_spaces, capture(id, name = "e50"), any_spaces, "/", any_spaces,
                   capture(cp, name = "cp"), any_spaces, "+", any_spaces, or("1", "1.0"), any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        c(rex::rex(capture(id, name = "emax"), any_spaces, "/", any_spaces,
                   "(", any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, "/", any_spaces,
                   capture(cp, name = "cp"), any_spaces, ")", any_spaces, "+",
                   any_spaces, or("1", "1.0"), any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),
        ## FIXME: -- log linerized forms?
        ## FIXME: Hill Variants
        ## Emax*C^g/(C^g+E50^g)
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## Emax*C^g/(E50^g+C^g)
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),

        ## C^g*Emax/(C^g+E50^g)
        c(rex::rex(cpg(TRUE), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "/",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g*Emax/(E50^g+C^g)
        c(rex::rex(cpg(TRUE), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")"), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g/(C^g+E50^g)*Emax
        c(rex::rex(cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g/(E50^g+C^g)*Emax
        c(rex::rex(cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")),
          "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),

        ## (C^g/(C^g+E50^g))*Emax
        c(rex::rex("(", any_spaces, cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")", any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")),
          "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),

        ## (C^g/(E50^g+C^g))*Emax
        c(rex::rex("(", any_spaces, cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"),
                   any_spaces, "+", any_spaces, cpg(),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")", any_spaces, ")", any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g/(C^g+E50^g)
        c(rex::rex(cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")"), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g/(E50^g+C^g)
        c(rex::rex(cpg(TRUE), any_spaces, "/",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"),
                   any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")"), "Hill(1, \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## Emax*C^g*(C^g+E50^g)^-1
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, cpg(TRUE), any_spaces, "*",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"1\")"),

        ## Emax*C^g*(E50^g+C^g)^-1
        c(rex::rex(capture(id, name = "emax"), any_spaces, "*", any_spaces, cpg(TRUE), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, or("^", "**"), any_spaces, capture_group("gamma"), any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),

        ## C^g*Emax*(C^g+E50^g)^-1
        c(rex::rex(cpg(TRUE), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "*",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"), any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g*Emax*(E50^g+C^g)^-1
        c(rex::rex(cpg(TRUE), any_spaces, "*", any_spaces, capture(id, name = "emax"), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, or("^", "**"), any_spaces, capture_group("gamma"), any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")", any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g*(C^g+E50^g)^-1*Emax
        c(rex::rex(cpg(TRUE), any_spaces, "*",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"), any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")",
                   any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"),
                   any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),
        ## C^g*(E50^g+C^g)^-1*Emax
        c(rex::rex(cpg(TRUE), any_spaces, "*",
                   any_spaces, "(", any_spaces, capture(id, name = "e50"), any_spaces, or("^", "**"), any_spaces, capture_group("gamma"), any_spaces, "+", any_spaces, cpg(),
                   any_spaces, ")",
                   any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"), any_spaces, "*",
                   any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")"),

        ## (C^g*(C^g+E50^g))*Emax
        c(rex::rex("(", any_spaces, cpg(TRUE), any_spaces, "*",
                   any_spaces, "(", any_spaces, cpg(), any_spaces, "+", any_spaces, capture(id, name = "e50"), any_spaces, or("^", "**"), any_spaces, capture_group("gamma"),
                   any_spaces, ")", any_spaces, ")",
                   any_spaces, or("**", "^"), any_spaces, or("-1", "(-1)", "-1.0", "(-1.0)"),
                   any_spaces, "*", any_spaces, capture(id, name = "emax")), "Hill(\"\\g{emax}\", \"\\g{cp}\", \"\\g{e50}\", \"\\g{gamma}\")")
        ##

        ## FIXME: Hodgkin?
        ## FIXME: Douglas?
        ## FIXME: Gompertz?
    );
    ret <- text
    for (i in 1:length(regs)){
        ret <- rxSub(regs[[i]][1], regs[[i]][2], ret);
    }
    return(ret);
}

nodeInfo <- function(x,       # RxODE normalized model
                     modVars, # model Vars
                     ...){
    ## nodeInfo returns a list containing edgeList, biList and nodes

    ##  The nodes object is the compartments in the model that will be
    ##  drawn.  Compartments that should be hidden are prefixed with a
    ##  "."

    ## The edgeList has a list of nodes and the directions they
    ## connect like c("A", "B") is an arrow from compartent A to B.

    ## The names of the edgeList is the label applied to the arrow.
    ##

    ## Move any of the equals logical operators to ~ operators so they
    ## wont be split...
    mod0 <- gsub(rex::rex(capture(one_of(">!<")), " = "),
                 "\\1~", x, perl = TRUE);
    mod0 <- gsub(rex::rex(" = ="), "~~", mod0, perl = TRUE);
    ## Make sure there are no spaces around remaining = sign
    mod0 <- gsub(rex::rex(any_spaces, " = ", any_spaces), " = ", mod0, perl = TRUE)
    ## remove the newlines
    mod0 <- gsub(rex::rex(newlines), "", mod0, perl = TRUE)
    ## Remove any if statements
    mod0 <- gsub(rex::rex(any_spaces, "if", any_spaces,
                          one_of("("), except_any_of(")"), one_of(")"),
                          any_spaces, one_of("{"),
                          except_any_of("}"), one_of("}")),
                 "", mod0, perl = TRUE);
    ## Remove any else statements
    mod0 <- gsub(rex::rex(any_spaces, "else", any_spaces, one_of("{"), except_any_of("}"), one_of("}")),
                 "", mod0, perl = TRUE);
    ## Split by semicolon, and =.
    mod0 <- strsplit(strsplit(mod0, rex::rex(some_of(";")))[[1]],
                     rex::rex(one_of(" = ")));
    mod <- eval(parse(text = sprintf("c(%s)", paste(unlist(lapply(mod0, function(x){
                                                  if (length(x) == 2){
                                                      x1 <- gsub(rex::rex(start, any_spaces, capture(anything), any_spaces, end),
                                                                 "\\1", x[1]);
                                                      x2 <- gsub(rex::rex(start, any_spaces, capture(anything), any_spaces, end),
                                                                 "\\1", x[2]);
                                                      return(sprintf("\"%s\" = \"%s\"", x1, x2));
                                                  } else {
                                                      return(NULL);
                                                  }
                                              })), collapse = ", "))))

    ## Expand  - abc*() and - ()*abc expressions
    expandFactor <- function(var = "F*KA*depot-C2*(CL+Q)+Q*C3" ){
        ## Expand any simple factored expressions.
        pr <- rex::rex(anything,
                       capture(one_of("-+"), any_spaces,
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
        neg <- rex::rex(start, any_spaces, one_of("-+"))
        ## Add initial + if needed
        if (regexpr(neg, var, perl = TRUE) == -1){
            var <- sprintf("+%s", var);
        }
        while(regexpr(pr, var, perl = TRUE) != -1){
            ## Get expression to expand
            totalExpr <- gsub(pr, "\\1", var, perl = TRUE);
            v <- gsub(rex::rex(any_of("-+"), any_spaces), "",
                      gsub(rex::rex(any_of("*"),
                                    any_spaces,
                                    one_of("("), except_any_of("()"), one_of(")"),
                                    any_spaces,
                                    any_of("*")),
                           "", totalExpr, perl = TRUE), perl = TRUE);
            ## Is this an negative expression?
            negVar <- FALSE
            if (regexpr(neg, totalExpr, perl = TRUE) != -1){
                negVar <- TRUE;
            }
            ## Capture the parenthetical expression
            p <- gsub(rex::rex(anything,
                               one_of("("), capture(except_any_of("()")), one_of(")"),
                               anything), "\\1", totalExpr, perl = TRUE)
            ## Add a positive to the beginning of the expression, if needed.
            if (regexpr(rex::rex(start, any_spaces, one_of("-")),
                        p, perl = TRUE) == -1){
                p <- paste0("+", p);
            }
            ## Protect negative expressions for later use
            p <- gsub(rex::rex(one_of("-")), "~~~~", p, perl = TRUE);
            ## Expand the positive expression
            p <- gsub(rex::rex(one_of("+")), sprintf("%s%s*", ifelse(negVar, "-", "+"), v), p, perl = TRUE)
            ## Expand the negative expression
            p <- gsub(rex::rex(n_times(one_of("~"), 4)), sprintf("%s%s*", ifelse(negVar, "+", "-"), v), p, perl = TRUE)
            ## replace in original
            var <- gsub(RxODE.quoteVar(totalExpr, FALSE), p, var, perl = TRUE);
        }
        ## Replace initial +
        var <- gsub(rex::rex(start, any_spaces, one_of("+"), any_spaces), "", var);
        return(var)
    }
    ## Expand factored expressions
    for (i in 1:length(mod)){
        mod[i] <- expandFactor(mod[i]);
    }
    print(mod);
    ## now replace defined variables in overall expressions.
    for (i in 1:length(mod)){
        if (i > 1){
            for (j in seq(1, i-1)){
                mod[i] <- gsub(RxODE.quoteVar(names(mod)[j]),
                               mod[j], mod[i], perl = TRUE);
            }
        }
    }
    ## Get the variables from the parser

    ## Subset to the ODE equations.
    de <- c()
    for (v in modVars$state){
        de <- c(de, mod[sprintf("d/dt(%s)", v)])
    }
    w <- which(regexpr(rex::rex(start, any_spaces, one_of("-")), de, perl = TRUE) == -1);
    de[w] <- paste0("+", de[w]);
    names(de) <- modVars$state;

    fullDe <- de;

    ## Currently only parses simple expressions.
    ## FIXME: Add michelis-menton / Hill kinetics
    de <- de[which(regexpr("[(]", de) == -1)];
    parDe <- de;

    sortDe <- function(x =c("centr/V2*CL", "centr/V2*Q")){
        ## sortDe returns
        ## First protect the / operators.
        x <- unlist(lapply(strsplit(gsub(rex::rex(any_spaces, one_of("/"), any_spaces),
                                         "/zzzzzzz", x[x != ""], perl = TRUE),
                                    rex::rex(any_spaces, one_of("/*"), any_spaces),
                                    perl = TRUE),
                           function(x){
            x <- paste(gsub(rex::rex(one_of("*"), n_times("z", 7)), "/", paste0("*", sort(x))), collapse = "");
            return(x)
        }));
        return(x);
    } # end function sortDe

    ## Find the negative exressions in the ODEs
    negDe <- gsub(rex::rex(one_of("+"), any_spaces, except_any_of("+-")), "", de, perl = TRUE);
    negDe <- negDe[negDe != ""];
    negDe <- lapply(strsplit(negDe, rex::rex(any_spaces, one_of("-"), any_spaces)), sortDe);


    ## Find the positive expressions in the ODEs
    de <- gsub(rex::rex(one_of("-"), any_spaces, except_any_of("+-")), "", de, perl = TRUE)
    nde <- names(de)
    de <- de[de != ""];
    de <- lapply(strsplit(de, rex::rex(any_spaces, one_of("+"), any_spaces)), sortDe);
    for (v in names(negDe)){
        negDe[[v]] <- gsub(rex::rex(any_spaces, one_of("*"), any_spaces, RxODE.quoteVar(v)), "", negDe[[v]])
    }
    ## Look for connections between compartments from the positive direction
    nodes <- c();
    edges <- NULL;
    edgeList <- list();
    biList <- list();
    for (n in names(de)){
        for (n2 in nde){
            reg <- rex::rex(any_spaces, one_of("*"), any_spaces, RxODE.quoteVar(n2))
            w <- which(regexpr(reg, de[[n]], perl = TRUE) != -1);
            if (length(w) == 1){
                nodes <- c(nodes, n2, n);
                var <- gsub(reg, "", de[[n]][w], perl = TRUE);
                for (e in edgeList){
                    if (e[1] == n & e[2] == n2){
                        biList[[length(biList) + 1]] <- c(n, n2);
                    }
                }
                edgeList[[var]] <- c(n2, n);
                ## If the corresponding negative direction exists,
                ## remove it from negDe.  The negDe should have overall body clearances.
                noFvar <- gsub(rex::rex(any_spaces, one_of("*"), any_spaces, one_of("Ff"), any_alnums), "", var, perl = TRUE)
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
                outNode <- sprintf(".out%s", outn)
                edgeList[[negDe[[i]]]] <- c(names(negDe)[i], outNode);
                nodes <- c(nodes, outNode);
                outn <- outn +1;
            }
        }
    }
    reg <- rex::rex(capture(anything), any_spaces, one_of("*"), any_spaces, capture(one_of("Ff"), any_alnums), capture(anything));
    w <- which(regexpr(reg, names(edgeList), perl = TRUE) != -1);
    if (length(w) > 0){
        names(edgeList)[w] <- gsub(reg, "\\1\\3\n\\2", names(edgeList)[w], perl = TRUE);
    }
    names(edgeList) <- gsub(rex::rex(any_spaces, one_of("*"), any_spaces), "", names(edgeList), perl = TRUE)
    ## Now see if we can recognize kin/kout from indirect response models
    idr <- sapply(1:length(fullDe), function(x){
        tmp <- idrInfo(names(fullDe)[x], fullDe[x]);
        if (regexpr(rex::rex(start, "IDR"), tmp) != -1){
            return(tmp)
        } else {
            return("")
        }
    })
    names(idr) <- names(fullDe);
    idr <- idr[idr != ""];
    w <- which(names(fullDe) %in% names(idr));
    fReg <- sprintf("^.*(%s).*$", paste(names(fullDe)[-w], collapse = "|"));
    IDR1 <- function(kin, kout, e50, ef, cp, edgeList, nodes){
        kout <- sprintf("%s\n", kout);
        kin  <- sprintf("%s\n\u25BC %s", kin, e50);
        idr <- "Indirect\nEffect (I)";
        idrFrom <- gsub(fReg, "\\1", cp);
        edgeList[[idr]] <- c(idrFrom, ef);
        edgeList[[kin]] <- c(".Kin", ef);
        edgeList[[kout]] <- c(ef, ".Kout");
        nodes <- c(nodes, ef, ".Kin", ".Kout");
        return(list(edgeList, nodes))
    }
    IDR2 <- function(kin, kout, e50, ef, cp, edgeList, nodes){
        kout <- sprintf("%s\n\u25BC %s", kout, e50);
        kin  <- sprintf("%s\n", kin);
        idr <- "Indirect\nEffect (II)";
        idrFrom <- gsub(fReg, "\\1", cp);
        edgeList[[idr]] <- c(idrFrom, ef);
        edgeList[[kin]] <- c(".Kin", ef);
        edgeList[[kout]] <- c(ef, ".Kout");
        nodes <- c(nodes, ef, ".Kin", ".Kout");
        return(list(edgeList, nodes))
    }
    IDR3 <- function(kin, kout, e50, ef, cp, edgeList, nodes){
        kout <- sprintf("%s\n", kout);
        kin  <- sprintf("%s\n\u25B3 %s", kin, e50);
        idr <- "Indirect\nEffect (III)";
        idrFrom <- gsub(fReg, "\\1", cp);
        edgeList[[idr]] <- c(idrFrom, ef);
        edgeList[[kin]] <- c(".Kin", ef);
        edgeList[[kout]] <- c(ef, ".Kout");
        nodes <- c(nodes, ef, ".Kin", ".Kout");
        return(list(edgeList, nodes))
    }
    IDR4 <- function(kin, kout, e50, ef, cp, edgeList, nodes){
        kout <- sprintf("%s\n\u25B3 %s", kout, e50);
        kin  <- sprintf("%s\n", kin);
        idr <- "Indirect\nEffect (IV)";
        idrFrom <- gsub(fReg, "\\1", cp);
        edgeList[[idr]] <- c(idrFrom, ef);
        edgeList[[kin]] <- c(".Kin", ef);
        edgeList[[kout]] <- c(ef, ".Kout");
        nodes <- c(nodes, ef, ".Kin", ".Kout");
        return(list(edgeList, nodes))
    }
    if (length(idr) > 0){
        tmp <- eval(parse(text = idr));
        edgeList <- tmp[[1]];
        nodes <- tmp[[2]];
    }
    return(list(nodes    = nodes,
                edgeList = edgeList,
                biList   = biList));
} # end function nodeInfo

##' @export
igraph.rxDll <- function(x,                                   #  object
                         shape      = c("square", "circle", "csquare", "rectangle", "crectangle", "vrectangle", "sphere", "none"),
                         size       = 30,                     # Size of square
                         colors     = c("accent1" = "#0460A9"), # Colors
                         fillColor  = "accent1",              # Color to fill
                         family     = "sans",                 # font Family
                         font       = 2,                      # Font (1: plain, 2: bold, 3: italic, 4: bold/italic, 5:symbol)
                         labelColor = "white",                # Label Color
                         lineColor  = "accent1",              # Line Color
                         shapeEnd   = c("none", "sphere", "circle", "square", "csquare", "rectangle", "crectangle", "vrectangle"),
                         sizeEnd    = 10,                     # Size of End
                         tk         = FALSE,                  # For tkplot; transparancy not supported.
                         ...){
    ## igraph.rxDll returns igraph object from rxDll
    if (!requireNamespace("igraph", quietly = TRUE)) {  # nocov start
        stop("Package igraph needed for this function to work. Please install it.",
             call. = FALSE)
    }  # nocov end
    with(nodeInfo(rxModelVars(x)$model["normModel"], rxModelVars(x)), {
        ret <- eval(parse(text = sprintf("igraph::graph_from_literal(%s);", paste(unlist(lapply(edgeList, function(x){sprintf("\"%s\" -+ \"%s\"", x[1], x[2])})), collapse = ", "))));
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
        ret <- igraph::set_vertex_attr(ret, "shape", value = shape);
        ret <- igraph::set_vertex_attr(ret, "size", value = size);
        ret <- igraph::set_vertex_attr(ret, "color", value = getColor(fillColor));
        ret <- igraph::set_vertex_attr(ret, "label.family", value = family)
        ret <- igraph::set_vertex_attr(ret, "label.font", value = font);
        ret <- igraph::set_vertex_attr(ret, "label.color", value = getColor(labelColor));
        for (n in nodes){
            if (substring(n, 0, 1) == "."){
                if (!tk){
                    ret <- igraph::set_vertex_attr(ret, "label.color", igraph::V(ret)[n], value = "transparent");
                    ret <- igraph::set_vertex_attr(ret, "shape", igraph::V(ret)[n], value = shapeEnd);
                    ret <- igraph::set_vertex_attr(ret, "size", igraph::V(ret)[n], value = sizeEnd);
                } else { # nocov start
                    ret <- igraph::set_vertex_attr(ret, "size", igraph::V(ret)[n], value = sizeEnd);
                } # nocov end
            }
        }
        ret <- igraph::set_edge_attr(ret, "color", value = getColor(lineColor))
        ret <- igraph::set_edge_attr(ret, "curved", value = FALSE)
        ret <- igraph::set_edge_attr(ret, "label.font", value = font)
        ## Set bidirectional to curved.
        eval(parse(text = paste(lapply(biList, function(x){sprintf("ret <- igraph::set_edge_attr(ret, \"curved\", igraph::E(ret)[\"%s\" %%--%% \"%s\"], value = TRUE)", x[1], x[2])}), collapse = ";")));
        ## add the labels
        for (l in names(edgeList)){
            v <- edgeList[[l]];
            eg <- igraph::E(ret)[v[1] %->% v[2]];
            ret <- igraph::set_edge_attr(ret, "label", eg, value = l);
        }
        return(ret);
    });
} # end function igraph.rxDll


##' Plot a model digram for simple ODE models
##'
##' This plots a line diagram for the current system of differential
##' equations.
##'
##' @param RxODEobj RxODE object
##' @param family Font family to use
##' @param interactive If true, use Tcl/tk to plot the diagram and
##'     allow movement of the compartments
##' @param shape Box shape, by default it is a square
##' @param size size of square for compartments
##' @param colors is a named list of colors.  If any of the names match
##'     the colors specified by this function, it will use the color in
##'     this list.  For example \code{accent1} would match the default
##'     color in this list.
##' @param fillColor Color of the fill
##' @param font Font (1: plain, 2: bold, 3: italic, 4: bold/italic,
##'     5:symbol)
##' @param labelColor Color of the label
##' @param lineColor Color of the line
##' @param shapeEnd End shape
##' @param sizeEnd Size of end
##' @param ... Other arguments passed to igraph's plot method.
##'
##' Currently this method supports simple elimination type models and
##' indirect-effect models.  Michelis Menton models are not yet
##' supported.
##'
##' @author Matthew L.Fidler
##' @export
rxPlot <- function(RxODEobj,
                   family      ="sans",
                   interactive = FALSE,
                   shape       = c("square","circle","csquare","rectangle","crectangle","vrectangle","sphere","none"),
                   size        = 30,                     # Size of square
                   colors      = c("accent1"="#0460A9"), # Colors
                   fillColor   = "accent1",              # Color to fill
                   font        = 2,                      # Font (1: plain, 2: bold, 3: italic, 4: bold/italic, 5:symbol)
                   labelColor  = "white",                # Label Color
                   lineColor   = "accent1",              # Line Color
                   shapeEnd    = c("none","sphere","circle","square","csquare","rectangle","crectangle","vrectangle"),
                   sizeEnd     = 10,
                   ...){ # nocov start
    ## rxPlot returns nothing, but plots a diagram
    if (class(RxODEobj) == "RxODE"){
        x <- RxODEobj$cmpMgr$rxDll()
    } else if (class(RxODEobj) == "RxCompilationManager") {
        x <- RxODEobj$rxDll()
    } else if (class(RxODEobj) == "rxDll"){
        x <- RxODEobj;
    } else {
        cat("The RxODEobj is not a supported RxODE object\n");
        return(invisible())
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
        cat("Package igraph needed for this function to work. Please install it.\n")
        return(invisible());
    }
    if (!rxDllLoaded(x)){
        cat("The shared RxODE library is not loaded, please load it by obj$dynLoad() or rxLoad(obj).\n");
        return(invisible());
    }
    ig <- igraph(x,family=family,tk=interactive,
                 shape = shape, size = size, colors = colors, fillColor = fillColor, font = font,
                 labelColor = labelColor, lineColor = lineColor, shapeEnd = shapeEnd, sizeEnd = sizeEnd);
    layout <- -igraph::layout.grid(ig);

    if (!interactive){
        op <- graphics::par()$mar
        graphics::par(mar=rep(0,4));
        plot(ig,edge.label.family=family, layout = layout, ...);
        suppressWarnings({graphics::par(mar=op)});
    } else {
        igraph::tkplot(ig,edge.label.family=family, layout = layout, ...);
    }
    invisible(NULL)
}# nocov end
## end function rxPlot

