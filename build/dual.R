fromC <- function(x) {
  gsub(" *[*]([^=\n]*)=","\\1=",gsub("double ","",gsub("[(][*]([^*)]+)[)]","\\1", x)), perl=TRUE)
}
toDual <- function(x="0.5*(s - sqrt(s*s - 4*(*k32)*(*k20)))"){
  is.num <- function(x){
    if (is.numeric(x)){
      return(TRUE)
    } else if (is.name(x)){
      if (any(as.character(x) == c("t", "r1", "r2", "b1", "b2", "tau", "tinf"))) {
        return(TRUE)
      }
    } else if (is.call(x)){
      if (identical(x[[1]], quote(`-`)) & length(x == 2)) {
        if (length(x[[2]]) == 1){
          if (any(as.character(x[[2]]) == c("t", "r1", "r2", "b1", "b2", "tau", "tinf"))) {
            return(TRUE)
          }
        }
      }
    }
    return(FALSE)
  }
  as.char <- function(x){
    if (is.numeric(x)) return(as.character(x))
    if (is.call(x)) {
      if (identical(x[[1]], quote(`-`)) & length(x == 2)) {
        if (length(x[[2]]) == 1){
          return(paste0("-SB_", x[[2]], "_SE"))}
        }
    }
    return(paste0("SB_", x, "_SE"))
  }
  f <- function(x) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      if (as.character(x) == "ka") assign(".hasKa2", TRUE, globalenv())
      return(x)
    } else if (is.pairlist(x)){
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`*`))) {
        if (is.num(x[[2]])) {
          return(eval(parse(text=paste0("quote(prodd2(", as.char(x[[2]]), ",", deparse1(f(x[[3]])), "))"))))
        } else if (is.num(x[[3]])){
          return(eval(parse(text=paste0("quote(prod2d(", deparse1(f(x[[2]])), ",", as.char(x[[3]]) , "))"))))
        }
        .x1 <- deparse1(f(x[[2]]))
        .x2 <- deparse1(f(x[[3]]))
        if (any(.x2 == c("(SB_tau_SE - SB_tinf_SE)", "(SB_tau_SE * SB_tinf_SE)"))){
          return(eval(parse(text=paste0("quote(prodd2(", .x2, ",", .x1, "))"))))
        } else {
          return(eval(parse(text=paste0("quote(prod2(", .x1, ",", .x2, "))"))))
        }
      } else if (identical(x[[1]], quote(`-`))) {
        if (length(x) == 2){
          if (is.num(x[[2]])){
            return(eval(parse(text=paste0("quote(-", deparse1(f(x[[2]])), ")"))))
          } else {
            return(eval(parse(text=paste0("quote(negD(", deparse1(f(x[[2]])), "))"))))
          }
        }
        if (is.num(x[[2]]) && is.num(x[[3]])) {
          return(eval(parse(text=paste0("quote(", as.char(x[[2]]), "-", as.char(x[[3]]), ")"))))
        }
        if (is.num(x[[2]])) {
          return(eval(parse(text=paste0("quote(subtrd2(", as.char(x[[2]]), ",", deparse1(f(x[[3]])), "))"))))
        } else if (is.num(x[[3]])){
          return(eval(parse(text=paste0("quote(subtr2d(", deparse1(f(x[[2]])), ",", as.char(x[[3]]), "))"))))
        }
        return(eval(parse(text=paste0("quote(subtr2(", deparse1(f(x[[2]])), ",", deparse1(f(x[[3]])), "))"))))
      } else if (identical(x[[1]], quote(`+`))) {
        if (length(x) == 2){
          return(x[[2]])
        }
        if (is.num(x[[2]]) & is.num(x[[3]])) {
          return(eval(parse(text=paste0("quote(", as.char(x[[2]]), "+", as.char(x[[3]]), ")"))))
        } else if (is.num(x[[2]])) {
          return(eval(parse(text=paste0("quote(addd2(", as.char(x[[2]]), ",", deparse1(f(x[[3]])), "))"))))
        } else if (is.num(x[[3]])){
          return(eval(parse(text=paste0("quote(add2d(", deparse1(f(x[[2]])), ",", as.char(x[[3]]), "))"))))
        }
        return(eval(parse(text=paste0("quote(add2(", deparse1(f(x[[2]])), ",", deparse1(f(x[[3]])), "))"))))
      } else if (identical(x[[1]], quote(`sqrt`))) {
        return(eval(parse(text=paste0("quote(sqrtD(", deparse1(f(x[[2]])), "))"))))
      } else if (identical(x[[1]], quote(`/`))){
        if (is.num(x[[2]])) {
          return(eval(parse(text=paste0("quote(divd2(", as.char(x[[2]]), ",", deparse1(f(x[[3]])), "))"))))
        } else if (is.num(x[[3]])){
          return(eval(parse(text=paste0("quote(div2d(", deparse1(f(x[[2]])), ",", as.char(x[[3]]), "))"))))
        }
        .x1 <- deparse1(f(x[[2]]))
        .x2 <- deparse1(f(x[[3]]))
        if (.x1 == "(SB_r2_SE + SB_r1_SE)"){
          return(eval(parse(text=paste0("quote(divd2(", .x1, ",", .x2, "))"))))
        } else {
          return(eval(parse(text=paste0("quote(div2(", .x1, ",", .x2, "))"))))
        }
      } else if (identical(x[[1]], quote(`exp`))) {
        return(eval(parse(text=paste0("quote(expD(", deparse1(f(x[[2]])), "))"))))
      } else if (identical(x[[1]], quote(`(`))) {
      } else {
        message("un-handled")
        print(x[[1]])
        print(x)
      }
      return(as.call(lapply(x, f)));
    }
  }
  return(gsub(" +", "", deparse1(f(eval(parse(text=paste0("quote(", fromC(x), ")")))))))
}

.lines <- readLines(devtools::package_file("src/lincmt.c"))

getFun <- function(x="oneCmtKaRateSSr1"){
  assign(".hasA3", FALSE, globalenv())
  assign(".hasA4", FALSE, globalenv())
  assign(".hasKa2", FALSE, globalenv())
  .w <- which(regexpr(paste0("void ", x, " *[(]"), .lines) != -1)[1]
  .l <- .lines[-seq(1, .w - 1)];
  .w <- which(regexpr("}", .l) != -1)[1]
  .l <- .l[seq(1, .w)]
  .l <- gsub("\t", " ", .l)
  while (length(grep("[{]", .l[1])) == 0L) {
    .l[2] <- paste(.l[1:2], collapse=" ")
    .l <- .l[-1]
  }
  .args <- strsplit(gsub(",", " ", gsub("double *[*]", " ", gsub(".*[(]([^)]*)[)].*", "\\1", .l[1]))), " +")[[1]]
  .args <- .args[.args != ""]
  .args <- .args[.args != "A"]
  .hasAlast <- any(.args == "Alast")
  .hasKa <- any(.args == "ka")
  .args <- .args[.args != "Alast"]
  .args2 <- .args[.args != "tinf"]
  .args2 <- .args2[.args2 != "tau"]
  .args2 <- .args2[.args2 != "r1"]
  .args2 <- .args2[.args2 != "r2"]
  .args2 <- .args2[.args2 != "b1"]
  .args2 <- .args2[.args2 != "b2"]
  .args2 <- .args2[.args2 != "t"]

  assign(".args", .args2, globalenv())
  .l <- .l[-1]
  .l <- paste(.l[-length(.l)], collapse="\n");
  .l <- fromC(.l)
  if (.hasAlast){
    .fargs <- c("A", "Alast", .args)
  } else {
    .fargs <- c("A", .args)
  }
  .fargs <- sapply(.fargs, function(x){
    if (any(x == c("ka", "k23", "k32", "k20", "k12", "k21", "k10"))) return(paste0(x, "d"))
    x
  })
  .fargs <- paste0("(double *", paste(.fargs, collapse=", double *"), ")")
  .fB <- paste0("static inline void ", x, "D", .fargs, " {")
  .l <- strsplit(.l, "\n")[[1]]
  .l <- strsplit(.l, "=")
  .l <- lapply(seq_along(.l), function(.i){
    .cur <- .l[[.i]];
    if (any(is.na(.cur))) return("");
    if (length(.cur) != 2) return("")
    .l1 <- gsub(" +", "", .cur[1])
    if (.l1 == "A3") assign(".hasA3", TRUE, globalenv())
    if (.l1 == "A3") assign(".hasA4", TRUE, globalenv())
    .l2 <- toDual(gsub(";", "", .cur[2]))
    if (.l1 == "A1" & .l2 == "0")
      return("dualN A1=iniD(0, -1);")
    return(paste0("dualN ", .l1, "=", .l2, ";"))
  })
  .l <- .l[.l != ""]
  .ret <- gsub("_SE", ")", gsub("SB_", "(*", paste(unlist(.l), collapse="\n")));

  if (.hasKa) {
    if (regexpr("oneCmt", x) != -1) {
      if (.hasAlast) {
        #define A1ka A[2]
        #define A2ka A[3]
        #define A2k20 A[4]
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "dualN ka  = iniD(*kad,  0);", "(void)(*kad);"),
                         "dualN k20 = iniD(*k20d, 1);",
                         "dualN A1last = iniD(Alast[0],-1);",
                         "A1last.grad[0] = Alast[2];",
                         "A1last.grad[1] = 0.0;", # doesn't depend on k20
                         "dualN A2last = iniD(Alast[1],-1);",
                         "A2last.grad[0] = Alast[3];", # ka
                         "A2last.grad[1] = Alast[4];", # k20
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         ## A1 derivatives
                         "A[2] = A1.grad[0];",
                         ## A2 derivatives
                         "A[3] = A2.grad[0];",
                         "A[4] = A2.grad[1];",
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "dualN ka  = iniD(*kad,  0);", "(void)(*kad);"),
                         "dualN k20 = iniD(*k20d, 1);",
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A1.grad[0];",
                         "A[3] = A2.grad[0];",
                         "A[4] = A2.grad[1];",
                         "}"),
                       collapse="\n")
      }
    } else if (regexpr("twoCmt", x) != -1) {
      #define A2ka A[6]
      #define A2k20 A[7]
      #define A2k23 A[8]
      #define A2k32 A[9]
      if (.hasAlast) {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "dualN ka  = iniD(*kad,  0);", "(void)(*kad);"),
                         "dualN k20 = iniD(*k20d, 1);",
                         "dualN k23 = iniD(*k23d, 2);",
                         "dualN k32 = iniD(*k32d, 3);",
                         "dualN A1last = iniD(Alast[0],-1);",
                         "A1last.grad[0] = Alast[3];",
                         "A1last.grad[1] = 0.0;", # doesn't depend on k20
                         "A1last.grad[2] = Alast[4];",
                         "A1last.grad[3] = Alast[5];",
                         "dualN A2last = iniD(Alast[1],-1);",
                         #define A2ka A[6]
                         #define A2k20 A[7]
                         #define A2k23 A[8]
                         #define A2k32 A[9]
                         "A2last.grad[0] = Alast[6];", # ka
                         "A2last.grad[1] = Alast[7];", # k20
                         "A2last.grad[2] = Alast[8];", # k23
                         "A1last.grad[3] = Alast[9];", # k32
                         "dualN A3last = iniD(Alast[2],-1);",
                         #define A3lastka Alast[10]
                         #define A3lastk20 Alast[11]
                         #define A3lastk23 Alast[12]
                         #define A3lastk32 Alast[13]
                         "A3last.grad[0] = Alast[10];", # ka
                         "A3last.grad[1] = Alast[11];", # k20
                         "A3last.grad[2] = Alast[12];", # k23
                         "A3last.grad[3] = Alast[13];", # k32
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         ## A1 derivatives
                         "A[3] = A1.grad[0];",
                         # doesn't depend on k20
                         "A[4] = A1.grad[2];",
                         "A[5] = A1.grad[3];",
                         ## A2 derivatives
                         "A[6] = A2.grad[0];", # ka
                         "A[7] = A2.grad[1];", # k20
                         "A[8] = A2.grad[2];", # k23
                         "A[9] = A2.grad[3];", # k32
                         ## A3 derivatives
                         "A[10] = A3.grad[0];", # ka
                         "A[11] = A3.grad[1];", # k20
                         "A[12] = A3.grad[2];", # k23
                         "A[13] = A3.grad[3];", # k32
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "dualN ka  = iniD(*kad,  0);", "(void)(*kad);"),
                         "dualN k20 = iniD(*k20d, 1);",
                         "dualN k23 = iniD(*k23d, 2);",
                         "dualN k32 = iniD(*k32d, 3);",
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         ## A1 derivatives
                         "A[3] = A1.grad[0];", # ka
                         "A[4] = A1.grad[2];", # k23
                         "A[5] = A1.grad[3];", # k32
                         ## A2 derivatives
                         "A[6] = A2.grad[0];", # ka
                         "A[7] = A2.grad[1];", # k20
                         "A[8] = A2.grad[2];", # k23
                         "A[9] = A2.grad[3];", # k32
                         ## A3 derivatives
                         "A[10] = A3.grad[0];", # ka
                         "A[11] = A3.grad[1];", # k20
                         "A[12] = A3.grad[2];", # k23
                         "A[13] = A3.grad[3];", # k32
                         "}"),
                       collapse="\n")
      }
    }
  } else {
    if (regexpr("oneCmt", x) != -1) {
      if (.hasAlast) {
        #define A1ka A[2]
        #define A2ka A[3]
        #define A2k20 A[4]
        .ret <- paste0(c(.fB,
                         "dualN k10 = iniD(*k10d, 0);",
                         "dualN A1last = iniD(Alast[0],-1);",
                         "A1last.grad[0] = Alast[1];",
                         .ret,
                         "A[0] = A1.f;",
                         ## A1 derivatives
                         "A[1] = A1.grad[0];",
                         ## A2 derivatives
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         "dualN k10 = iniD(*k10d, 0);",
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A1.grad[0];",
                         "}"),
                       collapse="\n")
      }
    } else if (regexpr("twoCmt", x) != -1) {
      if (.hasAlast) {
        #define A1k10 A[2]
        #define A1k12 A[3]
        #define A1k21 A[4]

        #define A2k10 A[5]
        #define A2k12 A[6]
        #define A2k21 A[7]
        .ret <- paste0(c(.fB,
                         "dualN k10 = iniD(*k10d, 0);",
                         "dualN k12 = iniD(*k12d, 1);",
                         "dualN k21 = iniD(*k21d, 2);",
                         "dualN A1last = iniD(Alast[0],-1);",
                         "A1last.grad[0] = Alast[2];",
                         "A1last.grad[1] = Alast[3];",
                         "A1last.grad[2] = Alast[4];",
                         "dualN A2last = iniD(Alast[1],-1);",
                         "A2last.grad[0] = Alast[5];", # k10
                         "A2last.grad[1] = Alast[6];", # k12
                         "A2last.grad[2] = Alast[7];", # k21
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         ## A1 derivatives
                         "A[2] = A1.grad[0];",
                         "A[3] = A1.grad[1];",
                         "A[4] = A1.grad[2];",
                         ## A2 derivatives
                         "A[5] = A2.grad[0];", # ka
                         "A[6] = A2.grad[1];", # k20
                         "A[7] = A2.grad[2];", # k23
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         "dualN k10 = iniD(*k10d, 0);",
                         "dualN k12 = iniD(*k12d, 1);",
                         "dualN k21 = iniD(*k21d, 2);",
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         ## A1 derivatives
                         "A[2] = A1.grad[0];",
                         "A[3] = A1.grad[1];",
                         "A[4] = A1.grad[2];",
                         ## A2 derivatives
                         "A[5] = A2.grad[0];", # ka
                         "A[6] = A2.grad[1];", # k20
                         "A[7] = A2.grad[2];", # k23
                         "}"),
                       collapse="\n")
      }
    }
  }
  return(.ret)
}

unlink(devtools::package_file("src/lincmtB1d.h"))
if (!file.exists(devtools::package_file("src/lincmtB1d.h"))){
  ## sink(devtools::package_file("src/lincmtB2d.h"));
  fs <- c("oneCmtRateSSr1", "oneCmtRateSS",
          "oneCmtRate", "oneCmtBolusSS",
          "oneCmtKaRateSSr1", "oneCmtKaRateSSr2", "oneCmtKaRateSStr1", "oneCmtKaRateSStr2",
          "oneCmtKaRate", "oneCmtKaSSb1", "oneCmtKaSSb2")
  .f <- paste(sapply(fs, getFun), collapse="\n");

  writeLines(c("
#ifndef linCmtB1_header
#define linCmtB1_header
#include \"dual.h\"
", .f, "#endif"), devtools::package_file("src/lincmtB1d.h"))
}

unlink(devtools::package_file("src/lincmtB2d.h"))
if (!file.exists(devtools::package_file("src/lincmtB2d.h"))){
  ## sink(devtools::package_file("src/lincmtB2d.h"));
  fs <- c("twoCmtRateSSr1", "twoCmtRateSS",
          "twoCmtRate", "twoCmtBolusSS",
          "twoCmtKaRateSSr1", "twoCmtKaRateSSr2", "twoCmtKaRateSStr1", "twoCmtKaRateSStr2",
          "twoCmtKaRate", "twoCmtKaSSb1", "twoCmtKaSSb2")
  .f <- paste(sapply(fs, getFun), collapse="\n");

  writeLines(c("
#ifndef linCmtB2_header
#define linCmtB2_header
#include \"dual.h\"
", .f, "#endif"), devtools::package_file("src/lincmtB2d.h"))
}

## toDual("exp(-(*ka)*(*t))")
