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
      if (identical(x[[1]], quote(`-`)) & length(x) == 2) {
        if (length(x[[2]]) == 1){
          if (is.numeric(x[[2]])) return(TRUE)
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
      if (identical(x[[1]], quote(`-`)) & length(x) == 2) {
        if (is.numeric(x[[2]])) return(as.character(paste0("-", x[[2]])))
        if (length(x[[2]]) == 1){
          return(paste0("-SB_", x[[2]], "_SE"))}
        }
    }
    return(paste0("SB_", x, "_SE"))
  }
  is.temp <- function(x){
    (regexpr("\\(", x) != -1)
  }
  f <- function(x) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      .x <- as.character(x)
      if (tolower(.x) == "ka") assign(".hasKa2", TRUE, globalenv())
      if (!any(.x == c("ka", "k23", "k32", "k20", "k12", "k21", "k10", "k13", "k31", "k24", "k42", "tau", "r1", "r2", "b1", "b2", "tinf"))) {
        if (regexpr("^[A-Za-z]", .x) != -1){
          return(eval(parse(text=paste0("quote(AND_", .x, ")"))))
        }
      }
      return(x)
    } else if (is.pairlist(x)){
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`*`))) {
        .x1 <- deparse1(f(x[[2]]))
        .x2 <- deparse1(f(x[[3]]))
        if (is.num(x[[2]])) {
          if (is.temp(.x2)) {
            return(eval(parse(text=paste0("quote(prodd2X(", as.char(x[[2]]), ",", .x2, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(prodd2(", as.char(x[[2]]), ",", .x2, ",", dn, "))"))))
          }
        } else if (is.num(x[[3]])){
          if (is.temp(.x1)) {
            return(eval(parse(text=paste0("quote(prod2dX(", .x1, ",", as.char(x[[3]]) , "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(prod2d(", .x1, ",", as.char(x[[3]]) , ",", dn, "))"))))
          }
        }
        if (any(.x2 == c("(SB_tau_SE - SB_tinf_SE)", "(SB_tau_SE * SB_tinf_SE)"))){
          if (is.temp(.x1)){
            return(eval(parse(text=paste0("quote(prodd2X(", .x2, ",", .x1, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(prodd2(", .x2, ",", .x1, ",", dn, "))"))))
          }
        } else {
          if (is.temp(.x1)) {
            return(eval(parse(text=paste0("quote(prod2X(", .x1, ",", .x2, "))"))))
          } else if (is.temp(.x2)){
            return(eval(parse(text=paste0("quote(prod2X2(", .x1, ",", .x2, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(prod2(", .x2, ",", .x1, ",", dn, "))"))))
          }
        }
      } else if (identical(x[[1]], quote(`-`))) {
        if (length(x) == 2){
          .x2 <- deparse1(f(x[[2]]));
          if (is.num(x[[2]])){
            return(eval(parse(text=paste0("quote(-", .x2, ")"))))
          } else {
            if (is.temp(.x2)) {
              return(eval(parse(text=paste0("quote(negDX(", deparse1(f(x[[2]])), "))"))))
            } else {
              .cur <- get(".curDual", globalenv())
              dn <- paste0("rx_dn", .cur)
              assign(".curDual", .cur + 1, globalenv())
              return(eval(parse(text=paste0("quote(negD(", deparse1(f(x[[2]])), ",", dn, "))"))))
            }
          }
        }
        if (is.num(x[[2]]) && is.num(x[[3]])) {
          return(eval(parse(text=paste0("quote(", as.char(x[[2]]), "-", as.char(x[[3]]), ")"))))
        }
        if (is.num(x[[2]])) {
          .x3 <- deparse1(f(x[[3]]))
          if (is.temp(.x3)){
            return(eval(parse(text=paste0("quote(subtrd2X(", as.char(x[[2]]), ",", .x3, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(subtrd2(", as.char(x[[2]]), ",", .x3, ",", dn, "))"))))
          }
        } else if (is.num(x[[3]])){
          .x2 <- deparse1(f(x[[2]]))
          if (is.temp(.x2)){
             return(eval(parse(text=paste0("quote(subtr2dX(", .x2, ",", as.char(x[[3]]), "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(subtr2d(", .x2, ",", as.char(x[[3]]),",", dn, "))"))))
          }
        }
        .x2 <- deparse1(f(x[[2]]))
        .x3 <- deparse1(f(x[[3]]))
        if (is.temp(.x2)){
          return(eval(parse(text=paste0("quote(subtr2X(", .x2, ",", .x3, "))"))))
        } else if (is.temp(.x3)) {
          return(eval(parse(text=paste0("quote(subtr2X2(", .x2, ",", .x3, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(subtr2(", .x2, ",", .x3, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`+`))) {
        if (length(x) == 2){
          return(x[[2]])
        }
        if (is.num(x[[2]]) & is.num(x[[3]])) {
          return(eval(parse(text=paste0("quote(", as.char(x[[2]]), "+", as.char(x[[3]]), ")"))))
        } else if (is.num(x[[2]])) {
          .x3 <- deparse1(f(x[[3]]))
          if (is.temp(.x3)) {
            return(eval(parse(text=paste0("quote(addd2X(", as.char(x[[2]]), ",", .x3, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(addd2(", as.char(x[[2]]), ",", .x3, ",", dn, "))"))))
          }
        } else if (is.num(x[[3]])){
          .x2 <- deparse1(f(x[[2]]))
          if (is.temp(.x2)) {
            return(eval(parse(text=paste0("quote(add2dX(", .x2, ",", as.char(x[[3]]), "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(add2d(", .x2, ",", as.char(x[[3]]), ",", dn, "))"))))
          }
        }
        .x2 <- deparse1(f(x[[2]]))
        .x3 <- deparse1(f(x[[3]]))
        if (is.temp(.x2)) {
          return(eval(parse(text=paste0("quote(add2X(", .x2, ",", .x3, "))"))))
        } else if (is.temp(.x3)) {
          return(eval(parse(text=paste0("quote(add2X2(", .x2, ",", .x3, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(add2(", .x2, ",", .x3, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`sqrt`))) {
        .x2 <- deparse1(f(x[[2]]))
        if (is.temp(.x2)){
          return(eval(parse(text=paste0("quote(sqrtDX(", .x2, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(sqrtD(", .x2, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`/`))){
        if (is.num(x[[2]])) {
          .x3 <- deparse1(f(x[[3]]))
          if (is.temp(.x3)) {
            return(eval(parse(text=paste0("quote(divd2X(", as.char(x[[2]]), ",", .x3, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(divd2(", as.char(x[[2]]), ",", deparse1(f(x[[3]])), ",", dn, "))"))))
          }
        } else if (is.num(x[[3]])){
          .x2 <- deparse1(f(x[[2]]))
          if (is.temp(.x2)) {
            return(eval(parse(text=paste0("quote(div2dX(", .x2, ",", as.char(x[[3]]), "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(div2d(", .x2, ",", as.char(x[[3]]), ",", dn, "))"))))
          }
        }
        .x1 <- deparse1(f(x[[2]]))
        .x2 <- deparse1(f(x[[3]]))
        if (.x1 == "(SB_r2_SE + SB_r1_SE)"){
          if (is.temp(.x2)) {
            return(eval(parse(text=paste0("quote(divd2X(", .x1, ",", .x2, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(divd2(", .x1, ",", .x2, ",", dn, "))"))))
          }
        } else {
          if (is.temp(.x1)){
            return(eval(parse(text=paste0("quote(div2X(", .x1, ",", .x2, "))"))))
          } else if (is.temp(.x2)) {
            return(eval(parse(text=paste0("quote(div2X2(", .x1, ",", .x2, "))"))))
          } else {
            .cur <- get(".curDual", globalenv())
            dn <- paste0("rx_dn", .cur)
            assign(".curDual", .cur + 1, globalenv())
            return(eval(parse(text=paste0("quote(div2(", .x1, ",", .x2, ",", dn, "))"))))
          }
        }
      } else if (identical(x[[1]], quote(`exp`))) {
        .x2 <- deparse1(f(x[[2]]))
        if (is.temp(.x2)) {
          return(eval(parse(text=paste0("quote(expDX(", .x2, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(expD(", .x2, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`sin`))) {
        .x2 <- deparse1(f(x[[2]]))
        if (is.temp(.x2)) {
          return(eval(parse(text=paste0("quote(sinDX(", .x2, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(sinD(", .x2, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`cos`))) {
        .x2 <- deparse1(f(x[[2]]))
        if (is.temp(.x2)) {
          return(eval(parse(text=paste0("quote(cosDX(", .x2, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(cosD(", .x2, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`R_pow`))) {
        .x2 <- deparse1(f(x[[2]]))
        .x3 <- deparse1(f(x[[3]]))
        if (is.temp(.x2)) {
          return(eval(parse(text=paste0("quote(pow2dX(", .x2, ",", .x3, "))"))))
        } else if (is.temp(.x3)) {
          return(eval(parse(text=paste0("quote(pow2dX2(", .x2, ",", .x3, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(pow2d(", .x2, ",", .x3, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`atan2`))) {
        .x2 <- deparse1(f(x[[2]]))
        .x3 <- deparse1(f(x[[3]]))
        if (is.temp(.x2)) {
          return(eval(parse(text=paste0("quote(atan2DX(", .x2, ",", .x3, "))"))))
        } else if (is.temp(.x3)) {
          return(eval(parse(text=paste0("quote(atan2DX2(", .x2, ",", .x3, "))"))))
        } else {
          .cur <- get(".curDual", globalenv())
          dn <- paste0("rx_dn", .cur)
          assign(".curDual", .cur + 1, globalenv())
          return(eval(parse(text=paste0("quote(atan2D(", .x2, ",", .x3, ",", dn, "))"))))
        }
      } else if (identical(x[[1]], quote(`(`))) {
      } else {
        message("un-handled")
        print(x[[1]])
        print(x)
      }
      return(as.call(lapply(x, f)));
    }
  }
  return(gsub("AND_", "&", gsub("rx_dn", "&rx_dn", gsub(" +", "", deparse1(f(eval(parse(text=paste0("quote(", fromC(x), ")")))))))))
}

.lines <- readLines(devtools::package_file("src/lincmt.c"))

doLines <- function(.l){
  .l <- strsplit(.l, "\n")[[1]]
  .l <- strsplit(.l, "=")
  .l <- lapply(seq_along(.l), function(.i){
    .cur <- .l[[.i]];
    if (any(is.na(.cur))) return("");
    if (length(.cur) != 2) return("")
    .l1 <- gsub(" +", "", .cur[1])
    if (.l1 == "A3") assign(".hasA3", TRUE, globalenv())
    if (.l1 == "A3") assign(".hasA4", TRUE, globalenv())
    assign(".curDual", 0, globalenv())
    .l2 <- toDual(gsub(";", "", .cur[2]))
    assign(".maxDual", max(get(".maxDual", globalenv()),
                           get(".curDual", globalenv())), globalenv())
    if (.l1 == "A1" & .l2 == "0")
      return("dualN A1;\niniD(0, -1, &A1);")
    return(paste0("dualN ", .l1, ";\nassignD(&", .l1, ",", .l2, ");"))
  })
  .l <- .l[.l != ""]
  .l
}

getFun <- function(x="oneCmtKaRateSSr1"){
  assign(".hasA3", FALSE, globalenv())
  assign(".hasA4", FALSE, globalenv())
  assign(".hasKa2", FALSE, globalenv())
  assign(".curDual", 0, globalenv())
  assign(".maxDual", 0, globalenv())
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
  .hasKa <- any(tolower(.args) == "ka")
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
    if (any(x == c("ka", "k23", "k32", "k20", "k12", "k21", "k10", "k13", "k31", "k24", "k42"))) return(paste0("dualN *", x))
    paste0("double *", x)
  })
  .fargs <- paste0("(", paste(.fargs, collapse=", "), ")")
  .fB <- paste0("static inline dualN ", x, "G", .fargs, " {")

  .l <- doLines(.l)
  .ret <- gsub("_SE", ")", gsub("SB_", "(*", paste(unlist(.l), collapse="\n")));
  .fB <- paste0(.fB, "\n", paste0("dualN ", paste(paste0("rx_dn", seq(0, get(".maxDual", globalenv()) - 1)), collapse=","), ";"))

  if (.hasKa) {
    if (regexpr("oneCmt", x) != -1) {
      if (.hasAlast) {
        #define A1ka A[2]
        #define A2ka A[3]
        #define A2k20 A[4]
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "// has ka", "(void)(ka);"),
                         "dualN A1last;\nA1last.f = Alast[0];",
                         "A1last.grad[dKa] = Alast[2];",
                         "A1last.grad[dP1] = Alast[3];",
                         "A1last.grad[dV1] = Alast[4];",
                         "dualN A2last;\nA2last.f = Alast[1];",
                         "A2last.grad[dKa] = Alast[5];", # ka
                         "A2last.grad[dP1] = Alast[6];", # k20
                         "A2last.grad[dV1] = Alast[7];", # dVp
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         ## A1 derivatives
                         "A[2] = A1.grad[dKa];",
                         "A[3] = A1.grad[dP1];",
                         "A[4] = A1.grad[dV1];",
                         ## A2 derivatives
                         "A[5] = A2.grad[dKa];",
                         "A[6] = A2.grad[dP1];",
                         "A[7] = A2.grad[dV1];",
                         "return A2;",
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "// has ka", "(void)(ka);"),
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A1.grad[dKa];",
                         "A[3] = A1.grad[dP1];",
                         "A[4] = A1.grad[dV1];",
                         "A[5] = A2.grad[dKa];",
                         "A[6] = A2.grad[dP1];",
                         "A[7] = A2.grad[dV1];",
                         "return A2;",
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
                         ifelse(.hasKa2, "// has ka", "(void)(ka);"),
                         "dualN A1last;\nA1last.f =Alast[0];",
                         "A1last.grad[dKa] = Alast[3];",
                         "A1last.grad[dP1] = Alast[4];", # doesn't depend on k20
                         "A1last.grad[dP2] = Alast[5];",
                         "A1last.grad[dP3] = Alast[6];",
                         "A1last.grad[dV1] = Alast[7];",
                         "dualN A2last;\nA2last.f = Alast[1];",
                         #define A2ka A[6]
                         #define A2k20 A[7]
                         #define A2k23 A[8]
                         #define A2k32 A[9]
                         "A2last.grad[dKa] = Alast[8];", # ka
                         "A2last.grad[dP1] = Alast[9];", # k20
                         "A2last.grad[dP2] = Alast[10];", # k23
                         "A1last.grad[dP3] = Alast[11];", # k32
                         "A1last.grad[dV1] = Alast[12];", # k32
                         "dualN A3last;\nA3last.f = Alast[2];",
                         #define A3lastka Alast[10]
                         #define A3lastk20 Alast[11]
                         #define A3lastk23 Alast[12]
                         #define A3lastk32 Alast[13]
                         "A3last.grad[dKa] = Alast[13];", # ka
                         "A3last.grad[dP1] = Alast[14];", # k20
                         "A3last.grad[dP2] = Alast[15];", # k23
                         "A3last.grad[dP3] = Alast[16];", # k32
                         "A3last.grad[dV1] = Alast[17];", # k32
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         ## A1 derivatives
                         "A[3] = A1.grad[dKa];",
                         "A[4] = A1.grad[dP1];",
                         "A[5] = A1.grad[dP2];",
                         "A[6] = A1.grad[dP3];",
                         "A[7] = A1.grad[dV1];",
                         # doesn't depend on k20
                         "A[8]  = A2.grad[dKa];",
                         "A[9]  = A2.grad[dP1];",
                         "A[10] = A2.grad[dP2];",
                         "A[11] = A2.grad[dP3];",
                         "A[12] = A2.grad[dV1];",
                         ## A2 derivatives
                         "A[13] = A3.grad[dKa];", # ka
                         "A[14] = A3.grad[dP1];", # k20
                         "A[15] = A3.grad[dP2];", # k23
                         "A[16] = A3.grad[dP3];", # k32
                         "A[17] = A3.grad[dV1];", # k32
                         "return A2;",
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "//has ka", "(void)(ka);"),
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         ## A1 derivatives
                         "A[3] = A1.grad[dKa];",
                         "A[4] = A1.grad[dP1];",
                         "A[5] = A1.grad[dP2];",
                         "A[6] = A1.grad[dP3];",
                         "A[7] = A1.grad[dV1];",
                         # A2
                         "A[8]  = A2.grad[dKa];",
                         "A[9]  = A2.grad[dP1];",
                         "A[10] = A2.grad[dP2];",
                         "A[11] = A2.grad[dP3];",
                         "A[12] = A2.grad[dV1];",
                         ## A3 derivatives
                         "A[13] = A3.grad[dKa];", # ka
                         "A[14] = A3.grad[dP1];", # k20
                         "A[15] = A3.grad[dP2];", # k23
                         "A[16] = A3.grad[dP3];", # k32
                         "A[17] = A3.grad[dV1];", # k32
                         "return A2;",
                         "}"),
                       collapse="\n")
      }
    } else if (regexpr("threeCmt", x) != -1) {
      if (.hasAlast) {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "// has ka", "(void)(ka);"),
                         "dualN A1last;\nA1last.f = Alast[0];",
                         "A1last.grad[dKa] = Alast[4];",
                         "A1last.grad[dP1] = Alast[5];", # doesn't depend on k20
                         "A1last.grad[dP2] = Alast[6];",
                         "A1last.grad[dP3] = Alast[7];",
                         "A1last.grad[dP4] = Alast[8];",
                         "A1last.grad[dP5] = Alast[9];",
                         "A1last.grad[dV1] = Alast[10];",
                         "dualN A2last;\nA2last.f = Alast[1];",
                         "A2last.grad[dKa] = Alast[11];",
                         "A2last.grad[dP1] = Alast[12];",
                         "A2last.grad[dP2] = Alast[13];",
                         "A1last.grad[dP3] = Alast[14];",
                         "A1last.grad[dP4] = Alast[15];",
                         "A1last.grad[dP5] = Alast[16];",
                         "A1last.grad[dV1] = Alast[17];",
                         "dualN A3last;\nA3last.f = Alast[2];",
                         "A3last.grad[dKa] = Alast[18];",
                         "A3last.grad[dP1] = Alast[19];",
                         "A3last.grad[dP2] = Alast[20];",
                         "A3last.grad[dP3] = Alast[21];",
                         "A3last.grad[dP4] = Alast[22];",
                         "A3last.grad[dP5] = Alast[23];",
                         "A3last.grad[dV1] = Alast[24];",
                         "dualN A4last;\nA4last.f = Alast[3];",
                         "A4last.grad[dKa] = Alast[25];",
                         "A4last.grad[dP1] = Alast[26];",
                         "A4last.grad[dP2] = Alast[27];",
                         "A4last.grad[dP3] = Alast[28];",
                         "A4last.grad[dP4] = Alast[29];",
                         "A4last.grad[dP5] = Alast[30];",
                         "A4last.grad[dV1] = Alast[31];",
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         "A[3] = A4.f;",
                         ## A1 derivatives
                         "A[4] = A1.grad[dKa];",
                         "A[5] = A1.grad[dP1];",
                         "A[6] = A1.grad[dP2];",
                         "A[7] = A1.grad[dP3];",
                         "A[8] = A1.grad[dP4];",
                         "A[9] = A1.grad[dP5];",
                         "A[10] = A1.grad[dV1];",
                         ## A2 derivatives
                         "A[11] = A2.grad[dKa];",
                         "A[12] = A2.grad[dP1];",
                         "A[13] = A2.grad[dP2];",
                         "A[14] = A2.grad[dP3];",
                         "A[15] = A2.grad[dP4];",
                         "A[16] = A2.grad[dP5];",
                         "A[17] = A2.grad[dV1];",
                         ## A3 derivatives
                         "A[18] = A3.grad[dKa];",
                         "A[19] = A3.grad[dP1];",
                         "A[20] = A3.grad[dP2];",
                         "A[21] = A3.grad[dP3];",
                         "A[22] = A3.grad[dP4];",
                         "A[23] = A3.grad[dP5];",
                         "A[24] = A3.grad[dV1];",
                         ## A4 derivatives
                         "A[25] = A4.grad[dKa];",
                         "A[26] = A4.grad[dP1];",
                         "A[27] = A4.grad[dP2];",
                         "A[28] = A4.grad[dP3];",
                         "A[29] = A4.grad[dP4];",
                         "A[30] = A4.grad[dP5];",
                         "A[31] = A4.grad[dV1];",
                         "return A2;",
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         ifelse(.hasKa2, "// has Ka", "(void)(ka);"),
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         "A[3] = A4.f;",
                         ## A1 derivatives
                         "A[4] = A1.grad[dKa];",
                         "A[5] = A1.grad[dP1];",
                         "A[6] = A1.grad[dP2];",
                         "A[7] = A1.grad[dP3];",
                         "A[8] = A1.grad[dP4];",
                         "A[9] = A1.grad[dP5];",
                         "A[10] = A1.grad[dV1];",
                         ## A2 derivatives
                         "A[11] = A2.grad[dKa];",
                         "A[12] = A2.grad[dP1];",
                         "A[13] = A2.grad[dP2];",
                         "A[14] = A2.grad[dP3];",
                         "A[15] = A2.grad[dP4];",
                         "A[16] = A2.grad[dP5];",
                         "A[17] = A2.grad[dV1];",
                         ## A3 derivatives
                         "A[18] = A3.grad[dKa];",
                         "A[19] = A3.grad[dP1];",
                         "A[20] = A3.grad[dP2];",
                         "A[21] = A3.grad[dP3];",
                         "A[22] = A3.grad[dP4];",
                         "A[23] = A3.grad[dP5];",
                         "A[24] = A3.grad[dV1];",
                         ## A4 derivatives
                         "A[25] = A4.grad[dKa];",
                         "A[26] = A4.grad[dP1];",
                         "A[27] = A4.grad[dP2];",
                         "A[28] = A4.grad[dP3];",
                         "A[29] = A4.grad[dP4];",
                         "A[30] = A4.grad[dP5];",
                         "A[31] = A4.grad[dV1];",
                         "return A2;",
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
                         "dualN A1last;\nA1last.f = Alast[0];",
                         "A1last.grad[dKa] = Alast[1];",
                         .ret,
                         "A[0] = A1.f;",
                         ## A1 derivatives
                         "A[1] = A1.grad[dP1];",
                         "A[2] = A1.grad[dV1];",
                         "return A1;",
                         ## A2 derivatives
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A1.grad[dP1];",
                         "A[2] = A1.grad[dV1];",
                         "return A1;",
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
                         "dualN A1last;\nA1last.f = Alast[0];",
                         "A1last.grad[dP1] = Alast[2];",
                         "A1last.grad[dP2] = Alast[3];",
                         "A1last.grad[dP3] = Alast[4];",
                         "A1last.grad[dV1] = Alast[5];",
                         "dualN A2last;\nA2last.f=Alast[1];",
                         "A2last.grad[dP1] = Alast[5];", # k10
                         "A2last.grad[dP2] = Alast[6];", # k12
                         "A2last.grad[dP3] = Alast[7];", # k21
                         "A2last.grad[dV1] = Alast[8];", # k21
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         ## A1 derivatives
                         "A[2] = A1.grad[dP1];",
                         "A[3] = A1.grad[dP2];",
                         "A[4] = A1.grad[dP3];",
                         "A[5] = A1.grad[dV1];",
                         ## A2 derivatives
                         "A[6] = A2.grad[dP1];", # ka
                         "A[7] = A2.grad[dP2];", # k20
                         "A[8] = A2.grad[dP3];", # k23
                         "A[9] = A2.grad[dV1];", # k23
                         "return A1;",
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         ## A1 derivatives
                         "A[2] = A1.grad[dP1];",
                         "A[3] = A1.grad[dP2];",
                         "A[4] = A1.grad[dP3];",
                         "A[5] = A1.grad[dV1];",

                         ## A2 derivatives
                         "A[6] = A2.grad[dP1];", # ka
                         "A[7] = A2.grad[dP2];", # k20
                         "A[8] = A2.grad[dP3];", # k23
                         "A[9] = A2.grad[dV1];", # k23
                         "return A1;",
                         "}"),
                       collapse="\n")
      }
    } else if (regexpr("threeCmt", x) != -1) {
      if (.hasAlast) {
        .ret <- paste0(c(.fB,
                         "dualN A1last;\nA1last.f = Alast[0];",
                         "A1last.grad[dP1] = Alast[3];",
                         "A1last.grad[dP2] = Alast[4];",
                         "A1last.grad[dP3] = Alast[5];",
                         "A1last.grad[dP4] = Alast[6];",
                         "A1last.grad[dP5] = Alast[7];",
                         "A1last.grad[dV1] = Alast[8];",
                         "dualN A2last;\nA2last.f = Alast[1];",
                         "A2last.grad[dP1] = Alast[9];",
                         "A2last.grad[dP2] = Alast[10];",
                         "A2last.grad[dP3] = Alast[11];",
                         "A2last.grad[dP4] = Alast[12];",
                         "A2last.grad[dP5] = Alast[13];",
                         "A2last.grad[dV1] = Alast[14];",
                         "dualN A3last;\nA3last.f = Alast[2];",
                         "A3last.grad[dP1] = Alast[15];",
                         "A3last.grad[dP2] = Alast[16];",
                         "A3last.grad[dP3] = Alast[17];",
                         "A3last.grad[dP4] = Alast[18];",
                         "A3last.grad[dP5] = Alast[19];",
                         "A3last.grad[dV1] = Alast[20];",
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         ## A1 derivatives
                         "A[3] = A1.grad[dP1];",
                         "A[4] = A1.grad[dP2];",
                         "A[5] = A1.grad[dP3];",
                         "A[6] = A1.grad[dP4];",
                         "A[7] = A1.grad[dP5];",
                         "A[8] = A1.grad[dV1];",
                         ## A2 derivatives
                         "A[9] = A2.grad[dP1];",
                         "A[10] = A2.grad[dP2];",
                         "A[11] = A2.grad[dP3];",
                         "A[12] = A2.grad[dP4];",
                         "A[13] = A2.grad[dP5];",
                         "A[14] = A2.grad[dV1];",
                         ## A3 derivatives
                         "A[15] = A3.grad[dP1];",
                         "A[16] = A3.grad[dP2];",
                         "A[17] = A3.grad[dP3];",
                         "A[18] = A3.grad[dP4];",
                         "A[19] = A3.grad[dP5];",
                         "A[20] = A3.grad[dV1];",
                         "return A1;",
                         "}"),
                       collapse="\n")
      } else {
        .ret <- paste0(c(.fB,
                         .ret,
                         "A[0] = A1.f;",
                         "A[1] = A2.f;",
                         "A[2] = A3.f;",
                         ## A1 derivatives
                         "A[3] = A1.grad[dP1];",
                         "A[4] = A1.grad[dP2];",
                         "A[5] = A1.grad[dP3];",
                         "A[6] = A1.grad[dP4];",
                         "A[7] = A1.grad[dP5];",
                         "A[8] = A1.grad[dV1];",
                         ## A2 derivatives
                         "A[9] = A2.grad[dP1];",
                         "A[10] = A2.grad[dP2];",
                         "A[11] = A2.grad[dP3];",
                         "A[12] = A2.grad[dP4];",
                         "A[13] = A2.grad[dP5];",
                         "A[14] = A2.grad[dV1];",
                         ## A3 derivatives
                         "A[15] = A3.grad[dP1];",
                         "A[16] = A3.grad[dP2];",
                         "A[17] = A3.grad[dP3];",
                         "A[18] = A3.grad[dP4];",
                         "A[19] = A3.grad[dP5];",
                         "A[20] = A3.grad[dV1];",
                         "return A1;",
                         "}"),
                       collapse="\n")
      }
    }
  }
  return(.ret)
}

unlink(devtools::package_file("src/lincmtB1g.h"))
if (!file.exists(devtools::package_file("src/lincmtB1g.h"))){
  ## sink(devtools::package_file("src/lincmtB2d.h"));
  fs <- c("oneCmtRateSSr1", "oneCmtRateSS",
          "oneCmtRate", "oneCmtBolusSS",
          "oneCmtKaRateSSr1", "oneCmtKaRateSSr2", "oneCmtKaRateSStr1", "oneCmtKaRateSStr2",
          "oneCmtKaRate", "oneCmtKaSSb1", "oneCmtKaSSb2")
  .f <- paste(sapply(fs, getFun), collapse="\n");

  writeLines(c("
#ifndef linCmtB1g_header
#define linCmtB1g_header
#include \"dual.h\"
", .f, "#endif"), devtools::package_file("src/lincmtB1g.h"))
}

unlink(devtools::package_file("src/lincmtB2g.h"))
if (!file.exists(devtools::package_file("src/lincmtB2g.h"))){
  ## sink(devtools::package_file("src/lincmtB2g.h"));
  fs <- c("twoCmtRateSSr1", "twoCmtRateSS",
          "twoCmtRate", "twoCmtBolusSS",
          "twoCmtKaRateSSr1", "twoCmtKaRateSSr2", "twoCmtKaRateSStr1", "twoCmtKaRateSStr2",
          "twoCmtKaRate", "twoCmtKaSSb1", "twoCmtKaSSb2")
  .f <- paste(sapply(fs, getFun), collapse="\n");
  writeLines(c("
#ifndef linCmtB2g_header
#define linCmtB2g_header
#include \"dual.h\"
", .f, "#endif"), devtools::package_file("src/lincmtB2g.h"))
}

unlink(devtools::package_file("src/lincmtB3g.h"))
if (!file.exists(devtools::package_file("src/lincmtB3g.h"))) {
  fs <- c("threeCmtRateSSr1", "threeCmtRateSS",
          "threeCmtRate", "threeCmtBolusSS",
          "threeCmtKaRateSSr1", "threeCmtKaRateSSr2", "threeCmtKaRateSStr1", "threeCmtKaRateSStr2",
          "threeCmtKaRate", "threeCmtKaSSb1", "threeCmtKaSSb2")
  .f <- paste(sapply(fs, getFun), collapse="\n");
  writeLines(c("
#ifndef linCmtB3g_header
#define linCmtB3g_header
#include \"dual.h\"
", .f, "#endif"), devtools::package_file("src/lincmtB3g.h"))
}

## ##toDual("(*b2)-a21+a22-a23-a24+a25")
