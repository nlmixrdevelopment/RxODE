rxodeTest({

  .rx <- loadNamespace("RxODE")

  test_that("comments are parsed correctly", {

    cmt <- c("function() {", "      ini({", "        ## You may label each parameter with a comment",
             "        tka <- 0.45 # Log Ka", "        tcl <- log(c(0, 2.7, 100)) # Log Cl",
             "        ## This works with interactive models", "        ## You may also label the preceding line with label(\"label text\")",
             "        tv <- 3.45; label(\"log V\")", "        ## the label(\"Label name\") works with all models",
             "        eta.ka ~ 0.6", "        eta.cl ~ 0.3", "        eta.v ~ 0.1",
             "        add.sd <- 0.7", "      })", "      model({", "        ka <- exp(tka + eta.ka)",
             "        cl <- exp(tcl + eta.cl)", "        v <- exp(tv + eta.v)",
             "        linCmt() ~ add(add.sd)", "      })", "    }")

    eq <- c("function () ", "{", "    ini({", "        tka <- 0.45", "        label(\"Log Ka\")",
            "        tcl <- log(c(0, 2.7, 100))", "        label(\"Log Cl\")",
            "        tv <- 3.45", "        label(\"log V\")", "        eta.ka ~ 0.6",
            "        eta.cl ~ 0.3", "        eta.v ~ 0.1", "        add.sd <- 0.7",
            "    })", "    model({", "        ka <- exp(tka + eta.ka)", "        cl <- exp(tcl + eta.cl)",
            "        v <- exp(tv + eta.v)", "        linCmt() ~ add(add.sd)",
            "    })", "}")


    expect_equal(.rx$.rxReplaceCommentWithLabel(cmt), eq)

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) | tmp
      })
    }

    str <- .rx$.rxFunction2string(one.cmt)

    if (!is.null(attr(one.cmt, "srcref"))) {
      expect_equal(str, eq)
      attr(one.cmt, "srcref") <- NULL
    }

    str <- .rx$.rxFunction2string(one.cmt)

    expect_equal(str,
                 c("function () ", "{", "    ini({", "        tka <- 0.45", "        label(\"Log Ka\")",
                   "        tcl <- log(c(0, 2.7, 100))", "        label(\"Log Cl\")",
                   "        tv <- 3.45", "        label(\"log V\")", "        eta.ka ~ 0.6",
                   "        eta.cl ~ 0.3", "        eta.v ~ 0.1", "        add.sd <- 0.7",
                   "    })", "    model({", "        ka <- exp(tka + eta.ka)", "        cl <- exp(tcl + eta.cl)",
                   "        v <- exp(tv + eta.v)", "        linCmt() ~ add(add.sd) | tmp",
                   "    })", "}"))

  })

  test_that("meta information parsing", {

    one.cmt <- function() {
      meta1 <- "meta"
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      meta2 <- "meta2"
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    tmp1 <- one.cmt()

    expect_equal(tmp1$meta$meta1, "meta")
    expect_equal(tmp1$meta$meta2, "meta2")


    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
        lambda <- c(-2, 1, 2)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) + boxCox(lambda) | tmp
      })
    }

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ lnorm(add.sd) | tmp
      })
    }

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
        bLambda <- c(0, 3)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ lnorm(add.sd) | tmp
        tmp2 ~ dpois(bLambda)
      })
    }


    cov <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        tvp <- 3.45; label("log V")
        cl.wt <- 0.1
        v.wt <- 0.1
        cl.sex <- 0.1
        v.sex <- 0.1
        cl.age <- 0.1
        v.age <- 0.1
        vp.wt <- 1
        vp.sex <- 1
        vp.age <- 1
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })

      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
        v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
        vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl/v * center
        cp = center/v
        cp ~ add(add.sd)
      })
    }




  })

}, test="lvl2")
