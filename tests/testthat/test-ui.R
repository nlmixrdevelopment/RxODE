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
        linCmt() ~ add(add.sd)
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
                   "        v <- exp(tv + eta.v)", "        linCmt() ~ add(add.sd)", 
                   "    })", "}"))


  })

}, test="lvl2")
