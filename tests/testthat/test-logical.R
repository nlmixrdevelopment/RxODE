context("Logical expressons test")
rxPermissive({
    rxClean()

    transTo <- function(model, syntax, match=TRUE){
        mod <- RxODE(model);
        if (match){
            test_that(sprintf("%s includes %s", model, syntax),
                      expect_true(regexpr(syntax,rxModelVars(mod)$model["parseModel"], fixed=TRUE) != -1));
        } else{
            test_that(sprintf("%s dose not include %s", model, syntax),
                      expect_false(regexpr(syntax,rxModelVars(mod)$model["parseModel"], fixed=TRUE) != -1));
        }
    }

    transTo("x=1;if (t <> 2){x=0}","(t !=2)")
    transTo("x=1;if (t ~= 2){x=0}","(t !=2)")
    transTo("x=1;if (t /= 2){x=0}","(t !=2)")
    transTo("x=1;if (t != 0 & t != 1){x=0}","&&")
    transTo("x=1;if ((t == 0) | (t == 1)){x=0}","||")
    transTo("x=1;if ((t == 0) & !(t == 1)){x=0}","&&")
    transTo("x=1;if ((t == 0) & !(t == 1)){x=0}","!(")

    rxClean()
}, silent=TRUE);
