context("Test Factorial operator");
rxPermissive({
    rxClean()

    transTo <- function(model, syntax, match=TRUE){
        mv <- rxModelVars(model);
        if (match){
            test_that(sprintf("%s includes %s", model, syntax),
                      expect_true(regexpr(syntax,RxODE:::.rxGetParseModel(), fixed=TRUE) != -1));
        } else{
            test_that(sprintf("%s dose not include %s", model, syntax),
                      expect_false(regexpr(syntax,RxODE:::.rxGetParseModel(), fixed=TRUE) != -1));
        }
    }

    transTo("d/dt(cmt)= factorial(1+fac)-cmt*ka", "factorial(1+fac)")
    transTo("d/dt(cmt)= lgamma(1+fac)-cmt*ka", "lgamma(1+fac)")
    transTo("d/dt(cmt)= gamma(1+fac)-cmt*ka", "lgammafn(1+fac)")
    transTo("d/dt(cmt)= lfactorial(1+fac)-cmt*ka", "lgamma1p(1+fac)")

}, silent=TRUE)
