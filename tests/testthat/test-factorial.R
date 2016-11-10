context("Test Factorial operator");

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

transTo("d/dt(cmt)= 3.0!-cmt*ka", "exp(lgamma1p(3.0))")
transTo("d/dt(cmt)= log(3.0!)-cmt*ka", "lgamma1p(3.0)")
transTo("d/dt(cmt)= log(3.0!)-cmt*ka", "exp(lgamma1p(3.0))", FALSE)

transTo("d/dt(cmt)= fac!-cmt*ka", "exp(lgamma1p(fac))")
transTo("d/dt(cmt)= log(fac!)-cmt*ka", "lgamma1p(fac)")
transTo("d/dt(cmt)= log(fac!)-cmt*ka", "exp(lgamma1p(fac))", FALSE)

transTo("d/dt(cmt)= (1+fac)!-cmt*ka", "exp(lgamma1p(1+fac))")
transTo("d/dt(cmt)= log((1+fac)!)-cmt*ka", "lgamma1p(1+fac)")
transTo("d/dt(cmt)= log((1+fac)!)-cmt*ka", "exp(lgamma1p(1+fac))", FALSE)
