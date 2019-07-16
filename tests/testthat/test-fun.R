context("User function tests")
rxPermissive({
    test_that("user function tests", {

        expect_error(RxODE("a=fun(d,b,c)"))
        expect_error(rxFromSE("Derivative(fun(a,b,c),a)"))

        fun <- "
 double fun(double a, double b, double c) {
   return a*a+b*a+c;
 }
"
        expect_error(rxFun(list("fun"), c("a", "b", "c"), fun))
        expect_error(rxFun("fun", c("a", "b", "c")))
        rxFun("fun", c("a", "b", "c"), fun);
        expect_error(rxFun("fun", c("a", "b", "c"), fun))

        tmp <- RxODE("a=fun(d,b,c)")

        tmp <- rxSolve(tmp, c(c=2, b=4, d=8), et(1))
        expect_equal(tmp$a, 8 ^ 2 + 4 * 8 + 2)

        expect_equal(rxToSE("fun(a,b,c)"), "fun(a,b,c)")
        expect_equal(rxFromSE("fun(a,b,c)"), "fun(a,b,c)")

        expect_error(rxToSE("fun(a)"))
        expect_error(rxFromSE("fun(a)"))

        expect_equal(rxFromSE("Derivative(fun(a,b,c),a)"),
                     "(fun(a-0.5*rx_c_delta,b,c)-fun(a+0.5*rx_c_delta,b,c))/rx_c_delta")
        expect_equal(rxFromSE("Derivative(fun(a,b,c),a)", unknownDerivatives = "forward"),
                     "(fun((a)+rx_f_delta,b,c)-fun(a,b,c))/rx_f_delta")
        expect_error(rxFromSE("Derivative(fun(a,b,c),a)", unknownDerivatives = "error"))

        ## now add derivative table
        rxD("fun", list(function(a,b,c){
                       paste0("2*",a,"+",b);
                   },
                   function(a,b,c){
                       return(a)
                   },
                   function(a,b,c){
                       return("0.0")
                   }))

        expect_equal(rxFromSE("Derivative(fun(a1,b1,c1),a1)"), "2*a1+b1")
        expect_equal(rxFromSE("Derivative(fun(a1,b1,c1),b1)"), "a1")
        expect_equal(rxFromSE("Derivative(fun(a1,b1,c1),c1)"), "0.0")

        expect_warning(rxD("fun", list(function(a,b,c){
                                      paste0("2*",a,"+",b);
                                  })))

        expect_equal(rxFromSE("Derivative(fun(a1,b1,c1),a1)"), "2*a1+b1")
        expect_equal(rxFromSE("Derivative(fun(a,b,c),b)"),
                     "(fun(a,b-0.5*rx_c_delta,c)-fun(a,b+0.5*rx_c_delta,c))/rx_c_delta")

        expect_warning(rxD("fun", list(NULL, function(a,b,c){
                                      paste0(a);
                                  })))

        expect_equal(rxFromSE("Derivative(fun(a,b,c),b)"),
                     "a")

        expect_equal(rxFromSE("Derivative(fun(a,b,c),a)"),
                     "(fun(a-0.5*rx_c_delta,b,c)-fun(a+0.5*rx_c_delta,b,c))/rx_c_delta")

        expect_error(rxD("fun", "matt"))
        expect_error(rxD("fun", list()))
        expect_error(rxD("fun", list(NULL, "a", function(x){x})))

        ## Errors do not replace derivative table
        expect_equal(rxFromSE("Derivative(fun(a,b,c),b)"),
                     "a")
        expect_equal(rxFromSE("Derivative(fun(a,b,c),a)"),
                     "(fun(a-0.5*rx_c_delta,b,c)-fun(a+0.5*rx_c_delta,b,c))/rx_c_delta")

    expect_error(rxRmFun(list("a")))
    expect_error(rxRmFun("a"))

    rxRmFun("fun")
    rxClean()
    expect_error(RxODE("a=fun(d,b,c)"))
    expect_error(rxFromSE("Derivative(fun(a,b,c),a)"))
})
})
