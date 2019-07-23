context("Conditional statements")
rxPermissive({

    ## else if is actually already supported...
    test_that("else if", {

        m <- RxODE({
            if (cnd <= 1){
                a = 1.0
            } else if (cnd <= 2){
                a = 2.0
            } else if (cnd <= 3){
                a = 3.
            } else {
                a = 100;
            }
            tmp = cnd
        })

        ## The prefered syntax is only if / else but it still works...
        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

    })

    test_that("ifelse", {

        m <- RxODE({
            a = ifelse(cnd <= 1, 1.0, ifelse(cnd <= 2, 2, ifelse(cnd <= 3, 3, 100)))
            tmp = cnd
        })

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)
    })

    test_that("embedded logical expressions", {

        m <- RxODE({
            a = (cnd == 1) * 1.0 + (cnd == 2) * 2 + (cnd == 3) * 3
            tmp = cnd
        })

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 0)

    })

    test_that("ifelse with assignments", {

        m <- RxODE({
            ifelse(cnd <= 1, a = 1.0, a = 2.0);
            tmp = cnd
        })

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        m <- RxODE({
            ifelse(cnd <= 1, a <- 1.0, a <- 2.0);
            tmp = cnd
        })

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)
    })
    ##
    context("Pruning checks");
    test_that("prune checks", {

        tmp <- "C2=centr/V;\nC3=peri/V2;\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nC4=CMT;\nif(CMT==1){\nprd=depot;\n}\nif(CMT==2){\nprd=centr;\n}\nif(CMT==3){\nprd=peri;\n}\n"
        expect_equal(rxPrune(tmp), "C2=centr/V\nC3=peri/V2\nd/dt(depot)=- KA*depot\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3\nd/dt(peri)=Q*C2-Q*C3\nC4=CMT\nprd=(CMT==1)*(depot)\nprd=(CMT==2)*(centr)+(1-((CMT==2)))*(prd)\nprd=(CMT==3)*(peri)+(1-((CMT==3)))*(prd)")

        ## Advanced context pruining:
        m <- RxODE({
            if (cnd <= 1){
                a = 1.0
            } else if (cnd <= 2){
                a = 2.0
            } else if (cnd <= 3){
                a = 3.
            } else {
                a = 100;
            }
            tmp = cnd
        })

        m <- RxODE(rxPrune(m));

        ## The prefered syntax is only if / else but it still works...
        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

        m <- RxODE(rxOptExpr(rxNorm(m)))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

        m <- RxODE({
            a = 100;
            if (cnd <= 1){
                a = 1.0
            }
            if (cnd > 1 && cnd <= 2){
                a = 2.0
            }
            if (cnd > 2 && cnd <= 3){
                a = 3.
            }
            tmp = cnd
        })

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

        m <- RxODE(rxPrune(m))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

        m <- RxODE(rxOptExpr(rxNorm(m)))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

        m <- RxODE({
            a = ifelse(cnd <= 1, 1.0, ifelse(cnd <= 2, 2, ifelse(cnd <= 3, 3, 100)))
            tmp = cnd
        })

        m <- RxODE(rxPrune(m))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

        m <- RxODE(rxOptExpr(rxNorm(m)))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)


        m <- RxODE({
            ifelse(cnd <= 1, a = 1.0, a = 2.0);
            tmp = cnd
        })

        m <- RxODE(rxPrune(m))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        m <- RxODE(rxOptExpr(rxNorm(m)))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        m <- RxODE({
            ifelse(cnd <= 1, a <- 1.0, a <- 2.0);
            tmp = cnd
        })

        m <- RxODE(rxPrune(m))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        m <- RxODE(rxOptExpr(rxNorm(m)))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        m <- RxODE({
            a = (cnd == 1) * 1.0 + (cnd == 2) * 2 + (cnd == 3) * 3
            tmp = cnd
        })

        m <- RxODE(rxPrune(m))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 0)

        m <- RxODE(rxOptExpr(rxNorm(m)))

        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 0)

    })

}, cran=TRUE)
