context("Test RxODE THETA/ETA support")
rxPermissive({

    rigid <- RxODE({
        y1(0)    = 1
        y2(0)    = 0
        y3(0)    = 0.9
        a1       = theta[1]
        a2       = theta[2]
        a3       = theta[3] + eta[1]
        d/dt(y1) = a1*y2*y3
        d/dt(y2) = a2*y1*y3
        d/dt(y3) = a3*y1*y2
    })

    et <- eventTable();
    et$add.sampling(seq(0,20,by=0.01))

    out <- solve(rigid,et, theta=c(-2, 1.25, -0.5), eta=c(0))

    test_that("Test rigid body example",{
        expect_equal(digest(round(as.data.frame(out),3)),
                     "14ac439006ceb32455ba00c7817f9163")
    })

    out <- solve(rigid,et, theta=c(-2, 1.25, -0.5), eta=c(1))

    test_that("Test rigid body example",{
        expect_equal(digest(round(as.data.frame(out),3)),
                     "c7cffaa650a47e2b28e4cba99c603dde")
    })

    m1 <- RxODE({
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V = exp(THETA[3] + ETA[2])
        d/dt(depot) = -KA*depot;
        d/dt(centr) = KA*depot - CL / V*centr;
    });

    et <- eventTable() %>% add.dosing(dose=30000) %>%
        add.sampling(c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 20, 24, 36, 48, 60, 71.99))

    pk <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V = exp(THETA[3] + ETA[2])
    }

    pred <- function(){
        return(cntr);
    }

    m2 <- RxODE({
        d/dt(depot) = -KA*depot;
        d/dt(centr) = KA*depot - CL / V*centr;
    })

    test_that("Error when pred dosen't depend on state varaibles", {
        expect_error(rxSymPySetupPred(m2, pred, pk));
    })

    pred <- function(){
        return(centr);
    }

    err <- function(f){
        return(f * theta[4]); ## Theta 4 is residual sd for proportional error.
    }

    m2a <- rxSymPySetupPred(m2, pred, pk, err)

    ## test_that("theta/eta solve works", {
    ##     expect_equal(suppressWarnings(digest(m2a %>% solve(et, theta=c(2, 1.6, 4.5, 0.1), eta=c(0.01, -0.01)) %>% as.data.frame %>% round(3))),
    ##                  "13baa6cabb14367e1328149d3c14ebf7");
    ## })

    ## m2a <- rxSymPySetupPred(m2, pred, pk, err)

    ## m2a %>% solve(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01))

}, silent=TRUE)
