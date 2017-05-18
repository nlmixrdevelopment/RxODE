context("Michelis Menton test (long lines)")
## Michelis Menton test
rxPermissive({

    mm <- RxODE({
        C2 = centr / V;
        d/dt(centr)  = -(VM*C2)/(KM+C2);
    })

    mypar3 <- function(lVM, lKM, lV )
    {
        VM = exp(theta[1] + eta[1])
        KM = exp(theta[2] + eta[2])
        V  = exp(theta[3] + eta[3])
    }

    pred <- function() C2

    ## test_that("Functions outside of RxODE/global/nlmixr raise errors", {
    ##     expect_error(rxSymPySetupPred(mm, pred, par, err=function(){err ~ prop(0.1)}, grad=TRUE));
    ## })

    focei.mm.mod2 <- rxSymPySetupPred(mm, pred, mypar3, err=function(){err ~ prop(0.1)}, grad=TRUE, logify=TRUE);

    ## Fixme: should the lines be split?  Logify is one approach, but perhaps just a split...?
    test_that("long lines are handled...",{
        expect_equal(class(focei.mm.mod2), "rxFocei");
    })

    ## Now test solved focei capability.
    sol.1c.ka <- RxODE({
        ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
        sf=solvedC(t, 1, 1, 1, V, CL, 0, 0, 0, 0, KA, 0);
    })

    prd <- function(){
        return(sf);
    }

    err <- function(f){
        return(f ^ 2* theta[4] ^ 2); ## Theta 4 is residual sd for proportional error.
    }

    mypar1 = function ()
    {
        CL = exp(THETA[1] + ETA[1])
        V = exp(THETA[2] + ETA[2])
    }

    rxSymPySetupPred(sol.1c.ka, prd, mypar1, err, run.internal=TRUE)
}, silent=TRUE)
