for (d in seq(1, ifelse(identical(Sys.getenv("RxODE_VALIDATION_FULL"), "true"), 4, rxSymInvCholN()))){
    if (identical(Sys.getenv("RxODE_VALIDATION_FULL"), "true")){
        dgs <- c("sqrt", "identity")
    } else {
        dgs <- c("log")
    }
    for (dg in dgs){
        test_that("omega chol", {
            context(sprintf("Omega Cholesky %sx%s, %s", d, d, dg));
            ## Creating covariance matrix
            tmp <- matrix(rnorm(d^2), d, d)
            mcov <- tcrossprod(tmp, tmp)
            v <- rxSymInvCholCreate(mcov, dg)
            expect_equal(v$omega,mcov)
            expect_equal(v$omegaInv, solve(mcov))
            expect_equal(v$chol.omegaInv, chol(solve(mcov)))
            expect_equal(v$chol.omega, chol(mcov))
            expect_equal(v$log.det.OMGAinv.5, 0.5 * log(det(solve(mcov))))
            expect_equal(length(v$d.omegaInv), v$ntheta)
            expect_equal(length(v$d.D.omegaInv), v$ntheta)
        })
    }
}
