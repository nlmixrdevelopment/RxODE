context("Removind dlls")
rxPermissive({
    options(RxODE.delete.unnamed=TRUE);
    rxSyncOptions();
    ode <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
    })
    dll <- rxDll(ode);
    test_that("dll exists", {
        expect_true(file.exists(dll));
    })
    rm(ode);
    gc();
    test_that("dll removed", {
        expect_false(file.exists(dll));
    })
    ## FIXME -- Parallelization removes dll?  I don't think so since it is loaded elsewhere.  Can we test this?
    options(RxODE.delete.unnamed=FALSE);
    rxSyncOptions();
    ode <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
    })
    dll <- rxDll(ode);
    test_that("dll exists #2", {
        expect_true(file.exists(dll));
    })
    rm(ode);
    gc();
    test_that("dll exists #3", {
        expect_true(file.exists(dll));
    })
    ## Note it won't delete the dll if not created when the RxODE.delete.unnamed=FALSE
    try(dyn.unload(dll), silent = TRUE)
    unlink(dll);
    options(RxODE.delete.unnamed=TRUE);
    rxSyncOptions();
}, silent=TRUE)
