require(RxODE)
context("Test loading, and unloading of models")
rxPermissive({

    rigid.txt <- "
y1(0)    = 1
y2(0)    = 0
y3(0)    = 0.9
a1       = -2
a2       = 1.25
a3       = -0.5
d/dt(y1) = a1*y2*y3
d/dt(y2) = a2*y1*y3
d/dt(y3) = a3*y1*y2
";

    rigid <- RxODE(rigid.txt);

    test_that("rxDll works for each object type",{
        expect_equal(rxDll(rigid.txt),rxDll(rigid));
    })

    test_that("loading and unloading works.",{
        dll <- rxDll(rigid);
        expect_true(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
        rxUnload(rigid)
        expect_false(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
        rxLoad(rigid);
        expect_true(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
        rxDelete(rigid);
        expect_false(rxDllLoaded(rigid));
        expect_false(file.exists(dll));
        options(RxODE.compile.on.load = TRUE)
        rxLoad(rigid);
        expect_true(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
        rxDelete(rigid);
        options(RxODE.compile.on.load = FALSE)
        expect_error(rxLoad(rigid),regexp="error loading dll file");
        options(RxODE.compile.on.load = TRUE)
        expect_false(rxDllLoaded(rigid));
        expect_false(file.exists(dll));
        ## Test $ syntax
        options(RxODE.compile.on.load = TRUE)
        rigid$compile();
        expect_true(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
        rigid$delete();
        expect_false(rxDllLoaded(rigid));
        expect_false(file.exists(dll));
        rigid$dynLoad();
        expect_true(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
        rigid$dynUnload();
        expect_false(rxDllLoaded(rigid));
        expect_true(file.exists(dll));
    })

}, silent=TRUE);
