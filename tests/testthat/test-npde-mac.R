rxodeTest({

  test_that("npde simulation works on mac nlmixr #460", {
    skip_if(!file.exists(test_path("si.qs")))
    si <- qs::qread(test_path("si.qs"))
    si$object <- RxODE(si$object)
    set.seed(1009)
    solve <- expect_error(suppressWarnings(do.call(RxODE::rxSolve, si)), NA)

  })


}, test="lvl2")
