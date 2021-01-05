rxodeTest({

  test_that("npde simulation works on mac nlmixr #460", {
    tempdir <- tempdir()
    unzip("si.zip", files="si.qs", exdir=tempdir)
    si <- qs::qread(file.path(tempdir, "si.qs"))
    si$object <- RxODE(si$object)
    set.seed(1009)
    solve <- expect_error(suppressWarnings(do.call(RxODE::rxSolve, si)), NA)

  })


}, test="lvl2")
