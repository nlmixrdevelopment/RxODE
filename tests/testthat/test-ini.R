library(RxODE);
context("Test Inis");

out <- RxODE("ini = 1; fun_ini = 2; fun = 4; addit = ini + fun_ini + pi")

test_that("Initial constants are correct",{
    expect_equal(rxInits(out)["pi"],c(pi=pi));
    expect_equal(rxInits(out)["ini"],c(ini=1));
    expect_equal(rxInits(out)["fun_ini"],c(fun_ini=2));
    expect_equal(rxInits(out)["fun"],c(fun=4));
})

rxClean();
