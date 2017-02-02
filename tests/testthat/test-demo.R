test_that("Demo runs successfully", {
    tmp <- demo("demo1", "RxODE", ask=FALSE, echo=FALSE);
    expect_true(tmp$value);
})
