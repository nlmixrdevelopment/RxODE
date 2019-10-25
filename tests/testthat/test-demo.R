rxPermissive({
test_that("Demo runs successfully", {
    tmp <- suppressWarnings({demo("demo1", "RxODE", ask=FALSE, echo=FALSE)});
    expect_true(any(names(tmp) == "value"));
})
})
