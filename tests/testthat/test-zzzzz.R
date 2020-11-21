rxPermissive({
  context("Cleanly unloads all dlls")
  test_that("Cleanly unloads all dlls", {
    expect_equal(rxUnloadAll(), list())
  })
}, test="focei")
