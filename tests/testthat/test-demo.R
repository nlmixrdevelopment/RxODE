rxodeTest(
  {
    test_that("Demo runs successfully", {
      expect_error(suppressWarnings({
        demo("demo1", "RxODE", ask = FALSE, echo = FALSE)
      }), NA)
    })
  },
  test = "demo"
)
