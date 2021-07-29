rxodeTest(
  {
    test_that("Demo runs successfully", {
      skip_if_not(dir.exists(file.path(system.file(package = "RxODE"), "demo")), "demo not installed")
      expect_error(suppressWarnings({
        demo("demo1", "RxODE", ask = FALSE, echo = FALSE)
      }), NA)
    })
  },
  test = "demo"
)
