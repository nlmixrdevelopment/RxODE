rxPermissive({
  test_that("run examples", {
    for (v in c("rxProgress", "rxExpandGrid", "rxFun", "genShinyApp.template", "forderForceBase")) {
      expect_error(eval(parse(text=paste0("example(", v, ")"))), NA)
    }
  })
},  test = "lvl2")
