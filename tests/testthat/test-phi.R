rxPermissive({

  context("phi")
  expect_equal(phi(1:3), pnorm(1:3))
  expect_equal(phi(as.double(1:3)), pnorm(as.double(1:3)))
  
}, test="cran")
