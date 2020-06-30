rxPermissive({

  f1 <- RxODE({
    a = 1
    if (id == "matt"){
      a = 2
    }
  })

  df <- data.frame(ID=c("matt", "not"),
                   TIME=c(0, 1))

  ref <- rxSolve(f1, df)

  test_that("id==matt #1",{
    expect_equal(ref$a, c(2, 1))
  })

  df <- data.frame(ID=c("not", "matt"),
                   TIME=c(0, 1))

  ref <- rxSolve(f1, df)

  test_that("id==matt #2",{
    expect_equal(ref$a, c(1, 2))
  })

  f2 <- RxODE({
    a = 1
    if ("matt" == ID){
      a = 2
    }
  })

  df <- data.frame(ID=c("matt", "not"),
                   TIME=c(0, 1))

  ref <- rxSolve(f2, df)

  test_that("id==matt #3",{
    expect_equal(ref$a, c(2, 1))
  })

  df <- data.frame(ID=c("not", "matt"),
                   TIME=c(0, 1))

  ref <- rxSolve(f2, df)

  test_that("id==mat #4",{
    expect_equal(ref$a, c(1, 2))
  })

  f1 <- RxODE({
    a = 1
    if (id != "matt"){
      a = 2
    }
  })

  df <- data.frame(ID=c("matt", "not"),
                   TIME=c(0, 1))

  ref <- rxSolve(f1, df)

  test_that("id!=matt #1",{
    expect_equal(ref$a, c(1, 2))
  })

  df <- data.frame(ID=c("not", "matt"),
                   TIME=c(0, 1))

  ref <- rxSolve(f1, df)

  test_that("id!=matt #2",{
    expect_equal(ref$a, c(2, 1))
  })

  f1 <- RxODE({
    a = 1
    if (id == "100"){
      a = 2
    }
  })

  df <- data.frame(ID=c(100, 200),
                   TIME=c(0, 1))

  ref <- rxSolve(f1, df)

  test_that("id2",{
    expect_equal(ref$a, c(2, 1))
  })

  expect_error(RxODE({
    a = 1
    if (id == 100){
      a = 2
    }
  }))


  f1 <- RxODE({
    a = 1 + (id == "100") * 2
  })

  ref <- rxSolve(f1, df)

  test_that("id3",{
    expect_equal(ref$a, c(3, 1))
  })

  f1 <- RxODE({
    a = 1 + (id == '100') * 2
  })

  f3 <- RxODE({
    if (cmt == "depot"){
      a = 2
    }
  })

  f4 <- RxODE({
    a = 1
    if ("depot" == cmt){
      a = 2
    }
  })

})
