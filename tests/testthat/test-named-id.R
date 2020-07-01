rxPermissive({

  context("test id==\"item\"")

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

  ref <- rxSolve(f1, df)

  test_that("id3",{
    expect_equal(ref$a, c(3, 1))
  })

  context("test cmt==\"item\"/no.")

  f3 <- RxODE({
    a = 3
    if (cmt == "depot"){
      a = 2
    }
  })

  df <- data.frame(ID=c("not", "matt"),
                   TIME=c(0, 1),
                   cmt=c("depot", "central"))

  tmp <- rxSolve(f3, df)

  test_that("cmt==\"depot\"", {
    expect_equal(tmp$a, c(2, 3))
  })

  f3 <- RxODE({
    a = 3
    if (cmt == 1){
      a = 2
    }
  })

  tmp <- rxSolve(f3, df)

  test_that("cmt==1", {
    expect_equal(tmp$a, c(2, 3))
  })

  f3 <- RxODE({
    a = 3
    if (cmt != "depot"){
      a = 2
    }
  })

  tmp <- rxSolve(f3, df)

  test_that("cmt!=\"depot\"", {
    expect_equal(tmp$a, c(3, 2))
  })

  f3 <- RxODE({
    a = 3
    if (cmt != 1){
      a = 2
    }
  })

  test_that("cmt!=\"depot\"", {
    expect_equal(tmp$a, c(3, 2))
  })

  f3 <- RxODE({
    a = 1
    if ("depot" == cmt){
      a = 2
    }
  })

  tmp <- rxSolve(f3, df)

  test_that("cmt==\"depot\"", {
    expect_equal(tmp$a, c(2, 1))
  })

  f3 <- RxODE({
    a = 1
    if ("depot" != cmt){
      a = 2
    }
  })

  tmp <- rxSolve(f3, df)

  test_that("cmt==\"depot\"", {
    expect_equal(tmp$a, c(1, 2))
  })

  test_that("translation to and from SE", {
    RxODE:::.clearSEstr()
    a <- rxToSE("id==\"matt\"")
    expect_equal(a, "rxEq(id,rxQ1)")
    b <- rxFromSE(a)
    expect_equal(b, "(id==\"matt\")")
  })

})
