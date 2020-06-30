rxPermissive({

  f1 <- RxODE({
    a = 1
    if (id == "matt"){
      a = 2
    }
  })

  f2 <- RxODE({
    a = 1
    if ("matt" == ID){
      a = 2
    }
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
