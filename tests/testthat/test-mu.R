rxMuRef(RxODE({
  ka <- exp(tka + eta.ka)
  cl <- exp(tcl + eta.cl)
  v <- exp(tv + eta.v)
  d/dt(depot) = -ka * depot
  d/dt(center) = ka * depot - cl/v * center
  cp = center/v
  ## cp ~ add(add.sd)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))


## Test a duplicated eta; It shouldn't be counted as mu-referenced
rxMuRef(RxODE({
  EmaxA <- exp(t.EmaxA + eta.emax)
  EmaxB <- exp(t.EmaxB + eta.emax)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))

## Composite expressions should be extracted to their own lines
rxMuRef(RxODE({
  ratio <- exp(t.EmaxA + eta.emaxA) / exp(t.EmaxB + eta.emaxB)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))

## This composite RxODE model should be extracted to something like above
rxMuRef(RxODE({
  d/dt(depot) = -exp(tka + eta.ka) * depot
  d/dt(center) = exp(tka + eta.ka) * depot - exp(tcl + eta.cl)/exp(tv + eta.v) * center
  cp = center/exp(tv + eta.v)
  #cp ~ add(add.sd)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))


## This should be expanded to eta.ka mu-referenced variables
rxMuRef(RxODE({
  ka <- tka * exp(eta.ka)
  cl <- tcl * exp(eta.cl)
  v <- tv * exp(eta.v)
  d/dt(depot) = -ka * depot
  d/dt(center) = ka * depot - cl/v * center
  cp = center/v
  ## cp ~ add(add.sd)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))
