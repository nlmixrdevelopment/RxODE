rxMuRef(quote({
  ka <- exp(tka + eta.ka)
  cl <- exp(tcl + eta.cl)
  v <- exp(tv + eta.v)
  d/dt(depot) = -ka * depot
  d/dt(center) = ka * depot - cl/v * center
  cp = center/v
  cp ~ add(add.sd)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))



rxMuRef(quote({
  d/dt(depot) = -exp(tka + eta.ka) * depot
  d/dt(center) = exp(tka + eta.ka) * depot - exp(tcl + eta.cl)/exp(tv + eta.v) * center
  cp = center/exp(tv + eta.v)
  cp ~ add(add.sd)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))


rxMuRef(quote({
  d/dt(depot) = -exp(tka + eta.ka) * depot
  d/dt(center) = exp(eta.ka + tka) * depot - exp(tcl + eta.cl)/exp(tv + v.wt * wt + eta.v) * center
  cp = center/exp(v.wt * wt + tv + eta.v)
  cp ~ add(add.sd)
}), theta=c("tka", "tcl", "tv", "add.sd"),
eta=c("eta.ka", "eta.cl", "eta.v"))
