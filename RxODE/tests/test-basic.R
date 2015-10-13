# Test basic RxODE() invocation (including errors)

library("RxODE")
library("tools")   # use assertError()

test.dir <- tempfile("Rx_base-")
dir.create(test.dir)

ode <- 'd/dt(y) = r * y * (1.0 - y/K);'

# ODE inside a string 
m1 <- RxODE(model = ode, modName = "m1", do.compile=FALSE)

# ODE in a text file
fn <- file.path(test.dir, "exam3.1.txt")
writeLines(ode, fn)
m2 <- RxODE(filename = fn, modName = "f1", do.compile=FALSE)

# Error: arguments model= and filename= are mutually exclusive.
tools::assertError(
   RxODE(model=ode, filename=fn, do.compile=FALSE)
)

unlink(test.dir, recursive = TRUE)
