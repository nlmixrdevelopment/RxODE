[![Build Status](https://travis-ci.org/mattfidler/RxODE.svg?branch=master)](https://travis-ci.org/mattfidler/RxODE)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mattfidler/RxODE?branch=master&svg=true)](https://ci.appveyor.com/project/mattfidler/RxODE)
[![codecov.io](https://codecov.io/github/mattfidler/RxODE/coverage.svg?branch=master)](https://codecov.io/github/mattfidler/RxODE?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/RxODE)](https://cran.r-project.org/package=RxODE)




## RxODE: A tool for performing simulations from Ordinary Differential Equation (ODE) models, with applications for pharmacometrics
***  

##### Authors: Melissa Hallow, Wenping Wang, and Matthew L. Fidler

***

`RxODE` installation under R for Windows
========================================

These notes briefly describe steps to properly install `RxODE` and to
ensure `Rtools` (https://cran.r-project.org/bin/windows/Rtools/) are properly 
configured to avoid compilation issues during the use of `RxODE`. 

In a nutshell, installing `RxODE` is very straight forwad, but installing
and configuring `Rtools` is a bit more delicate and you need to 
carefully follow the instructions in the "R Installation and Adminstration" 
manual, in particular Section 6.3, and Appendix D "The Windows Toolset". 
We point out a couple of details worth extra attention.  Please read on.

Steps:
------

1. Install the appropriate `Rtools` for your R for Windows version,
   e.g., `Rtools` 3.2 for R versions 3.1.x through 3.2.x (for full details
   see http://cran.r-project.org/bin/windows/Rtools/). A couple of 
   important details:

   * When installing `Rtools`, in the "Select Components" dialog box, 
     you may select the default "Package authoring installation".

   * In the "Select Additional Tasks" dialog window, check the
     option "Edit the system PATH".  This is important to be able to
     locate the C, Fortran compilers and other tools needed during 
     the use of `RxODE`.  

   * A simple way to test whether `Rtools` was properly installed is
     to compile the `hello.c` program.  Simply open a new MSDOS command 
     window, create a text file `hello.c` and compile it as follows:
   
     ```
     C:\hello> type hello.c
     #include<stdio.h>
     
     void main(int argc, char **argv)
     {
         printf("Hello World!\n");
     }

     C:\hello> gcc -o hello hello.c

     C:\hello> .\hello
     Hello World!
     ```

     If you get the error `gcc: error: CreateProcess: No such file or directory`     then you know `Rtools` was not properly installed, in particular,
     it did not update your system `PATH` variable.

2.  Obtain the `RxODE` package, either from github or CRAN.  The 
    installation requires use of the gcc compiler, so you'll know if Step 1 
    was successfully executed.

    * CRAN. Use the usual method for installing pacakges from CRAN.

    * GitHub. First install the `devtools` package (if needed) and 
      then `RxODE` from GitHub.  You may want to avoid using a library 
      folder that has spaces in its name (see question 4.1 in the 
      "R for Windows FAQ" and the pointers therein).  As of `RxODE`
      version 0.5-1, we've been able to test installations on folder with 
      spaces in their name, but you may want to be on the safe side.
      
      ``` 
      install.packages("devtools")
      library("devtools", lib = "C:/Rlib")
      install_github("hallowkm/RxODE/RxODE")
      ```

3. Test the `RxODE` installation:

    ``` 
    library("RxODE", lib = "C:/Rlib")
    demo("demo1","RxODE")
    ```

If the demo runs without error, click on the plot window and see if a 
new plot comes up each time. If so, `RxODE` has been installed correctly.

See `browseVignettes("RxODE")` for an extended example on using 
`RxODE` for simulations.

#### Introduction
`RxODE` is an R package that facilitates simulation with ODE models in
R. It is designed with pharmacometrics models in mind, but can be
applied more generally to any ODE model.

***
#### Description of RxODE illustrated through an example
The model equations are specified through a text string in R. Both
differential and algebraic equations are permitted. Differential
equations are specified by `d/dt(var_name) = `. Each
equation is separated by a semicolon.

```r
ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"
```
To load `RxODE` package and compile the model: 

```r
library(RxODE)
work <- tempfile("Rx_intro-")
mod1 <- RxODE(model = ode, modName = "mod1", wd = work)
```

A typical pharmacokinetics-pharmacodynamics (PKPD) model can be
plotted in `RxODE`. This model, as shown in the figure below:

```r
plot(mod1);
```

```
##                              C2                              C3 
##                      "centr/V2"                       "peri/V3" 
##                     d/dt(depot)                     d/dt(centr) 
##                     "-KA*depot"      "KA*depot-CL*C2-Q*C2+Q*C3" 
##                      d/dt(peri)                       d/dt(eff) 
##                     "Q*C2-Q*C3" "Kin-Kout*(1-C2/(EC50+C2))*eff"
```

![plot of chunk unnamed-chunk-3](vignettes/figure/unnamed-chunk-3-1.png)

Sometimes the size of the boxes may need to be adjusted, you can do
this by adjusting the `size` argument:

```r
plot(mod1,size=40);
```

```
##                              C2                              C3 
##                      "centr/V2"                       "peri/V3" 
##                     d/dt(depot)                     d/dt(centr) 
##                     "-KA*depot"      "KA*depot-CL*C2-Q*C2+Q*C3" 
##                      d/dt(peri)                       d/dt(eff) 
##                     "Q*C2-Q*C3" "Kin-Kout*(1-C2/(EC50+C2))*eff"
```

![plot of chunk unnamed-chunk-4](vignettes/figure/unnamed-chunk-4-1.png)

Model parameters can be defined as named vectors. Names of parameters in
the vector must be a superset of parameters in the ODE model, and the
order of parameters within the vector is not important. 

```r
theta <- 
   c(KA=2.94E-01, CL=1.86E+01, V2=4.02E+01, # central 
     Q=1.05E+01,  V3=2.97E+02,              # peripheral
     Kin=1, Kout=1, EC50=200)               # effects  
```

Initial conditions (ICs) are defined through a vector as well. If the
vector is not named, the number of ICs must equal exactly the number of
ODEs in the model, and the order must be the same as the order in
which the ODEs are listed in the model. 

```r
inits <- c(0, 0, 0, 1)    
```

When elements are named, missing elements are added and set to
zero. Also when named, the order of initilizations does not matter.
Therefore the following code is an equvalent initialization as the
code above.

```r
inits <- c(eff=1);
```


`RxODE` provides a simple and very flexible way to specify dosing and
sampling through functions that generate an event table. First, an
empty event table is generated through the "eventTable()" function:

```r
ev <- eventTable(amount.units='mg', time.units='hours')
```

Next, use the `add.dosing()` and `add.sampling()` functions of the
`EventTable` object to specify the dosing (amounts, frequency and/or
times, etc.) and observation times at which to sample the state of the
system.  These functions can be called multiple times to specify more
complex dosing or sampling regiments.  Here, these functions are used
to specify 10mg BID dosing for 5 days, followed by 20mg QD dosing for
5 days:

```r
ev$add.dosing(dose=10000, nbr.doses=10, dosing.interval=12)
ev$add.dosing(dose=20000, nbr.doses=5, start.time=120, dosing.interval=24)
ev$add.sampling(0:240)
```

The functions `get.dosing()` and `get.sampling()` can be used to
retrieve information from the event table.  

```r
head(ev$get.dosing())
```

```
##   time evid   amt
## 1    0  101 10000
## 2   12  101 10000
## 3   24  101 10000
## 4   36  101 10000
## 5   48  101 10000
## 6   60  101 10000
```

```r
head(ev$get.sampling())
```

```
##    time evid amt
## 16    0    0  NA
## 17    1    0  NA
## 18    2    0  NA
## 19    3    0  NA
## 20    4    0  NA
## 21    5    0  NA
```

The simulation can now be run by calling the model object's run
function. Simulation results for all variables in the model are stored
in the output matrix x. 

```r
x <- mod1$solve(theta, ev, inits)
```

```
## Warning in rxInits(object, inits, rxState(object), 0): Assiged depot,
## centr, peri to 0.
```

```r
head(x)
```

```
##      time     depot    centr      peri      eff       C2        C3
## [1,]    0 10000.000    0.000    0.0000 1.000000  0.00000 0.0000000
## [2,]    1  7452.765 1783.897  273.1895 1.084664 44.37555 0.9198298
## [3,]    2  5554.370 2206.295  793.8758 1.180825 54.88296 2.6729825
## [4,]    3  4139.542 2086.518 1323.5783 1.228914 51.90343 4.4564927
## [5,]    4  3085.103 1788.795 1776.2702 1.234610 44.49738 5.9807076
## [6,]    5  2299.255 1466.670 2131.7169 1.214742 36.48434 7.1774981
```

This can also be solved by the `predict()` or `solve()` methods:

```r
x <- predict(mod1,theta, ev, inits)
```

```
## Warning in rxInits(object, inits, rxState(object), 0): Assiged depot,
## centr, peri to 0.
```

```r
head(x)
```

```
##      time     depot    centr      peri      eff       C2        C3
## [1,]    0 10000.000    0.000    0.0000 1.000000  0.00000 0.0000000
## [2,]    1  7452.765 1783.897  273.1895 1.084664 44.37555 0.9198298
## [3,]    2  5554.370 2206.295  793.8758 1.180825 54.88296 2.6729825
## [4,]    3  4139.542 2086.518 1323.5783 1.228914 51.90343 4.4564927
## [5,]    4  3085.103 1788.795 1776.2702 1.234610 44.49738 5.9807076
## [6,]    5  2299.255 1466.670 2131.7169 1.214742 36.48434 7.1774981
```
or

```r
x <- solve(mod1,theta, ev, inits)
```

```
## Warning in rxInits(object, inits, rxState(object), 0): Assiged depot,
## centr, peri to 0.
```

```r
head(x)
```

```
##      time     depot    centr      peri      eff       C2        C3
## [1,]    0 10000.000    0.000    0.0000 1.000000  0.00000 0.0000000
## [2,]    1  7452.765 1783.897  273.1895 1.084664 44.37555 0.9198298
## [3,]    2  5554.370 2206.295  793.8758 1.180825 54.88296 2.6729825
## [4,]    3  4139.542 2086.518 1323.5783 1.228914 51.90343 4.4564927
## [5,]    4  3085.103 1788.795 1776.2702 1.234610 44.49738 5.9807076
## [6,]    5  2299.255 1466.670 2131.7169 1.214742 36.48434 7.1774981
```

```r
par(mfrow=c(1,2))
matplot(x[,"C2"], type="l", ylab="Central Concentration")
matplot(x[,"eff"], type="l", ylab = "Effect")
```

![plot of chunk unnamed-chunk-15](vignettes/figure/unnamed-chunk-15-1.png)

#### Simulation of Variability with RxODE
Variability in model parameters can be simulated by creating a matrix
of parameter values for use in the simulation. In the example below,
40% variability in clearance is simulated. 

```r
options(RxODE.warn.on.assign = FALSE); ## Turn of warning on assignment of initial conditions
nsub <- 100						  #number of subproblems
CL <- 1.86E+01*exp(rnorm(nsub,0,.4^2))
theta.all <- 
	cbind(KA=2.94E-01, CL=CL, V2=4.02E+01,  # central 
	Q=1.05E+01, V3=2.97E+02,                # peripheral
	Kin=1, Kout=1, EC50=200)                # effects  
head(theta.all)
```

```
##         KA       CL   V2    Q  V3 Kin Kout EC50
## [1,] 0.294 20.28104 40.2 10.5 297   1    1  200
## [2,] 0.294 21.35152 40.2 10.5 297   1    1  200
## [3,] 0.294 21.62325 40.2 10.5 297   1    1  200
## [4,] 0.294 15.21664 40.2 10.5 297   1    1  200
## [5,] 0.294 15.93629 40.2 10.5 297   1    1  200
## [6,] 0.294 24.54465 40.2 10.5 297   1    1  200
```

Each subproblem can be simulated by using an explicit loop (or the `apply()`
function) to run the simulation for each set of parameters of in the parameter
matrix. 

```r
nobs <- ev$get.nobs()
set.seed(1)
cp.all <- matrix(NA, nobs, nsub)
for (i in 1:nsub)
{
	theta <- theta.all[i,]
	x <- mod1$solve(theta, ev, inits=inits)
	cp.all[, i] <- x[, "C2"]
}

matplot(cp.all, type="l", ylab="Central Concentration")
```

![plot of chunk unnamed-chunk-17](vignettes/figure/unnamed-chunk-17-1.png)

It is now straightforward to perform calculations and generate plots
with the simulated data. Below,  the 5th, 50th, and 95th percentiles
of the simulated data are plotted. 

```r
cp.q <- apply(cp.all, 1, quantile, prob = c(0.05, 0.50, 0.95))
matplot(t(cp.q), type="l", lty=c(2,1,2), col=c(2,1,2), ylab="Central Concentration")
```

![plot of chunk unnamed-chunk-18](vignettes/figure/unnamed-chunk-18-1.png)

#### Facilities for generating R shiny applications

An example of creating an R [shiny application](http://shiny.rstudio.com) to
interactively explore responses of various complex dosing regimens is available
at http://qsp.engr.uga.edu:3838/RxODE/RegimenSimulator.  Shiny applications
like this one may be programmatically created with the experimental function
`genShinyApp.template()`.

The above application includes widgets for varying the dose, dosing
regimen, dose cycle, and number of cycles.

```

genShinyApp.template(appDir = "shinyExample", verbose=TRUE)

library(shiny)
runApp("shinyExample")

```

### [Click here to go to the Shiny App](http://qsp.engr.uga.edu:3838/RxODE/RegimenSimulator)

