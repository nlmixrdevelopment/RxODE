<img src="https://raw.githubusercontent.com/nlmixrdevelopment/RxODE/master/vignettes/logo.png" width="120" alt="RxODE logo (by Justin Wilkins)" style="position:absolute; top:0; right:0; padding:10px; border: 0;">

# RxODE

[![Build Status](https://travis-ci.org/nlmixrdevelopment/RxODE.svg?branch=master)](https://travis-ci.org/nlmixrdevelopment/RxODE)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/8vv1e3hncve9tnva?svg=true)](https://ci.appveyor.com/project/mattfidler/rxode)
[![codecov.io](https://codecov.io/github/nlmixrdevelopment/RxODE/coverage.svg?branch=master)](https://codecov.io/github/nlmixrdevelopment/RxODE?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/RxODE)](https://cran.r-project.org/package=RxODE)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/grand-total/RxODE)](https://cran.r-project.org/package=RxODE)
[![Dependency Status](https://dependencyci.com/github/nlmixrdevelopment/RxODE/badge)](https://dependencyci.com/github/nlmixrdevelopment/RxODE)

**RxODE** is an R package for solving and simulating from ode-based
models. These models are convert the RxODE mini-language to C and create
a compiled dll for fast solving.

## Basic Installation 

To install from source, you will have to have a working c and c++
compiler (either gcc, icc, clang) as well as a fortran compiler
(gfortran). 

To run RxODE, you need a working c compiler.  To use parallel threaded
solving in RxODE, this c compiler needs to support open-mp.

## Windows -- Setup the c/fortran compilers

These notes briefly describe steps to properly install `RxODE` and to
ensure `Rtools` (https://cran.r-project.org/bin/windows/Rtools/) are properly 
configured to avoid compilation issues during the use of `RxODE`. 

In a nutshell, installing `RxODE` is very straight forward, but
installing and configuring `Rtools` is a bit more delicate and you
need to carefully follow the instructions in the "R Installation and
Administration" manual, in particular Section 6.3, and Appendix D "The
Windows Toolset".  We point out a couple of details worth extra
attention.  Please read on.

### Steps:

1. Install the appropriate `Rtools` for your R for Windows version,
   (see http://cran.r-project.org/bin/windows/Rtools/). A couple of 
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

     If you get the error `gcc: error: CreateProcess: No such file or
     directory` then you know `Rtools` was not properly installed, in
     particular, it did not update your system `PATH` variable.

2.  Obtain the `RxODE` package, either from github or CRAN.  The 
    installation requires use of the gcc compiler, so you'll know if Step 1 
    was successfully executed.

    * CRAN. Use the usual method for installing packages from CRAN.

    * GitHub. First install the `devtools` package (if needed) and 
      then `RxODE` from GitHub.  You may want to avoid using a library 
      folder that has spaces in its name (see question 4.1 in the 
      "R for Windows FAQ" and the pointers therein).  As of `RxODE`
      version 0.5-1, we've been able to test installations on folder with 
      spaces in their name, but you may want to be on the safe side.
      
      ``` 
      install.packages("devtools")
      library("devtools", lib = "C:/Rlib")
      install_github("nlmixrdevelopment/RxODE")
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

## Mac OSX

Installation on a mac is much similar to RxODE installation under
windows.  To enable open mp on R and RxODE, you will need to install
the gfortran and clang compilers located at
https://cran.r-project.org/bin/macosx/tools/

## Linux

To install on linux make sure you install gcc (with openmp support)
and gfortran using your distribution's package manager.

