# test ODE parsing for syntax errors
library("RxODE")
library("tools") 

tmp <- tempdir()

# list of model specs with errors (test description and code)
errs <- list()

errs[[1]] <- 
   c(desc = 'incorrect d/dt operator',
     code = 'd/dt(y = 1);'
   )

errs[[2]] <- 
   c(desc = 'comments must be outside statements',
     code = 'd/dt(y) = 1   # bad comment;' 
   )

errs[[3]] <- 
   c(desc = 'missing end of statement ";"',
     code = paste(sep = "\n",
        'd/dt(depot) = -ka * depot',
        'd/dt(centr) =  ka * depot - kout * centr;')
   )

errs[[4]] <- 
   c(desc = 'arithmetic syntax error',
     code = paste(sep = "\n",
        '# comment, just to show error in line 3',
        'd/dt(y) = -ka;',
        'C1 = /y;')
   )

errs[[5]] <- 
   c(desc = 'unexistent operator **',
     code = paste(sep = "\n",
        'd/dt(y) = -ka;',
        'C1 = ka *  y**2;')
   )

errs[[6]] <- 
   c(desc = 'unexistent operator %',
     code = paste(sep = "\n",
        'remainder = 4 % 3;',
        'd/dt(y) = -ka;',
        'C1 = ka * y;')
   )

errs[[7]] <-
   c(desc = 'incorrect "if" statement',
     code = paste(sep = "\n",
        'if(comed==0){',
        '   F = 1.0;',
        'else {',         # missing "}"' 
        '   F = 0.75;',
        '};',
        'd/dt(y) = F * y;')
   )

errs[[8]] <-
   c(desc = 'illegal variable name (starting w. a digit)',
     code = paste(sep = "\n",
        'F = 0.75;',
        '12foo_bar = 1.0/2.0;',
        'd/dt(y) = F * y;')
   )

errs[[9]] <-
   c(desc = 'illegal variable name (illegal ".")',
     code = paste(sep = "\n",
        'F = 0.75;',
        'foo.bar = 1.0/2.0;',
        'd/dt(y) = F * y;')
   )

errs[[10]] <-
   c(desc = 'illegal variable name in d/dt()',
     code = paste(sep = "\n",
        'd/dt(y_1) = F * y;',   # okay
        'd/dt(y.1) = F * y;')   # not okay
   )

N <- length(errs)
for(i in 1:N){

   desc <- errs[[i]]["desc"]
   code <- errs[[i]]["code"]

   cat(sprintf('Syntax test %d of %d (%s)\n', i, N, desc))
   cat("==========================================================\n")
   cat("Input:\n", code, "\n", sep="")
   cat("\nRxODE message is:\n")

   assertError(RxODE(model = code, wd = tmp, modName=paste0("err",i)))

   cat("\n")
}

unlink(tmp, recursive = TRUE)
