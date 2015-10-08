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
   c(desc = 'unexistent operator %%',
     code = paste(sep = "\n",
        'd/dt(y) = -ka;',
        'C1 = ka ** y;')
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
