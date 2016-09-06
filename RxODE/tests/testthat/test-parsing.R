## test ODE parsing for syntax errors
library("RxODE")
context("Test Parsing of models")

badParse <- function(desc,code){
    test_that(desc,{
        expect_error(RxODE(code));
    })
}

goodParse <- function(desc,code){
    test_that(desc,{
        rx <- RxODE(code);
        expect_equal(class(rx),"RxODE");
    })
}

equivSyntax <- function(desc,code1,code2){
    test_that(desc,{
        rx1 <- RxODE(code1);
        rx2 <- RxODE(code2);
        expect_equal(rxMd5(rx1)["parsed_md5"],rxMd5(rx2)["parsed_md5"])
    })
}

badParse('incorrect d/dt operator','d/dt(y = 1);')

## Statements don't require ; now.
goodParse('comments must be outside statements','d/dt(y) = 1   # bad comment;')
goodParse('missing end of statement ";" dosen\'t cause errors',
          paste(sep = "\n",
                'd/dt(depot) = -ka * depot',
                'd/dt(centr) =  ka * depot - kout * centr;'))

badParse('arithmetic syntax error',
         paste(sep = "\n",
               '# comment, just to show error in line 3',
               'd/dt(y) = -ka;',
               'C1 = /y;'))
## added ** operator
goodParse('existing operator **',
          code = paste(sep = "\n",
                       'd/dt(y) = -ka;',
                       'C1 = ka *  y**2;'))

badParse('unexistent operator %',
         code = paste(sep = "\n",
                      'remainder = 4 % 3;',
                      'd/dt(y) = -ka;',
                      'C1 = ka * y;')
         );

badParse(desc = 'incorrect "if" statement',
         code = paste(sep = "\n",
                      'if(comed==0){',
                      '   F = 1.0;',
                      'else {',         # missing "}"' 
                      '   F = 0.75;',
                      '};',
                      'd/dt(y) = F * y;')
         )

badParse(desc = 'illegal variable name (starting w. a digit)',
         code = paste(sep = "\n",
                      'F = 0.75;',
                      '12foo_bar = 1.0/2.0;',
                      'd/dt(y) = F * y;')
         )

badParse(desc = 'illegal variable name (illegal ".")',
         code = paste(sep = "\n",
                      'F = 0.75;',
                      'foo.bar = 1.0/2.0;',
                      'd/dt(y) = F * y;')
         )

badParse(desc = 'illegal variable name in d/dt()',
         code = paste(sep = "\n",
                      'd/dt(y_1) = F * y;',   # okay
                      'd/dt(y.1) = F * y;')   # not okay
         )

badParse(desc = 'Assignment with <<- not supported',
         'd/dt(y_1) <<- F*y')

goodParse(desc = 'Assignment with <- supported',
          'd/dt(y_1) <- F*y')

equivSyntax(desc = "time and t are equivalent",
            "d/dt(depot) = time^2","d/dt(depot) = t^2")

rxClean();
